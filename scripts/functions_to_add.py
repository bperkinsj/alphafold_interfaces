def mean_plddt(df, path, fnt='af_filename'):
    # Calculate mean plDDT for our region of interest.

    print('Calculating mean plDDT...')

    # Turn region ranges into list of residues
    df = utils.region_search_range(df).reset_index(drop=True)

    for i in range(len(df)):
        fn = df.loc[i, fnt] # Either af_filename or cf_filename
        if fnt == 'cf_filename':
            uniprot = df.loc[i, 'uniprot']
            fp = join(path, uniprot, fn)
        else:
            fp = join(path, fn)

        region_1_range = df.loc[i, 'region_1 search']
        region_2_range = df.loc[i, 'region_2 search']

        if '.cif' in fn:
            fp = utils.cif_to_pdb(fp)

        # Convert to pandas pdb object
        ppdb = PandasPdb().read_pdb(fp)
        protein = ppdb.df['ATOM']

        # Get average pLDDT for entire protein
        complex_mean = protein['b_factor'].mean()

        # Get average pLDDT for regions 1 and 2
        r1 = protein[protein['residue_number'].isin(region_1_range)]
        r2 = protein[protein['residue_number'].isin(region_2_range)]
        r1_mean = r1['b_factor'].mean()
        r2_mean = r2['b_factor'].mean()

        df.loc[i, 'complex_mean_plddt'] = round(complex_mean, 3)
        df.loc[i, 'r1_mean_plddt'] = round(r1_mean, 3)
        df.loc[i, 'r2_mean_plddt'] = round(r2_mean, 3)

    df.drop(columns=['region_1 search', 'region_2 search'], inplace=True)
    return df


def mean_paes(df, path, affix, suffix):
    # Calculate the average pae for region 1 to region 1, region 2 to region 2, and region 1 to region 2

    print('Calculating mean pae...')

    for i in range(len(df)):
        uniprot = df.loc[i, 'uniprot']
        fn = affix + uniprot + suffix
        region_1 = df.loc[i, 'region_1']
        region_2 = df.loc[i, 'region_2']

        # Region bounds are in the format [start, end] for each region. Regions with multiple sections look like [[start, end], [start, end], ...]
        reg1_bounds = utils.region_bounds(region_1)
        reg2_bounds = utils.region_bounds(region_2)

        reg1_array = np.array(reg1_bounds)
        reg2_array = np.array(reg2_bounds)

        # Read in json file
        prot_array = utils.pae_from_json(path, fn)

        '''
        We want means of reg1 compared against reg1, reg1 compared against reg2, and reg2 compared against reg2.
        '''

        if prot_array.any() == np.nan:
            mean11 = 0
            mean12 = 0
            mean22 = 0

        else:
            mean11 = utils.calculate_pae_mean(prot_array, reg1_array, reg1_array)
            mean12 = utils.calculate_pae_mean(prot_array, reg1_array, reg2_array)
            mean22 = utils.calculate_pae_mean(prot_array, reg2_array, reg2_array)
        

        df.loc[i, 'mean_pae_1_1'] = round(mean11, 3)
        df.loc[i, 'mean_pae_1_2'] = round(mean12, 3)
        df.loc[i, 'mean_pae_2_2'] = round(mean22, 3)
    
    return df

def region_bounds(x):
    """
    This function takes in a `string` representing a region of interest in a
    protein. The region of interest can be a single region or multiple regions
    of a protein. Returns the bounds of the regions in list form.
    
    Parameters:
    
        x (string): String containing a region or several regions of interest in a 
            protein.
            Format of x: single region -> 'start-end'
                         multiple regions -> 'start1-end1,start2-end2'
                     
    Returns:
    
        region boundaries in list form

            Format: single region -> [start, end+1]
                    multiple region -> [[start1, end1+1], [start2, end2+1]]
    """
    # Handle instances with more than one range
    if ',' in x:
        list_temp = x.split(sep = ',') #list_temp = ['123-456,' '789-1111']
        for y in range(len(list_temp)): 
            list_temp[y] = list_temp[y].split(sep = '-') #list_temp[y] = [['123', '456'], ['789', '1111']]
        for y in range(len(list_temp)): 
            for x in range(len(list_temp[y])):
                list_temp[y][x] = int(list_temp[y][x]) #turns each list item into an integer

        # Make a range object with the bounds of the range. Note to the 
        # end a 1 has to be added in order to include the last position in the range
        for y in range(len(list_temp)): #[1, 2] where 1=[123, 456] and 2=[789, 1111]
            for x in range(len(list_temp[y])): #[123, 456]       
                list_temp[y] = [list_temp[y][x], list_temp[y][x+1]+1] #list_temp[0][0] = [123], list_temp[0][0+1]+1 or [456] + 1 = [457]
                break

        return list_temp

    # Handle instances with only one range
    else:
        list_temp = x.split(sep = '-')
        for y in range(len(list_temp)):
            list_temp[y] = int(list_temp[y]) #

        # Make a range object with the bounds of the region. Note to the 
        # end a 1 has to be added in order to include the last position in the range
        return [[list_temp[0], list_temp[1]+1]]

def pae_from_json(path, fn):
    '''Read in the json file, which is in the format:
    [{"predicted_aligned_error":[[0, 1, 3, 5, 19, ...], [0, 4, 12, 38, ...], ...]}]
    '''

    try:
        f = open(path + fn)
        data = json.load(f)
        data = data[0]
        pae = data['predicted_aligned_error']
        array = np.array(pae)
    except FileNotFoundError:
        print(f'File {fn} not found')
        array = np.empty((1,1,))
        array[:] = np.nan

    return array

def calculate_pae_mean(prot_array, reg_a, reg_b):
    '''
    Gives the mean pae for all regions of interest compared against all regions of interest (reg1 to reg1, reg1 to reg2, reg2 to reg2)
    Reg_a and Reg_b are given as arrays.
    '''

    '''
    Method proceeds like this:
    Let's say we're handed reg1 from a protein, which is "1-2, 3-4". It's given to us in the form:
    reg_a = [[1, 2], [3, 4]], reg_b = [[1, 2], [3, 4]]]
    So we perform:
    mean(prot_array[1:3, 1:3]) + mean(prot_array[1:3, 3:5]) + mean(prot_array[3:5, 1:3]) + mean(prot_array[3:5, 3:5]). Then take mean of all that.
    '''

    # First is reg1 on reg1
    means = []
    for i in range(len(reg_a)):
        for n in range(len(reg_b)):
            a_start = reg_a[i][0]
            a_end = reg_a[i][1]
            b_start = reg_b[n][0]
            b_end = reg_b[n][1]
            sub_array = prot_array[a_start:a_end, b_start:b_end]
            sub_mean = np.mean(sub_array)
            means.append(sub_mean)

    mean = np.mean(means)

    return (mean)