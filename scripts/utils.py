import os
from alphafetcher import AlphaFetcher

def download_alphafold_structures(uniprot_accs, output_dir):
    '''
    Download the AlphaFold2 structures from the AlphaFold Protein Structure Database for
    the given UniProt accessions.
    '''
    fetcher = AlphaFetcher(base_savedir=output_dir)

    # Add UniProt access codes
    fetcher.add_proteins(uniprot_accs)

    # Retrieve metadata
    fetcher.fetch_metadata(multithread=True, workers=4)

    # Download the structures
    fetcher.download_all_files(pdb=True, cif=False, multithread=True, workers=4)

def determine_missing(uniprots, output_dir):
    '''
    Determine which structures were not downloaded. Takes a list of UniProt accessions and the output directory
    and returns a list of UniProt accessions that were not downloaded.
    '''

    struct_dir = os.path.join(output_dir, 'pdb_files')
    structures = os.listdir(struct_dir)
    downloaded_uniprots = [s.split('.')[0] for s in structures]
    missing_uniprots = [u for u in uniprots if u not in downloaded_uniprots]

    return missing_uniprots

def string2range(x):
    
    """
    This function takes in a `string` representing a region of interest in a
    protein. The region of interest can be a single region or multiple regions
    of a protein. Returns a range for single regions or a list of ranges for
    multiple regions.
    
    Parameters:
    
        x (string): String containing a region or several regions of interest in a 
            protein.
            Format of x: single region -> 'start-end'
                         multiple regions -> 'start1-end1,start2-end2'
                     
    Returns:
    
        range or list of ranges: For single region proteins a range is returned. For 
            multiple region proteins a list of ranges is returned

            Format: single region -> range(start, end+1)
                    multiple region -> [range(start1, end1+1), range(start2, end2+1)]
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
                list_temp[y] = list(range(list_temp[y][x], list_temp[y][x+1]+1)) #list_temp[0][0] = [123], list_temp[0][0+1]+1 or [456] + 1 = [457]
                break

        return list(set([item for sublist in list_temp for item in sublist]))

    # Handle instances with only one range
    else:
        list_temp = x.split(sep = '-')
        for y in range(len(list_temp)):
            list_temp[y] = int(list_temp[y]) #

        # Make a range object with the bounds of the region. Note to the 
        # end a 1 has to be added in order to include the last position in the range
        return list(range(list_temp[0], list_temp[1]+1))

def region_search_range(df):
    # Convert the domain region strings to ranges or lists of ranges
    df['region_1 search'] = df['region_1'].apply(lambda x: string2range(x))
    df['region_2 search'] = df['region_2'].apply(lambda x: string2range(x))

    return df

def get_pdb_struct_dict(name, fn, path):
    # Join the path and the file name
    full_path = os.path.join(path, fn)

    # To load a PDB file make a parser object
    parser = PDBParser(QUIET=True)
            
    # Then make a structure object
    structure = parser.get_structure(name, full_path)

    return structure

def get_domain_residues(region1, region2, structure, label_model, label_chain):
    # Iterate through all the models in the structure
    for model in structure:
        
        # Analyze only the model that corresponds to the current row in df_prot
        if model.get_id() == label_model:
        
            for chain in model:
                
                # Analyze only the chain that corresponds to the current row in df_prot
                if chain.get_id() == label_chain:

                    print(f'Looking at chain {label_chain} for interface')
                    # Get all the atoms in the chain
                    atom_list = Selection.unfold_entities(chain, "A")
                    
                    # Make an empty list to store the atoms in the region_1 and in region_2
                    atoms_ns = []
                    
                    # Iterate through all the atoms in the chain and append the ones that occur inside
                    # the region_1 or the region_2 to the list
                    for atom in atom_list:
                        
                        # Get the parent residue for the atom
                        res = atom.get_parent()
                        
                        # Make sure the residue is an amino acid
                        if res.get_id()[0] == ' ':
                            
                            # Check whether the residue lies inside the region_1 or the region_2 and append
                            # the atoms of these residues into a list
                            if res.get_id()[1] in region1:
                                atoms_ns.append(atom)
                                
                            elif res.get_id()[1] in region2:
                                atoms_ns.append(atom)

                    return atoms_ns

def domain_neighborsearch(region1, region2, atoms):
    # Make an NeighborSearch object with all the atoms inside the region_1 and the region_2
    ns = NeighborSearch(atoms)
      
    # Search for all the interacting residues in the region_1 and in the region_2
    # with atoms that are within a 5 A radius 
    ns_all = ns.search_all(5, 'R')
    
    # Make a set to store the residues at the interface
    interface_res = set()
    
    # Save the interacting residue pairs as a list of tuples
    interacting_pairs = []
    
    # Iterate thorugh all the interacting residue pairs and append those that have a residue
    # in the region_1 and another in the region_2 to a list. Save the residue positions in a set
    for pairs in ns_all:
        
        res_0 = pairs[0].get_id()[1]
        res_1 = pairs[1].get_id()[1]
        
        if res_0 in region1 and res_1 in region2:
            interface_res.add(res_0)
            interface_res.add(res_1)
            interacting_pairs.append((res_0, res_1))
            
        elif res_1 in region1 and res_0 in region2:
            interface_res.add(res_0)
            interface_res.add(res_1)
            interacting_pairs.append((res_0, res_1))
            
    # Save the results in the appropriate columns of df_prot
    if len(interface_res) > 0 and len(interacting_pairs) > 0:
        return str(interacting_pairs), str(interface_res), len(interface_res)  
        
    else: 
        return np.nan, np.nan, np.nan