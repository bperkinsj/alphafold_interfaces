import argparse
import utils
import os
import pandas as pd

def get_filenames(args):
    '''
    Add the filenames to the input dataframe
    '''

    filenames = os.listdir(args.structures)

    # Add the filenames to the dataframe
    df = pd.read_csv(args.file)

    df['filename'] = ''

    for index, row in df.iterrows():
        uniprot = row['uniprot']
        for fn in filenames:
            if uniprot in fn:
                df.at[index, 'filename'] = fn
    
    df.to_csv(args.file, index=False)


def get_interfaces(args):
    '''
    Take a list of AlphaFold2 structures and return a file detailing the interface contacts between
    the functional domain and the inhibitory module.
    '''

    # Make the output directory
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Read the input file
    df = pd.read_csv(args.file)

    # Define columns for new stats
    df['pdb_mutations'] = ''
    df['interacting_residue_pairs'] = ''
    df['interface_residues'] = ''
    df['number_interface_residues'] = ''

    df = utils.region_search_range(df)

    for index, row in df.iterrows():

        # Define values for retrieval 
        region_1_res = row['region_1 search']
        region_2_res = row['region_2 search']
        uniprot = row['uniprot']
        chain = 'A'
        model = 0
        fn = row['filename']

    # Get structure
    structure = utils.get_pdb_struct_dict(uniprot, fn, args.structures)

    atoms_ns = utils.get_domain_residues(region_1_res, region_2_res, structure, model, chain)

    # Get interacting residues
    interacting_pairs, interface_res, len_interface_res = utils.domain_neighborsearch(region_1_res, region_2_res, atoms_ns)

    df.loc[i, 'interacting_residue_pairs'] = interacting_pairs
    df.loc[i, 'interface_residues'] = interface_res
    df.loc[i, 'number_interface_residues'] = len_interface_res

    # Save the dataframe
    df.to_csv(os.path.join(args.output, args.interface_file), index=False)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--file', required=True, type=str, help='Path to the input file')
    parser.add_argument('-s', '--structures', required=True, type=str, help='Path to the directory containing the AlphaFold2 structures',
                        default='data/pdb')
    parser.add_argument('-o', '--output', required=True, type=str, help='Path to the output directory', default='data/output')
    parser.add_argument('-i', '--interface_file', required=True, type=str, help='Interface filename',
                        default='interface_analysis.csv')
    args = parser.parse_args()

    get_filenames(args)


if __name__ == '__main__':
    main()