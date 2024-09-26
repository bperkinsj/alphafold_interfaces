import argparse
import utils
import os
import pandas as pd

def download_alphafold(args):
    '''
    Take a list of UniProt accessions and download the AlphaFold2 structures from the AlphaFold Protein Structure Database.
    Then return a list of UniProt accessions that were not downloaded.
    '''

    # Read the input file
    df = pd.read_csv(args.file)

    print(df.head())

    # Get the UniProt accessions
    uniprots = df['uniprot'].tolist()

    # Download the structures
    utils.download_alphafold_structures(uniprots, args.output)

    # Determine what structures were not downloaded
    missing = utils.determine_missing(uniprots, args.output)

    # Save the missing structures
    missing_df = df[df['uniprot'].isin(missing)]
    missing_df.to_csv(os.path.join(args.output, 'missing_structures.csv'), index=False)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--file', required=True, type=str, help='Path to the input file')
    parser.add_argument('-o', '--output', type=str, help='Path to the output directory', default='data')

    args = parser.parse_args()

    download_alphafold(args)

if __name__ == '__main__':
    main()