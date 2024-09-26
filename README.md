# alphafold_interfaces

A small repository for analyzing the interface residues of AlphaFold2-predicted structures.

This document contains the guide to collecting, creating, and analyzing the interfaces of AlphaFold2 structural predictions.

## Cloning this repository
```
cd ~/Documents/  # or any other base directory you want
git clone https://github.com/bperkinsj/alphafold_interfaces.git
cd alphafold_interfaces/
```

## Setting up your virtual environment
First, you'll have to (install Anaconda)[https://conda.io/projects/conda/en/latest/user-guide/install/index.html] if you don't already have it.
Navigate to the folder containing this repository.

```
$cd ~/some/directory/path/alphafold_interfaces
```

Create a new environment using the environment.yml file.

```
$ conda env create -f environment.yml
```

Activate the environment

```
$ conda activate interfaces
```

## Downloading AlphaFold2 structures

The list of proteins should be in the format:

uniprot,region_1,region_2
P28482,1-41,223-409
P62826,47-79,234-798,
...

where region_1 represents the inhibitory module and region_2 represents the functional domain.

Once you have the list of proteins of interest, run scripts/download_alphafold.py, making sure to pass it the file.

```
$ python scripts/download_alphafold.py -f data/proteins_of_interest.csv
```

It will download the appropriate AlphaFold2 structures into the data folder and create a file with the proteins it did not find an AlphaFold2 structure for.

## Creating AlphaFold2 structures
For each protein that requires a novel AlphaFold2 structure, you will first need to generate a multiple sequence alignment (MSA) for it in the format of an .a3m file, which you will do using the (ColabFold notebook)[https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=G4yBrceuFbf3].

First, navigate to [uniprot.org] and search using the UniProt accession ID of your protein. Copy the sequence and paste it into the ```query_sequence``` field in the ColabFold notebook (I recommend also changing the ```jobname``` to the UniProt ID). Run the notebook. Click the folder icon to the left. A folder should appear with the jobname. As the MSA is one of the first steps that the notebook runs, the .a3m file should appear relatively quickly and you're free to halt the notebook once it does. I recommend starting the notebook and then opening a new window for the next protein so you don't have to wait for each one to finish. Repeat this process for all your proteins.

You will now be using these a3m files to create the structures.

## Creating the Predicted Structures with AlphaFold2
On Alliance Canada, follow the instructions [here](https://github.com/alirezaomidi/colabfold-scripts) to set up ColabFold. Place all of your MSAs into a single folder and provide the path to the folder when you run single-gpu.sh instead of to the individual files.

## Analyzing the Interfaces
