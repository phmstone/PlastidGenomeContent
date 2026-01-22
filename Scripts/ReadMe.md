# Scripts

## fetchGenbank.py
This script takes an input of a list of genbank IDs in a text file and makes an output file of one big genbank file. At the moment that big genbank file goes into the `Outputs` directory within `Scripts`. The file name of the output file has to be given as an argument on the command line. An email address associated with an account on genbank/NCBI must also be given to run the script.

## presenceAbsence.py
This script converts the annotations present in genbank files to a TSV with a column for each gene and a row for each taxon. The value entered represents the status of the gene in that organism (0 = present, 1 = missing, 2 = pseudogene)

