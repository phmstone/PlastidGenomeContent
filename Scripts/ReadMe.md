# Scripts

## fetchGenbank.py
This script takes an input of a list of genbank IDs in a text file and makes an output file of one big genbank file. At the moment that big genbank file goes into the `Outputs` directory within `Scripts`.
The file name of the output file has to be given as an argument on the command line. An email address associated with an account on genbank/NCBI must also be given to run the script.

## presenceAbsence.py
This script converts the annotations present in genbank files to a TSV with a column for each gene and a row for each taxon. The value entered represents the status of the gene in that organism (0 = present, 1 = missing, 2 = pseudogene)

## presentGeneMultiFasta.py
This script makes multi fastas by gene for all chloroplast genes (could be manually changed for desired genes of interest) where the gene is annotated as being present in the genbank file.
It creates a directory called `PresentGeneMultiFastas` where these multifasta files are kept, but the name of this directory can be chosen using `--outdir` on the command line

## presentGeneMultiFasta.py
This script uses the .TSV generated from presentGeneMultiFasta.py and a list genbank IDs to use as plastid gene reference sequences input by the user.
The presence/absence profiles from the TSV are read in by gene for each species. Species from this TSV that have all genes present will be added to the reference species list.
The reference species list should ideally contain genbank IDs of whole plastid sequences of species closely to the taxa of interest that contain all plastid genes of interest.
For each species in the TSV, blast searches are performed on genes that are categorised as either missing or pseudogenised. The blast searches are done with a reference sequence from the user's input reference species list.
The script outputs a directory `Blast` with four subdirectories; `Databases`, `PlastidSequences`, `ReferenceGeneSequences`, `ReferenceGenomes`, and `Results`.
* `Databases` database files for every genbank accession with missing genes or pseudogenes that are needed for blast
* `PlastidSequences` fasta files for every genbank accession with missing or pseudogenes
* `ReferenceGeneSequences` A fasta file for every gene of interest created from the first reference sequence
* `ReferenceGenomes` Genbank files for all reference sequences (including those added after checking presence/absence profiles)
* `Results` Text files from the output of each blast search performed to look for missing or pseudo genes using a reference sequence