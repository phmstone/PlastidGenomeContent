import os
import re
import argparse
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

###############################################################################
# This program:
# 1. Adds reference gene sequences to gene-specific alignment FASTA files
# 2. Parses BLAST tabular output files
# 3. Extracts matching gene regions from plastid genomes
# 4. Handles strand orientation correctly
# 5. Appends extracted sequences to alignment files
###############################################################################


# ---------------------------------------------------------------------------------------------------
# Command line inputs
# ---------------------------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract gene sequences from plastid genomes using BLAST results "
            "and append them to gene-specific alignment FASTA files."))

    parser.add_argument(
        "--blast-dir",
        required=True,
        help="Directory containing BLAST result subfolders")

    parser.add_argument(
        "--reference-dir",
        required=True,
        help="Directory containing reference gene FASTA files")

    parser.add_argument(
        "--genome-dir",
        required=True,
        help="Directory containing plastid genome FASTA files")

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for output alignment FASTA files")

    parser.add_argument(
        "--flank",
        type=int,
        default=0,
        help="Number of base pairs to extract on each side of BLAST hit (default: 0)")

    parser.add_argument(
        "--ir-cutoff",
        type=int,
        default=5000,
        help="Maximum expected gene length; larger hits are flagged (default: 5000)")

    return parser.parse_args()


args = parse_args()

# ---------------------------------------------------------------------------------------------------
# Gene list
# ---------------------------------------------------------------------------------------------------

# seperate the genes into an actual list and make everything lowercase for matching ease
geneList = "ndhA	ndhB	ndhC	ndhD	ndhE	ndhF	ndhG	ndhH	ndhI	ndhJ	ndhK	ccsA	cemA	petA	petB	petD	petG	petL	petN	psaA	psaB	psaC	psaI	psaJ	psbA	psbB	psbC	psbD	psbE	psbF	psbH	psbI	psbJ	psbK	psbL	psbM	psbN	psbT	psbZ	rbcL	ycf3	ycf4	rpoA	rpoB	rpoC1	rpoC2	atpA	atpB	atpE	atpF	atpH	atpI	infA	rpl2	rpl14	rpl16	rpl20	rpl22	rpl23	rpl32	rpl33	rpl36	rps2	rps3	rps4	rps7	rps8	rps11	rps12	rps14	rps15	rps16	rps18	rps19	accD	clpP	matK	ycf1	ycf2	rrn4.5	rrn5	rrn16	rrn23	trnA-UGC	trnC-GCA	trnD-GUC	trnE-UUC	trnF-GAA	trnfM-CAU	trnG-GCC	trnG-UCC	trnH-GUG	trnI-CAU	trnI-GAU	trnK-UUU	trnL-CAA	trnL-UAA	trnL-UAG	trnM-CAU	trnN-GUU	trnP-UGG	trnQ-UUG	trnR-ACG	trnR-UCU	trnS-GCU	trnS-GGA	trnS-UGA	trnT-GGU	trnT-UGU	trnV-GAC	trnV-UAC	trnW-CCA	trnY-GUA"
geneList = geneList.lower()
geneList = geneList.split('\t')


# -----------------------------------------------------------------------------------------------------------------------------------
# build a mapping of genome accessions to FASTA file paths because decimals are sometimes stripped from genbank IDs
# -----------------------------------------------------------------------------------------------------------------------------------
genome_files = {}
for f in os.listdir(args.genome_dir):
    if f.endswith(".fasta"):
        name_no_ext = os.path.splitext(f)[0]        # removes .fasta
        genome_files[name_no_ext] = os.path.join(args.genome_dir, f)


# -----------------------------------------------------------------------------------------------------------------------------------
# Add one copy of each reference sequence to the alignment file
# -----------------------------------------------------------------------------------------------------------------------------------

# make the output directory if it does not already exist
os.makedirs(args.output_dir, exist_ok=True)

# iterate over the genes
for gene in geneList:
    alignment_path = os.path.join(
        args.output_dir, f"{gene}-alignment-unaligned.fasta")
    
    # path to the reference multifasta for this gene
    reference_fasta = os.path.join(args.reference_dir, f"{gene}.fasta")
 
    # write out the sequences from the reference multifasta to the new file   
    # only add reference sequences if the alignment file does not already exist
    if os.path.exists(reference_fasta) and not os.path.exists(alignment_path):
        with open(alignment_path, "w") as alignment_file:
            for record in SeqIO.parse(reference_fasta, "fasta"):
                SeqIO.write(record, alignment_file, "fasta")


# -----------------------------------------------------------------------------------------------------------------------------------
# Process blast results
# -----------------------------------------------------------------------------------------------------------------------------------

# keep track of files with no hits or issues (likely IR problems)
bigNoHitsList = []
bigProblemFilesList = []

# list of folders inside blast results directory
blast_folders = [
    os.path.join(args.blast_dir, d)
    for d in os.listdir(args.blast_dir)
    if os.path.isdir(os.path.join(args.blast_dir, d))
]

# loop through each plastid genome BLAST folder
for directory in blast_folders:

    somethingWrongWithTheseFiles = []
    noHitsFiles = []

    # loop through each gene BLAST file
    for filename in os.listdir(directory):
        if not filename.endswith(".txt"):
            continue

        blastResults = []

        # read BLAST output
        with open(os.path.join(directory, filename)) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                blastResults.append(line)

        # no hits in this file
        if not blastResults:
            noHitsFiles.append(f"{genome_files}/{filename}")
            continue

        # split blast hits into columns
        tabularResults = [line.split("\t") for line in blastResults]

        # extract gene name from filename
        geneName = os.path.splitext(filename)[0].lower()

        # extract IDs from first hit
        refSeqID = tabularResults[0][0]   # qseqid
        genomeID = tabularResults[0][1]   # sseqid

        # clean subject IDs from blast outfmt6 (gb|XXXX.1|)
        if "|" in genomeID:
            genomeID = genomeID.split("|")[-2]
        # strip version number if present
        genomeID = genomeID.split(".")[0]

        # lists to store coordinates and strand info
        geneBounds = []
        genomeBounds = []
        strands = []

        # parse each BLAST hit
        for row in tabularResults:
            q_start = int(row[3])
            q_end   = int(row[4])
            s_start = int(row[5])
            s_end   = int(row[6])

            geneBounds.append([q_start, q_end])
            genomeBounds.append([s_start, s_end])

            strand = "+" if s_start < s_end else "-"
            strands.append(strand)

        # -------------------------------------------------
        # deal with duplicated hits
        # -------------------------------------------------

        sorter = {}
        for q, s, strand in zip(geneBounds, genomeBounds, strands):
            q_tuple = tuple(q)
            if q_tuple not in sorter:
                sorter[q_tuple] = [(s, strand)]
            else:
                sorter[q_tuple].append((s, strand))

        duplicateGeneBounds = [[]]
        duplicateGenomeBounds = [[]]

        for k, v in sorter.items():
            for i, (coords, strand) in enumerate(v):
                if i >= len(duplicateGeneBounds):
                    duplicateGeneBounds.append([])
                    duplicateGenomeBounds.append([])
                duplicateGeneBounds[i].append(list(k))
                duplicateGenomeBounds[i].append((coords, strand))

        genomeBoundList = []

        for numberSet in range(len(duplicateGeneBounds)):
            gene_array = np.array(duplicateGeneBounds[numberSet])
            genome_array = np.array([x[0] for x in duplicateGenomeBounds[numberSet]])
            strand_list = [x[1] for x in duplicateGenomeBounds[numberSet]]

            min_idx = np.unravel_index(np.argmin(gene_array), gene_array.shape)
            max_idx = np.unravel_index(np.argmax(gene_array), gene_array.shape)

            s1 = genome_array[min_idx]
            s2 = genome_array[max_idx]

            start = min(s1, s2)
            end = max(s1, s2)
            strand = strand_list[min_idx[0]]

            genomeBoundList.append((start, end, strand))

            if end - start > args.ir_cutoff:
                somethingWrongWithTheseFiles.append(f"{genome_files}/{filename}")

        # -------------------------------------------------
        # read test plastid genome FASTA
        # -------------------------------------------------

        if genomeID not in genome_files:
            print(f"WARNING: Genome file for {genomeID} not found in {args.genome_dir}. Skipping.")
            continue

        genome_record = SeqIO.read(genome_files[genomeID], "fasta")
        genome_seq = genome_record.seq

        # -------------------------------------------------
        # extract sequences and append to alignment files
        # -------------------------------------------------

        for start, end, strand in genomeBoundList:
            start0 = max(start - args.flank - 1, 0)
            end0 = min(end + args.flank, len(genome_seq))

            subseq = genome_seq[start0:end0]

            if strand == "-":
                subseq = subseq.reverse_complement()

            header = (
                f"{genome_record.id}|{geneName}|"
                f"{start}-{end}|refSeq:{refSeqID}"
            )

            new_record = SeqRecord(subseq, id=header, description="")

            aln_path = os.path.join(
                args.output_dir, f"{geneName}-alignment-unaligned.fasta"
            )

            with open(aln_path, "a") as out:
                SeqIO.write(new_record, out, "fasta")

    bigProblemFilesList.append(somethingWrongWithTheseFiles)
    bigNoHitsList.append(noHitsFiles)


# -----------------------------------------------------------------------------------------------------------------------
# write summary files for manual inspection
# -----------------------------------------------------------------------------------------------------------------------
problem_file_path = os.path.join(args.output_dir, "FilesToCheckAgain.txt")
nohits_file_path = os.path.join(args.output_dir, "NoHitsFiles.txt")

with open(problem_file_path, "w") as out:
    for f in sorted(set(sum(bigProblemFilesList, []))):
        out.write(f + "\n")

with open(nohits_file_path, "w") as out:
    for f in sorted(set(sum(bigNoHitsList, []))):
        out.write(f + "\n")

