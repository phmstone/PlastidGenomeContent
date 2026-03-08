
# ----------------------------------------------------------------------------------------------------------------
# import packages
# ----------------------------------------------------------------------------------------------------------------
import argparse
from pathlib import Path
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ----------------------------------------------------------------------------------------------------------------
# External tools
# ----------------------------------------------------------------------------------------------------------------
MAFFT = "mafft"  # must be in PATH

# ----------------------------------------------------------------------------------------------------------------
# Command line arguments
# ----------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Align reference and test sequences")
parser.add_argument("--input", required=True, help="Input directory with unaligned multifasta files ending in '-unaligned.fa'")
parser.add_argument("--output", required=True, help="Output directory for alignments")
args = parser.parse_args()

input_dir = Path(args.input)
output_dir = Path(args.output)
output_dir.mkdir(exist_ok=True)

# ----------------------------------------------------------------------------------------------------------------
# Process each unaligned fasta
# ----------------------------------------------------------------------------------------------------------------
for fasta_file in input_dir.glob("*-unaligned.fasta"):
    # derive gene/sample name by stripping "-unaligned"
    gene_name = fasta_file.stem.replace("-unaligned", "")
    print(f"[INFO] Processing {gene_name} from {fasta_file}")

    # create subdirectory for this sample/gene
    gene_dir = output_dir / gene_name
    gene_dir.mkdir(exist_ok=True)

    # read sequences
    records = list(SeqIO.parse(fasta_file, "fasta"))
    ref_records = [r for r in records if "refSeq" not in r.id]
    test_records = [r for r in records if "refSeq" in r.id]
    print(f"[INFO] Found {len(ref_records)} reference sequences and {len(test_records)} test sequences")

    # ----------------------------------------------------------------------------------------------------------------
    # Skip protein coding alignment steps for tRNA or rRNA genes
    # ----------------------------------------------------------------------------------------------------------------
    is_non_coding = gene_name.startswith("trn") or gene_name.startswith("rrn")
    if not is_non_coding:
        # ----------------------------------------------------------------------------------------------------------------
        # Translate sequences to proteins
        # ----------------------------------------------------------------------------------------------------------------
        protein_records = []
        for rec in records:
            seq = rec.seq
            protein_seq = seq.translate(to_stop=False)
            protein_records.append(SeqRecord(protein_seq, id=rec.id, description=""))

        # write unaligned proteins
        prot_fa = gene_dir / f"{gene_name}.proteins.fa"
        SeqIO.write(protein_records, prot_fa, "fasta")
        print(f"[INFO] Written unaligned protein sequences to {prot_fa}")

        # ----------------------------------------------------------------------------------------------------------------
        # Align proteins with MAFFT
        # ----------------------------------------------------------------------------------------------------------------
        prot_aln = gene_dir / f"{gene_name}.proteins.aln.fa"
        print(f"[INFO] Running MAFFT for protein alignment...")
        with prot_aln.open("w") as aln_out:
            subprocess.run([MAFFT, "--auto", str(prot_fa)], stdout=aln_out, check=True)
        print(f"[INFO] Protein alignment saved to {prot_aln}")

        # ----------------------------------------------------------------------------------------------------------------
        # Back-translate to nucleotide codon alignment
        # ----------------------------------------------------------------------------------------------------------------
        aln_records = list(SeqIO.parse(prot_aln, "fasta"))
        codon_records = []

        for rec in aln_records:
            orig_seq = next(r.seq for r in records if r.id == rec.id)
            codons = [orig_seq[i:i+3] for i in range(0, len(orig_seq), 3)]
            codon_i = 0
            nt_aln = ""
            for aa in rec.seq:
                if aa == "-":
                    nt_aln += "---"
                else:
                    if codon_i < len(codons):
                        nt_aln += codons[codon_i]
                        codon_i += 1
                    else:
                        nt_aln += "---"
            codon_records.append(SeqRecord(Seq(nt_aln), id=rec.id, description=""))

        codon_aln_file = gene_dir / f"{gene_name}.codon.aln.fa"
        SeqIO.write(codon_records, codon_aln_file, "fasta")
        print(f"[INFO] Codon-aware nucleotide alignment saved to {codon_aln_file}")

        # ----------------------------------------------------------------------------------------------------------------
        # Reference-only alignments
        # ----------------------------------------------------------------------------------------------------------------
        ref_prot_aln_file = gene_dir / f"{gene_name}.ref_only.proteins.aln.fa"
        SeqIO.write([r for r in aln_records if r.id in [x.id for x in ref_records]],
                    ref_prot_aln_file, "fasta")
        print(f"[INFO] Reference-only protein alignment saved to {ref_prot_aln_file}")

        ref_codon_aln_file = gene_dir / f"{gene_name}.ref_only.codon.aln.fa"
        SeqIO.write([r for r in codon_records if r.id in [x.id for x in ref_records]],
                    ref_codon_aln_file, "fasta")
        print(f"[INFO] Reference-only codon alignment saved to {ref_codon_aln_file}")
    else:
        print(f"[INFO] Skipping protein/codon alignment for non-coding gene {gene_name}")

    # ----------------------------
    # Always create plain nucleotide alignment with MAFFT
    # ----------------------------
    nt_aln_input = fasta_file
    nt_aln_output = gene_dir / f"{gene_name}.nucleotide.aln.fa"
    print(f"[INFO] Running MAFFT for plain nucleotide alignment for {gene_name}...")
    with nt_aln_output.open("w") as out_handle:
        subprocess.run([MAFFT, "--auto", str(nt_aln_input)], stdout=out_handle, check=True)
    print(f"[INFO] Nucleotide alignment saved to {nt_aln_output}")

print("[INFO] All alignments completed successfully!")