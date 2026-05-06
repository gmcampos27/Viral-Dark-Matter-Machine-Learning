import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os
from collections import Counter
from itertools import product


base_path = "results_annotation"
fasta_path = "final_contigs_filtered"
pfam_path = os.path.join(base_path, "pfam")
output_final = "ML_features_matrix.csv"

# 256 Tetranucleotides
kmer_list = [''.join(i) for i in product('ACGT', repeat=4)]

def fast_tetranucleotide(sequence):
    seq = str(sequence).upper()
    kmers = [seq[i:i+4] for i in range(len(seq) - 3)]
    counts = Counter(kmers)
    total = sum(counts.values())
    return {k: counts.get(k, 0) / total for k in kmer_list} if total > 0 else {k: 0 for k in kmer_list}

all_samples_data = []

# HMMER
colunas_hmmer = [
    "Query_ID", "t_accession", "tlen", "HMM_ID", "q_accession", "qlen", 
    "Full_Evalue", "Full_Score", "Full_Bias", "Domain_Num", "Total_Domains", 
    "c_Evalue", "i_Evalue", "Domain_Score", "Domain_Bias", 
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc", "description"
]

print("Features Extraction...")
fasta_files = [f for f in os.listdir(fasta_path) if f.endswith(".fasta")]

for f_name in fasta_files:
    sample_name = f_name.replace("_scaffolds_min500.fasta", "")
    print(f"Sample: {sample_name}")

    # GC + K-mers
    features_list = []
    for record in SeqIO.parse(os.path.join(fasta_path, f_name), "fasta"):
        row = {
            "Sample": sample_name,
            "Node": record.id,
            "Length": len(record.seq),
            "GC_Content": round(gc_fraction(record.seq) * 100, 2)
        }
        row.update(fast_tetranucleotide(record.seq))
        features_list.append(row)
    
    df_genomic = pd.DataFrame(features_list)

    # 2. Pfam
    pfam_file = os.path.join(pfam_path, f"{sample_name}_pfam.domtbl")
    
    if os.path.exists(pfam_file):
        try:
            pfam_df = pd.read_csv(pfam_file, sep=r'\s+', comment='#', names=colunas_hmmer, usecols=range(23))
            
            pfam_pivot = pfam_df.groupby(["Query_ID", "HMM_ID"])["Full_Score"].max().unstack(fill_value=0).reset_index()
            
            df_final = pd.merge(df_genomic, pfam_pivot, left_on="Node", right_on="Query_ID", how="left").fillna(0)
            if "Query_ID" in df_final.columns: df_final.drop(columns=["Query_ID"], inplace=True)
            
            all_samples_data.append(df_final)
        except Exception as e:
            print(f"ERROR in {sample_name}: {e}")
    else:
        print(f"Pfam not found for sample {sample_name}")

if all_samples_data:
    master_matrix = pd.concat(all_samples_data, axis=0, ignore_index=True).fillna(0)
    master_matrix.to_csv(output_final, index=False)
    print(f"\nFinal matrix: {output_final} ({master_matrix.shape})")