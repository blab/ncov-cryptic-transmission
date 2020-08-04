"""
Take a sequences FASTA and metadata TSV and output filtered metadata TSV
with only strains present in the the sequences FASTA
"""
import argparse
import Bio
import Bio.SeqIO
from Bio.Seq import Seq
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter metadata TSV based on sequences FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequences", required=True, help="sequences FASTA")
    parser.add_argument("--metadata", required=True, help="metadata TSV")
    parser.add_argument("--output", required=True, help="filtered metadata TSV")
    args = parser.parse_args()

    strains = set()
    for record in Bio.SeqIO.parse(args.sequences, 'fasta'):
        strains.add(record.id)

    df = pd.read_csv(args.metadata, sep='\t')
    print(len(df))
    df = df.loc[df['strain'].isin(strains)]
    print(len(df))
    df.to_csv(args.output, sep='\t', index=False)
