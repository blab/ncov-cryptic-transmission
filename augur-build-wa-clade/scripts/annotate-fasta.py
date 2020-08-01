"""
Annotate FASTA with metadata for BEAST
"""
import argparse
import Bio
import Bio.SeqIO
from Bio.Seq import Seq
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Annotate FASTA with metadata for BEAST",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="FASTA file of alignment")
    parser.add_argument("--metadata", required=True, help="metadata tsv")
    parser.add_argument("--output", required=True, help="FASTA file of output alignment")
    args = parser.parse_args()

    df = pd.read_csv(args.metadata, sep='\t')
    strain_to_date = {}
    strain_to_region = {}
    for index, row in df.iterrows():
        strain_to_date[row['strain']] = row['date']
        strain_to_region[row['strain']] = row['region']

    with open(args.output, 'w') as outfile:
        for record in Bio.SeqIO.parse(args.alignment, 'fasta'):
            print(record)
            if strain_to_region[record.id] == "Washington":
                new_id = record.id + "|" + strain_to_date[record.id] + "|" + strain_to_region[record.id]
                record.id = new_id
                record.name = ""
                record.description = ""
                Bio.SeqIO.write(record, outfile, 'fasta')
