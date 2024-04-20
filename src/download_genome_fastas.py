# from Bio import Entrez
import argparse
import requests
import os
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="get_ftp_links_for_genomes",
        description="Get ftp links from tsv summary of the sequences",
    )
    parser.add_argument("source_summary_tsv")
    parser.add_argument("genome_data_dir")

    args = parser.parse_args()
    if not args.source_summary_tsv.endswith(".tsv"):
        raise Exception("Expected tsv file here")



    os.makedirs(args.genome_data_dir, exist_ok=True)

    df = pd.read_csv(args.source_summary_tsv, sep="\t")

    for idx, row in df.iterrows():
        genome_accession_id = row["Assembly"]

        request_url = (r"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/" +
            r"{}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&".format(genome_accession_id) +
            r"include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&"+
            r"include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED")


        r = requests.get(request_url)

        break

