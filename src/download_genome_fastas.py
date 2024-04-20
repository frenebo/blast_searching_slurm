# from Bio import Entrez
import argparse
import zipfile
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


        zip_file_path = os.path.join(args.genome_data_dir, genome_accession_id + ".zip")
        accession_data_dir = os.path.join(args.genome_data_dir, genome_accession_id + "_genomedata")
        os.makedirs(accession_data_dir, exist_ok=True)

        with requests.get(request_url, stream=True) as r:
            r.raise_for_status()
            with open(zip_file_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192): 
                    # If you have chunk encoded response uncomment if
                    # and set chunk_size parameter to None.
                    #if chunk: 
                    f.write(chunk)

            
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(accession_data_dir)
        
        # break
        # if idx > 20:
        #     break

