# from Bio import Entrez
import argparse
import zipfile
import requests
import os
import pandas as pd

def download_genomes(source_summary_tsv, genome_data_dir, record_fhandle):
    os.makedirs(args.genome_data_dir, exist_ok=True)
    df = pd.read_csv(args.source_summary_tsv, sep="\t")

    for idx, row in df.iterrows():
        genome_accession_id = str(row["Assembly"])
        
        print("{idx}/{rowcount} ".format(idx=idx+1,rowcount=len(df.index)), end="",flush=True)
        if str(genome_accession_id) == 'nan':
            print(" (Skipping row, accession id missing)")
            continue

        request_url = (r"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/" +
            r"{}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&".format(genome_accession_id) +
            r"include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&"+
            r"include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED")


        zip_file_path = os.path.join(genome_data_dir, genome_accession_id + ".zip")
        accession_data_dir = os.path.join(genome_data_dir, genome_accession_id + "_genomedata")
        os.makedirs(accession_data_dir, exist_ok=True)
        
        with requests.get(request_url, stream=True) as r:
            r.raise_for_status()
            with open(zip_file_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192): 
                    f.write(chunk)

        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(accession_data_dir)
        
        catalog_json_fp = os.path.join(accession_data_dir, "ncbi_dataset/data/dataset_catalog.json")
        assembly_data_json_fp = os.path.join(accession_data_dir, "ncbi_dataset/data/assembly_data_report.jsonl")

        record_f.write("{}\t{}\t{}\n".format(genome_accession_id, catalog_json_fp, assembly_data_json_fp))
        record_f.flush()
        os.fsync()




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="get_ftp_links_for_genomes",
        description="Get ftp links from tsv summary of the sequences",
    )
    parser.add_argument("source_summary_tsv")
    parser.add_argument("genome_data_dir")
    parser.add_argument("download_record_tsv")

    args = parser.parse_args()
    if not args.source_summary_tsv.endswith(".tsv"):
        raise Exception("Expected tsv file here")
    if not args.download_record_tsv.endswith(".tsv"):
        raise Exception("Expected tsv file here")
    

    # df = 

    with open(args.download_record_tsv, "w") as record_f:
        record_f.write("genome_accession\tcatalog_json_fp\tassembly_data_json_fp\n")
    
        download_genomes(args.source_summary_tsv, args.genome_data_dir, record_f)

