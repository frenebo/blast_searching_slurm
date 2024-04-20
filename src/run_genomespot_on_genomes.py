import argparse
import os
import pandas as pd

def run_genomes_and_save_preds_to_file(genome_info_tsv, results_filehandle):
    genome_info_df = pd.read_csv(genome_info_tsv, sep="\t")

    # genome_info_df
    for row, idx in genome_info_df.iterrows():
        print(row)
        break





def main():
    parser = argparse.ArgumentParser(
        prog="Run genome spot from given tsv of genome data stuff",
    )
    parser.add_argument("source_genomepath_tsv", help="Path to the tab separated value file containing genome json info from ncbi")
    parser.add_argument("result_save_tsv", help="Tab separated file containing the output of GenomeSPOT")

    args = parser.parse_args()

    if not args.source_genomepath_tsv.endswith(".tsv"):
        raise Exception("Expected source genome path tsv file to end in .tsv")
    if not args.result_save_tsv.endswith(".tsv"):
        raise Exception("Expected result save tsv filepath to end in .tsv")
    
    with open(args.result_save_tsv, "w") as result_file:
        run_genomes_and_save_preds_to_file(args.source_genomepath_tsv, result_file)
    


if __name__ == "__main__":
    main()