import argparse
import os

def read_genomespot_resultfile(filepath):
    pass

def read_directory_files_and_save_results(source_dirpath, result_tsv_fp):
    all_files = os.listdir(source_dirpath)
    pred_suffix = ".predictions.tsv"
    tsv_files = [ f for f in all_files if f.endswith(pred_suffix)]

    for fname in tsv_files:
        fpath = os.path.join(source_dirpath, fname)
        accession = [:-len(".predictions.tsv")]
        print(fpath, accession)
        # pred_suffix

    # tsv


    # pass

def main():
    parser = argparse.ArgumentParser(
        help="Takes a directory of genome spot results and collates them into a table",
    )
    parser.add_argument("source_directory_genomespot_results")
    parser.add_argument("output_tsv_file")

    args = parser.parse_args()

    read_directory_files_and_save_results(
        source_dirpath=args.source_directory_genomespot_results,
        result_tsv_fp=args.output_tsv_file)

    # if args
    # os


if __name__ == "__main__":
    main()