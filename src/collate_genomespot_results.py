import argparse
import os
import pandas as pd

def read_genomespot_resultfile(filepath):
    pass

def flatten_prediction_df_to_dict(pred_df, genome_accession):
    pred_dict = {
        "genome_accession": genome_accession
    }
    for idx, row in pred_df.iterrows():
        target = row["target"]
        value = row["value"]
        error = row["error"]
        # units = row[""]
        is_novel = row["is_novel"]
        warning = row["warning"]

        pred_dict[target + "_value"] = value
        pred_dict[target + "_error"] = error
        pred_dict[target + "_is_novel"] = is_novel
        pred_dict[target + "_warning"] = warning



    """
    Example:

    target  value   error   units   is_novel        warning
    oxygen  tolerant        0.9164311302519951      probability     False   None
    ph_max  9.17161751827928        1.3156701965912105      pH      False   None
    ph_min  5.190154346282054       0.8895952514190193      pH      False   None
    ph_optimum      7.085979582004057       0.9085824364534854      pH      False   None
    salinity_max    6.739979245954464       3.97333676534293        % w/v NaCl      False   None
    salinity_min    0       1.1827437075980203      % w/v NaCl      False   min_exceeded
    salinity_optimum        1.7614621239525108      2.048273557162601       % w/v NaCl      False   None
    temperature_max 39.70402546795431       5.487577411119473       C       False   None
    temperature_min 11.337712558738229      6.5935753513485995      C       False   None
    temperature_optimum     29.90546059538
    """


def read_directory_files_and_save_results(source_dirpath, result_tsv_fp):
    all_files = os.listdir(source_dirpath)
    pred_suffix = ".predictions.tsv"
    tsv_files = [ f for f in all_files if f.endswith(pred_suffix)]

    genome_pred_entries = []

    for fname in tsv_files:
        fpath = os.path.join(source_dirpath, fname)
        accession = fname[:-len(".predictions.tsv")]
        print(fpath, accession)
        genome_preds_df = pd.read_csv(fpath, sep="\t")
        print(genome_preds_df)

        genome_flattened_preds = flatten_prediction_df_to_dict(
            pred_df=genome_preds_df,
            genome_accession=accession,
        )
        genome_pred_entries.append(genome_flattened_preds)
    
    result_df = pd.DataFrame(genome_pred_entries)

    result_df.to_csv(result_tsv_fp, sep="\t")

    # tsv


    # pass

def main():
    parser = argparse.ArgumentParser(
        prog="Takes a directory of genome spot results and collates them into a table",
    )
    parser.add_argument("source_directory_genomespot_results")
    parser.add_argument("output_tsv_file")

    args = parser.parse_args()

    if not args.output_tsv_file.endswith(".tsv"):
        raise Exception("Expected '{}' to end with .tsv".format(args.output_tsv_file))

    read_directory_files_and_save_results(
        source_dirpath=args.source_directory_genomespot_results,
        result_tsv_fp=args.output_tsv_file)

    # if args
    # os


if __name__ == "__main__":
    main()