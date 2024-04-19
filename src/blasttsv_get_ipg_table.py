import argparse
# import pandas as pd
# "/home/pk5192/Documents/blast_searching_slurm/data"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='blasttsv_get_ipg_table',
                    description="Take the tsv output from the blast and use eutils to get source info"
    )
    parser.add_argument("source_blast_tsv")
    parser.add_argument("output_efetch_tsv")
    args = parser.parse_args()
    if args.source_blast_tsv[-4:] != ".tsv":
        raise Exception("Expected tsv file!")
    elif args.output_efetch_tsv[-4:] != ".tsv":
        raise Exception("expected tsv file")
    
    with open(args.source_blast_tsv, "r") as f:
        for line in f:
            line_vals = line.split("\t")
            if len(line_vals) <= 1:
                continue
            
            refs_for_protseq_match = line_vals[3]
            for protref in refs_for_protseq_match.split(";"):
                print(protref.split("|")[1])
            # print(refs_for_protseq_match.split(";"))
            break
