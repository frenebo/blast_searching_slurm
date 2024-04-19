import argparse
# import pandas as pd
import subprocess
# "/home/pk5192/Documents/blast_searching_slurm/data"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='blasttsv_get_ipg_table',
                    description="Take the tsv output from the blast and use eutils to get source info"
    )
    parser.add_argument("source_blast_tsv")
    parser.add_argument("intermediate_prot_accession_file")
    parser.add_argument("output_efetch_tsv")
    args = parser.parse_args()
    if args.source_blast_tsv[-4:] != ".tsv":
        raise Exception("Expected tsv file!")
    if args.intermediate_prot_accession_file[-4:] != ".txt":
        raise Exception("expected txt file!")
    elif args.output_efetch_tsv[-4:] != ".tsv":
        raise Exception("expected tsv file")
    
    prot_accessions_to_search = []
    
    with open(args.source_blast_tsv, "r") as f:
        for line in f:
            line_vals = line.split("\t")
            if len(line_vals) <= 1:
                continue
            
            refs_for_protseq_match = line_vals[3]
            for protref in refs_for_protseq_match.split(";"):
                prot_accession = protref.split("|")[1]
                # print(prot_accession)
                prot_accessions_to_search.append(prot_accession)
                break
            # print(refs_for_protseq_match.split(";"))
            break
    with open(args.intermediate_prot_accession_file, "w") as f:
        f.write("\t".join(prot_accessions_to_search))
    
    comm = "cat '{accession_text_list}' | epost -db protein | efetch -format ipg  > {output_efetch_tsv}".format(
        accession_text_list=args.intermediate_prot_accession_file,
        output_efetch_tsv=args.output_efetch_tsv,
        )

    process = subprocess.Popen(comm,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    output, error=process.communicate()
    errcode = process.returncode
    
    output = output.decode("utf-8")
    error = error.decode("utf-8")

    if errcode != 0:
        print(output)
        print(error)
        raise Exception("Got error code {}".format(errcode))
    
    print(output)
    print(error)