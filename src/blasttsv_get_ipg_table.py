import argparse
from Bio import Entrez
# import pandas as pd
import xml.etree.ElementTree
import subprocess
import xml
import math

def search_proteins_in_entrez(all_prot_accession_ids):
    proteins_not_found_in_database = []
    # Split requests into groups of 1000
    group_search_size = 1000
    for i in range(math.ceil(len(all_prot_accession_ids) / group_search_size)):
        start_idx = i * group_search_size
        end_idx =  (i + 1) * group_search_size
        if end_idx >= len(all_prot_accession_ids):
            end_idx = len(all_prot_accession_ids)
        print("Searching entrez for protein indices {}-{}".format(start_idx,end_idx-1))
        
        prot_search_group = all_prot_accession_ids[start_idx : end_idx]
        http_response = Entrez.esearch("protein", ",".join(prot_search_group))
        text_response = http_response.read().decode('utf-8')
        # print(prot_search_group)
        # print(text_response)

        xml_resp = xml.etree.ElementTree.fromstring(text_response)
        # xml_resp.getroot().
        xml_errlist = xml_resp.find("ErrorList")
        if xml_errlist is None:
            continue
        else:
            for xml_notfoundphrase in xml_errlist.findall("PhraseNotFound"):
                missing_protein_accession = xml_notfoundphrase.text
                proteins_not_found_in_database.append(missing_protein_accession)
                # print(xml_notfoundphrase.text)
    
    proteins_present_in_db = [acc for acc in all_prot_accession_ids if (acc not in proteins_not_found_in_database)]

    return proteins_present_in_db

        # break


def get_protein_info_from_entrez(prot_accessions):
    group_search_size = 500
    output_tsv_string = ""
    for i in range(math.ceil(len(prot_accessions) / group_search_size)):
        start_idx = i * group_search_size
        end_idx =  (i + 1) * group_search_size
        if end_idx >= len(prot_accessions):
            end_idx = len(prot_accessions)
        
        print("Getting info for proteins {}-{}".format(start_idx,end_idx-1))
        search_prots = prot_accessions[start_idx:end_idx]
        html_response = Entrez.efetch(db="protein", id=",".join(search_prots), rettype='ipg', retmode='text')
        text_response = html_response.read().decode("utf-8")
        
        
        lines = text_response.split("\n")
        
        # Remove extra header if this is the first request
        if i != 0:
            lines = lines[1:]
        
        without_errlines = []
        for l in lines:
            if len(l.split("\t")) != 11:
                print("   line err: {}".format(l))
            else:
                without_errlines.append(l)
        lines = without_errlines

            
        output_tsv_string += "\n".join(lines)
    
    return output_tsv_string
        # print(text_response)
        # break
        # text_respo
        # print("Searching entrez for protein indices {}-{}".format(start_idx,end_idx-1))
        

    # Entrez.efetch(db="protein", id=prot_accessions.join(","))
# "/home/pk5192/Documents/blast_searching_slurm/data"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='blasttsv_get_ipg_table',
                    description="Take the tsv output from the blast and use eutils to get source info"
    )
    parser.add_argument("source_blast_tsv")
    # parser.add_argument("intermediate_prot_accession_file")
    parser.add_argument("output_efetch_tsv")
    parser.add_argument("entrez_email")
    args = parser.parse_args()
    if args.source_blast_tsv[-4:] != ".tsv":
        raise Exception("Expected tsv file!")
    # if args.intermediate_prot_accession_file[-4:] != ".txt":
    #     raise Exception("expected txt file!")
    elif args.output_efetch_tsv[-4:] != ".tsv":
        raise Exception("expected tsv file")
    
    Entrez.email = args.entrez_email
    
    prot_accessions_to_search = []
    
    with open(args.source_blast_tsv, "r") as f:
        for line in f:
            line_vals = line.split("\t")
            if len(line_vals) <= 1:
                continue
            
            refs_for_protseq_match = line_vals[3]
            cnt = 0
            for protref in refs_for_protseq_match.split(";"):
                prot_accession = protref.split("|")[1]
                # print(prot_accession)
                prot_accessions_to_search.append(prot_accession)
                # break
                # cnt += 1
                # if cnt > 2:
                #     break
            # prot_accessions_to_search
            # print(refs_for_protseq_match.split(";"))
            # break
    # prot_accessions_to_search = prot_accessions_to_search[0:120]
    print("Searching entrez for {} protein accessions".format(len(prot_accessions_to_search)))
    prot_accession_presentindb = search_proteins_in_entrez(prot_accessions_to_search)
    print("Found {} accessions through entrez, continuing with those".format(len(prot_accession_presentindb)))

    restext = get_protein_info_from_entrez(prot_accession_presentindb)
    with open(args.output_efetch_tsv, "w") as f:
        f.write(restext)