import argparse
import os
import pandas as pd

def get_genomesubdir_names(genomes_dirpath, expected_suffix):
    all_dircontents = os.listdir(genomes_dirpath)
    genomesubdir_names = []
    for name in all_dircontents:
        itempath = os.path.join(genomes_dirpath, name)
        if not os.path.isdir(itempath):
            continue
        if not itempath.endswith(expected_suffix):
            continue
        
        genomesubdir_names.append(name)
    
    return genomesubdir_names

def make_flattened_obj(source_obj, keyprefix=""):
    flattened = {

    }
    for key, val in source_obj.items():
        if isinstance(val, dict):
            subobj_flattened = make_flattened_obj(val, keyprefix=keyprefix + key + "_")

            flattened  = flattened | subobj_flattened
        elif isinstance(val, list):
            pass # Unimplemented
        else:
            flattened[keyprefix + key] = val
        

    return flattened

# def populate_flattened_data_report(flattened_rep, json_data):
#     for key in flattened_rep:

#         if "_" in key:
#             pass
#         else:
#             if key in json_data:

# def

        

def get_info_from_data_report(ass_data_rep_path):
    expected_structure = {
        # "accession": None,
        "annotationInfo": {
            "stats": {
                "geneCounts": {
                    "nonCoding": None,
                    "proteinCoding": None,
                    "pseudogene": None,
                    "total": None,
                }
            }
        },
        "assemblyInfo": {
            "assemblyLevel": None,
            "assemblyLevel": None,
            "assemblyName": None,
            "assemblyStatus": None,
            "assemblyType": None,
            "bioprojectAccession": None,
        },
        "assemblyStats": {
            "contigL50": None,
            "contigN50": None,
            "gcCount": None,
            "gcPercent": None,
            "numberOfComponentSequences": None,
            "numberOfContigs": None,
            "numberOfScaffolds": None,
            "scaffoldL50": None,
            "scaffoldN50": None,
            "totalNumberOfChromosomes": None,
            "totalSequenceLength": None,
            "totalUngappedLength": None,
        },
        "averageNucleotideIdentity": {
            "bestAniMatch": {
                "ani": None,
                "assembly": None,
                "assemblyCoverage": None,
                "category": None,
                "organismName": None,
                "typeAssemblyCoverage": None,
            },
            "category": None,
            "comment": None,
            "matchStatus": None,
            "submittedOrganism": None,
            "submittedSpecies": None,
            "taxonomyCheckStatus": None,
        },
        "checkmInfo": {
            "checkmMarkerSet": None,
            "checkmMarkerSetRank": None,
            "checkmSpeciesTaxId": None,
            "checkmVersion": None,
            "completeness": None,
            "completenessPercentile": None,
            "contamination": None,
        },
        "organism": {
            "organismName": None,
            "taxId": None,
        },
        "sourceDatabase": None,
    }


    genome_data = make_flattened_obj(expected_structure)
    genome_data["datafilemissing"] = None
    genome_data["datafileunreadable"] = None

    if not os.path.exists(ass_data_rep_path):
        genome_data["datafilemissing"] = True
        return genome_data
    genome_data["datafilemissing"] = False
    
    
    # else:
    json_contents = None
    try:
        with open(ass_data_rep_path, "r") as f:
            json_contents = json.load(f)
    except:
        genome_data["datafileunreadable"] = True
        return genome_data

    genome_data["datafileunreadable"] = False
    
    flattened_data_report = make_flattened_obj(json_contents)

    # Merge the actual data into our expected data structure (leaves out unwanted features)
    for k in genome_data:
        if k in flattened_data_report:
            genome_data[k] = flattened_data_report[k]
    
    return genome_data




def get_info_about_genomes(genomes_dirpath):
    genomesubdir_expected_suffix = "_genomedata"
    genome_subdir_names = get_genomesubdir_names(genomes_dirpath, expected_suffix=genomesubdir_expected_suffix)

    # genome_info = {
    #     "could_not_access_data": [],
    #     ""
    # }
    all_genomes_info = []

    for gsubdir_name in genome_subdir_names:
        genome_dirpath = os.path.join(genomes_dirpath, gsubdir_name)
        genome_accession = gsubdir_name[:-len(genomesubdir_expected_suffix)]
        
        data_dirpath = os.path.join(genome_dirpath, "ncbi_dataset/data")
        ass_data_rep_path = os.path.join(data_dirpath, "assembly_data_report.jsonl")

        
        flattened_data = get_info_from_data_report(ass_data_rep_path)
        flattened_data["Genome Accession"] = genome_accession

        all_genomes_info.append(flattened_data)
    
    return pd.DataFrame(all_genomes_info)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Get genome sequence quality data from genome files"
    )

    parser.add_argument("source_genomes_directory", help="Directory containing downloaded genomes")
    parser.add_argument("results_stats_tsv", help="TSV file to store resulting statistics inside")

    args = parser.parse_args()
    
    if not os.path.exists(args.source_genomes_directory):
        raise Exception("Directory {} does not exist".format(args.source_genomes_directory))
    if not os.path.isdir(args.source_genomes_directory):
        raise Exception("{} is not a directory".format(args.source_genomes_directory))
    if not args.results_stats_tsv.endswith(".tsv"):
        raise Exception("Expected results_stats_tsv to end in .tsv")
    
    
    genome_stats = get_info_about_genomes(args.source_genomes_directory)

    genome_stats.to_csv(args.results_stats_tsv, sep="\t")
