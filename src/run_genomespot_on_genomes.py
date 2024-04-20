import argparse
import os
import math
import json
import pandas as pd
import subprocess

def build_slurm_job(genomes_and_protein_info_for_job, genomespot_models_path, jobname, timestring, ntasks, memorystring):
    slurm_setup = "#!/bin/bash\n" +\
        "#SBATCH --job-name={jobname}\n".format(jobname=jobname) +\
        "#SBATCH --nodes=1\n" +\
        "#SBATCH --ntasks=1\n" +\
        "#SBATCH --mem={memorystring}\n".format(memorystring) +\
        "#SBATCH --time={timestring}\n".format(timestring=timestring) +\
        "#SBATCH --mail-type=begin\n" +\
        "#SBATCH --mail-type=end\n" +\
        "\n"
    
    tellbash_echo_commands = "set -x\n"
    
    module_setup = "module purge\n" +\
        "module load anaconda3/2024.2\n" +\
        "conda activate genomespotstuff\n"
    
    start_timestamp = "echo 'start time'\n" +\
        "date\n"
    
    do_genomespot = ""
    
    for idx, infoline in enumerate(genomes_and_protein_info_for_job):
        genome_fasta_fp = infoline["genomic_nucleotide_fasta_fp"]
        protein_fasta_fp = infoline["protein_fasta_fp"]
        output_fp_without_predictions_tsv = infoline["save_res_fp_prefix"]

        do_genomespot +=(
            "python -m genome_spot.genome_spot --models {genomespot_models_path} ".format(genomespot_models_path=genomespot_models_path) +
            " --contigs {genome_fasta_fp} ".format(genome_fasta_fp=genome_fasta_fp) + 
            " --proteins {protein_fasta_path} ".format(protein_fasta_path=protein_fasta_fp) + 
            " --output {outputname} ".format(outputname=output_fp_without_predictions_tsv)
        )
        do_genomespot += " & \n"

        if idx % ntasks == 0:
            do_genomespot += "wait\n"

    do_genomespot += "wait\n"
    
    # do_genomespot += "\n"

    end_timestamp = "echo 'end time'\n" +\
        "date\n"
    
    return slurm_setup + tellbash_echo_commands + module_setup + start_timestamp + do_genomespot + end_timestamp
    

def start_slurm_job(slurmfile_path):
    comm = "sbatch {slurmfile_path}".format(slurmfile_path = slurmfile_path)
    print(comm)
    process = subprocess.Popen(comm,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    # process = ubpr
    output, error=process.communicate()
    errcode = process.returncode
    if errcode != 0:
        print(output.decode("utf-8"))
        print(error.decode("utf-8"))
        raise Exception("Got error code {}".format(errcode))
    output = output.decode("utf-8")
    expected_prefix = "Submitted batch job "
    if output[:len(expected_prefix)] != expected_prefix:
        print("Unexpected output!")
        print(output)
    else:
        print("slurm job path: {}".format(slurmfile_path))
        print("job number:")
        print("    ",end="")
        print(output[len(expected_prefix):])

def run_genomespot_slurmjobs(genome_and_proteins_and_save_infos, genomespot_models_path, slurmjobs_dir):
    # Actual time is about 5 seconds per genome - give it maybe 20 per genome to be on the safe side?
    # 10 minute jobs with 30 genomes each

    genomes_per_job = 100
    nthreads = 10
    # each batch will take about 6 seconds, so give some leeway
    # (genomes_per_job/nthreads) * 6 seconds * safety factor
    timestr_per_job = "00:3:00"
    jobmemory="10G" # Need about 500 mb per thread

    for i in range(math.ceil(len(genome_and_proteins_and_save_infos) / genomes_per_job)):
        start_idx = i * genomes_per_job
        end_idx = (i + 1) * genomes_per_job

        jobname = "gsp{}".format(i)

        job_genome_infos = genome_and_proteins_and_save_infos[start_idx:end_idx]

        jobfile_string = build_slurm_job(
            job_genome_infos,
            genomespot_models_path,
            jobname,
            timestr_per_job,
            nthreads,
            jobmemory,
            )
        jobslurm_path = os.path.join(slurmjobs_dir, "{}.slurm".format(jobname))

        with open(jobslurm_path, "w") as f:
            f.write(jobfile_string)

        
        # start_slurm_job(jobslurm_path)



def run_genomes_and_save_preds(
    genome_info_tsv,
    genomedata_absolute_dirpath,
    results_save_dir,
    genomespot_models_path,
    slurmjobs_dir):
    genome_info_df = pd.read_csv(genome_info_tsv, sep="\t")

    os.makedirs(results_save_dir, exist_ok=True)
    os.makedirs(slurmjobs_dir, exist_ok=True)

    collected_info = []
    # genome_info_df
    for idx,row in genome_info_df.iterrows():
        print("{}/{} ".format(idx+1,len(genome_info_df.index)),end="",flush=True)

        genome_accession_id = str(row["genome_accession"])

        with open(row["catalog_json_fp"], "r") as cat_f:
            catalog_info = json.load(cat_f)
        
        with open(row["assembly_data_json_fp"]) as asm_f:
            assembly_data = json.load(asm_f)
        
        assert assembly_data["accession"] == genome_accession_id

        genomic_nucleotide_fasta_fp = None
        protein_fasta_fp = None
        for assem in catalog_info["assemblies"]:
            if "accession" not in assem:
                continue
            
            assert assem["accession"] == genome_accession_id
            for fileentry in assem["files"]:
                if fileentry["fileType"] == "GENOMIC_NUCLEOTIDE_FASTA":
                    genomic_nucleotide_fasta_fp = fileentry["filePath"]
                
                if fileentry["fileType"] == "PROTEIN_FASTA":
                    protein_fasta_fp = fileentry["filePath"]
            
            break
        
        if genomic_nucleotide_fasta_fp is None or protein_fasta_fp is None:
            raise Exception("Could not find protein and genome nucleotide fasta data in {}".format(row["catalog_json_fp"]))
        
        # Get absolute paths where the filepaths will be inside
        # For example /scratch/gpfs/pk5192/April20_MtrCgenomes/GCF_030008505.1_genomedata/ncbi_dataset/data/
        absolute_data_parent_dirpath = os.path.join(
            genomedata_absolute_dirpath,
            genome_accession_id + "_genomedata",
            "ncbi_dataset",
            "data",
            )
        
        genomic_nucleotide_fasta_fp = os.path.join(absolute_data_parent_dirpath, genomic_nucleotide_fasta_fp)
        protein_fasta_fp = os.path.join(absolute_data_parent_dirpath, protein_fasta_fp)
        
        save_res_fp_prefix = os.path.join(results_save_dir, genome_accession_id)

        collected_info.append({
            "genomic_nucleotide_fasta_fp": genomic_nucleotide_fasta_fp,
            "protein_fasta_fp": protein_fasta_fp,
            "save_res_fp_prefix": save_res_fp_prefix,
        })

        if idx > 100:
            break
    
    run_genomespot_slurmjobs(collected_info, genomespot_models_path, slurmjobs_dir)



def main():
    parser = argparse.ArgumentParser(
        prog="Run genome spot from given tsv of genome data stuff",
    )
    parser.add_argument("source_genomepath_tsv", help="Path to the tab separated value file containing genome json info from ncbi")
    parser.add_argument("pathto_genomedata_rootdir", help="Path to the directory with all of the genomes and their downloaded ncbi sequences")
    parser.add_argument("result_save_dir", help="Directory to put results for the genome predictions")
    parser.add_argument("genomespot_models_path", help="Path to GenomeSPOT's models directory")
    parser.add_argument("slurmjobs_dir", help="Directory to store the slurm job files while running")

    args = parser.parse_args()

    if not args.source_genomepath_tsv.endswith(".tsv"):
        raise Exception("Expected source genome path tsv file to end in .tsv")
    # if not args.result_save_tsv.endswith(".tsv"):
    #     raise Exception("Expect?ed result save tsv filepath to end in .tsv")
    
    run_genomes_and_save_preds(
        args.source_genomepath_tsv,
        args.pathto_genomedata_rootdir,
        args.result_save_dir,
        args.genomespot_models_path,
        args.slurmjobs_dir)
    # with open(args.result_save_tsv, "w") as result_file:
    #     run_genomes_and_save_preds_to_file(args.source_genomepath_tsv, result_file)
    


if __name__ == "__main__":
    main()