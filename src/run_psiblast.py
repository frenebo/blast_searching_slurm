import os
import subprocess


def make_data_dir(dirpath):
    os.makedirs(dirpath,exist_ok=True)

def prepare_blast_slurmjob_text(
    jobname,
    query_fasta_path,
    blastdb_path,
    working_dirpath,
    outputfilename,
    thread_count,
    max_target_seqs,
    evalue,
    word_size,
    num_iterations,
    timestring,
    ):
    slurm_setup = "#!/bin/bash\n" +\
        "#SBATCH --job-name={}\n".format(jobname) +\
        "#SBATCH --nodes=1\n" +\
        "#SBATCH --ntasks={}\n".format(thread_count) +\
        "#SBATCH --mem=10G\n" +\
        "#SBATCH --time={}\n".format(timestring) +\
        "#SBATCH --mail-type=begin\n" +\
        "#SBATCH --mail-type=end\n" +\
        "\n"
    
    tellbash_echo_commands = "set -x\n"
    
    module_setup = "module purge\n" +\
        "module load anaconda3/2024.2\n" +\
        "conda activate mtrtest\n"
    
    goto_workingpath = "cd {dirpath}\n".format(dirpath=working_dirpath)
    
    start_timestamp = "echo 'start time'\n" +\
        "date\n"
    
    psiblast_command = "export BLASTDB=\"{blastdb_path}\"\n".format(blastdb_path=blastdb_path) +\
        "psiblast -db nr -query {query_fasta_path} ".format(query_fasta_path=query_fasta_path)+\
        "-num_iterations {num_iterations} ".format(num_iterations=num_iterations)+\
        "-out_pssm pssm_stuff.smp "+\
        "-save_each_pssm "+\
        "-evalue {evalue} ".format(evalue=evalue)+\
        "-word_size {word_size} ".format(word_size=word_size)+\
        "-max_target_seqs {max_target_seqs} ".format(max_target_seqs=max_target_seqs)+\
        "-num_threads {thread_count} ".format(thread_count=thread_count)+\
        "-out {outputfilename} ".format(outputfilename=outputfilename)+\
        "-outfmt \"6 qseqid sseqid sgi sallseqid sallgi sacc sallacc staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand mismatch gapopen qstart qend sstart send evalue pident length bitscore\" \n"
        
    end_timestamp = "echo 'end time'\n" +\
        "date\n"
    
    return slurm_setup + tellbash_echo_commands + module_setup +  goto_workingpath +  start_timestamp + psiblast_command + end_timestamp
    
    # "\n" +
    # "python3 myscript.py""
    
#     export BLASTDB=/scratch/gpfs/pk5192/ncbi_blastdatabase_downloads/nrstuff/

# psiblast -db nr -query /home/pk5192/Documents/seqtest/doing_with_slurm/MtrC_WP_164927685.1.fasta -num_iterations 1  -out_pssm /home/pk5192/Documents/seqtest/doing_with_slurm/mypssm.smp -save_each_pssm   -evalue 0.05 -word_size 3 -outfmt 10 -max_target_seqs 500 -num_threads 4 -out myresultfile.txt




def start_slurm_job(slurmjob_path):
    comm = r"sbatch {}".format(slurmjob_path)
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
        print("slurm job path: {}".format(slurmjob_path))
        print("job number:")
        print("    ",end="")
        print(output[len(expected_prefix):])
    # print(out)

def run_slrm_blast(query_fasta_path, eval_thresh, num_iterations, stringname, timestring):
    job1_working_dir = "/home/pk5192/Documents/blast_searching_slurm/data/MR1queries/{}_test/".format(stringname)
    
    make_data_dir(job1_working_dir)
    
    job1_text = prepare_blast_slurmjob_text(
        jobname="{}".format(stringname),
        query_fasta_path=query_fasta_path,
        blastdb_path="/scratch/gpfs/pk5192/ncbi_blastdatabase_downloads/nrstuff/",
        working_dirpath=job1_working_dir,
        outputfilename="psiblastoutput_{}.tsv".format(stringname),
        thread_count=4,
        max_target_seqs=100000,
        evalue=eval_thresh,
        word_size=3,
        num_iterations=num_iterations,
        timestring=timestring,
    )
    job1_slurmpath = os.path.join(job1_working_dir, "job.slurm")
    
    with open(job1_slurmpath, "w") as f:
        f.write(job1_text)
    
    start_slurm_job(job1_slurmpath)
    

if __name__ == "__main__":
    # try_amt(100,"1hd")
    # try_amt(500,"5hd")
    # try_amt(1000,"1k")
    # try_amt(1000,"5k")
    # run_slrm_blast(eval_thresh=0.05, num_iterations=1, stringname="pk05", timestring="04:00:00")
    # run_slrm_blast(eval_thresh=0.01, num_iterations=1, stringname="pk01", timestring="04:00:00")
    # run_slrm_blast(eval_thresh=0.005, num_iterations=1, stringname="pk005", timestring="04:00:00")
    # run_slrm_blast(eval_thresh=0.001, num_iterations=1, stringname="pk001", timestring="04:00:00")
    run_slrm_blast(query_fasta_path="/home/pk5192/Documents/blast_searching_slurm/query_seqs/mr1_mtrproteins/MtrC_WP_011071901.1.fasta",
        eval_thresh=0.005,
        num_iterations=1,
        stringname="MtrC_search",
        timestring="02:30:00")
    run_slrm_blast(query_fasta_path="/home/pk5192/Documents/blast_searching_slurm/query_seqs/mr1_mtrproteins/MtrA_WP_011071900.1.fasta",
        eval_thresh=0.005,
        num_iterations=1,
        stringname="MtrA_search",
        timestring="02:30:00")
    run_slrm_blast(query_fasta_path="/home/pk5192/Documents/blast_searching_slurm/query_seqs/mr1_mtrproteins/MtrB_WP_011071899.1.fasta",
        eval_thresh=0.005,
        num_iterations=1,
        stringname="MtrB_search",
        timestring="02:30:00")
    # try_amt(10, "10its", "40:00:00")