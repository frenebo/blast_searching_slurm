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
    ):
    slurm_setup = "#!/bin/bash\n" +\
        "#SBATCH --job-name={}\n".format(jobname) +\
        "#SBATCH --nodes=1\n" +\
        "#SBATCH --ntasks=4\n" +\
        "#SBATCH --mem=16G\n" +\
        "#SBATCH --time=01:00:00\n" +\
        "#SBATCH --mail-type=begin\n" +\
        "#SBATCH --mail-type=end\n" +\
        "#SBATCH --mail-user=pk5192@princeton.edu\n" +\
        "\n"
    
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
        "-outfmt 10 "+\
        "-max_target_seqs {max_target_seqs} ".format(max_target_seqs=max_target_seqs)+\
        "-num_threads {thread_count} ".format(thread_count=thread_count)+\
        "-out {outputfilename}\n".format(outputfile=outputfile)
        
    end_timestamp = "echo 'end time'\n" +\
        "date\n"
    
    return slurm_setup + module_setup +  goto_workingpath +  start_timestamp + psiblast_command + end_timestamp
    
    # "\n" +
    # "python3 myscript.py""
    
#     export BLASTDB=/scratch/gpfs/pk5192/ncbi_blastdatabase_downloads/nrstuff/

# psiblast -db nr -query /home/pk5192/Documents/seqtest/doing_with_slurm/MtrC_WP_164927685.1.fasta -num_iterations 1  -out_pssm /home/pk5192/Documents/seqtest/doing_with_slurm/mypssm.smp -save_each_pssm   -evalue 0.05 -word_size 3 -outfmt 10 -max_target_seqs 500 -num_threads 4 -out myresultfile.txt




def start_slurm_job(slurmjob_path):
    comm = r"sbatch {}".format(slurmjob_path)
    print(comm)
    process = subprocess.Popen(comm,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    # process = ubpr
    output, error=proces.communicate()
    errcode = proces.returncode
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

if __name__ == "__main__":
    job1_working_dir = "/home/pk5192/Documents/blast_searching_slurm/data/first_1000_test/"
    make_data_dir(job1_working_dir)
    
    job1_text = prepare_blast_slurmjob_text(
        jobname="job_top1000_blast",
        query_fasta_path="/home/pk5192/Documents/seqtest/doing_with_slurm/MtrC_WP_164927685.1.fasta",
        blastdb_path="/scratch/gpfs/pk5192/ncbi_blastdatabase_downloads/nrstuff/",
        working_dirpath=job1_working_dir,
        outputfilename="top_1000.csv",
        thread_count=4,
        max_target_seqs=500,
        evalue=0.05,
        word_size=3,
        num_iterations=2,
    )
    job1_slurmpath = os.path.join(job1_working_dir, "job.slurm"
    
    with open(job1_slurmpath, "w") as f:
        f.write(job1_text)
    
    start_slurm_job(job1_slurmpath)
    