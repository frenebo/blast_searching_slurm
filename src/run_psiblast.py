import os
import subprocess


def make_data_dir(dirpath):
    os.makedirs(dirpath,exist_ok=True)

def prepare_blast_slurmjob_text(
    jobname,
    query_fasta_path,
    blastdb_path,
    working_dirpath,
    output_tsv,
    nthreads,
    max_target_seqs,
    evalue_threshold,
    word_size,
    num_iterations,
    timestring,
    database_name,
    ):
    slurm_setup = "#!/bin/bash\n" +\
        "#SBATCH --job-name={}\n".format(jobname) +\
        "#SBATCH --nodes=1\n" +\
        "#SBATCH --ntasks={}\n".format(nthreads) +\
        "#SBATCH --mem=10G\n" +\
        "#SBATCH --time={}\n".format(timestring) +\
        "#SBATCH --mail-type=begin\n" +\
        "#SBATCH --mail-type=end\n" +\
        "\n"
    
    tellbash_echo_commands = "set -x\n"
    
    module_setup = "module purge\n" +\
        "module load anaconda3/2024.2\n" +\
        "conda activate mtrtest\n"
    
    # goto_workingpath = "cd {dirpath}\n".format(dirpath=working_dirpath)
    
    start_timestamp = "echo 'start time'\n" +\
        "date\n"
    
    psiblast_command = "export BLASTDB=\"{blastdb_path}\"\n".format(blastdb_path=blastdb_path) +\
        "psiblast -db {database_name} -query {query_fasta_path} ".format(database_name=database_name, query_fasta_path=query_fasta_path)+\
        "-num_iterations {num_iterations} ".format(num_iterations=num_iterations)+\
        "-evalue {evalue} ".format(evalue=evalue_threshold)+\
        "-word_size {word_size} ".format(word_size=word_size)+\
        "-max_target_seqs {max_target_seqs} ".format(max_target_seqs=max_target_seqs)+\
        "-num_threads {thread_count} ".format(thread_count=nthreads)+\
        "-out {outputfilename} ".format(outputfilename=output_tsv)+\
        "-outfmt \"6 qseqid sseqid sgi sallseqid sallgi sacc sallacc staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand mismatch gapopen qstart qend sstart send evalue pident length bitscore\" \n"
        # "-out_pssm pssm_stuff.smp "+\
        # "-save_each_pssm "+\
        
    end_timestamp = "echo 'end time'\n" +\
        "date\n"
    
    return slurm_setup + tellbash_echo_commands + module_setup + start_timestamp + psiblast_command + end_timestamp
    
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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run local PSI blast on a specified local downloaded database, with a given query file and BLAST settings"
    )

        
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    have_defaults = parser.add_argument_group('arguments with default values')

    required.add_argument('--queryfasta', required=True, help="Path to the fasta file with the query sequence")
    required.add_argument('--evalue_threshold', required=True, type=float, help="The maximum evalue for BLAST to consider a protein a match")
    required.add_argument('--blastdb_path', required=True, help="Dirpath for the downloaded NCBI database files")
    required.add_argument('--output_tsv', required=True, help="Output .tsv file to write results to")
    required.add_argument('--db_name', required=True, help="Name of the database to search for matches. Should probably be nr.")
    required.add_argument('--slurmfile_path', required=True, help="Filepath for where to put the slurm file")
    
    have_defaults.add_argument('--jobname', required=True, default="psiblast", help="Name for the slurm job")
    have_defaults.add_argument('--nthreads', required=True, default=4, type=int, help="How many threads to use in blast search. Any increase above 6 likely will not help.")
    have_defaults.add_argument('--max_target_seqs', required=True, default=100000, type=int, help="Max number of matches before BLAST cuts off the output file. Make this bigger than it needs to be!")
    have_defaults.add_argument('--word_size', required=True, default=3, help="Blast algorithm word size")
    have_defaults.add_argument('--timestring', required=True, default="04:00:00", help="HH:MM:SS formatted time to allow slurm job to run.")


    args = parser.parse_args()

    if not os.path.isfile(args.queryfasta):
        raise Exception("Invalid query fasta, is not a file: '{}'".format(args.queryfasta))
    
    if not os.path.isdir(args.blastdb_path):
        raise Exception("Invalid blastdb_path, is not a directory: '{}'".format(args.blastdb_path))
    
    if args.evalue_threshold < 0:
        raise Exception("Invalid evalue threshold, must be positive: {}".format(args.evalue_threshold))
    if args.nthreads < 1 or args.nthreads > 20:
        raise Exception("Invalid nthreads, must be between 1 and 20: {}".format(args.nthreads))
    if args.word_size < 2:
        raise Exception("Invalid word size")

    slurmjob_text = prepare_blast_slurmjob_text(
        query_fasta_path=args.queryfasta,
        evalue_threshold=args.evalue_threshold,
        num_iterations=1,
        jobname=args.jobname,
        timestring=args.timestring,
        blastdb_path=args.blastdb_path,
        database_name=args.db_name,
        word_size=args.word_size,
        max_target_seqs=args.max_target_seqs,
        nthreads=args.nthreads,
        output_tsv=args.output_tsv,
    )

    with open(args.slurmfile_path, "w") as f:
        f.write(job_text)
    
    start_slurm_job(slurmfile_path)