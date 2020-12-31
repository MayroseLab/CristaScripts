__author__ = 'Shiran'


def create_job_file(job_name, command, file_name, error_files_path, job_files_path, python_version="3.3"):
	with open(job_files_path + "/" + file_name, "w") as handle:
		handle.write("#!/bin/bash\n\n")  # Vladi said it should be tcsh!
		handle.write("#PBS -N " + job_name + "\n")
		handle.write("#PBS -r y\n")
		handle.write("#PBS -q lifesciweb\n")
		handle.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n")
		handle.write("#PBS -e " + error_files_path + "\n")
		handle.write("#PBS -o " + error_files_path + "\n")
		handle.write("module load python/python-anaconda3.6.5\n")
		handle.write("module load bwa/bwa-0.7.17\n")
		handle.write("module load samtools/samtools-1.6\n")
		handle.write(command + "\n")
	return job_files_path + "/" + file_name

