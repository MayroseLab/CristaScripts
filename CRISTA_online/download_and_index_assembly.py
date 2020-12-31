import argparse
import re
import sys
sys.path.append("/groups/itay_mayrose/shiranabad/CRISPR/CRISTA_online/")
import pandas as pd
import shutil
import time
from datetime import datetime, timedelta
import logging
from utils import init_commandline_logger, change_path_permissions_to_777
import urllib.request, urllib.error, os
from definitions import *

UCSC_DB_PAGE = "http://hgdownload.cse.ucsc.edu/goldenPath/{0}/bigZips/{0}.fa.gz"
SPECIES_FTP_PAGE = "ftp://ftp.ensemblgenomes.org/pub/plants/release-44/fasta/{0}/dna/"


def write_current_timestamp_for_assembly(assembly_dirpath):
	with open(assembly_dirpath + SEP + "timestamp.txt", "w") as fpw:
		fpw.write(time.strftime("%x")) # to be read with time.strptime(str, "%x")


def read_current_timestamp_for_assembly(assembly_dirpath):
	"""
	:param assembly_dirpath: where timestamp.txt is
	:return: time_struct object
	"""
	with open(assembly_dirpath + "timestamp.txt") as fpr:
		time_str = fpr.read().strip()
	date = time.strptime(time_str, "%x")
	return date


def is_time_older_in_months(assembly_dirpath, months=6):
	#read date and transform time_struct to datetime object
	modification_date = read_current_timestamp_for_assembly(assembly_dirpath)
	modification_date = datetime.fromtimestamp(time.mktime(modification_date))
	x_months_ago = datetime.now() - timedelta(days=months*30)
	return modification_date < x_months_ago


def download_error(assembly_dirpath, debugging_msg="Can't download assembly."):
	set_error_message("Can't download assembly.")
	if os.path.exists(assembly_dirpath + "currently_processing"):
		os.remove(assembly_dirpath + "currently_processing")
	raise RuntimeError(debugging_msg)


def get_download_url_from_ensemble_plants(assembly_name, assembly_dirpath):
	ftp_download_path = assembly_dirpath + "/" + assembly_name + ".txt"
	ftp_url = SPECIES_FTP_PAGE.format(assembly_name)

	# first get ftp page
	for i in range(5):
		try:
			urllib.request.urlretrieve(ftp_url, ftp_download_path)
		except Exception as err:
			pass

	if not os.path.exists(ftp_download_path):
		download_error(assembly_dirpath)

	# parse ftp page to get correct url
	with open(ftp_download_path) as fpr:
		for line in fpr:
			lookup = re.search("\.dna\.toplevel\.fa\.gz", line)
			if lookup:
				download_dest_filename = re.search(assembly_name + ".*gz", line, re.I)
				if not download_dest_filename:
					download_error(assembly_dirpath, "Can't find dna.toplevel.gz file for plant.")
				else:
					return ftp_url + download_dest_filename.group()


def get_download_url_from_ucsc(assembly_name):
	return UCSC_DB_PAGE.format(assembly_name)


def plant_assembly_name(assembly_name):
	"""
	checks if assembly is PLANTS. if so, returns the assembly nanem, else returns False
	:param assembly_name: the value name (the name of the assembly version,
	different from the organism name that appears on the url
	:return: organism name, if not plant returns None
	"""
	df = pd.read_csv(PLANTS_MAPPING_FILE, sep=",")
	print("assembly_name:", assembly_name)
	if assembly_name in df["value"].values:
		print(df.loc[df["value"]==assembly_name, "assembly_name"])
		return df.loc[df["value"]==assembly_name, "assembly_name"].values[0]


def download_assembly(assembly_name, assembly_dirpath):
	"""
	downloads assembly genome assembly, extracts fasta from it (converts 2bit or extracts gz)
	:return: returns the extracted fasta file
	"""
	tries = 1

	plant_name = plant_assembly_name(assembly_name)

	if plant_name is None:
		url_for_assembly_file = get_download_url_from_ucsc(assembly_name)
	else:
		url_for_assembly_file = get_download_url_from_ensemble_plants(plant_name, assembly_dirpath)

	print(url_for_assembly_file)
	file_suffix = re.search("\.(fa.*)|(2bit)$", url_for_assembly_file.split(SEP)[-1]).group()
	fa_file = assembly_dirpath + assembly_name + ".fa"
	downloaded_file = assembly_dirpath + assembly_name
	max_tries = 10
	while not os.path.exists(downloaded_file + file_suffix) and tries <= max_tries:
		try:
			urllib.request.urlretrieve(url_for_assembly_file, downloaded_file + file_suffix)
		except Exception as err:
			if url_for_assembly_file.endswith("fa.gz"):
				file_suffix = ".2bit"
				try:
					urllib.request.urlretrieve(url_for_assembly_file[:-6] + file_suffix, downloaded_file + file_suffix)
				except:
					file_suffix = ".fa.gz"
					time.sleep(60)
					tries += 1

	if not os.path.exists(downloaded_file + file_suffix):  # loop has ended with no results
		download_error(assembly_dirpath)

	try:
		if file_suffix == ".2bit":
			# transform twobit to fa.gz
			os.system(TWOBITTOFA_EXE + " " + downloaded_file + file_suffix + " " + fa_file)
			os.remove(downloaded_file + ".2bit")
		elif file_suffix == ".fa.gz":
			# unzip gzip (to rezip with bgzip)
			os.system("gunzip " + downloaded_file + file_suffix)

	except:
		set_error_message("Failed to extract assembly fasta file.")
		raise RuntimeError("Failed to extract assembly fasta file.")

	return fa_file #returns the path of the downloaded fa file


def process_assembly(assembly_name, logger=None):
	assembly_dirpath = GENOMES_PATH + assembly_name + SEP
	fa_gz_file = assembly_dirpath + "/" + assembly_name + ".fa.gz"

	if not os.path.exists(assembly_dirpath):
		os.mkdir(assembly_dirpath)

	# if another process is currently processing the file
	tries = 1
	max_tries = 300
	while os.path.exists(assembly_dirpath + "currently_processing") and tries < max_tries:
		if logger:
			logger.info("Waiting for another process to download file...")
		time.sleep(120)
		tries += 1

	if tries == max_tries:
		download_error(assembly_dirpath)
	elif os.path.exists(fa_gz_file + ".bwt"):
		return

	open(assembly_dirpath + "currently_processing", "w").close()

	fasta_path = download_assembly(assembly_name, assembly_dirpath)

	try:
		# zip to a fa.gz file
		os.system("/groups/itay_mayrose/shiranabad/applications/htslib-1.3.2/bgzip " + fasta_path)

		# create fa.gz.fai file
		if not os.path.exists(fa_gz_file + ".fai"):
			os.system(SAMTOOLS_EXE + " faidx " + fa_gz_file)

		# create fa.gz.bwt file
		if not os.path.exists(fa_gz_file + ".bwt"):
			os.system(BWA_EXE + " index " + fa_gz_file)

		write_current_timestamp_for_assembly(assembly_dirpath)

	except Exception as err:
		with open(assembly_dirpath + "processing error message.txt", "w") as fpw:
			fpw.write(err)

		set_error_message("Assembly file processing failed.")
		raise RuntimeError("Assembly file processing failed.")
	finally:
		if os.path.exists(assembly_dirpath + "currently_processing"):
			os.remove(assembly_dirpath + "currently_processing")
		change_path_permissions_to_777(assembly_dirpath)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Looks for off-targets for a single sgRNA')
	parser.add_argument('--genome_database', '-g', default=None)

	args = parser.parse_args()
	logger = logging.getLogger('CRISTA indexing bwa')
	init_commandline_logger(logger)
	process_assembly(args.genome_database, logger)
