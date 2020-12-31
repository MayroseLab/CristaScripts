from __future__ import division
import math, argparse, os, logging, sys, urllib.request, pickle, csv, traceback
import numpy as np
import scipy
import scipy.stats
from Bio import Entrez, SeqIO, Seq
import regex as re


def get_gene_from_ncbi(gene_accession, ncbi_db="nucleotide", rettype="gb"):
	"""
	:param gene_accession:
	:param db: "nucleotide" | "gene" | "genome" ...
	:return:
	"""
	Entrez.email = "A.N.Other@example.com"
	handle = Entrez.efetch(db=ncbi_db, id=gene_accession, rettype=rettype, retmode="text")
	record = SeqIO.read(handle, "genbank")
	handle.close()

	return record



def fasta_to_dictionary(file_path):
	"""Returns a dictionary that match each protein name to sequence in the input fasta file"""
	fasta_dic={}
	content = open(file_path, 'r').read().splitlines()

	for line in content:
		identifier = re.search("(?=\> ?).*(?=chr)", line)
		if identifier:
			identifier = identifier.group()
		else:
			fasta_dic[identifier] = line

	return fasta_dic


def change_path_permissions_to_777(path):
	os.chmod(path, 0o777)
	for root, dirs, files in os.walk(path):
		for dir in dirs:
			try:
				os.chmod(os.path.join(root, dir), 0o777)
			except:
				pass
		for file in files:
			try:
				os.chmod(os.path.join(root, file), 0o777)
			except:
				pass


def init_commandline_logger(logger):
	logger.setLevel(logging.DEBUG)
	# create console handler and set level to debug
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.DEBUG)
	# create formatter
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	# add formatter to ch
	ch.setFormatter(formatter)
	# add ch to logger
	logger.addHandler(ch)


def download_file(url, filename):
	if not os.path.exists(filename):
		f = urllib.request.urlretrieve(url, filename)
	return filename


def RMSE(observed, predictions):
	observed = np.asarray(observed, np.float)
	predictions = np.asarray(predictions, np.float)
	return np.sqrt(((predictions - observed) ** 2).mean())


def squared_pearson(observed, predicted):
	return scipy.stats.pearsonr(observed, predicted)[0]**2


def squared_spearman(observed, predicted):
	return scipy.stats.spearmanr(observed, predicted)[0]**2


def gc_content(seq):
	useq = seq.upper()
	return (useq.count("G") + useq.count("C"))/float(len(seq))


def get_avg(l):
	return sum(l)/float(len(l))


def get_var(l):
	avg = get_avg(l)
	dist_from_avg = list(map(lambda x: (x - avg)**2, l))
	var = sum(dist_from_avg)/len(dist_from_avg)
	return var


def get_std(l):
	var = get_var(l)
	std = math.sqrt(var)
	return std


def sum_of_squares(l):
	return sum([c**2 for c in l])


def interval_defaultdict():
	return intervalmap()


def remove_gaps_from_sequence(seq):
	return re.sub("[^a^g^c^t^A^G^C^T]+", "", seq, re.I)


def calculate_column_entropy(col_string):
	# column_entropy = - sum(for every nucleotide x) {count(x)*log2(Prob(nuc x in col i))}
	col_gapless = remove_gaps_from_sequence(col_string).upper()
	col_entropy = 0
	for x in ['A', 'G', 'C', 'T']:
		count_x = str.count(col_gapless, x)
		if count_x == 0:
			entropy_x = 0
		else:
			prob_x = count_x/len(col_gapless)
			entropy_x = count_x*math.log2(prob_x)
		col_entropy += entropy_x

	return -col_entropy


def get_msa_avg_entropy(msa):
	#msa is a Bio.AlignIO object
	msa_length = msa.get_alignment_length()
	sum_entropy = 0
	for col_i in range(0, msa_length):
		sum_entropy += calculate_column_entropy(msa.get_column(col_i))

	return sum_entropy/msa_length


def get_num_of_seqs_in_alignment(msa):
	return len(msa.get_column(0))


def readlines(f, newline):
	"""
	:param f: file or filepath
	:param newline: delimiter
	:return: a generator that reads the file as if line by line according to the delimiter

	:usage:
		with open('file') as f:
			for line in readlines(f, "."):
				print line
	"""
	if isinstance(f, str):
		f = open(f, "r")
	buf = ""
	while True:
		while newline in buf:
			pos = buf.index(newline)
			yield buf[:pos]
			buf = buf[pos + len(newline):]
		chunk = f.read(4096)
		if not chunk:
			yield buf
			break
		buf += chunk


def reverse_complement(seq):
	return Seq.reverse_complement(seq).upper()


def change_path_permissions_to_777(path):
	os.chmod(path, 0o777)
	for root, dirs, files in os.walk(path):
		for dir in dirs:
			os.chmod(os.path.join(root, dir), 0o777)
		for file in files:
			os.chmod(os.path.join(root, file), 0o777)