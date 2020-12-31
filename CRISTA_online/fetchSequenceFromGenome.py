import glob
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import download_and_index_assembly
from definitions import *
from xml.dom.minidom import parseString
from utils import download_file
from subprocess import Popen, PIPE


############################# Fetch from online UCSC and save to cache ###############################################
def get_xml_filename(genome, chrom, startpos, endpos, cache_dir):
	current_filename = "{0}{1}_{2}_{3}_{4}.xml".format(cache_dir, genome, chrom, str(startpos), str(endpos))
	return current_filename


def get_dna_coordinates_xmlfile(genome, chromosome, startpos, endpos, cache_dir):
	url = "http://genome.ucsc.edu/cgi-bin/das/{0}/dna?segment={1}:{2},{3}" #chr15:65637530,65637553"
	current_filename = get_xml_filename(genome, chromosome, startpos, endpos, cache_dir)
	current_url = url.format(genome, chromosome, startpos, endpos)
	download_file(current_url, current_filename)
	return current_filename


def fetch_dna_coordinates(genome, chrom, startpos, endpos, cache_dir):
	current_filename = get_dna_coordinates_xmlfile(genome, chrom, startpos, endpos, cache_dir)
	with open(current_filename) as fp:
		xmldata = parseString(fp.read())
	seq = re.sub("\s", "", xmldata.childNodes[1].childNodes[1].childNodes[1].childNodes[0].data.upper())
	return seq
######################################################################################################################


def get_chop_filename(genome_assembly, chromosome, position):
	"""
	:param chromosome: relevant chromosome
	:param position: position in the chromosome
	:return: the chop filename and relative index in it
	"""
	### NOTE TO MYSELF: here the sequences begin at 0, in the online uscs browser at 1
	position -= 1
	chops_path = LOCAL_DB_DIRS[genome_assembly] + "/chops/"
	idx = position % CHOP_LENGTH
	chop_number = (position // CHOP_LENGTH) * CHOP_LENGTH
	filename = str.format(chops_path + "{0}.{1}", chromosome, str(chop_number))

	return filename, idx


def fetch_dna_coordinates_offline(genome_assembly, chrom, startpos, endpos, target_site=None):
	"""
	reads the sequence from the local copy of the genome assembly. assumes that the sequence does not spread over more than 2 chop files
	:return:
	"""
	if not platform.system() == 'Linux':
		return fetch_dna_coordinates(genome_assembly, target_site, chrom, startpos, endpos)
	seq = fetch_dna_coordinates_faidx(genome_assembly, chrom=chrom, start=startpos, end=endpos)
	if len(seq) == 0:
		# depends on how the samtools query file was created - if it was without "chr", then the query should be without it too
		seq = fetch_dna_coordinates_faidx(genome_assembly, chrom=chrom[3:], start=startpos, end=endpos)
	if len(seq) == 0:
		add_warning_message("Failed to fetch DNA coordinates for chromosome {}: {}-{}".format(chrom, startpos, endpos))
	if 1 < len(seq) < endpos-startpos+1:
		add_warning_message("Genomic coordinates exceed chromosome ends")
	return seq

	start_chop_filename, start_idx = get_chop_filename(genome_assembly, chrom, startpos)
	end_chop_filename, end_idx = get_chop_filename(genome_assembly, chrom, endpos)
	seq_len = endpos - startpos + 1

	if start_chop_filename == end_chop_filename:
		with open(start_chop_filename) as fp:
			fp.seek(start_idx, 0)
			seq = fp.read(seq_len)
	else:
		with open(start_chop_filename) as fp:
			fp.seek(start_idx, 0)
			seq = fp.read(CHOP_LENGTH - start_idx)  # the end of the file contains extra few nucleotides, so I have to calculate exactly how many

		with open(end_chop_filename) as fp:
			seq += (fp.read(seq_len - len(seq)))  # read the rest

	return seq.upper()


def get_bwa_genome_path(genome_assembly):
	assembly_dir = GENOMES_PATH + genome_assembly
	genome_path = glob.glob(assembly_dir + "/*.fa*.gz")
	if len(genome_path) == 0:
		download_and_index_assembly.process_assembly(genome_assembly)
		genome_path = glob.glob(assembly_dir + "/*.fa*.gz")
	else:
		download_and_index_assembly.write_current_timestamp_for_assembly(GENOMES_PATH + genome_assembly + SEP)

	return genome_path[0]


def fetch_dna_coordinates_faidx(genome_assembly, chrom, start, end):
	"""Call tabix and generate an array of strings for each line it returns."""
	filename = get_bwa_genome_path(genome_assembly)
	if not os.path.exists(filename):
		raise ValueError("Genome assembly {} is not valid. Please contact us for assistance.".format(genome_assembly))
	query = '{}:{}-{}'.format(chrom, start, end)
	#print(SAMTOOLS_EXE + ' faidx ', filename, " ", query)
	process = Popen([SAMTOOLS_EXE + ' faidx ' + filename + " " + query], stdout=PIPE, shell=True)
	seq = ""
	lines = process.stdout.readlines()
	if len(lines) > 1:
		for line in lines[1:]:
			seq += line.decode('ascii').strip()

	return seq.upper()


def fetch_many_coordinates(genome_assembly, list_of_regions):
	#todo: not verified
	"""
	:param genome_assembly:
	:param list_of_regions: of items: (chromosome, start, end)
	:return:
	"""
	regions_str = ""
	for (chrom, start, end) in list_of_regions:
		regions_str += '{}:{}-{} '.format(chrom, start, end)
	filename = "/groups/itay_mayrose/shiranabad/CRISPR/bwa_genomes/{0}/{0}.fa.gz".format(genome_assembly)

	process = Popen([SAMTOOLS_EXE + ' faidx ' + filename + " " + regions_str], stdout=PIPE, shell=True)

	lines = process.stdout.readlines()
	res = {}
	input_counter = -1
	if len(lines) > 1:
		for line in lines:
			if line.startswith(">"):
				input_counter += 1
				res[list_of_regions[input_counter]] = ""
			else:
				res[list_of_regions[input_counter]] += line.decode('ascii').strip().upper()
	return res


def get_seq_by_orientation(seq, strand):
	"""
		returns the original sequence if it's the plus strand, o/w returns the reverse-complement
	"""
	if strand == "-":
		seqobject = Seq(seq, generic_dna)
		seq = str(Seq.reverse_complement(seqobject))
	return seq


def find_true_coordinates(seq, chromosome, genome_db, optional_start_coordinate, optional_end_coordinate=0, strand=None, offset=100):
	if optional_end_coordinate == 0:
		optional_end_coordinate = optional_start_coordinate
	if type(optional_start_coordinate) is str:
		optional_start_coordinate = int(optional_start_coordinate)
	if type(optional_end_coordinate) is str:
		optional_end_coordinate = int(optional_end_coordinate)
	start_lookup_position = max(optional_start_coordinate - offset, 0)
	end_lookup_position = optional_end_coordinate + offset
	is_strand_known = strand is not None

	# look in the forward and reverse orientation if needed
	candidate_fwd_match = candidate_rev_match = None
	if strand == "+" or not is_strand_known:
		seq = get_seq_by_orientation(seq, "+")
		fwd_surrounding_seq = fetch_dna_coordinates_offline(genome_db, chromosome, start_lookup_position, end_lookup_position)
		candidate_fwd_match = re.search(seq, fwd_surrounding_seq, re.I)
	if strand == "-" or not is_strand_known:
		seq = get_seq_by_orientation(seq, "-")
		rev_surrounding_seq = fetch_dna_coordinates_offline(genome_db, chromosome, start_lookup_position, end_lookup_position)
		candidate_rev_match = re.search(seq, rev_surrounding_seq, re.I)

	if candidate_fwd_match is None and candidate_rev_match is None:
		return None
	elif candidate_fwd_match is not None:
		return (start_lookup_position + candidate_fwd_match.start(), start_lookup_position + candidate_fwd_match.end()-1, "+", fwd_surrounding_seq)
	elif candidate_rev_match is not None:
		return (start_lookup_position + candidate_rev_match.start(), start_lookup_position + candidate_rev_match.end()-1, "-", rev_surrounding_seq)


def fetch_extended_5_prime_end(genome, chromosome, start_position, end_position, strand, extension, target_site, offtarget_seq):

	rev_strand = strand == "-"
	if rev_strand:
		end_position += extension
	else:
		start_position -= extension

	seq = fetch_dna_coordinates_offline(genome, chromosome, start_position, end_position)

	#if it's the minus strand, fetch the reverse complement
	seq = get_seq_by_orientation(seq, strand)

	if seq[3:].upper() != offtarget_seq.upper():
		print("Inconsistency:", target_site, "myDB:", offtarget_seq, "UCSC:", seq[3:].upper())

	return seq
