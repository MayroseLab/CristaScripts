import re
import PA_limitedIndel_unite_mm as PA_script
import fetchSequenceFromGenome
from definitions import *

MAX_ALLOWED_GAPS = 3


def fasta_to_dictionary(file_path):
	"""Returns a dictionary that match each protein name to sequence in the input fasta file"""
	fasta_dic={}
	content = open(file_path, 'r').read().splitlines()

	for line in content:
		identifier = re.search("(?=\> ?).*(?=chr)", line)
		if identifier:
			iden = identifier.group()
		elif re.search("[actg]+", line, re.IGNORECASE):
			fasta_dic[iden] = line

	return fasta_dic


def best_pa(offtarget, match=MATCH_SCORE, mm=MISMATCH_PENALTY, gap=GAP_PENALTY, sgRNA=None):
	"""
	:param offtarget: offtarget class object
	:param offtarget if iit's a string, then contains PAM and 3 nucleotides to the left
	:return:
	"""
	if type(offtarget) is not str:

		dna_5prime_extension = fetchSequenceFromGenome.fetch_extended_5_prime_end(offtarget.genome_database,
								offtarget.chromosome, offtarget.start_position, offtarget.end_position,
								offtarget.strand, 3, offtarget.ontarget_name, offtarget.offtarget_seq)[:3]
		extended_offtarget_seq = dna_5prime_extension + offtarget.offtarget_seq
		sgRNA = offtarget.sgRNA_seq
	else:
		extended_offtarget_seq = offtarget

	original_sg_length = len(sgRNA)-3
	max_score = float("-inf")

	for i in [0, 6, 1, 5, 2, 4, 3]: #because starting at 3 is the original - we'd prefer that
		current_dna = extended_offtarget_seq[i:-3]
		(alnA, alnB, score) = PA_script.align_pair(seqA=sgRNA[:-3], seqB=current_dna, match_score=match, mismatch_score=mm, gap_score=gap, gaps_allowed=MAX_ALLOWED_GAPS) #regular pa

		if re.search("^\-", alnA) is None and score >= max_score:
			# the target can begin a '-', it means that the last nt of the sg is not paired
			# the sg cannot begin with a '-' - in that case, a better alignment would be found (shorter DNA).
			#   However, if we first found this target and then another with a different score- we'd prefer the other
			(alignmentA, alignmentB, max_score) = (alnA, alnB, score)

	aligned_sgRNA = alignmentA + sgRNA[-3:]
	aligned_offtarget = alignmentB + extended_offtarget_seq[-3:]

	return aligned_sgRNA, aligned_offtarget, max_score


def calculate_pa_score(seq1, seq2, match=MATCH_SCORE, mm=MISMATCH_PENALTY, gap=GAP_PENALTY):

	return PA_script.calculate_pa_score(seq1, seq2, match, mm, gap)


def count_matches(aligned_sgRNA, aligned_offtarget, ignore_PAM = False):
	cnt = 0
	if ignore_PAM:
		ending = len(aligned_offtarget) - 3
	else:
		ending = len(aligned_offtarget)

	for i in range(ending):
		cnt += int(aligned_offtarget[i] == aligned_sgRNA[i])
	return cnt


def count_mismatches(aligned_sgRNA, aligned_offtarget, ignore_PAM = False):
	cnt = 0
	if ignore_PAM:
		ending = len(aligned_offtarget) - 3
	else:
		ending = len(aligned_offtarget)

	for i in range(ending):
		cnt += int(aligned_offtarget[i] != aligned_sgRNA[i] and aligned_offtarget[i] != "-" and aligned_sgRNA[i] != "-")
	return cnt


def count_bulges_in_seq(aln):
	return aln.count("-")


def cnt_RNA_bulge(aligned_sgRNA, aligned_offtarget):
	return count_bulges_in_seq(aligned_offtarget)


def cnt_DNA_bulge(aligned_sgRNA, aligned_offtarget):
	return count_bulges_in_seq(aligned_sgRNA)

	
def count_bulges(seq1, seq2):
	cnt = 0
	ending = len(seq1) - 3
	for i in range(ending):
		cnt += seq1[i] == "-" or seq2[i] == "-"
	return cnt

	
def count_consecutive_inconsistencies(aligned_sgRNA, aligned_offtarget):
	cnt = 0
	current_cnt = 0

	for i in range(len(aligned_offtarget) - 3):
		if aligned_offtarget[i] != aligned_sgRNA[i]:
			current_cnt += 1
		else:
			cnt += current_cnt > 0
			current_cnt = 0
	return cnt