import sys
sys.path.append("/bioseq/crista/CRISTA_online/")
from main import check_input_vars, get_features, predict_on_df
import numpy as np

#######################################################################
#
#       CRISTA: for run from any other location on the cluster
#
#######################################################################


def predict_cleavage_score(sgseq_lst, extended29_genomic_seq_lst):
	"""
	:param sgseq_lst: list of 20 nucleotides sgRNA sequence
	:param extended29_genomic_seq_lst: list of complementary 3nt upstream + 23-nt target site + 3-nt downstream
	:return: cleavage score by CRISTA
	"""
	assert len(sgseq_lst) == len(extended29_genomic_seq_lst), "input lists are not of the same length!"
	features_mat = None

	for sgseq, extended29_genomic_seq, in sgseq_lst, extended29_genomic_seq_lst:
		rna_seq, dna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome, \
		strand, start_position, end_position, include_genomic_features, w_flanking = \
			check_input_vars(rna_seq=sgseq, extended100_genomic_seq=extended29_genomic_seq,
							 genome_database=None, cell_type=None,
							 chromosome=None, strand=None,
							 start_position=None)

		features = get_features(rna_seq, dna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome, strand,
								start_position, end_position, include_genomic_features, w_flanking=w_flanking)
		if features_mat is None:
			features_mat = np.asmatrix(features)
		else:
			features_mat = np.concatenate((features_mat, np.asmatrix(features)))

	score_df = predict_on_df(features_mat, include_genomic_features, w_flanking)
	return score_df["CRISTA score"]