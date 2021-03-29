import pandas as pd

import del_unused_genome_assemblies
import main
import web_utils
from main import EXTENSION
from definitions import *
from format_data import format_csv_allfeatures_table, sort_dataframe_by_score, convert_pos_cols_to_int
from utils import *
import regex as re

def search_sites_in_seq(ref_seq, is_fwd):
	"""
	:param ref_seq: reference dna
	:param sgseq: SG sequence without NGG
	:param is_fwd: indicates if ref_seq is on the + or - strand
	:return: list of potential hits that end with GG or AG PAMs - the start_pos and end_pos are if the ref-seq is 1-indexed
	"""
	SG_LEN = 20
	pam_iter = re.finditer(r"[AGCT]{20}.[AG]G", ref_seq, re.I, overlapped=True)
	out_list = []  # each item is a tuple of: (dna_seq, relative_start_pos, relative_end_pos)

	for hit in pam_iter:
		hit_end = hit.end(0)
		if hit_end < EXTENSION + SG_LEN + 3:
			continue
		if hit_end > len(ref_seq) - EXTENSION:
			break

		dna_seq = ref_seq[hit_end - (SG_LEN + 3):hit_end]

		if is_fwd:
			relative_start_pos, relative_end_pos = hit_end - (SG_LEN + 2), hit_end
		else:
			relative_start_pos, relative_end_pos = len(ref_seq) - hit_end + 1, len(ref_seq) - hit_end + SG_LEN + 3
		out_list.append((dna_seq, relative_start_pos, relative_end_pos))
	#print(out_list)
	return out_list


def get_targets_features(ref_seq): #genome_db, chromosome_num=None, strand=None, ref_seq_start_ingenome=0
	plus_sites_list = search_sites_in_seq(ref_seq, True)
	minus_sites_list = search_sites_in_seq(reverse_complement(ref_seq), False)
	sites_list = list(zip("+"*len(plus_sites_list), plus_sites_list)) + \
				 list(zip("-"*len(minus_sites_list), minus_sites_list))

	if len(sites_list) == 0:
		set_error_message("No available targets were found! Please check the given input.")
		raise(Exception)

	features_mat = None
	data_df = pd.DataFrame(columns=[SGRNA_HEADING, OFFTARGET_HEADING, STRAND_HEADING,
									STARTPOS_HEADING, ENDPOS_HEADING])

	for strand, site_record in sites_list:
		dna_seq, relative_start_pos, relative_end_pos = site_record
		sgrna_seq = dna_seq[:-3] + "NGG"
		start_pos = relative_start_pos  #ref_seq_start_ingenome+relative_start_pos
		end_pos = relative_end_pos  #ref_seq_start_ingenome+relative_end_pos
		strand = strand # need to correct when input strand is given
		features = main.get_features(rna_seq=sgrna_seq, dna_seq=dna_seq,
									 extended100_genomic_seq=ref_seq[relative_start_pos-1 - EXTENSION:relative_end_pos + EXTENSION],
									 genome_database="hg19", cell_type='negatives',
									 chromosome=None if True else None, strand=strand,
									 start_position=start_pos, end_position=end_pos,
									 include_genomic_features=False, w_flanking=False,
									 aligned_sgRNA=sgrna_seq, aligned_offtarget=dna_seq, pa_score=len(dna_seq)-3)

		if features_mat is None:
			features_mat = np.asmatrix(features)
		else:
			features_mat = np.concatenate((features_mat, np.asmatrix(features)))

		data_df.loc[data_df.shape[0]] = [sgrna_seq, dna_seq, strand, start_pos, end_pos]

	return data_df, features_mat


def get_seq_from_gene(gene_annotation):
	# todo: get sequence from gene annotation, find in exons
	pass


def update_output_html(output_df, results_df, logger, isDailyTest=False):
	filename = "CRISTA_offtargets_scores{0}.csv".format(RUN_NUMBER)
	output_df.to_csv(RESULTS_PATH + "/" + filename)
	html_addition = "<a href=\"" + filename + "\" download>Download results with full features values</a>\n"  # ref to file
	html_addition += "</p>\n" \
					"<br><div class=\"tooltip\">The score scale is in the interval [0,1].<span class=\"tooltiptext\">" \
					"The score represents the frequency of indels relative to a cleavage at the on-target of a highly " \
					"efficient sgRNA. </span></div><br>" \
					 "\n"
	html_addition += "<font face=Verdana size=2>\n" \
					 "<br><b>Click the relevant sgRNAs for finding off-targets.<br></b>\n" \
					 "In the next window you will be able to edit your selection and select the reference genome.</b></font><br>\n" \
					 "<ul id=\"dynamic-list\"></ul>\n" \
					 "<div id=\"search_button\" class=\"hidden\">\n" \
					 "<input type=submit value=\"Take me to search off-targets for my selection!\"\n" \
					 "style=\"font-weight:bold;height: 35px;font-size:14\"\n" \
					 "onclick=\"javascript:send_sgrna_lst_to_multiofftargets_form()\">\n" \
					 "</div>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n" \
					 "<br><br>\n"
	web_utils.add_to_results_html(RESULTS_PATH, html_addition)
	web_utils.add_table_to_html_results_page(RESULTS_PATH, results_df, True)
	
	#send email
	if not isDailyTest:
		web_utils.send_last_email(RECIPIENT, RUN_NUMBER, logger=logger)
	else:
		web_utils.write_daily_test('crista_in_sequence', RUN_NUMBER, 'PASS')
		


if __name__ == '__main__':
	if False:
		global EXTENSION
		EXTENSION = 3
		ref_seq = "ATGGCTTGTTGGCCTCAGCTGAGGTTGCTGCTGTGGAAGAACCTCACTTTCAGAAGAAGACAAACAGTAAGCTTGGGTTTTTCAGCAGCGGGGGGTTCTCTCAT" \
				  "TTTTTCTTTGTGGTTTTGAGTTGGGGATTGGAGGAGGGAGGGAGGGAAGGAAGCTGTGTTGGTTTTCACACAGGGATTGATGGAATCTGGCTCTTATGGACACA" \
				  "GGACTGTGTGGTCCGGATATGGCATGTGGCTTATCATAGAGGGCAGATTTGCAGCCAGGTAGAAATAGTAGCTTTGGTTTGTGCTACTGCCCAGGCATGAGTTC" \
				  "TGATCCCTAGGACCTGGCTCCGAATCGCCCCTGAGCACCCCACTTTTTCCTTTTGCTGCAGCCCTGGGAGCCACCTGGCTCTCCAAAAGCCCCTAATGGGCCCC" \
				  "TGTATTTCTGGAAGCTGTGGGTGAAGTGAGTTAGTGGCCCCACTCTTAGAGATCAATACTGGGTATCTTGGTGTCAATCTGGATTCTTTCCTTCAGGCCTGGAG" \
				  "GAATATAATAACTGAGACTTGTTTTATTTCTGCAGAGGGTTCTAAGCCATTCACTTCCCAGATGGGCCAATAATGCTTTGAGTAATCTGGAGATCATCTTTAAT" \
				  "GCGCAGGTGAATGGAACTCTTCCACAGAGGGATGTGAGGGCTGTAGAGCAGAGTGAACTCCCTGAAACTCAGACGTCAGCTCTTTGTCTCTCTATCTCTGAACA" \
				  "CCCTTCCTTAGAGATCCCATCTCTAGGATGCATTTCTCTGTAGTTAGTTTCTAAGTCTCTTGTTCCTGTTCTGCCTTTATTTTTTTTTCCTGGATTCTAAGCCA" \
				  "GTATCCCCACTTGGCTGTCTTAATGTAGCTTAACATGTCTGTAATCAAAATGATCATCTTTCTGAGATTCAAAGGGCTATAAGGGACTTTGGAGAGAATTTCAT" \
				  "TCAGTTTTCCTCAAACTAGAATAATGCTTGCACTGTCTGTAAAAGAACAAAAGTGTCAAAGCATCCTTTTGTTCACTAAATTTCCTTTTTTATTATAGTGTTAC" \
				  "TTAAATATTAGGAAGTAAAAGTAGGTATAAACTTCTTATAGGCTGTTATTATACAACTATATGACCCATACATATTTACAAATTAAGTGCAGCCAAAATTGCAA" \
				  "AATCAATACCATTCAAATTAATACCTTAAATGTGGTGAGGCAGCTGTTGTTCAACTGAAACCAAATTATAAGTTGCATGGCAGTAAATGCTATCATGCTGATCA" \
				  "TTTTGAGTTTGGCCAGTCTATATTATCATGTGCTAATGATTTGTGGTCCGGATATACCCATGAAGGAATTCTCCACCCATTTTTCTACTTGTATGACCTTAATT" \
				  "GAGTTTGCTACAATTATACTGGTGCCAACACAATCATAAACACAAATATAAACTTGGGCTTTGAAATCTTGTGCCAGAACTTGGCTTTAAAGTAAGCATTTAAA" \
				  "AAATCCATATGTGTTTATTAGACTTTGTTTAGATGACTGTTGAAATGAAAACAAAGTGTTTAAAATCCTCTTAGAGAACTTAAATATAATCCCTCAGCAATATG" \
				  "TATACAGATCTTCCTTTGAGAAAAACTGATTGTGTTCAGCCTCTCATGTTACAAATGGGGAACCTGAATTCTGAGGTCTCTAGTGAGAGAACAGGGACTGGAAT" \
				  "CTGTGGATCCTATCTGTTTTAATAATAATTGTAAAGTATAATAGATAATATTATATTAATAAAATAAAAGCAAACACTTAGAATGAGCTTCCATGTGTGAGGCA" \
				  "CTAACTGATTAGGCATTATTAACTAGATTTATTCCTTTTAAGGCCCCGCGATGTACTGTTATTTCCACATGTTGTAGCTGGGGAACGTGCTACTCAGAGAGGTT" \
				  "AAGTAACTTGTCTGAGGTCCACACCACTAACAAGGAGCACAGGTATGTGGTCCGGATATGAGCATGTGGGGGTTCAAATCCAGATAATCTGACTTTGGAGCTGG" \
				  "GCTTTTCAGTGGTGTCATTATTTTGCCTATTCTCCATCTGAGAATATTGAAGTTTCTGACTCCTTCCTTGCCTTTCTCCCTGCCTCCCGTGGTTATCCCCAGGT" \
				  "CTTGGTGTTCCAGTCCTCTATGTCCGTCCTTACTCTTATTCCTTTGCTACAGTGTGATCCAGGGCTCCTGCCCCTTCTTATCCTGGTAGAGGGGGCCCACTTGC" \
				  "TGGGAAATTGTCTCCGCCATGGTTTATCCATGTTGTGTGTCCATTAGTGAGTAGTG"
		print(get_targets_features(ref_seq))
		#print(search_sites_in_seq(reverse_complement(ref_seq), False))
		exit()

	logger = logging.getLogger('CRISTA main script')
	init_commandline_logger(logger)

	parser = argparse.ArgumentParser(description='find optional targets in gene or sequence and run feature collection for every sgrna')
	parser.add_argument('--genomic_seq_file', '-f', default=None)
	parser.add_argument('--genome_database', '-d', default=None)
	parser.add_argument('--gene_annotation', '-r', default=None)
	#parser.add_argument('--cell_type', '-e', default=None)
	parser.add_argument('--path', '-p', default=None)
	parser.add_argument('--run_number', '-n', default=None)
	parser.add_argument('--daily_test', action='store_true')

	args = parser.parse_args()
	global RESULTS_PATH
	global RUN_NUMBER
	EXTENSION = 3
	RESULTS_PATH = args.path
	RUN_NUMBER = args.run_number
	RECIPIENT = web_utils.get_email_recipient(RESULTS_PATH)

	try:
		if not args.daily_test:
			web_utils.send_first_email(RECIPIENT, RUN_NUMBER, logger=logger)
		with open(args.genomic_seq_file) as fpr:
			genomic_seq = fpr.read().strip().upper()
		data_df, features_mat = get_targets_features(genomic_seq)
		data_df.to_csv(RESULTS_PATH + "data.csv")

		score_df = main.predict_on_df(features_mat, include_genomic_features=False, w_flanking=False, logger=logger)
		all_df = pd.concat((data_df, pd.DataFrame(features_mat[:,2:], columns=DATASET_HEADER_NONGENOMIC_MAT), score_df), axis=1)
		results_df = pd.concat((data_df, score_df), axis=1)
		format_csv_allfeatures_table(all_df, logger)
		sort_dataframe_by_score(results_df)
		convert_pos_cols_to_int(results_df)

		#main.update_output_html(output_df=all_df)
		update_output_html(all_df, results_df, logger, isDailyTest=args.daily_test)

	except:
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)
		if not args.daily_test: 
			web_utils.send_error_email(RECIPIENT, RUN_NUMBER, logger=logger)
		else: 
			web_utils.write_daily_test('crista_in_sequence', RUN_NUMBER, 'FAIL')
		web_utils.send_error_email("shiranos@gmail.com", RUN_NUMBER, logger)
		web_utils.load_error_page(RESULTS_PATH, get_error_message())
	web_utils.stop_refreshing(RESULTS_PATH)
	del_unused_genome_assemblies.delete_unused_genome_assemblies()