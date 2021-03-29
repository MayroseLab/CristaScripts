import del_unused_genome_assemblies
from utils import *
from format_data import *
import argparse, glob
import regex as re
import time
import numpy as np
import web_utils
from format_data import format_csv_allfeatures_table
from utils import readlines, reverse_complement
import fetchSequenceFromGenome
from definitions import *
import pandas as pd
import PA_limitedIndel_unite_mm as PA_script
import main
from main import EXTENSION
import download_and_index_assembly

N_GAPS = 3
N_DISTEDIT = 4


def parse_sam_out_table(samtable_filepath, sgrna_seq, genomic_db="hg19"):

	sam_gen = readlines(samtable_filepath, ";")

	first_line = sam_gen.__next__()
	firsts = first_line.split("XA:Z:")
	(first_line, alter) = firsts[0], None if len(firsts)==1 else firsts[1]

	data_df = pd.DataFrame(columns=[CHROMOSOME_HEADING, STRAND_HEADING, STARTPOS_HEADING,
									ENDPOS_HEADING, SGRNA_HEADING, OFFTARGET_HEADING])
	#parse best match (first line)
	first_line_tokens = first_line.split("\t")
	chromosome_num = "chr" + first_line_tokens[2]
	start_pos = int(first_line_tokens[3])
	cigar_str = first_line_tokens[5]
	dna_seq = first_line_tokens[9]

	if chromosome_num == "chr*" and start_pos == 0:
		return None, None

	features_mat = None
	is_first_line = True

	while alter and not alter.isspace() or is_first_line:
		#parse token
		if not is_first_line:
			alter_tokens = alter.split(",")
			chromosome_num = alter_tokens[0]

			chromosome_num = "chr" + chromosome_num
			strand = alter_tokens[1][0]
			start_pos = int(alter_tokens[1][1:])
			cigar_str = alter_tokens[2]

		# fetch sequence
		cigar_tokens = re.findall("([0-9]+[MIDNSHPX\=])", cigar_str)
		ot_pattern = ""
		for token in cigar_tokens:
			ot_pattern += (token[-1]*int(token[:-1]))
			ot_len = ot_pattern.count("M") + ot_pattern.count("S") + ot_pattern.count("D")

		end_pos = start_pos + ot_len - 1

		extended100_genomic_seq = fetchSequenceFromGenome.fetch_dna_coordinates_offline(genomic_db, chromosome_num,
																						start_pos - EXTENSION,
																						end_pos + EXTENSION)

		if re.search("[^ACGT]", extended100_genomic_seq) or len(extended100_genomic_seq) < end_pos-start_pos+1+2*EXTENSION:
			if not is_first_line:
				alter = sam_gen.__next__()
			is_first_line = False
			continue

		if is_first_line:
			fwd_fuzzy_match = re.search("(" + dna_seq + "){i<=3,d<=3,s<=5,e<8}", sgrna_seq, re.I|re.BESTMATCH)
			rev_fuzzy_match = re.search("(" + reverse_complement(dna_seq) + "){i<=3,d<=3,s<=5,e<8}", sgrna_seq, re.I|re.BESTMATCH)
			fwd_fuzzy_cnt = float("inf") if fwd_fuzzy_match is None else sum(fwd_fuzzy_match.fuzzy_counts)
			rev_fuzzy_cnt = float("inf") if rev_fuzzy_match is None else sum(rev_fuzzy_match.fuzzy_counts)

			strand = "+" if fwd_fuzzy_cnt < rev_fuzzy_cnt else "-"

		extended100_genomic_seq = fetchSequenceFromGenome.get_seq_by_orientation(extended100_genomic_seq, strand)
		dna_seq = extended100_genomic_seq[EXTENSION:-1 * EXTENSION]

		PAM_chunk = dna_seq[-2:]
		if (PAM_chunk.count("G") + PAM_chunk.count("A")) < 2 or "random" in chromosome_num:
			if not is_first_line:
				alter = sam_gen.__next__()
			is_first_line = False
			continue

		alignmentA, alignmentB, score = PA_script.align_pair(seqA=sgrna_seq[:-3], seqB=dna_seq[:-3], match_score=MATCH_SCORE,
															 mismatch_score=MISMATCH_PENALTY, gap_score=GAP_PENALTY,
															 gaps_allowed=3)
		alignmentA += sgrna_seq[-3:]
		alignmentB += dna_seq[-3:]

		features = main.get_features(rna_seq=sgrna_seq, dna_seq=dna_seq, extended100_genomic_seq=extended100_genomic_seq,
					 genome_database=genomic_db, cell_type='negatives', chromosome=None if True else chromosome_num, strand=strand,
					 start_position=start_pos, end_position=end_pos, include_genomic_features=False, w_flanking=True,
					 aligned_sgRNA=alignmentA, aligned_offtarget=alignmentB, pa_score=score)

		if features_mat is None:
			features_mat = np.asmatrix(features)
		else:
			features_mat = np.concatenate((features_mat, np.asmatrix(features)))

		data_df.loc[data_df.shape[0]] = [chromosome_num[3:], strand, start_pos, end_pos, alignmentA, alignmentB]

		if not is_first_line:
			alter = sam_gen.__next__()

		is_first_line = False

	return data_df, features_mat


def run_bwa(sgRNA_seq, logger, genomic_db="hg19", n_gaps=N_GAPS,
			mm_penalty=MISMATCH_PENALTY*(-1), gap_open_penalty=GAP_PENALTY*(-1), gap_extension_penalty=0,
			exec_path=None):
	"""
	:param sgRNA_seq: 23nt
	:param ref_genome_path:
	:param n_mms:
	:param n_gaps:
	:return: bwa output table for
	"""
	if exec_path is None:
		exec_path = RESULTS_PATH
	ref_genome_path = fetchSequenceFromGenome.get_bwa_genome_path(genomic_db)
	sgrna_file = exec_path + "sgrna.fa"
	sai_output_file = sgrna_file + ".out.sai"
	sam_output_file = sgrna_file + ".out.sam"
	output_table_file = sgrna_file + ".out.table"
	# BWA_COMMAND_ON_SERVER = "/share/apps/bwa-0.7.15/bwa" #lecs2


	with open(sgrna_file, "w") as fpw:
		fpw.write(">seq\n" + sgRNA_seq)
	cmd = BWA_EXE + " aln "
	cmd += "-N " #don't stop once reaches the optimal match
	cmd += "-l 20 " #seed length
	cmd += "-i 0 " #Disallow an indel within INT bp towards the ends
	#cmd += "-e 3 " #Maximum number of gap extensions
	cmd += "-n {0} ".format(str(N_DISTEDIT+1)) #max number of mms in the whole region
	cmd += "-o {0} ".format(str(n_gaps)) #max number of gaps in the whole region
	cmd += "-d 3 " #Disallow a long deletion within INT bp towards the 3'-end
	cmd += "-k {0} ".format(str(N_DISTEDIT)) #Maximum edit distance in the seed
	cmd += "-M {0} ".format(str(0))
	cmd += "-O {0} ".format(str(1))
	cmd += "-E {0} ".format(str(gap_extension_penalty))
	cmd += ref_genome_path + " " + sgrna_file + " > " + sai_output_file
	os.system(cmd)
	if not os.path.exists(sai_output_file):
		raise Exception(IOError, "BWA samse failed, sai file not created")
	time.sleep(1)
	logger.info("BWA align - execution of bwa completed")

	# convert sai to sam
	os.system(BWA_EXE + " samse -n 1000000 {0} {1} {2} > {3}".format(ref_genome_path, sai_output_file, sgrna_file, sam_output_file))
	if not os.path.exists(sam_output_file):
		raise Exception(IOError, "BWA samse failed, sam file not created")
	time.sleep(1)
	logger.info("BWA samse - conversion from sai to sam completed!")

	# read sam file with samtools
	os.system(SAMTOOLS_EXE + " view -S " + sam_output_file + " > " + output_table_file)
	if not os.path.exists(output_table_file):
		raise Exception(IOError, "BWA samse failed, output table file not created")
	time.sleep(1)
	logger.info(SAMTOOLS_EXE + " view - readable sam file conversion completed")

	data_df, features_mat = parse_sam_out_table(output_table_file, sgRNA_seq, genomic_db=genomic_db)
	if features_mat is None:
		return None, None
	data_df.to_csv(exec_path + "data.csv")

	score_df = main.predict_on_df(features_mat, include_genomic_features=False, w_flanking=True, logger=logger)
	logger.info("CRISTA predictions completed")
	all_df = pd.concat((data_df, pd.DataFrame(features_mat[:,2:], columns=DATASET_HEADER_NONGENOMIC_MAT), score_df), axis=1)
	results_df = pd.concat((data_df, score_df), axis=1)

	remove_duplicates_from_df(results_df)
	remove_duplicates_from_df(all_df)

	logger.info("Formatting results")
	format_csv_allfeatures_table(all_df, logger)
	logger.info("Done")

	#main.update_output_html(output_df=all_df)
	#all_df.to_csv(TEMP_EXECUTION_DIR+"outscores.csv")
	return all_df, results_df


def update_output_html(sgrna=None, output_df=None, logger=None, isDailyTest=False):

	filename = sgrna + "_CRISTA_offtargets_scores.csv"
	output_df.to_csv(RESULTS_PATH + "/" + filename)
	html_addition = "<a href=\"" + filename + "\" download>Download results with full features values</a>\n"  # ref to file
	html_addition += "</p>\n" \
					"<br><div class=\"tooltip\">The score scale is in the interval [0,1].<span class=\"tooltiptext\">" \
					"The score represents the frequency of indels relative to a cleavage at the on-target of a highly " \
					"efficient sgRNA. </span></div><br>" \
					 "\n"\
					"" \
					"</font></p>\n"
	web_utils.add_to_results_html(RESULTS_PATH, html_addition)
	
	#send email
	if not isDailyTest:
		web_utils.send_last_email(RECIPIENT, RUN_NUMBER, logger=logger)
	else:
		web_utils.write_daily_test('crista_find_offtargets', RUN_NUMBER, 'PASS')
		


if __name__ == '__main__':
	logger = logging.getLogger('CRISTA main script')
	init_commandline_logger(logger)

	parser = argparse.ArgumentParser(description='Looks for off-targets for a single sgRNA')
	parser.add_argument('--sgseq', '-s', default=None, help="20bp long")
	parser.add_argument('--genome_database', '-g', default=None)
	parser.add_argument('--path', '-p', default=None)
	parser.add_argument('--run_number', '-n', default=None)
	parser.add_argument('--daily_test', action='store_true')
		
	args = parser.parse_args()
	global RESULTS_PATH
	global RUN_NUMBER
	RESULTS_PATH = args.path
	RUN_NUMBER = args.run_number
	RECIPIENT = web_utils.get_email_recipient(RESULTS_PATH)
	genome_assembly = args.genome_database

	try:
		if not args.daily_test:
			web_utils.send_first_email(RECIPIENT, RUN_NUMBER, logger=logger)
		sgrna_seq = re.sub("U", "T", args.sgseq.upper()) + "NGG"
		all_df, results_df = run_bwa(sgrna_seq, logger=logger, genomic_db=genome_assembly)
		if all_df is None:
			set_error_message("Matches were not found in the selected genome reference.")
			raise ValueError("Matches were not found in the selected genome reference.")
		sort_dataframe_by_score(results_df)
		update_output_html(sgrna_seq, all_df, logger, isDailyTest=args.daily_test)
		convert_pos_cols_to_int(results_df)
		web_utils.add_table_to_html_results_page(RESULTS_PATH, results_df)

	except:
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)
		if not args.daily_test: 
			web_utils.send_error_email(RECIPIENT, RUN_NUMBER, logger)
		else: 
			web_utils.write_daily_test('crista_find_offtargets', RUN_NUMBER, 'FAIL')
		web_utils.send_error_email("shiranos@gmail.com", RUN_NUMBER, logger)
		web_utils.load_error_page(RESULTS_PATH, get_error_message())
	web_utils.stop_refreshing(RESULTS_PATH)

	# parse_sam_out_table(r"D:\Dropbox\lab\CRISPR\data\results_files\sgrna.fa.out.table", sgrna_seq="GAGTTAGAGCAGAAGAAGAAAGG")
	del_unused_genome_assemblies.delete_unused_genome_assemblies()