import bwa
import del_unused_genome_assemblies
from createJobFile import create_job_file
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
import shutil

GENOMES_PATH = "/groups/itay_mayrose/shiranabad/CRISPR/bwa_genomes/"
HG19_PATH = GENOMES_PATH + "hg19/human_g1k_v37.fasta.gz"

N_GAPS = 3
N_DISTEDIT = 4

TEMP_LINE_PATTERN = "<b>{}:</b> <font style=\"color:#808080\"><i>the results file will be available for download once the execution is accomplished.</i></font>"


def make_html(sgseq_lst, create_html="0", results_path="", run_number="", running_mode="", running_params=""):
	create_html = bool(int(create_html))
	if create_html:
		if not os.path.exists(results_path):
			global RESULTS_PATH
			RESULTS_PATH = "/bioseq/data/results/crista/" + run_number + "/"
			results_path = RESULTS_PATH
			if not os.path.exists(RESULTS_PATH):
				os.mkdir(RESULTS_PATH)
		web_utils.make_new_results_page(results_path, run_number, running_mode, running_params)

	sgline_pattern = TEMP_LINE_PATTERN + "\n<br><br>\n"
	html_addition = "\n\n"

	for sgrna in sgseq_lst:
		html_addition += sgline_pattern.format(sgrna)

	html_addition += "</p>\n" \
					 "<br><div class=\"tooltip\">The score scale is in the interval [0,1].<span class=\"tooltiptext\">" \
					 "The score represents the frequency of indels relative to a cleavage at the on-target of a highly " \
					 "efficient sgRNA. </span></div><br>" \
					 "\n" \
					 "" \
					 "</font></p>\n"
	web_utils.add_to_results_html(RESULTS_PATH, html_addition)
	# web_utils.send_last_email(RECIPIENT, RUN_NUMBER, logger=logger)


def run_single_sgRNA_job(sgrna_seq, genome_assembly, logger, run_number):
	current_path = RESULTS_PATH + sgrna_seq + "/"
	if not os.path.exists(current_path):
		os.mkdir(current_path)

	all_df, results_df = bwa.run_bwa(sgrna_seq, logger=logger, genomic_db=genome_assembly, exec_path=current_path)
	if all_df is None:
		link_text = "Matches were not found in the selected genome reference."
	else:
		sort_dataframe_by_score(results_df)
		filename = sgrna_seq + "_CRISTA_offtargets_scores.csv"
		all_df.to_csv(RESULTS_PATH + "/" + filename)
		web_utils.make_new_results_page(current_path, run_number, running_mode="Multiple sgRNAs - a single run", running_params="")
		convert_pos_cols_to_int(results_df)
		html_out = web_utils.add_table_to_html_results_page(current_path, results_df)
		shutil.copyfile(html_out, RESULTS_PATH + sgrna_seq + "_output.php")
		os.chmod(html_out, 777)
		html_url = 'http://crista.tau.ac.il/results/' + str(run_number) + "/" + sgrna_seq + "_output.php"
		link_text = "<a href=\"" + filename + "\" download>Download results with full features values</a>; " \
					"<a href=\"" + html_url + "\">or click here to view them on web</a>.\n"
	link_text = "<b>{}:</b> " + link_text
	sgseq = sgrna_seq[:-3]
	web_utils.replace_line_in_html(RESULTS_PATH, TEMP_LINE_PATTERN.format(sgseq), link_text.format(sgseq))

	#todo if all files appear (the temp line doesn't) send email that it's done


def qsub_jobs_for_single_sgrnas(sgseq_lst, genome_assembly, RESULTS_PATH, RUN_NUMBER):
	errors_path = RESULTS_PATH + "/jobs_errors/"
	if not os.path.exists(errors_path):
		os.mkdir(errors_path)
	change_path_permissions_to_777(RESULTS_PATH)

	for sgrna in sgseq_lst:
		job_file = create_job_file("CRISTA_sub_"+RUN_NUMBER + "_" + sgrna,
								   "python /bioseq/crista/CRISTA_online/multiple_sgRNAs_offtargets_search.py "
								   "-r {0} -g {1} -p {2} -n {3}".format(sgrna, genome_assembly, RESULTS_PATH, RUN_NUMBER),
								   "CRISTA_sub_" + RUN_NUMBER + "_" + sgrna + ".sh",
								   errors_path, errors_path)
		os.system("qsub " + job_file)


if __name__ == '__main__':

	logger = logging.getLogger('CRISTA main script')
	init_commandline_logger(logger)

	parser = argparse.ArgumentParser(description='Looks for off-targets for a single sgRNA')
	parser.add_argument('--sgseq_lst', '-s', default=None, help="20bp long sequences, joined by a comma/whitespace")
	parser.add_argument('--sgseq', '-r', default=None, help="20bp long") #for a single sgRNA run
	parser.add_argument('--genome_database', '-g', default=None)
	parser.add_argument('--create_html', '-m', default="0")
	parser.add_argument('--path', '-p', default="")
	parser.add_argument('--run_number', '-n', default=None)

	args = parser.parse_args()
	global RESULTS_PATH
	global RUN_NUMBER
	RESULTS_PATH = args.path
	RUN_NUMBER = args.run_number
	RECIPIENT = web_utils.get_email_recipient(RESULTS_PATH)
	genome_assembly = args.genome_database
	create_html = args.create_html

	try:
		web_utils.send_first_email(RECIPIENT, RUN_NUMBER, logger=logger)
		if args.sgseq is None:
			sgseq_lst = re.sub(",+", " ", args.sgseq_lst.strip())
			sgseq_lst = re.sub("\s+", " ", sgseq_lst.strip())
			sgseq_lst = set(re.split(" ", sgseq_lst.strip()))
			make_html(sgseq_lst, create_html, results_path=RESULTS_PATH, run_number=RUN_NUMBER,
					  running_mode="Detect possible off-targets throughout the genome for multiple sgRNAs",
					  running_params="sgRNAs: " + " ".join(sgseq_lst) + "<br><br>Genome assembly: " + genome_assembly)
			qsub_jobs_for_single_sgrnas(sgseq_lst, genome_assembly, RESULTS_PATH, RUN_NUMBER)
		else:
			sgrna_seq = args.sgseq.upper() + "NGG"
			run_single_sgRNA_job(sgrna_seq, genome_assembly, logger, RUN_NUMBER)

	except:
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)
		web_utils.send_error_email(RECIPIENT, RUN_NUMBER, logger=logger)
		web_utils.send_error_email("shiranos@gmail.com", RUN_NUMBER, logger)
		web_utils.load_error_page(RESULTS_PATH, get_error_message())

	web_utils.stop_refreshing(RESULTS_PATH)
	# parse_sam_out_table(r"D:\Dropbox\lab\CRISPR\data\results_files\sgrna.fa.out.table", sgrna_seq="GAGTTAGAGCAGAAGAAGAAAGG")
	del_unused_genome_assemblies.delete_unused_genome_assemblies()