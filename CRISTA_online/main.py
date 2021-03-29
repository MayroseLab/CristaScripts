import random
import sys
import argparse
import getDNAshape, get_genomic_features
import web_utils
from definitions import *
import targets_pa
import fetchSequenceFromGenome
import del_unused_genome_assemblies
from utils import *
import numpy as np
import pandas as pd
from format_data import *


global EXTENSION
EXTENSION = 100
SEND_EMAIL_SCRIPT = "/bioseq/pupkoSVN/trunk/www/bioSequence_scripts_and_constants/sendEmail.pl"
SERVER_CONSTANTS_FILE = "/bioseq/pupkoSVN/trunk/www/bioSequence_scripts_and_constants/GENERAL_CONSTANTS.pm"
global ERROR_MESSAGE_ON_EXCEPTION
ERROR_MESSAGE_ON_EXCEPTION = ""
global RESULTS_PATH
global RUN_NUMBER
global RECIPIENT


RF_PREDICTORS = None


def get_DNA_enthalpy(dna_seq):
	seq = re.sub("[^ACGT]", lambda x: random.choice(["A", "C", "G", "T"]), dna_seq)
	return sum([DNA_PAIRS_THERMODYNAMICS[seq[i-1:i+1]] for i in range(1, len(seq))])


def get_features(rna_seq, dna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome,
				 strand, start_position, end_position, include_genomic_features,
				 w_flanking=True,
				 aligned_sgRNA=None, aligned_offtarget=None, pa_score=None):
	"""
	correct extended100_genomic_seq if not in the right orientation
	:return: list of collecting_features
	"""
	global EXTENSION
	if not w_flanking:
		EXTENSION = 3

	############ DONE assertions and definitions

	## align sequences
	site_length = 23

	if not dna_seq in extended100_genomic_seq:
		extended100_genomic_seq = reverse_complement(extended100_genomic_seq)

	if aligned_sgRNA is None:
		aligned_sgRNA, aligned_offtarget, pa_score = \
			targets_pa.best_pa(extended100_genomic_seq[EXTENSION-3:EXTENSION+site_length], sgRNA=rna_seq)

	# get alignment features
	pa_score_17bp = targets_pa.calculate_pa_score(aligned_offtarget[-20:], aligned_sgRNA[-20:],
	                                              MATCH_SCORE, MISMATCH_PENALTY, GAP_PENALTY)
	mms_cnt = targets_pa.count_mismatches(aligned_sgRNA, aligned_offtarget, True)
	rna_bulges = targets_pa.cnt_RNA_bulge(aligned_sgRNA, aligned_offtarget)
	dna_bulges = targets_pa.cnt_DNA_bulge(aligned_sgRNA, aligned_offtarget)
	gapless_dnaseq = re.sub("-", "", aligned_offtarget)

	# quartets mismatches counts
	rev_rna = (aligned_sgRNA[::-1])[3:]
	rev_dna = (aligned_offtarget[::-1])[3:]
	mismatches_1_4 = targets_pa.count_mismatches(rev_rna[:4], rev_dna[:4], False)
	mismatches_5_8 = targets_pa.count_mismatches(rev_rna[4:8], rev_dna[4:8], False)
	mismatches_9_12 = targets_pa.count_mismatches(rev_rna[8:12], rev_dna[8:12], False)
	mismatches_13_16 = targets_pa.count_mismatches(rev_rna[12:16], rev_dna[12:16], False)
	mismatches_17_end = targets_pa.count_mismatches(rev_rna[16:], rev_dna[16:], False)

	# get from alignment mismatches per position
	mismatches = [-2] * 23  # 5' -> 3', without PAM                                             # undefined
	offset = 26 - len(aligned_offtarget)
	for i in range(len(aligned_offtarget) - 3):
		rna_base = aligned_sgRNA[i]
		dna_base = aligned_offtarget[i]
		# Categorization of mismatch type
		if rna_base == dna_base:                                                                # match
			mismatches[offset + i] = 0
		elif rna_base == "-" or dna_base == "-":                                                # indel
			mismatches[offset + i] = -1
		elif (rna_base == "T" and dna_base == "C") or (rna_base == "G" and dna_base == "A"):    # wobble: rG:dA, rT:dC
			mismatches[offset + i] = 1
		elif rna_base in ["G", "A"] and dna_base in ["T", "C"]:                                 # R-R pairing
			mismatches[offset + i] = 2
		elif rna_base in ["C", "T"] and dna_base in ["G", "A"]:                                 # Y-Y pairing
			mismatches[offset + i] = 3
		elif (rna_base == "A" and dna_base == "G") or (rna_base == "C" and dna_base == "T"):    # other transversion
			mismatches[offset + i] = 4

	# total types mismatches
	wobble_total = mismatches.count(1)
	RR_total = mismatches.count(2)
	YY_total = mismatches.count(3)
	Tv_total = mismatches.count(4)

	# pairs of nucleotides in positions 1-5 upstream to PAM (1-2, 2-3, 3-4, 4-5)
	seed_couples = []
	for i in range(4):
		seed_couples.append(agct2numerals(dna_seq[-8 + i:-6 + i]))

	# PAM and 5'-end nucleotides
	PAM_2_first = dna_seq[-2:]
	PAM_N_id = dna_seq[-3]
	last_pos_nucleotide = agct2numerals(gapless_dnaseq[0])

	# mismatches and bulges - linked
	consecutive_inconsistencies_cnt = targets_pa.count_consecutive_inconsistencies(aligned_sgRNA, aligned_offtarget)
	avg_inconsistency_length = (mms_cnt + rna_bulges + dna_bulges) / float(
		consecutive_inconsistencies_cnt) if consecutive_inconsistencies_cnt > 0 else 0

	# nucleotides occupancies in DNA target sequence
	nA = gapless_dnaseq.count("A")
	nC = gapless_dnaseq.count("C")
	nG = gapless_dnaseq.count("G")
	nT = gapless_dnaseq.count("T")

	# GC content
	extended_genomic_GC_content = (extended100_genomic_seq.count("C") + extended100_genomic_seq.count("G")) / float(
		len(extended100_genomic_seq))

	# five nucleotides downstream to PAM
	if w_flanking:
		ext = 5
	else:
		ext = 3
	nucleotides_down_pam = [agct2numerals(c) for c in extended100_genomic_seq[EXTENSION+site_length:EXTENSION + site_length + ext]]
	for i in range(5-ext):
		nucleotides_down_pam.append(0)

	# geometry features: dna_enthalpy
	extended100_dna_enthalpy = get_DNA_enthalpy(extended100_genomic_seq)
	dna_enthalpy = get_DNA_enthalpy(dna_seq)

	# geometry features: DNA shape per pentamer
	#print(extended100_genomic_seq[EXTENSION - 3:EXTENSION + site_length + 3])
	#print(extended100_genomic_seq[EXTENSION - 3:EXTENSION + site_length + 6])
	dna_shape_features = getDNAshape.get_DNAshape_features(
		extended100_genomic_seq[EXTENSION - 3:EXTENSION + site_length + 3], True)

	# genomic features
	genomic_features = []
	if include_genomic_features:
		dist_from_centromere = get_genomic_features.get_distance_from_centromere(genome_database, chromosome,
		                                                                         start_position, end_position)
		dist_from_telomere = get_genomic_features.get_distance_from_telomere(genome_database, chromosome,
		                                                                     start_position, end_position)

		DHS_signalValue, distance_from_DHS = get_genomic_features.get_DHS_signalVal(genome_database, chromosome,
		                                                                            start_position, end_position,
		                                                                            cell_type)
		genomic_subcompartment = get_genomic_features.get_subcompartments(genome_database, chromosome,
		                                                                  start_position, end_position)

		in_tx, in_cds = get_genomic_features.get_gene_properties(genome_database, chromosome,
		                                                         start_position, end_position)
		NGG_exon_signalValue, opposite_exon_signalValue = \
			get_genomic_features.get_exon_signalVal(genome_database, chromosome, strand,
			                                        start_position, end_position, cell_type)
		in_NGG_exon, in_opposite_exon = \
			get_genomic_features.get_exon_signalVal(genome_database, chromosome, strand,
			                                        start_position, end_position)

		nucleosome_mapping, length_from_nucleosome = get_genomic_features.is_located_in_nucleosome(
			genome_database, chromosome,
			start_position, end_position, cell_type)

		genomic_features = 	[dist_from_centromere, dist_from_telomere, DHS_signalValue, distance_from_DHS,
		                       nucleosome_mapping, length_from_nucleosome, genomic_subcompartment,
		                       NGG_exon_signalValue, opposite_exon_signalValue,
		                       in_NGG_exon, in_opposite_exon, in_tx, in_cds]

	#[db_source, ontarget_name, kind, rna_seq, dna_seq,
	 #strand, str(start_position), str(end_position),

	if chromosome is None:
		chromosome_num = []
	else:
		chromosome_num = [get_chromosome_number(chromosome)]


	sources = [aligned_sgRNA, aligned_offtarget] + chromosome_num
	features = [pa_score, rna_bulges + dna_bulges, rna_bulges, dna_bulges, agct2numerals(PAM_2_first),
	        agct2numerals(PAM_N_id), last_pos_nucleotide, mms_cnt,
	        consecutive_inconsistencies_cnt, avg_inconsistency_length] \
	       + [mismatches_1_4, mismatches_5_8, mismatches_9_12, mismatches_13_16, mismatches_17_end] + \
	       [wobble_total, YY_total, RR_total, Tv_total] + \
	       seed_couples + [agct2numerals(dna_seq[i]) for i in range(-8, -3)] + \
	       [extended_genomic_GC_content, #upstream_50_extension_gc, downstream_50_extension_gc,
	        dna_enthalpy, extended100_dna_enthalpy,
	        nA, nC, nT, nG] + nucleotides_down_pam + \
	       [min(dna_shape_features["MGW"])] + \
	       [get_avg(dna_shape_features["HelT"])] + \
	       [get_avg(dna_shape_features["Roll"])] + \
	       [get_avg(dna_shape_features["ProT"])] + \
	       [dna_shape_features["MGW"][-3], dna_shape_features["HelT"][-3],
	        dna_shape_features["Roll"][-3], dna_shape_features["ProT"][-3]] + \
	       genomic_features
	return sources+features


def predict_on_df(features_mat, include_genomic_features, w_flanking, logger=None):
	"""
	:param features_df: dataframe: first col: rna, second: dna, the rest are features
	mode: either full, nogenomic or noflanking
	:return: features df + prediction col
	"""
	n_predictors = 10 #100
	global RF_PREDICTORS

	if RF_PREDICTORS is None:
		if include_genomic_features:
			print("genomic full")
			path = RF_MODEL_DIR
		elif w_flanking:
			print("with flanking")
			path = RF_MODEL_DIR_WOGENOMIC
		else:
			print("wo flanking")
			path = RF_MODEL_DIR_WOFLANKING

		predictors = [path + str(i)+ "/RFpredictor.pkl" for i in range(n_predictors)]
		RF_PREDICTORS = [pickle.load(open(predictors[i], "rb")) for i in range(n_predictors)]

	prediction_df = pd.DataFrame()
	features_mat = np.delete(features_mat, [0,1], axis=1)

	for i in range(n_predictors):
		rf_predictor = RF_PREDICTORS[i]
		prediction_df[i] = rf_predictor.predict(features_mat)
		if i % 10 == 9 and logger:
			logger.info(str(i + 1) + "% are done")
	final_pred_df = pd.DataFrame(np.divide(prediction_df.mean(axis=1), 8.22), columns=["CRISTA score"])
	final_pred_df.loc[final_pred_df["CRISTA score"]>1, "CRISTA score"] = 1.0

	return final_pred_df


def update_output_html(sgrna=None, dna=None, score=None, output_df=None, logger=None, isDailyTest=False):

	html_path = RESULTS_PATH + "/output.php"

	with open(html_path) as fpr:
		html_content = fpr.read()

	html_content += "\n<p><font face=Verdana size=4>\n" \
	                "<br>-------------------------------------------------------------------------------------------------\n" \
	                "<br><br>\n" \
	                "Results:<br>\n" \
	                "</font></p>\n\n" \
	                "<p><font face=Verdana size=2>\n"
	header_message = "Completed "
	# Check if there are warnings:
	if len(WARNING_MESSAGES) > 0:
		html_content += "<p><font color=\"red\">\n" \
		                "Warnings found during the execution:<br>" + get_warning_message() +\
		                "</font></p>\n\n"
		header_message += "with warnings."
	else:
		header_message += "successfully."

	if output_df is None:
		html_content += "Genomic target site: " + dna + "<br>\n" \
	                "sgRNA: " + sgrna + "<br>\n" \
	                "Cleavage score: <font color='red'>{0:.2f}</font><br>\n".format(score)
	else:
		output_df.to_csv(RESULTS_PATH + "/CRISTA_scores.csv", index=False)
		html_content += "<a href=\"CRISTA_scores.csv\" download>Download results with full features values</a>\n" #ref to file

	html_content += "</p>\n" \
					"<br><div class=\"tooltip\">The score scale is in the interval [0,1].<span class=\"tooltiptext\">" \
					"The score represents the frequency of indels relative to a cleavage at the on-target of a highly " \
					"efficient sgRNA. </span></div><br>" \
					"\n" \
	                "" \
	                "</font></p>\n"
	html_content = re.sub("(?<=Your job status is.*?)Running", header_message, html_content)

	with open(html_path, "w") as fpw:
		fpw.write(html_content)
		fpw.flush()

	#send email
	if not isDailyTest:
		web_utils.send_last_email(RECIPIENT, RUN_NUMBER, logger=logger)
	else:
		web_utils.write_daily_test('crista_pair_score', RUN_NUMBER, 'PASS')



def check_input_vars(rna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome, strand, start_position, index_name=""):

	global EXTENSION, ERROR_MESSAGE_ON_EXCEPTION
	prefix = ""
	if index_name != "":
		prefix = "In line indexed as {}: ".format(index_name)

	if len(rna_seq) == 23 and rna_seq.upper().endswith("GG"):
		rna_seq = rna_seq[:-3]
		add_warning_message(prefix + "The RNA seq was truncated to 20nt (the NGG end was deleted)")
	rna_seq = re.sub("U", "T", rna_seq.upper())
	rna_seq = rna_seq.upper() + "NGG"


	end_position = include_genomic_features = w_flanking = None

	if extended100_genomic_seq is not None:
		extended100_genomic_seq = extended100_genomic_seq.upper()
		include_genomic_features = False
		assert isinstance(extended100_genomic_seq, str), "genomic argument is not sequence"
		len_extended = len(extended100_genomic_seq)
		w_flanking = len_extended==223
		if w_flanking:
			if not (len(re.search("[AGCT]+", extended100_genomic_seq).group()) == 223):
				add_warning_message(prefix + "In line indexed {}: DNA sequence must be ACGT only of length 223")
				return #error
		else:
			if not (len(re.search("[AGCT]+", extended100_genomic_seq).group()) == 29):
				add_warning_message(prefix + "In line indexed {}: DNA sequence must be ACGT only of length 29")
				return #error
			EXTENSION = 3


	else:
		include_genomic_features = True
		start_position = int(start_position)
		end_position = start_position + 22

		if chromosome.upper().startswith("CHR"):
			chromosome = chromosome[3:]
		if chromosome.isnumeric():
			chromosome = "chr" + chromosome
		else:
			chromosome = "chr" + chromosome.upper()

		if strand == "plusStrand":
			strand = "+"
		elif strand == "minusStrand":
			strand = "-"

		extended100_genomic_seq = \
			fetchSequenceFromGenome.fetch_dna_coordinates_offline(genome_database,
																  chromosome,
																  start_position - EXTENSION,
																  end_position + EXTENSION)
		if len(extended100_genomic_seq) < end_position-start_position+1+2*EXTENSION:
			return # error on this sample. if in file it remains a warning, else would be an error.

		extended100_genomic_seq = fetchSequenceFromGenome.get_seq_by_orientation(extended100_genomic_seq, strand)

		if genome_database == "hg19":
			cell_type = CELLS_DICT[cell_type]
		else:
			include_genomic_features = False  # it's only the sequence, not other genomic parameters
			chromosome = None #because we use the model without genomic features

	if not isinstance(rna_seq, str):
		add_warning_message(prefix + "In line indexed {}: RNA input is not a sequence")
		return #error
	rna_seq = rna_seq.upper()
	if not len(re.search("[AGCTU]+", rna_seq).group()) == 20:
		add_warning_message(prefix + "In line indexed {}: sgRNA sequence must be ACGTU only")
		return #error
	if 'U' in rna_seq:
		rna_seq = rna_seq.replace('U', 'T')
	dna_seq = extended100_genomic_seq[EXTENSION:-1*EXTENSION]

	#print(rna_seq)
	#print(dna_seq)

	return rna_seq, dna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome, \
		   strand, start_position, end_position, include_genomic_features, w_flanking


def parse_input_file(file):
	#input_df = pd.read_csv(file, skiprows=1, index_col="target_id")
	if not web_utils.is_csv(file):
		set_error_message("The file is not in csv format. Please save in excel as comma-delimited file. For instructions, see "
		                  "<a href=https://support.microsoft.com/en-us/office/save-a-workbook-to-text-format-txt-or-csv-3e9a9d6c-70da-4255-aa28-fcacf1f081e6>Microsoft Support website</a>.")

		raise ValueError(get_error_message())
	fpr = open(file)
	csvr = csv.reader(fpr)

	begin_flag = 0
	ids_list = []
	features_mat = None

	for row in csvr:
		if begin_flag == 0:
			if row[0]=="target_id":
				begin_flag = 1
				is_genomic_seqs = row[2]=="dna_seq_29nt"
				file_headers = row
		else:
			if is_genomic_seqs:
				res = check_input_vars(rna_seq=row[1], extended100_genomic_seq=row[2],
					                 genome_database=None, cell_type=None,
					                 chromosome=None, strand=None,
					                 start_position=None, index_name=row[0])
				if res is not None:
					rna_seq, dna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome, \
					strand, start_position, end_position, include_genomic_features, w_flanking =\
						res
				else:
					continue
			else:
				res = check_input_vars(rna_seq=row[1], extended100_genomic_seq=None,
					                 genome_database=row[2], cell_type=row[3],
					                 chromosome=row[4], strand=row[5],
					                 start_position=row[6], index_name=row[0])

				if res is not None:
					rna_seq, dna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome, \
					strand, start_position, end_position, include_genomic_features, w_flanking = res
				else:
					continue

			ids_list.append(row)
			features = get_features(rna_seq, dna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome,
							strand, start_position, end_position, include_genomic_features,
							w_flanking=w_flanking)

			if features_mat is None:
				features_mat = np.asmatrix(features)
			else:
				features_mat = np.concatenate((features_mat, np.array(features).reshape((1,len(features)))))

	return pd.DataFrame(ids_list, columns=file_headers), features_mat, include_genomic_features, w_flanking


if __name__ == '__main__':
	logger = logging.getLogger('CRISTA main script')
	init_commandline_logger(logger)
	parser = argparse.ArgumentParser(description='run feature collection for every sgrna')
	parser.add_argument('--file', '-f', default=None)
	parser.add_argument('--sgseq', '-s', default=None)
	parser.add_argument('--extended100_genomic_seq', '-d', default=None)
	parser.add_argument('--genome_database', '-g', default=None)
	parser.add_argument('--cell_type', '-e', default=None)
	parser.add_argument('--chromosome', '-c', default=None)
	parser.add_argument('--strand', '-w', default=None)
	parser.add_argument('--start_position', '-a', default=None)
	parser.add_argument('--path', '-p', default=None)
	parser.add_argument('--run_number', '-n', default=None)
	parser.add_argument('--daily_test', action='store_true')

	args = parser.parse_args()

	global RESULTS_PATH
	global RUN_NUMBER
	global RECIPIENT
	RESULTS_PATH = args.path
	RUN_NUMBER = args.run_number
	RECIPIENT = web_utils.get_email_recipient(RESULTS_PATH)

	try:
		if not args.daily_test:
			web_utils.send_first_email(RECIPIENT, RUN_NUMBER, logger=logger)
		input_file = args.file
		if input_file is not None:
			input_ids, features_mat, include_genomic_features, w_flanking = parse_input_file(input_file)
			score_df = predict_on_df(features_mat, include_genomic_features, w_flanking)
			input_ids["CRISTA prediction"] = score_df
			update_output_html(output_df=input_ids, logger=logger, isDailyTest=args.daily_test)

		else:
			res = check_input_vars(rna_seq=args.sgseq, extended100_genomic_seq=args.extended100_genomic_seq,
								 genome_database=args.genome_database, cell_type=args.cell_type,
								 chromosome=args.chromosome, strand=args.strand,
								 start_position=args.start_position)
			if res is None:
				set_error_message(get_warning_message())
				raise ValueError(get_error_message())
			rna_seq, dna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome, \
			strand, start_position, end_position, include_genomic_features, w_flanking = res

			logger.info([rna_seq, dna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome, strand, str(start_position), str(end_position), str(include_genomic_features)])
			features = get_features(rna_seq, dna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome, strand, start_position, end_position, include_genomic_features, w_flanking=w_flanking)
			features_mat = np.asmatrix(features)
			score_df = predict_on_df(features_mat, include_genomic_features, w_flanking, logger)
			update_output_html(sgrna=rna_seq, dna=dna_seq, score=score_df.iloc[0]["CRISTA score"], logger=logger, isDailyTest=args.daily_test)
	except:
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)
		if not args.daily_test: 
			web_utils.send_error_email(RECIPIENT, RUN_NUMBER, logger)
		else: 
			web_utils.write_daily_test('crista_pair_score', RUN_NUMBER, 'FAIL')
		web_utils.send_error_email("shiranos@gmail.com", RUN_NUMBER, logger)
		web_utils.load_error_page(RESULTS_PATH, get_error_message())

	web_utils.stop_refreshing(RESULTS_PATH)

	del_unused_genome_assemblies.delete_unused_genome_assemblies()


