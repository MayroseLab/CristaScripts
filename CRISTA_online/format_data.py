import traceback

import sys

from definitions import *
import pandas as pd
import numpy as np
import regex as re


ACGT_REPLACEMENT = {"A": '1', "C": '2', "G": '3', "T": '4', 'N': '0'}
REV_ACGT_REPLACEMENT = {v: k for k, v in ACGT_REPLACEMENT.items()}
ACGT_COUPLES_REP = {}

for i in range(1,5):
	for j in range(1,5):
		ACGT_COUPLES_REP[str(i*10+j)] = REV_ACGT_REPLACEMENT[str(i)] + REV_ACGT_REPLACEMENT[str(j)]

MMS_TYPE_REPLACEMENT = {'0': "match", '-1': "indel", '1': "wobble", '2': "RR transition", '3': "YY transition", '4': "transversion"}


def agct2numerals(st):
	try:
		new_st = ""
		for x in st:
			#print(x)
			try:
				new_st += ACGT_REPLACEMENT[x]
			except KeyError:
				new_st += "0"
		return new_st
	except:
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)
		set_error_message("Invalid positions in genomic sequence. The genomic site you selected might contain unannotated positions.")
		raise


def get_chromosome_number(chrom):
	n = re.search("(?<=chr)[XYM0-9]+", chrom).group()
	if n == "X":
		return '23'
	if n == "Y":
		return '24'
	if n == "M":
		return '25'
	return n


def convert_column_data(df, column_names, conversion_dict):
	"""inplace"""
	for col in column_names:
		if col in df.columns:
			df[col] = df[col].map(conversion_dict)


def convert_nucleotides_cols_from_numbers(df):
	"""inplace"""

	columns_of_nucleotides = ["PAM_2_first", "PAM_N_id", "last_nucleotide", "PamProximalCouple1", "PamProximalCouple2",
							  "PamProximalCouple3", "PamProximalCouple4", "id_pos5", "id_pos4", "id_pos3", "id_pos2", "id_pos1",
							  "nuc_id_afterPAM1", "nuc_id_afterPAM2", "nuc_id_afterPAM3", "nuc_id_afterPAM4", "nuc_id_afterPAM5"]
	one_dict = dict(ACGT_COUPLES_REP)
	one_dict.update(REV_ACGT_REPLACEMENT)
	convert_column_data(df, columns_of_nucleotides, one_dict)


def convert_mms_cols_from_numbers(df):
	"""inplace"""
	columns_of_mms = ["MM_in_site_" + str(i) for i in range(20, 0, -1)]
	convert_column_data(df, columns_of_mms, MMS_TYPE_REPLACEMENT)


def remove_duplicates_from_df(df):
	sort_dataframe_by_score(df)
	df.drop_duplicates(subset=[CHROMOSOME_HEADING, STRAND_HEADING, STARTPOS_HEADING], inplace=True)
	df.drop_duplicates(subset=[CHROMOSOME_HEADING, STRAND_HEADING, ENDPOS_HEADING], inplace=True)


def sort_dataframe_by_score(df):
	df.sort_values(by="CRISTA score", axis=0, ascending=False, inplace=True)
	df.reset_index(inplace=True, drop=True)
	df.index += 1 #start with 1


def convert_pos_cols_to_int(df):
	# make positions cols int
	df["start position"] = df["start position"].astype(int)
	df["end position"] = df["end position"].astype(int)
	if 'index' in df.columns:
		df.drop('index', axis=1, inplace=True)


def format_csv_allfeatures_table(df, logger):
	"""inplace"""

	#remove index col
	sort_dataframe_by_score(df)
	# sort according to score
	logger.info("step 1")

	# change nucleotides columns from numbers to nucleotides
	convert_nucleotides_cols_from_numbers(df)
	logger.info("step 2")

	# remove second chromosome column #todo

	# change mms to Tv, ts, wobble, YY, RR, match, gap
	convert_mms_cols_from_numbers(df)
	convert_pos_cols_to_int(df)
	logger.info("step 3")

	# change heading
	features_map = get_features_names_mapping()
	logger.info("step 4")
	heading = pd.DataFrame(df.columns.values, columns=["feature name"]) #wo index col
	logger.info("step 5")
	heading_map = pd.merge(heading, features_map, how="left", right_index=True, left_on=["feature name"])
	logger.info("step 6")
	heading_map["display name"].fillna(heading_map["feature name"], inplace=True)
	logger.info("step 7")
	df.columns = heading_map["display name"]
	logger.info("step 8")

