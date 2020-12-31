import sys
sys.path.append("/groups/itay_mayrose/shiranabad/CRISPR/code/")
from collections import defaultdict
from definitions import *
from utils import *
from parse_u2os_gene_expression import build_u2os_ge_interval_map



def get_ucsc_table_file(genome_db, hgta_group, hgta_track, hgta_table, hgta_outFileName, hgta_regionType, hgta_outputType, position=None):
	"""
	if hgta_regionType == 'range' then position != None
	:return:
	"""
	return hgta_outFileName


def get_centromere_positions_dict(genome_db):
	"""
	:param genome_db: "hg19"/"mm9"/"hg38"/...
	:return: dictionary - for each chromosome with a tuple- endpos and startpos of the centromere
	"""
	if CENTROMERE_DICT_BY_GENOME[genome_db] is not None:
		return CENTROMERE_DICT_BY_GENOME[genome_db]
	#else:
	fp = open(get_ucsc_table_file(genome_db=genome_db, hgta_group="allTables", hgta_track=genome_db,
								  hgta_table="gap", hgta_outFileName=GAP_TABLE_FILENAME_BY_GENOME[genome_db],
								  hgta_regionType="genome", hgta_outputType="primaryTable"))
	d = {}
	for line in fp:
		if re.search("centromere", line):
			row = line.strip().split("\t")
			d[row[1]] = (int(row[2]), int(row[3]))    # d[chr#] = (startpos, endpos)
	fp.close()

	CENTROMERE_DICT_BY_GENOME[genome_db] = d
	return d


def get_distance_from_centromere(genome_db, chrom, offtarget_startpos, offtarget_endpos):
	d = get_centromere_positions_dict(genome_db)
	cent_pos = d[chrom]

	if cent_pos[0] <= offtarget_startpos <= cent_pos[1] or cent_pos[0] <= offtarget_endpos <= cent_pos[1]:
		return 0
	return min(abs(offtarget_startpos - cent_pos[1]), abs(cent_pos[0] - offtarget_endpos))


def get_chromosome_legnths_dict(genome_db):
	"""
	:param genome_db: "hg19"/"mm9"/"hg38"/...
	:return: dictionary - for each chromosome with a list- first position is the ending position of the first telomere,
							and second is the starting position of the last telomere
	"""
	if CHROM_LENGTHS_DICT_BY_GENOME[genome_db] is not None:
		return CHROM_LENGTHS_DICT_BY_GENOME[genome_db]
	#else:
	f = download_file(USCS_CHROM_LENGTHS_URL[genome_db], CHROM_LENGTHS_FILENAME_BY_GENOME[genome_db])
	d = {}
	fp = open(f)
	for line in fp:
		if re.match("chr[XYM0-9]+\t[0-9]+", line):
			row = line.strip().split("\t")
			d[row[0]] = int(row[1])
	CHROM_LENGTHS_DICT_BY_GENOME[genome_db] = d
	fp.close()
	return d


def get_telomere_positions_dict(genome_db):
	"""
	:param genome_db: "hg19"/"mm9"/"hg38"/...
	:return: dictionary - for each chromosome with a list- first position is the ending position of the first telomere,
							and second is the starting position of the last telomere
	"""
	if TELOMERE_DICT_BY_GENOME[genome_db] is not None:
		return TELOMERE_DICT_BY_GENOME[genome_db]
	#else:
	fp = open(get_ucsc_table_file(genome_db=genome_db, hgta_group="allTables", hgta_track=genome_db,
								  hgta_table="gap", hgta_outFileName=GAP_TABLE_FILENAME_BY_GENOME[genome_db],
								  hgta_regionType="genome", hgta_outputType="primaryTable"))

	d = defaultdict(lambda: [None, None])
	for line in fp:
		if re.search("telomere", line):
			row = line.strip().split("\t")
			if int(row[2]) == 0: #first telomere
				d[row[1]][0] = int(row[3])   # d[chr#] = [first telomere length, last telomere position]
			else:
				d[row[1]][1] = int(row[2])
	fp.close()
	chrom_lengths_dict = get_chromosome_legnths_dict(genome_db)
	for k in chrom_lengths_dict.keys():
		if k not in d.keys():
			if "chr1" in d.keys(): #not empty
				d[k] = [d["chr1"][0], chrom_lengths_dict[k] - d["chr1"][0]] #taking chr1's telomere as the default size for this genome
			else:
				d[k] = [10000, chrom_lengths_dict[k] - 10000] #taking 10000 as default size
	TELOMERE_DICT_BY_GENOME[genome_db] = d
	return d


def get_distance_from_telomere(genome_db, chrom, offtarget_startpos, offtarget_endpos):
	d = get_telomere_positions_dict(genome_db)
	tel_pos = d[chrom]

	if offtarget_endpos <= tel_pos[0] or offtarget_startpos >= tel_pos[1]:
		return 0
	return min(offtarget_startpos - tel_pos[0], tel_pos[1] - offtarget_endpos)


def map_CpG_islands(genome_db):
	"""
	:param genome_db:
	:return:
	"""
	if CpG_ISLANDS_MAPPING_DICT_BY_GENOME[genome_db] is not None:
		return CpG_ISLANDS_MAPPING_DICT_BY_GENOME[genome_db]

	if os.path.exists(CpG_ISLANDS_MAPPING_DICT_BY_GENOME_PICKLE[genome_db]):
		with open(CpG_ISLANDS_MAPPING_DICT_BY_GENOME_PICKLE[genome_db], "rb") as fpr:
			CpG_ISLANDS_MAPPING_DICT_BY_GENOME[genome_db] = pickle.load(fpr)
		return CpG_ISLANDS_MAPPING_DICT_BY_GENOME[genome_db]

	d = defaultdict(interval_defaultdict)
	with open(get_ucsc_table_file(genome_db=genome_db, hgta_group="allTables", hgta_track=genome_db,
								  hgta_table="cpgIslandExt", hgta_outFileName=CpG_ISLANDS_MAPPING_FILENAME_BY_GENOME[genome_db],
								  hgta_regionType="genome", hgta_outputType="primaryTable")) as fp:
		headerline = next(fp)
		for line in fp:
			row = line.strip().split("\t")
			if re.match("chr[XYM0-9]+$", row[1]):
				d[row[1]][int(row[2]):int(row[3])+1] = row[-1] #I took the field "obsExp"
				#http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=cpgIslandExt&hgta_table=cpgIslandExt&hgta_doSchema=describe+table+schema

	CpG_ISLANDS_MAPPING_DICT_BY_GENOME[genome_db] = d
	pickle.dump(CpG_ISLANDS_MAPPING_DICT_BY_GENOME[genome_db],
	            open(CpG_ISLANDS_MAPPING_DICT_BY_GENOME_PICKLE[genome_db], "wb"))
	return d


def is_in_CpG_island(genome_db, chrom, offtarget_startpos, offtarget_endpos):
	d = map_CpG_islands(genome_db)

	return (d[chrom][offtarget_startpos] or d[chrom][offtarget_endpos] or "0")


def map_methylation(genome_db, cell_line):
	"""
	:param genome_db:
	:return:
	"""

	if METHYL_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] is not None: # previously loaded
		return METHYL_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	if os.path.exists(METHYL_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE[genome_db][cell_line]):
		METHYL_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] = pickle.load(open(
			METHYL_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE[genome_db][cell_line], "rb"))
		return METHYL_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	tables_dicts = []
	for tablename in METHYL_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE[genome_db][cell_line]:
		d_plus = defaultdict(interval_defaultdict)
		d_minus = defaultdict(interval_defaultdict)
		with open(get_ucsc_table_file(genome_db=genome_db, hgta_group="Regulation", hgta_track=METHYL_ENCODE_TRACKNAME_BY_GENOME[genome_db],
								  hgta_table=tablename,
								  hgta_outFileName=get_ucsc_filename(tablename, cell_line, genome_db),
								  hgta_regionType="genome", hgta_outputType="primaryTable")) as fp:
			headerline = next(fp)
			for line in fp:
				row = line.strip().split("\t")
				if re.match("chr[XYM0-9]+$", row[1]):
					if row[-4] == "+":
						d_plus[row[1]][int(row[2]):int(row[3]) + 1] = row[5] #(strand, score)
					else:
						d_minus[row[1]][int(row[2]):int(row[3]) + 1] = row[5] #(strand, score)
			tables_dicts.append({"+": d_plus, "-": d_minus})

	METHYL_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] = tables_dicts
	pickle.dump(METHYL_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line],
	            open(METHYL_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE[genome_db][cell_line], "wb"))
	return METHYL_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]


def is_methylated(genome_db, chrom, strand, offtarget_startpos, offtarget_endpos, cell_line):
	if cell_line in [mESCs_CELLINE]:
		#print("no data for:" + cell_line)
		return None

	tables_dicts = map_methylation(genome_db, cell_line)
	val = 0
	for d in tables_dicts:
		val += float(d[strand][chrom][offtarget_startpos] or d[strand][chrom][offtarget_endpos] or '0')

	return str(val/float(len(tables_dicts)))


def map_Dnase1_hsSites(genome_db, cell_line):
	"""
	:param genome_db:
	:return:
	"""

	if DHS_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] is not None: # previously loaded
		return DHS_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	if os.path.exists(DHS_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE[genome_db][cell_line]):
		DHS_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] = pickle.load(
			open(DHS_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE[genome_db][cell_line], "rb"))
		return DHS_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	d = defaultdict(interval_defaultdict)
	tablename = DHS_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE[genome_db][cell_line]
	with open(get_ucsc_table_file(genome_db=genome_db, hgta_group="Regulation", hgta_track=DHS_ENCODE_TRACKNAME_BY_GENOME[genome_db],
								  hgta_table=tablename,
								  hgta_outFileName=get_ucsc_filename(tablename, cell_line, genome_db),
								  hgta_regionType="genome", hgta_outputType="primaryTable")) as fp:
		headerline = next(fp)
		for line in fp:
			row = line.strip().split("\t")
			if re.match("chr[XYM0-9]+$", row[1]):
				d[row[1]][int(row[2]):int(row[3]) + 1] = row[-4] #signalValue

	DHS_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] = d
	pickle.dump(DHS_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line],
	            open(DHS_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE[genome_db][cell_line], "wb"))
	return d


def get_DHS_signalVal(genome_db, chrom, offtarget_startpos, offtarget_endpos, cell_line):
	d = map_Dnase1_hsSites(genome_db, cell_line)

	in_DHS = d[chrom][offtarget_startpos] or d[chrom][offtarget_endpos] or "0"

	# distance from dhs
	dist = 0
	found = 1 if in_DHS != "0" else None
	while dist < 5000 and found is None:
		found = d[chrom][offtarget_startpos-dist] or d[chrom][offtarget_endpos+dist]
		dist += 50

	return in_DHS, (str(dist) if found is not None else "-1")


def map_nucleosomes(genome_db, cell_line):
	"""
	:param genome_db:
	:return:
	"""

	if NUCLEOSOME_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] is not None: # previously loaded
		return NUCLEOSOME_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	if os.path.exists(NUCLEOSOME_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE[genome_db][cell_line]):
		NUCLEOSOME_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] = pickle.load(
			open(NUCLEOSOME_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE[genome_db][cell_line], "rb"))
		return NUCLEOSOME_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	d = defaultdict(interval_defaultdict)
	if cell_line == K562_CELLINE:
		filename = get_ucsc_filename("wgEncodeSydhNsomeK562Sig", K562_CELLINE, genome_db)
		get_ucsc_table_file(genome_db=genome_db, hgta_group="Regulation", hgta_track="wgEncodeSydhNsome",
									  hgta_table="wgEncodeSydhNsomeK562Sig",
									  hgta_outFileName=filename,
									  hgta_regionType="genome", hgta_outputType="primaryTable")
	if cell_line == HEK293_CELLINE:
		filename = NUCLEOSOME_DATA_FILE[genome_db][cell_line]

	with open(filename) as fp:
		for line in fp:
			row = line.strip().split("\t")
			if re.match("chr[XYM0-9]+$", row[0]):
				d[row[0]][int(row[1]):int(row[2]) + 1] = float(row[3])

	NUCLEOSOME_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] = d
	pickle.dump(NUCLEOSOME_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line],
	            open(NUCLEOSOME_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE[genome_db][cell_line], "wb"))
	return d


def is_located_in_nucleosome(genome_db, chrom, offtarget_startpos, offtarget_endpos, cell_line):
	d = map_nucleosomes(genome_db, HEK293_CELLINE) #TODO: add queries to K562 Nsome data (bigwig file with bx-python)
	in_nucleosome = d[chrom][offtarget_startpos] or d[chrom][offtarget_endpos] or "0"

	# distance from nucleosome
	dist = 0
	found = 1 if in_nucleosome != "0" else None
	while dist < 1000 and found is None:
		found = d[chrom][offtarget_startpos-dist] or d[chrom][offtarget_endpos+dist]
		dist += 1

	return in_nucleosome, (str(dist) if found is not None else "-1")


def map_repeats(genome_db, w_repeatFamily=False):
	"""
	:param genome_db:
	:param w_repeatFamily: if false gets the class only, ow the type
	:return:
	"""

	mapping_dict = REPEATS_MAPPING_DICT[genome_db]
	mapping_dict_pkl = REPEATS_MAPPING_DICT_PICKLE[genome_db]

	if w_repeatFamily:
		mapping_dict = REPEATS_MAPPING_DICT_wFAMILY[genome_db]
		mapping_dict_pkl = REPEATS_MAPPING_DICT_PICKLE_wFAMILY[genome_db]

	if mapping_dict is not None: # previously loaded
		return mapping_dict

	if os.path.exists(mapping_dict_pkl):
		mapping_dict = pickle.load(
			open(mapping_dict_pkl, "rb"))
		return mapping_dict

	d = defaultdict(interval_defaultdict)

	with open(get_ucsc_table_file(genome_db=genome_db, hgta_group="Repeats", hgta_track="rmsk",
									  hgta_table="rmsk",
									  hgta_outFileName=get_ucsc_filename("rmsk", "", genome_db),
									  hgta_regionType="genome", hgta_outputType="primaryTable")) as fp:

		for line in fp:
			row = line.strip().split("\t")
			if re.match("chr[XYM0-9]+$", row[5]):
				repeat_class = row[11].replace("?","")
				repeat_family = row[12].replace("?","")
				if (not w_repeatFamily):
					val = repclass_dict[repeat_class]
				else:
					val = repfamily_dict[repeat_class + "/" + repeat_family]
				d[row[5]][int(row[6]):int(row[7]) + 1] = val

				#print(inv_repclass_dict[val])

	mapping_dict = d
	pickle.dump(mapping_dict,
	            open(mapping_dict_pkl, "wb"))
	return d


def is_repetitive_element(genome_db, chrom, offtarget_startpos, offtarget_endpos, w_repeatFamily=False):
	d = map_repeats(genome_db, w_repeatFamily)
	repeat_type = d[chrom][offtarget_startpos] or d[chrom][offtarget_endpos] or 0

	return str(repeat_type)



def map_Histone_h3k4me_Modifications(genome_db, cell_line):
	"""
	:param genome_db:
	:return:
	"""
	if cell_line in [U2OS_CELLINE]:
		HISTONE_H3k4me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line] = None
		return None

	if HISTONE_H3k4me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line] is not None: #previously loaded
		return HISTONE_H3k4me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	#elif os.path.exists(get_ucsc_picle_filename("HistoneModifications", cell_line, genome_db)): #load from pickle
	#	HISTONE_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] = pickle.load(open(get_ucsc_picle_filename("histone", cell_line, genome_db), "rb"))
	#	return HISTONE_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	#non existant - create it

	tables_dicts = []
	for tablename in HISTONE_H3k4me3_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE[genome_db][cell_line]:
		d = {}
		d = defaultdict(interval_defaultdict, d)
		with open(get_ucsc_table_file(genome_db=genome_db, hgta_group="Regulation", hgta_track=HISTONE_ENCODE_TRACKNAME_BY_CELLINE[genome_db][cell_line],
									  hgta_table=tablename,
									  hgta_outFileName=get_ucsc_filename(tablename, cell_line, genome_db),
									  hgta_regionType="genome", hgta_outputType="primaryTable")) as fp:
			headerline = next(fp)
			for line in fp:
				row = line.strip().split("\t")
				if re.match("chr[XYM0-9]+$", row[1]):
					d[row[1]][int(row[2]):int(row[3]) + 1] = row[-4] #signalValue
			tables_dicts.append(d)

	#pickle.dump(d, open(get_ucsc_pickle_filename("HistoneModifications", cell_line, genome_db), "wb"))
	HISTONE_H3k4me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line] = tables_dicts
	return tables_dicts


def map_Histone_h3k27me_Modifications(genome_db, cell_line):
	"""
	:param genome_db:
	:return:
	"""
	if cell_line in [HEK293_CELLINE, U2OS_CELLINE, mESCs_CELLINE]:
		HISTONE_H3k27me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line] = None
		return None

	if HISTONE_H3k27me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line] is not None: #previously loaded
		return HISTONE_H3k27me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	#elif os.path.exists(get_ucsc_picle_filename("HistoneModifications", cell_line, genome_db)): #load from pickle
	#	HISTONE_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] = pickle.load(open(get_ucsc_picle_filename("histone", cell_line, genome_db), "rb"))
	#	return HISTONE_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	#non existant - create it

	tables_dicts = []
	for tablename in HISTONE_H3k27me3_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE[genome_db][cell_line]:
		d = {}
		d = defaultdict(interval_defaultdict, d)
		with open(get_ucsc_table_file(genome_db=genome_db, hgta_group="Regulation", hgta_track=HISTONE_ENCODE_TRACKNAME_BY_CELLINE[genome_db][cell_line],
									  hgta_table=tablename,
									  hgta_outFileName=get_ucsc_filename(tablename, cell_line, genome_db),
									  hgta_regionType="genome", hgta_outputType="primaryTable")) as fp:
			headerline = next(fp)
			for line in fp:
				row = line.strip().split("\t")
				if re.match("chr[XYM0-9]+$", row[1]):
					d[row[1]][int(row[2]):int(row[3]) + 1] = row[-4] #signalValue
			tables_dicts.append(d)

	#pickle.dump(d, open(get_ucsc_pickle_filename("HistoneModifications", cell_line, genome_db), "wb"))
	HISTONE_H3k27me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line] = tables_dicts
	return tables_dicts


def map_Histone_h3k36me_Modifications(genome_db, cell_line):
	"""
	:param genome_db:
	:return:
	"""
	if cell_line in [HEK293_CELLINE, mESCs_CELLINE]:
		HISTONE_H3k36me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line] = None
		return None

	if HISTONE_H3k36me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line] is not None: #previously loaded
		return HISTONE_H3k36me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	#elif os.path.exists(get_ucsc_picle_filename("HistoneModifications", cell_line, genome_db)): #load from pickle
	#	HISTONE_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] = pickle.load(open(get_ucsc_picle_filename("histone", cell_line, genome_db), "rb"))
	#	return HISTONE_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	#non existant - create it

	tables_dicts = []
	for tablename in HISTONE_H3k36me3_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE[genome_db][cell_line]:
		d = {}
		d = defaultdict(interval_defaultdict, d)
		with open(get_ucsc_table_file(genome_db=genome_db, hgta_group="Regulation", hgta_track=HISTONE_ENCODE_TRACKNAME_BY_CELLINE[genome_db][cell_line],
									  hgta_table=tablename,
									  hgta_outFileName=get_ucsc_filename(tablename, cell_line, genome_db),
									  hgta_regionType="genome", hgta_outputType="primaryTable")) as fp:
			headerline = next(fp)
			for line in fp:
				row = line.strip().split("\t")
				if re.match("chr[XYM0-9]+$", row[1]):
					d[row[1]][int(row[2]):int(row[3]) + 1] = row[-4] #signalValue
			tables_dicts.append(d)

	#pickle.dump(d, open(get_ucsc_pickle_filename("HistoneModifications", cell_line, genome_db), "wb"))
	HISTONE_H3k36me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line] = tables_dicts
	return tables_dicts


def map_Histone_h3k9me_Modifications(genome_db, cell_line):
	"""
	:param genome_db:
	:return:
	"""
	if cell_line in [HEK293_CELLINE, hIPS_CELLINE, mESCs_CELLINE]:
		HISTONE_H3k9me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line] = None
		return None

	if HISTONE_H3k9me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line] is not None: #previously loaded
		return HISTONE_H3k9me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line]


	tables_dicts = []
	for tablename in HISTONE_H3k9me3_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE[genome_db][cell_line]:
		d = {}
		d = defaultdict(interval_defaultdict, d)
		with open(get_ucsc_table_file(genome_db=genome_db, hgta_group="Regulation",
		                              hgta_track=HISTONE_ENCODE_TRACKNAME_BY_CELLINE[genome_db][cell_line],
									  hgta_table=tablename,
									  hgta_outFileName=get_ucsc_filename(tablename, cell_line, genome_db),
									  hgta_regionType="genome", hgta_outputType="primaryTable")) as fp:
			headerline = next(fp)
			for line in fp:
				row = line.strip().split("\t")
				if re.match("chr[XYM0-9]+$", row[1]):
					sv = row[-4]
					if sv == ".":
						sv = "0"
					d[row[1]][int(row[2]):int(row[3]) + 1] = sv #signalValue
			tables_dicts.append(d)

	#pickle.dump(d, open(get_ucsc_pickle_filename("HistoneModifications", cell_line, genome_db), "wb"))
	HISTONE_H3k9me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line] = tables_dicts
	return HISTONE_H3k9me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE[genome_db][cell_line]


def get_HistoneModifications_signalVal(genome_db, chrom, offtarget_startpos, offtarget_endpos, cell_line, histone_type):
	"""
	:param genome_db:
	:param chrom:
	:param offtarget_startpos:
	:param offtarget_endpos:
	:param cell_line:
	:param histone_type:
	:return:
	"""

	if ALL_HISTONES_DICT[cell_line] is None:
		ALL_HISTONES_DICT[cell_line] = {"H3k4me": map_Histone_h3k4me_Modifications(genome_db, cell_line),
							 "H3k9me": map_Histone_h3k9me_Modifications(genome_db, cell_line),
							 "H3k27me": map_Histone_h3k27me_Modifications(genome_db, cell_line),
							 "H3k36me": map_Histone_h3k36me_Modifications(genome_db, cell_line)}
	tables_dicts = ALL_HISTONES_DICT[cell_line][histone_type]
	if tables_dicts == None:
		return None
	val = 0
	for d in tables_dicts:
		val += float(d[chrom][offtarget_startpos] or d[chrom][offtarget_endpos] or '0')

	return str(val/float(len(tables_dicts)))


def map_exon_sites(genome_db, cell_line):
	"""
	:param genome_db:
	:return:
	"""
	if cell_line in [mESCs_CELLINE]: #U2OS- took Hela instead...
		#print("no data for:" + cell_line)
		return

	if EXON_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] is not None:
		return EXON_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	if os.path.exists(EXON_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE[genome_db][cell_line]):
		EXON_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] = pickle.load(
			open(EXON_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE[genome_db][cell_line], "rb"))
		return EXON_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]

	if cell_line == U2OS_CELLINE:
		tables_dicts = [build_u2os_ge_interval_map()]

	else:
		tables_dicts = []
		for tablename in EXON_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE[genome_db][cell_line]:
			d_plus = defaultdict(interval_defaultdict)
			d_minus = defaultdict(interval_defaultdict)
			with open(get_ucsc_table_file(genome_db=genome_db, hgta_group="Expression",
			                              hgta_track=EXON_ENCODE_TRACKNAME_BY_GENOME[genome_db],
										  hgta_table=tablename,
										  hgta_outFileName=get_ucsc_filename(tablename, cell_line, genome_db),
										  hgta_regionType="genome", hgta_outputType="primaryTable")) as fp:
				headerline = next(fp)
				for line in fp:
					row = line.strip().split("\t")

					if re.match("chr[XYM0-9]+$", row[0]):
						if row[5] == "+":
							d_plus[row[0]][int(row[1]):int(row[2]) + 1] = row[-3] #(strand, signalValue)
						else:
							d_minus[row[0]][int(row[1]):int(row[2]) + 1] = row[-3] #(strand, signalValue)
				tables_dicts.append({"+": d_plus, "-": d_minus})

	EXON_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line] = tables_dicts

	pickle.dump(EXON_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line],
	            open(EXON_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE[genome_db][cell_line], "wb"))

	return EXON_MAPPING_DICT_BY_GENOME_BY_CELLINE[genome_db][cell_line]


def is_mapped_to_exon(chrom, strand, offtarget_startpos, offtarget_endpos, tables_dicts):

	exon_val = 0
	opposite_strand = "-" if strand=="+" else "+"
	opposite_exon_val = 0

	for d in tables_dicts:
		start_val = d[strand][chrom][offtarget_startpos-1]
		end_val = d[strand][chrom][offtarget_endpos+1]
		vals = []
		if start_val is not None:
			vals.append(float(start_val))
		if end_val is not None:
			vals.append(float(end_val))
		if len(vals)>0:
			exon_val += max(vals)

		start_val = d[opposite_strand][chrom][offtarget_startpos-1]
		end_val = d[opposite_strand][chrom][offtarget_endpos+1]
		vals = []
		if start_val is not None:
			vals.append(float(start_val))
		if end_val is not None:
			vals.append(float(end_val))
		if len(vals)>0:
			opposite_exon_val += max(vals)

	return str(exon_val / float(len(tables_dicts))), str(opposite_exon_val / float(len(tables_dicts)))


def get_exon_signalVal(genome_db, chrom, strand, offtarget_startpos, offtarget_endpos, cell_line=None):
	if cell_line in [mESCs_CELLINE]:
		#print("no data for:" + cell_line)
		return None, None
	if cell_line is not None:
		return is_mapped_to_exon(chrom, strand, offtarget_startpos, offtarget_endpos, map_exon_sites(genome_db, cell_line))
	else:
		return is_mapped_to_exon(chrom, strand, offtarget_startpos, offtarget_endpos, [get_exon_mapping(genome_db)])


def map_genes(genome_db):
	"""
	:param genome_db:
	:return:
	"""

	if GENES_MAPPING_DICT_BY_GENOME[genome_db] is not None:
		return GENES_MAPPING_DICT_BY_GENOME[genome_db]

	if os.path.exists(GENES_MAPPING_DICT_BY_GENOME_PICKLE[genome_db]) and os.path.exists(EXONS_MAPPING_DICT_BY_GENOME_PICKLE[genome_db]):
		GENES_MAPPING_DICT_BY_GENOME[genome_db] = pickle.load(open(GENES_MAPPING_DICT_BY_GENOME_PICKLE[genome_db], "rb"))
		EXONS_MAPPING_DICT_BY_GENOME[genome_db] = pickle.load(open(EXONS_MAPPING_DICT_BY_GENOME_PICKLE[genome_db], "rb"))
		return GENES_MAPPING_DICT_BY_GENOME[genome_db]

	tablename = GENES_TABLENAME_BY_GENOME[genome_db]

	d_tx = {"+": defaultdict(interval_defaultdict), "-": defaultdict(interval_defaultdict)}
	d_cds = {"+": defaultdict(interval_defaultdict), "-": defaultdict(interval_defaultdict)}
	d_exons = {"+": defaultdict(interval_defaultdict), "-": defaultdict(interval_defaultdict)}
	with open(get_ucsc_table_file(genome_db=genome_db, hgta_group="genes", hgta_track="refGene",
								hgta_table=tablename, hgta_outFileName=get_ucsc_filename(tablename, "all", genome_db),
								hgta_regionType="genome", hgta_outputType="primaryTable")) as fp:
		headerline = next(fp)
		for line in fp:
			row = line.strip().split("\t")
			chromosome = row[2]
			tx_start = int(row[4])
			tx_end = int(row[5])
			cds_start = int(row[6])
			cds_end = int(row[7])
			strand = row[3]
			exons_start = [int(x) for x in row[9].strip(",").split(",")]
			exons_end = [int(x) for x in row[10].strip(",").split(",")]
			exons_cnt = int(row[8])

			if re.match("chr[XYM0-9]+$", chromosome):
				d_tx[strand][chromosome][tx_start:tx_end + 1] = "1"
				d_cds[strand][chromosome][cds_start:cds_end + 1] = "1"

				for j in range(exons_cnt):
					d_exons[strand][chromosome][exons_start[j]: exons_end[j]+1] = '1' #exon


	GENES_MAPPING_DICT_BY_GENOME[genome_db] = {"tx": d_tx, "cds": d_cds}
	EXONS_MAPPING_DICT_BY_GENOME[genome_db] = d_exons

	pickle.dump(GENES_MAPPING_DICT_BY_GENOME[genome_db],
	            open(GENES_MAPPING_DICT_BY_GENOME_PICKLE[genome_db], "wb"))
	pickle.dump(EXONS_MAPPING_DICT_BY_GENOME[genome_db],
	            open(EXONS_MAPPING_DICT_BY_GENOME_PICKLE[genome_db], "wb"))

	return GENES_MAPPING_DICT_BY_GENOME[genome_db]


def get_gene_properties(genome_db, chrom, offtarget_startpos, offtarget_endpos):

	tables_dict = map_genes(genome_db)
	in_tx = tables_dict["tx"]["+"][chrom][offtarget_startpos] or tables_dict["tx"]["+"][chrom][offtarget_endpos] or \
	        tables_dict["tx"]["-"][chrom][offtarget_startpos] or tables_dict["tx"]["-"][chrom][offtarget_endpos] or "0"

	in_cds = tables_dict["cds"]["+"][chrom][offtarget_startpos+6] or tables_dict["cds"]["+"][chrom][offtarget_endpos-6] or \
	         tables_dict["cds"]["-"][chrom][offtarget_startpos+6] or tables_dict["cds"]["-"][chrom][offtarget_endpos-6] or "0"
	#according to Doench, 6nt further are possible

	return in_tx, in_cds


def get_exon_mapping(genome_db):
	map_genes(genome_db)
	return EXONS_MAPPING_DICT_BY_GENOME[genome_db]


def map_HK_genes(genome_db):
	"""
	:param genome_db:
	:return:
	"""
	if HK_GENES_MAPPING_DICT_BY_GENOME[genome_db] is not None:
		return HK_GENES_MAPPING_DICT_BY_GENOME[genome_db]
	#else:
	d = {}
	d = defaultdict(interval_defaultdict)

	if not os.path.exists(HK_GENES):
		tablename = GENES_TABLENAME_BY_GENOME[genome_db]

		d_tx = {"+": defaultdict(interval_defaultdict), "-": defaultdict(interval_defaultdict)}
		with open(get_ucsc_table_file(genome_db=genome_db, hgta_group="genes", hgta_track=tablename,
		                              hgta_table=tablename, hgta_outFileName=get_ucsc_filename(tablename, "all", genome_db),
		                              hgta_regionType='genome', hgta_outputType="primaryTable")) as fp:
			headerline = next(fp)
			for line in fp:
				row = line.strip().split("\t")
				chromosome = row[2]
				tx_start = int(row[4])
				tx_end = int(row[5])
				strand = row[3]

				if re.match("chr[XYM0-9]+$", chromosome) and (d_tx[strand][chromosome][tx_start] is None or d_tx[strand][chromosome][tx_end] is None):
					d_tx[strand][chromosome][tx_start:tx_end + 1] = [row[-4], row[1], chromosome, strand, tx_start, tx_end]


		arrange_data.format_hk_genes_file(HK_GENES_ORIG, PARSED_HK_GENES)
		with open(PARSED_HK_GENES) as fp:
			headerline = next(fp)
			csvw = open(HK_GENES, "w", newline='')
			csv_writer = csv.writer(csvw, dialect='excel')
			csv_writer.writerow(["gene name", "refseq code", "chromosome", "strand", "txstart", "txend"])

			for line in fp:
				row = line.strip().split(",")

				if (d_tx[row[3]][row[2]][int(row[4])]):
					csv_writer.writerow(d_tx[row[3]][row[2]][int(row[4])])
				else:
					csv_writer.writerow(d_tx["-" if row[3]=="+" else "+"][row[2]][int(row[4])])

			csvw.close()

	with open(HK_GENES) as fp:
		headerline = next(fp)
		for line in fp:

			row = line.strip().split(",")
			d[row[2]][int(row[4]):int(row[5])+1] = 1 # added 2000 bp upstream to account for the promoter

	HK_GENES_MAPPING_DICT_BY_GENOME[genome_db] = d
	return d


def is_in_HK_gene(genome_db, chrom, offtarget_startpos, offtarget_endpos):
	d = map_HK_genes(genome_db)

	return int((d[chrom][offtarget_startpos] or d[chrom][offtarget_endpos]) != None)


def map_nucleosome_occupancy_prediction(genome_db, chromosome):
	"""
	:param genome_db:
	:return:
	"""

	if NUC_OCC_PREDICTIONS_MAPPING_DICT_BY_GENOME[genome_db] is not None and \
					chromosome in NUC_OCC_PREDICTIONS_MAPPING_DICT_BY_GENOME[genome_db].keys():
		return NUC_OCC_PREDICTIONS_MAPPING_DICT_BY_GENOME[genome_db][chromosome]

	if NUC_OCC_PREDICTIONS_MAPPING_DICT_BY_GENOME[genome_db] is None:
		NUC_OCC_PREDICTIONS_MAPPING_DICT_BY_GENOME[genome_db] = defaultdict(interval_defaultdict)

	NUC_OCC_PREDICTIONS_MAPPING_DICT_BY_GENOME[genome_db][chromosome] = \
		pickle.load(open(NUC_OCC_PREDICTIONS_DIR[genome_db] + chromosome + ".fa.tab.pkl", "rb"))

	return NUC_OCC_PREDICTIONS_MAPPING_DICT_BY_GENOME[genome_db][chromosome]


def get_nucleosome_occupancy_prediction(genome_db, chrom, offtarget_startpos, offtarget_endpos):
	d = map_nucleosome_occupancy_prediction(genome_db, chrom)
	if d[offtarget_endpos-3] is not None:
		return d[offtarget_endpos-3]
	else:
		return 0


def map_subcompartments(genome_db):
	"""
	:param genome_db:
	:return:
	"""

	if THREE_D_SUBCOMPARTMENTS_DICT[genome_db] is not None: # previously loaded
		return THREE_D_SUBCOMPARTMENTS_DICT[genome_db]

	if os.path.exists(THREE_D_SUBCOMPARTMENTS_DICT_PICKLE[genome_db]):
		THREE_D_SUBCOMPARTMENTS_DICT[genome_db] = pickle.load(
			open(THREE_D_SUBCOMPARTMENTS_DICT_PICKLE[genome_db], "rb"))
		return THREE_D_SUBCOMPARTMENTS_DICT[genome_db]

	d = defaultdict(interval_defaultdict)
	filename = THREE_D_SUBCOMPARTMENTS[genome_db]
	with open(filename) as fp:
		headerline = next(fp)
		for line in fp:
			row = line.strip().split("\t")
			if re.match("chr[XYM0-9]+$", row[0]):
				d[row[0]][int(row[1]):int(row[2]) + 1] = SUBCOMPARTMENTS_2_CATEGORIES_MAP[row[3]] #subcompartment

	THREE_D_SUBCOMPARTMENTS_DICT[genome_db] = d
	pickle.dump(THREE_D_SUBCOMPARTMENTS_DICT[genome_db],
	            open(THREE_D_SUBCOMPARTMENTS_DICT_PICKLE[genome_db], "wb"))
	return d


def get_subcompartments(genome_db, chrom, offtarget_startpos, offtarget_endpos):
	d = map_subcompartments(genome_db)

	subcompartment = d[chrom][offtarget_startpos] or d[chrom][offtarget_endpos] or "0"

	return str(subcompartment)


if __name__ == '__main__':
	map_repeats("hg38")
	exit()
	print(get_subcompartments("hg19", 'chr13', 100741280, 100741300))
	print (is_methylated("hg19", 'chr13', "+", 100741280, 100741300, K562_CELLINE))

	print (get_exon_signalVal("hg19", 'chr13', "+", 100741280, 100741300, K562_CELLINE))
	print (get_exon_signalVal("hg19", 'chr13', "+", 100741280, 100741300, U2OS_CELLINE))
	print (get_exon_signalVal("hg19", 'chr13', "+", 100741280, 100741300, HEK293_CELLINE))
