import sys
from collections import defaultdict

sys.path.append("/groups/itay_mayrose/shiranabad/CRISPR/code/")
sys.path.append("D:\\Dropbox\\lab\\CRISPR\\crisprProject")
from Bio import Entrez, SeqIO
from utils import *
from definitions import *

DELIMITER = "\t"
EVENT_TYPE_COL = 12
EVENT_LOCATION_COL = 20
REF_ACCESSION_COL = 23
LOCATION_IN_REF_COL = 21
ENTREZ_GENE_ID_COL = 11
U2OS_DIR = SEP.join([DATA_DIR, "u2os_geneExpression", ""])
BLAT_DATA = SEP.join([U2OS_DIR, "blat", ""])
GE_FILEPATH = SEP.join([U2OS_DIR, "GSE22149_Gene_level.csv"])


def parse_gene_expression_file():
	"""
	parses "GSE22149_Gene_level.csv" so that we have the mean of control GE measurements (NS) for every entrez id
	:return: dict: entrez_gene_id -> avg ge
	"""
	ge_dict = {}
	with open(GE_FILEPATH) as ge_fpr:
		next(ge_fpr) # skip header
		for line in ge_fpr:
			ge_row = line.strip().split(",")
			ge_dict[ge_row[0]] = get_avg([float(ge_row[4]), float(ge_row[5]), float(ge_row[6])])
	return ge_dict


def get_gene_from_ncbi(gene_accession):
	Entrez.email = "A.N.Other@example.com"
	handle = Entrez.efetch(db="nucleotide", id=gene_accession, rettype="gb", retmode="text")
	record = SeqIO.read(handle, "genbank")
	handle.close()

	return record


def get_cds_length_from_record(record):
	for item in record.features:
		if item.type == "CDS":
			return item.location


def get_exon_locations_from_record(record):
	exons_locations = []
	for item in record.features:
		if item.type == "exon":
			exons_locations.append(item.location)
	return exons_locations


def is_in_cds(location_in_ref, event_location):
	"""
	:param location_in_ref: a string: CDS / 3 UTR / 5 UTR
	:return: True if in cds, else false
	"""
	if location_in_ref == 'CDS' and not(event_location.startswith("<<<") or event_location.endswith(">>>")):
		return True
	else:
		return False


def include_intron(introns_events, event_location):
	"""
	:param introns_events: a list of introns
	:param event_location: list of indices [x, x+1]
	:return: updates intron_events s.t the intron between exon x-1 and x is included
	"""
	introns_events[event_location[0]-1] = True


def parse_location(location_event_str):
	"""
	:param location_event_str: strings of types: x, x<>y, x><y
	 excluded cases: <<1, n<< - because I ignore UTRs
	:return:
	"""
	if re.match("^[0-9]+$", location_event_str):
		location_event = [int(location_event_str)]
	else:
		first = int(re.search("^[0-9]+", location_event_str).group())
		sec = int(re.search("[0-9]+$", location_event_str).group())
		first, sec = min([first, sec]), max([first, sec])
		location_event = [i for i in range(first, sec+1)]
	return location_event


def remove_exons(exons_events, event_location):
	for i in event_location:
		exons_events[i-1] = False


def add_exon(introns_events, event_location):
	# since I don't know the coordinates, I'll include the whole intron
	include_intron(introns_events, event_location)


def update_exon_start(exons_locations_in_genome):
	pass # I have no way to deal with it


def update_exon_end(exons_locations_in_genome):
	pass # I have no way to deal with it


def get_exon_locations_from_blat(blockSizes, tStarts):
	exons_locations = []
	start_poses = tStarts.strip().strip(",").split(",")
	exons_lengths = blockSizes.strip().strip(",").split(",")
	for i in range(len(start_poses)):
		exons_locations.append((int(start_poses[i]) + 1, int(start_poses[i]) + int(exons_lengths[i])))
	return exons_locations


def parse_source_file(in_filename):

	logger = logging.getLogger('Map U2OS GE')
	init_commandline_logger(logger)

	#PROTEINS_MAPPING = {"XM_936315.1": "NM_001164462.1", "NM_005795.2": "NM_005795.5", "NM_170605.2": "NM_176877.2", "NM_006091.1": "NM_006091.4"}
	#fpw_unreferenced_cdnas = open("/groups/itay_mayrose/shiranabad/CRISPR/data/u2os_geneExpression/unreferenced_cdnas.fa", "a")
	with open(in_filename) as fpr:
		last_loc = ""
		last_event = ""
		last_gene_ref = ""
		last_entrez_gid = ""
		skip_gene = False

		fpw = open(in_filename + "_parsed_included_regions.csv", "w", newline='')
		csv_writer = csv.writer(fpw, dialect='excel')
		csv_writer.writerow(["accesion", "type", "chromosome", "start_location_in_cds", "end_location_in_cds", "strand", "entrez_gid", "ge_val"])

		for line in fpr:

			row = line.strip().split(DELIMITER)
			row = [x.strip('"') for x in row]

			if not is_in_cds(row[LOCATION_IN_REF_COL], row[EVENT_LOCATION_COL]):
				continue

			event_type = row[EVENT_TYPE_COL]
			event_location = parse_location(row[EVENT_LOCATION_COL])
			gene_accession = row[REF_ACCESSION_COL]
			entrez_gid = row[ENTREZ_GENE_ID_COL]

			#if gene_accession in PROTEINS_MAPPING.keys():
			#	gene_accession = PROTEINS_MAPPING[gene_accession]

			if last_loc == event_location and event_type == last_event and last_gene_ref == gene_accession:
				# if it's the same as the previous event
				continue

			if last_gene_ref != gene_accession:

				# save previous arrays to fpw

				if last_gene_ref != "" and not skip_gene:
					for i in range(len(exons_locations_in_genome)):
						start, end = exons_locations_in_genome[i][0], exons_locations_in_genome[i][1]
						if exons_events[i]:
							csv_writer.writerow([last_gene_ref, "exon", chromosome, start, end, strand, last_entrez_gid, ge_val])
						if i!=len(exons_locations_in_genome)-1 and introns_events[i]:
							csv_writer.writerow([last_gene_ref, "intron", chromosome, end, exons_locations_in_genome[i+1][1], strand, last_entrez_gid, ge_val])

					fpw.flush()
				skip_gene = False


				# retrieve the accesion CDS with the exons from ncbi
				"""
				cds_record = get_gene_from_ncbi(gene_accession)
				exons_locations_in_genome = get_exon_locations_from_record(cds_record)
				corrected_gene_accession = gene_accession
				try:
					while len(exons_locations_in_genome) == 0:
						corrected_gene_accession = corrected_gene_accession[:-1] + str(int(corrected_gene_accession[-1])+1)
						cds_record = get_gene_from_ncbi(corrected_gene_accession)
						#gene_accession = corrected_gene_accession
						exons_locations_in_genome = get_exon_locations_from_record(cds_record)
				except:
					logger.info(corrected_gene_accession + " is non-existant")
					skip_gene = True
					exons_locations_in_genome = [get_cds_length_from_record(cds_record)]
					last_loc = row[EVENT_LOCATION_COL]
					last_event = event_type
					last_gene_ref = gene_accession
				"""

				current_dir = BLAT_DATA + gene_accession + SEP
				outputfile = current_dir + "out.psl"


				if not os.path.exists(outputfile):
					try:
						cds_record = get_gene_from_ncbi(gene_accession)
						cds_seq = cds_record.seq._data
						chromosome = "chr" + cds_record.features[0].qualifiers['chromosome'][0]
						chromosome_ref = LOCAL_DB_DIRS["hg19"] + chromosome + ".fa"

						if not os.path.exists(current_dir):
							os.mkdir(current_dir)
						cds_file = current_dir + "seq.txt"

						with open(cds_file, "w") as fpw:
							fpw.write(">" + gene_accession + "\n" + cds_seq)

						os.system(" ".join([BLAT_EXE, chromosome_ref, cds_file, outputfile]))
						"""
						job_filename = \
							createJobFile.create_job_file(gene_accession, " ".join([BLAT_EXE, chromosome_ref, cds_file, outputfile]),
					                                      gene_accession + ".sh", "/groups/itay_mayrose/shiranabad/CRISPR/data/u2os_geneExpression/error_files/",
					                                        "/groups/itay_mayrose/shiranabad/CRISPR/data/u2os_geneExpression/job_files/")

						os.system('qsub -p -1 ' + job_filename)
						"""
					except:
						#fpw_unreferenced_cdnas.write(">" + gene_accession + "\n" + cds_seq + "\n")
						print(gene_accession, "- can't load ncbi record")
						last_loc = row[EVENT_LOCATION_COL]
						last_event = event_type
						last_gene_ref = gene_accession
						last_entrez_gid = entrez_gid
						skip_gene = True
						continue


				with open(outputfile) as blat_fp:
					min_row = None
					max_matches = -1
					blat_lines = blat_fp.readlines()
					if len(blat_lines) == 5:
						skip_gene = True
						last_loc = row[EVENT_LOCATION_COL]
						last_event = event_type
						last_gene_ref = gene_accession
						last_entrez_gid = entrez_gid
						continue

					for blat_line in blat_lines[5:]:
						blat_row = blat_line.split("\t")
						if len(blat_row) > 20 and int(blat_row[0]) > max_matches:
							max_matches = int(blat_row[0])
							min_row = blat_row
					num_exons = int(min_row[17])
					exons_locations_in_genome = get_exon_locations_from_blat(min_row[18], min_row[20])
					chromosome = min_row[13]
					strand = min_row[8]

				exons_events = [True]*num_exons #True means that the exon is included
				introns_events = [False]*(num_exons-1) #"included" in position 0 means that the intron between exon 0 and exon 1 is included (retention)
				ge_val = get_avg([float(row[4]), float(row[5]), float(row[6])])

			if skip_gene:
				continue

			if max(event_location) > len(exons_events):
				logger.info("in " + gene_accession + " there is a larger exon than in ncbi")
				continue

			logger.info("taking care of event " + event_type + " of " + gene_accession + "at locations" + str(event_location))
			if event_type == "intron retention":
				# like exon skipping but the exon is not flanked by introns (so it's like an intron - might be within an exon)
				# x<>x+1 means that the intron between exons x and exon x+1 is not removed
				include_intron(introns_events, event_location)
			elif event_type == "exon skipped":
				remove_exons(exons_events, event_location)
			elif event_type == "alternative splice acceptor":
				# the 3' exon start boundary is more downstream - in the middle of the exon
				# i.e., a part of the following exon is alternative and the continuance is constitutive)
				update_exon_start(exons_locations_in_genome, )
			elif event_type == "alternative splice donor":
				# the 5' exon end boundary is not consistent (sometimes more upstream - in the middle of the exon)
				update_exon_end(exons_locations_in_genome)
			elif event_type == "novel exon":
				add_exon(introns_events, event_location)
			elif event_type == "exons skipped":
				remove_exons(exons_events, event_location)
			elif event_type == "novel exons":
				add_exon(introns_events, event_location)
			elif event_type == "novel intron":
				update_exon_start(exons_locations_in_genome, )
				update_exon_end(exons_locations_in_genome)
			elif event_type == "---":
				continue
			else:
				print("Event not included:", event_type)


			last_loc = row[EVENT_LOCATION_COL]
			last_event = event_type
			last_gene_ref = gene_accession
			last_entrez_gid = entrez_gid

	#fpw_unreferenced_cdnas.close()


def build_u2os_ge_interval_map():
	ge_dict = parse_gene_expression_file()
	d_plus = defaultdict(interval_defaultdict)
	d_minus = defaultdict(interval_defaultdict)
	for file in U2OS_GE_FILES:
		filename = file + "_parsed_included_regions.csv"
		with open(filename) as fpr:
			next(fpr)
			for line in fpr:
				row = line.strip().split(",")
				if not row[6] in ge_dict.keys():
					d_minus[row[2]][int(row[3]):int(row[4]) + 1] = row[7]
				elif row[5] == "-":
					d_minus[row[2]][int(row[3]):int(row[4]) + 1] = str(ge_dict[row[6]])
				else:
					d_plus[row[2]][int(row[3]):int(row[4]) + 1] = str(ge_dict[row[6]])
	return {"+": d_plus, "-": d_minus}

