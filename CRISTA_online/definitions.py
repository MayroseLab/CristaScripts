import platform, socket, os, subprocess

hostname = socket.gethostname()
print(hostname)

if platform.system() == 'Linux':
	SEP = "/"
	CODE_DIR = "/bioseq/crista/CRISTA_online/"
	PLANTS_MAPPING_FILE = "/bioseq/crista/CRISTA_online/plants_names_mapping.csv"
	DATA_DIR = "/groups/itay_mayrose/shiranabad/CRISPR/data/"
	OUTPUT_DIR = "/groups/itay_mayrose/shiranabad/CRISPR/output/"
	CACHE_DIR = "/groups/itay_mayrose/shiranabad/CRISPR/data/URLcache/"
	ROOT_CRISPR_DIR = "/groups/itay_mayrose/shiranabad/CRISPR/"
	# SAMTOOLS_EXE = "/groups/itay_mayrose/shiranabad/applications/samtools-1.2/samtools "
	GENOMES_PATH = "/groups/itay_mayrose/shiranabad/CRISPR/bwa_genomes/"
	SAMTOOLS_EXE = "samtools "
	BWA_EXE = "bwa"

	TWOBITTOFA_EXE = "/share/apps/bin/twoBitToFa"
else:
	SEP = "\\"
	CODE_DIR = "D:\\Dropbox\\lab\\CRISPR\\crisprProject\\"
	OUTPUT_DIR = "D:\\Dropbox\\lab\\CRISPR\\output\\"
	DATA_DIR = "D:\\Dropbox\\lab\\CRISPR\\data\\"
	ROOT_CRISPR_DIR = "D:\\Dropbox\\lab\\CRISPR\\"
	CACHE_DIR = "D:\\urlcache\\"


HEK293_CELLINE = "HEK293"
K562_CELLINE = "K56"
mESCs_CELLINE = "mESCs"
U2OS_CELLINE = "U2OS"
hIPS_CELLINE = "hIPS"
EL4_CELLINE = "EL4"
CELLS_DICT = {"293": HEK293_CELLINE,
				"293T": HEK293_CELLINE,
				"hIPS": hIPS_CELLINE,
				"HEK293": HEK293_CELLINE,
				"Hek293": HEK293_CELLINE,
				"HEK": HEK293_CELLINE,
				"U2OS": U2OS_CELLINE,
				"K56": K562_CELLINE,
				"k56": K562_CELLINE,
				"K562": K562_CELLINE,
				"k562": K562_CELLINE,
				"invitro": "invitro",
				"A549": "A549",
				"mESCs": mESCs_CELLINE,
				"negatives": "negatives",
                "EL4": EL4_CELLINE}

RNAFOLD_DIR = SEP.join([DATA_DIR, "RNAfold", ""])
HK_GENES_ORIG = SEP.join([DATA_DIR, "HKgenes", "HK_exons_hg19.csv"])
PARSED_HK_GENES = SEP.join([DATA_DIR, "HKgenes", "HK_exons_hg19_parsed.csv"])
HK_GENES = SEP.join([DATA_DIR, "HKgenes", "HK_exons_hg19_final.csv"])

SUBCOMPARTMENTS_2_CATEGORIES_MAP = {"A1": 1, "A2": 2, "B1": 3, "B2": 4, "B3": 5, "B4": 6, "NA": 0}
THREE_D_SUBCOMPARTMENTS = {"hg19": SEP.join([DATA_DIR, "3dmaps", "GSE63525_GM12878_subcompartments.bed"]), "mm9": None}
THREE_D_SUBCOMPARTMENTS_DICT = {"hg19": None, "mm9": None}
THREE_D_SUBCOMPARTMENTS_DICT_PICKLE = {"hg19": SEP.join([DATA_DIR, "3dmaps", "GSE63525_GM12878_subcompartments.pkl"]), "mm9": None}

#UCSC downloaded files
UCSC_FILES_DIR = SEP.join([DATA_DIR, "UCSC_datafiles", ""])
UCSC_PICKES_DIR = SEP.join([UCSC_FILES_DIR, "pickles", ""])
GAP_TABLE_FILENAME_BY_GENOME = {"hg19": UCSC_FILES_DIR+"UCSC_gapTable_hg19.txt",
							"mm9": UCSC_FILES_DIR+"UCSC_gapTable_mm9.txt"}
TELOMERE_DICT_BY_GENOME = {"hg19": None, "mm9": None}
CENTROMERE_DICT_BY_GENOME = {"hg19": None, "mm9": None}
ORGANISM_BY_DB = {"hg19": "Human", "hg38": "Human", "mm9": "Mouse"}
USCS_CHROM_LENGTHS_URL = {"hg19": "http://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes",
						  "mm9": "http://hgdownload-test.cse.ucsc.edu/goldenPath/mm9/encodeDCC/referenceSequences/male.mm9.chrom.sizes"}
CHROM_LENGTHS_DICT_BY_GENOME = {"hg19": None, "mm9": None}
CHROM_LENGTHS_FILENAME_BY_GENOME = {"hg19": UCSC_FILES_DIR + "UCSC_chromosomesLengths_hg19.txt",
									 "mm9": UCSC_FILES_DIR + "UCSC_chromosomesLengths_mm9.txt"}

# CpG islands
CpG_ISLANDS_MAPPING_FILENAME_BY_GENOME = {"hg19": UCSC_FILES_DIR + "UCSC_cpgIslands_hg19.txt",
										  "mm9": UCSC_FILES_DIR + "UCSC_cpgIslands_mm9.txt"}
CpG_ISLANDS_MAPPING_DICT_BY_GENOME = {"hg19": None, "mm9": None}
CpG_ISLANDS_MAPPING_DICT_BY_GENOME_PICKLE = {"hg19": UCSC_PICKES_DIR + "UCSC_cpgIslands_hg19.pkl",
										     "mm9": UCSC_PICKES_DIR + "UCSC_cpgIslands_mm9.pkl"}

# DHS sites
DHS_MAPPING_DICT_BY_GENOME_BY_CELLINE = {"hg19": {HEK293_CELLINE: None, K562_CELLINE: None, hIPS_CELLINE: None, U2OS_CELLINE: None}, "mm9": {mESCs_CELLINE: None}}
DHS_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE = {"hg19": {HEK293_CELLINE: "wgEncodeOpenChromDnaseHek293tPk",
													  hIPS_CELLINE: "wgEncodeOpenChromDnaseIpscwru1Pk", #ips nihi7 is female, 11 is male - the paper doesn't indicate which
													  K562_CELLINE: "wgEncodeOpenChromDnaseK562PkV2",
													  U2OS_CELLINE: "wgEncodeOpenChromDnaseOsteoblPk"}, #the paper http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1356136/ suggests that different cell-types are correlated
											 "mm9": {mESCs_CELLINE: "wgEncodeUwDnaseEse14129olaME0PkRep1"}}
DHS_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE = {"hg19": {HEK293_CELLINE: UCSC_PICKES_DIR + "dhs_hg19_hek293.pkl",
													  hIPS_CELLINE: UCSC_PICKES_DIR + "dhs_hg19_hips.pkl",
													  K562_CELLINE: UCSC_PICKES_DIR + "dhs_hg19_k562.pkl",
													  U2OS_CELLINE: UCSC_PICKES_DIR + "dhs_hg19_u2os.pkl"},
											 "mm9": {mESCs_CELLINE: UCSC_PICKES_DIR + "dhs_mm9_mesc.pkl"}}
DHS_ENCODE_TRACKNAME_BY_GENOME = {"hg19": "wgEncodeOpenChromDnase", "mm9": "wgEncodeUwDnase"}


# nucleosomes
#ATACseq files - from Neville
ATACseq_NUCLEOSOME_MAPPING_DICTS_BY_GENOME = {"hg19": None, "mm9": None}
ATACseq_NUCLEOSOME_MAPPING_DICTS_BY_GENOME_PICKLE = {"hg19": UCSC_PICKES_DIR + "atacseq_nucleosome_occ_hg19.pkl", "mm9": None}
ATACseq_NUCLEOSOME_DATA_FILES_DIR = {"hg19": DATA_DIR + "nucleosome_from_Neville/", "mm9": None}


# from ENCODE
NUCLEOSOME_MAPPING_DICT_BY_GENOME_BY_CELLINE = {"hg19": {HEK293_CELLINE: None, K562_CELLINE: None, hIPS_CELLINE: None, U2OS_CELLINE: None}, "mm9": {mESCs_CELLINE: None}}
NUCLEOSOME_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE = {"hg19": {HEK293_CELLINE: UCSC_PICKES_DIR + "nucleosome_occ_hg19_hek293.pkl",
                                                                K562_CELLINE: UCSC_PICKES_DIR + "nucleosome_occ_hg19_k562.pkl",
                                                                hIPS_CELLINE: UCSC_PICKES_DIR + "nucleosome_occ_hg19_hips.pkl",
                                                                U2OS_CELLINE: UCSC_PICKES_DIR + "nucleosome_occ_hg19_u2os.pkl"},
                                                       "mm9": {mESCs_CELLINE: UCSC_PICKES_DIR + "nucleosome_occ_mm9_mesc.pkl"}}
NUCLEOSOME_DATA_FILE = {"hg19": {HEK293_CELLINE: UCSC_FILES_DIR + "HEK293_hglft_genome_3d13_1580c0.bed",
                                 K562_CELLINE: None, hIPS_CELLINE: None, U2OS_CELLINE: None}, "mm9": {mESCs_CELLINE: None}}
NUCLEOSOME_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE = {"hg19": {HEK293_CELLINE: None,
													  hIPS_CELLINE: None, #ips nihi7 is female, 11 is male - the paper doesn't indicate which
													  K562_CELLINE: "wgEncodeSydhNsomeK562Sig",
													  U2OS_CELLINE: None}, #the paper http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1356136/ suggests that different cell-types are correlated
											 "mm9": {mESCs_CELLINE: "wgEncodeUwDnaseEse14129olaME0PkRep1"}}

#repeats datafile
REPEATS_MAPPING_DICT = {"hg19": None, "hg38": None} #Table Browser: Repeats > RepeatMasker > rmsk
REPEATS_MAPPING_DICT_PICKLE = {"hg19": UCSC_PICKES_DIR + "rmsk_hg19.pkl",
                               "hg38": UCSC_PICKES_DIR + "rmsk_hg38.pkl"}
REPEATS_MAPPING_DICT_wFAMILY = {"hg19": None, "hg38": None}
REPEATS_MAPPING_DICT_PICKLE_wFAMILY = {"hg19": UCSC_PICKES_DIR + "rmsk_hg19_wFam.pkl",
                                     "hg38": UCSC_PICKES_DIR + "rmsk_hg38_wFam.pkl"}

# Methylated RRBS sites
METHYL_MAPPING_DICT_BY_GENOME_BY_CELLINE = {"hg19": {HEK293_CELLINE: None, K562_CELLINE: None, hIPS_CELLINE: None, U2OS_CELLINE: None}, "mm9": {mESCs_CELLINE: None}}
METHYL_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE = {"hg19": {HEK293_CELLINE: UCSC_PICKES_DIR + "methyl_hg19_hek293.pkl",
                                                            K562_CELLINE: UCSC_PICKES_DIR + "methyl_hg19_k562.pkl",
	                                                        hIPS_CELLINE: UCSC_PICKES_DIR + "methyl_hg19_hips.pkl",
                                                            U2OS_CELLINE: UCSC_PICKES_DIR + "methyl_hg19_u2os.pkl"},
                                                   "mm9": {mESCs_CELLINE: UCSC_PICKES_DIR + "methyl_mm9_mesc.pkl",}}
METHYL_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE = {"hg19": {HEK293_CELLINE: ["wgEncodeHaibMethylRrbsHek293StanfordSitesRep1", "wgEncodeHaibMethylRrbsHek293StanfordSitesRep2"],
													  K562_CELLINE: ["wgEncodeHaibMethylRrbsK562HaibSitesRep1", "wgEncodeHaibMethylRrbsK562HaibSitesRep2"],
													  hIPS_CELLINE: ["wgEncodeHaibMethylRrbsH1hescHaibSitesRep1", "wgEncodeHaibMethylRrbsH1hescHaibSitesRep1"],
													  U2OS_CELLINE: ["wgEncodeHaibMethylRrbsHelas3HaibSitesRep1", "wgEncodeHaibMethylRrbsHelas3HaibSitesRep2"]},
											    "mm9": {mESCs_CELLINE: ""}}
METHYL_ENCODE_TRACKNAME_BY_GENOME = {"hg19": "wgEncodeHaibMethylRrbs", "mm9": None}


# Histone modifications
#active marks: h3k4me*, h3k9ac, h3k27ac
#repressive: h3k9me3, h3k27me
ALL_HISTONES_DICT = {HEK293_CELLINE: None, K562_CELLINE: None, U2OS_CELLINE: None, hIPS_CELLINE: None}

HISTONE_H3k4me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE = {"hg19": {HEK293_CELLINE: None, K562_CELLINE: None, U2OS_CELLINE: None, hIPS_CELLINE: None},
											 "mm9": {mESCs_CELLINE: None}}
HISTONE_H3k4me3_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE = {"hg19": {HEK293_CELLINE: ["wgEncodeUwHistoneHek293H3k4me3StdPkRep1", "wgEncodeUwHistoneHek293H3k4me3StdPkRep2"],
													  K562_CELLINE: ["wgEncodeUwHistoneK562H3k4me3StdPkRep1", "wgEncodeUwHistoneK562H3k4me3StdPkRep2"],
													  hIPS_CELLINE: ["wgEncodeUwHistoneH7esH3k4me3StdPkRep1", "wgEncodeUwHistoneH7esH3k4me3StdPkRep2"],
													  U2OS_CELLINE: ["wgEncodeUwHistoneGm12878H3k4me3StdPkRep1", "wgEncodeUwHistoneGm12878H3k4me3StdPkRep2"]},
												 "mm9": {mESCs_CELLINE: ["wgEncodeSydhHistEse14H3k04me1StdPk", "wgEncodeSydhHistEse14H3k04me3StdPk"]}}

HISTONE_H3k27me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE = {"hg19": {HEK293_CELLINE: None, K562_CELLINE: None, U2OS_CELLINE: None, hIPS_CELLINE: None},
											 "mm9": {mESCs_CELLINE: None}}
HISTONE_H3k27me3_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE = {"hg19": {K562_CELLINE: ["wgEncodeUwHistoneK562H3k27me3StdPkRep1", "wgEncodeUwHistoneK562H3k27me3StdPkRep2"],
													  hIPS_CELLINE: ["wgEncodeUwHistoneH7esH3k27me3StdPkRep1", "wgEncodeUwHistoneH7esH3k27me3StdPkRep2"],
													  U2OS_CELLINE: []}}

HISTONE_H3k36me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE = {"hg19": {HEK293_CELLINE: None, K562_CELLINE: None, U2OS_CELLINE: None, hIPS_CELLINE: None},
											 "mm9": {mESCs_CELLINE: None}}
HISTONE_H3k36me3_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE = {"hg19": {K562_CELLINE: ["wgEncodeUwHistoneK562H3k36me3StdPkRep1", "wgEncodeUwHistoneK562H3k36me3StdPkRep2"],
													  hIPS_CELLINE: ["wgEncodeUwHistoneH7esH3k36me3StdPkRep1", "wgEncodeUwHistoneH7esH3k36me3StdPkRep2"],
													  U2OS_CELLINE: ["wgEncodeSydhHistoneU2osH3k36me3bUcdPk"]}}

HISTONE_H3k9me3_MAPPING_DICTS_BY_GENOME_BY_CELLINE = {"hg19": {HEK293_CELLINE: None, K562_CELLINE: None, U2OS_CELLINE: None, hIPS_CELLINE: None},
											 "mm9": {mESCs_CELLINE: None}}
HISTONE_H3k9me3_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE = {"hg19": {K562_CELLINE: ["wgEncodeBroadHistoneK562H3k9me3StdPk", ],
													  U2OS_CELLINE: ["wgEncodeSydhHistoneU2osH3k9me3UcdPk", "wgEncodeSydhHistoneU2osH3k36me3bUcdPk"]},
												 "mm9": {mESCs_CELLINE: ["wgEncodeSydhHistEse14H3k09me3StdPk"]}}

HISTONE_ENCODE_TRACKNAME_BY_CELLINE = {"hg19": {HEK293_CELLINE: "wgEncodeUwHistone", K562_CELLINE: "wgEncodeUwHistone",
												hIPS_CELLINE: "wgEncodeUwHistone", U2OS_CELLINE: "wgEncodeSydhHistone"},
									   "mm9": {mESCs_CELLINE: "wgEncodeSydhHist"}}


#Nucleosome occupancy
# GEO source: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM971947
# cite: Leroy, G., Chepelev, I., Dimaggio, P.A., Blanco, M.A., Zee, B.M., Zhao, K., and Garcia, B.A.
#(2012). Proteogenomic characterization and mapping of nucleosomes decoded by Brd and HP1 proteins. Genome Biol. 13, R68.
# converted with UCSC liftOver to hg19


DNA_PAIRS_THERMODYNAMICS = {"AA": 9.1, "AT": 8.6, "TA": 6.0, "CA": 5.8, "GT": 6.5, "CT": 7.8, "GA": 5.6, "CG": 11.9,
							"GC": 11.1, "GG": 11.0, "TT": 9.1, "TG": 5.8, "AC": 6.5, "AG": 7.8, "TC": 5.6, "CC": 11.0} #Breslauer et al.
#DNA_PAIRS_THERMODYNAMICS = {"AA": 1.0, "AT": 0.88, "TA": 0.58, "CA": 1.45, "GT": 1.44, "CT": 1.28, "GA": 1.3, "CG": 2.17,
#							"GC": 2.24, "GG": 1.84, "TT": 1.0, "TG": 1.45, "AC": 1.44, "AG": 1.28, "TC": 1.3, "CC": 1.84} #Santalucia

# These are calculated by the "harmonic" mean (1) or mean (2) according to the paper http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2932579/#SEC3
# The calculated means are not identical to what they have in the paper (Table S3) because they used another parameters and estimated those accordingly
DNA_PAIRS_PERSISTANCE_LEN = {"AA": 50.4, "AT": 40.9, "TA": 44.7, "CA": 46.7, "GT": 55.4, "CT": 51.0, "GA": 54.4, "CG": 56.0,
							"GC": 44.6, "GG": 41.7, "TT": 50.4, "TG": 46.7, "AC": 55.4, "AG": 51.0, "TC": 54.4, "CC": 41.7}
DNA_HELICAL_REPEAT = {"AA": 10.27, "AT": 10.59, "TA": 10.65, "CA": 10.6, "GT": 10.58, "CT": 10.42, "GA": 10.36, "CG": 10.35,
					  "GC": 10.43, "GG": 10.76, "TT": 10.27, "TG": 10.6, "AC": 10.58, "AG": 10.42, "TC": 10.36, "CC": 10.76}



# DGF sites
EXON_MAPPING_DICT_BY_GENOME_BY_CELLINE = {"hg19": {HEK293_CELLINE: None, K562_CELLINE: None, hIPS_CELLINE: None, U2OS_CELLINE: None}}
EXON_MAPPING_DICT_BY_GENOME_BY_CELLINE_PICKLE = {"hg19": {HEK293_CELLINE: UCSC_PICKES_DIR + "exon_hg19_hek293.pkl",
                                                           K562_CELLINE: UCSC_PICKES_DIR + "exon_hg19_k562.pkl",
                                                           U2OS_CELLINE: UCSC_PICKES_DIR + "exon_hg19_u2os.pkl",
                                                           hIPS_CELLINE: UCSC_PICKES_DIR + "exon_hg19_hips.pkl"}}
EXON_ENCODE_TABLENAME_BY_GENOME_BY_CELLINE = \
	{"hg19": {K562_CELLINE: ["wgEncodeDukeAffyExonK562SimpleSignalRep1V2", "wgEncodeDukeAffyExonK562SimpleSignalRep2V2",
							   "wgEncodeDukeAffyExonK562SimpleSignalRep3V2", "wgEncodeDukeAffyExonK562SimpleSignalRep4"],
			  HEK293_CELLINE: ["wgEncodeDukeAffyExonHek293tSimpleSignalRep1", "wgEncodeDukeAffyExonHek293tSimpleSignalRep2"],
			  hIPS_CELLINE: ["wgEncodeDukeAffyExonH7esSimpleSignalRep1", "wgEncodeDukeAffyExonH7esSimpleSignalRep2",
							 "wgEncodeDukeAffyExonH7esSimpleSignalRep3"]},
	 "mm9": {mESCs_CELLINE: None}}
#	         U2OS_CELLINE: ["wgEncodeDukeAffyExonHelas3SimpleSignalRep1V2", "wgEncodeDukeAffyExonHelas3SimpleSignalRep2V2", "wgEncodeDukeAffyExonHelas3SimpleSignalRep3V2"],
EXON_ENCODE_TRACKNAME_BY_GENOME = {"hg19": "wgEncodeDukeAffyExon", "mm9": None}

GENES_MAPPING_DICT_BY_GENOME = {"hg19": None}
GENES_MAPPING_DICT_BY_GENOME_PICKLE = {"hg19": UCSC_PICKES_DIR + "genes_mapping_hg19.pkl"}
GENES_TABLENAME_BY_GENOME = {"hg19": "refGene"}
EXONS_MAPPING_DICT_BY_GENOME = {"hg19": None}
EXONS_MAPPING_DICT_BY_GENOME_PICKLE = {"hg19": UCSC_PICKES_DIR + "exons_mapping_hg19.pkl"}

NUC_OCC_PREDICTIONS_DIR = {"hg19": SEP.join([DATA_DIR, "nucleosome_occ_pred", "fasta_inputs", ""])}
NUC_OCC_PREDICTIONS_MAPPING_DICT_BY_GENOME = {"hg19": None}
NUC_OCC_PREDICTIONS_MAPPING_DICT_BY_GENOME_PICKLE = {"hg19": UCSC_PICKES_DIR + "nuc_occ_prediction_mapping_hg19.pkl"}

#U2OS gene expression files
# GEO: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22149
# paper: https://nar.oxfordjournals.org/content/early/2015/12/30/nar.gkv1500.full
U2OS_GE_FILES = [SEP.join([DATA_DIR, "u2os_geneExpression", "GSE22149_SpliceArray_Analysis_Results_ALL_Evidenced_Probesets.txt"]),
                 SEP.join([DATA_DIR, "u2os_geneExpression", "GSE22149_SpliceArray_Analysis_Results_ALL_Discovery_Probesets.txt"])]


# HK genes
HK_GENES_MAPPING_FILENAME_BY_GENOME = {"hg19": PARSED_HK_GENES,
										  "mm9": ""}
HK_GENES_MAPPING_DICT_BY_GENOME = {"hg19": None, "mm9": None}
HK_GENES_MAPPING_DICT_BY_GENOME_PICKLE = {"hg19": UCSC_PICKES_DIR + "hk_hg19.pkl",
                                          "mm9": UCSC_PICKES_DIR + "hk_mm9.pkl"}

#sources db files
source_path = DATA_DIR
GUIDESEQ_SOURCE = SEP.join([source_path, "guideSeq", "nbt.3117-S2.csv"])
GUIDESEQ_HF1_SOURCE = SEP.join([source_path, "guideSeq", "nature16526-s5.csv"])
#IDLV_SOURCE = SEP.join([source_path, "IDLV", "suppTable5-CR-4.csv"])
IDLV_SOURCES = [SEP.join([source_path, "IDLV", "suppTable2-CR-4_offtargets.csv"]),
				SEP.join([source_path, "IDLV", "suppTable3-CR-5_offtargets.csv"])]
CHO_SOURCE = SEP.join([source_path, "Cho", "Supplemental_Table_2.csv"])
HTGTS_SOURCE = SEP.join([source_path, "HTGTS", "htgts_ots.csv"])
MULTIPLEX_DIGENOMESEQ_SOURCE = SEP.join([source_path, "multiplex_digenomeseq", "corrected_muldiseq_data.csv"])
DIGENOMESEQ_SOURCE = SEP.join([source_path, "digenomeSeq", "digenomeSeq_ots2.csv"])
DIGENOMESEQ_VALIDATED_SOURCE = SEP.join([source_path, "digenomeSeq", "digenomeSeq_ots2_validated.csv"])
KUSCO_SOURCE = SEP.join([source_path, "Kusco", "OTs.csv"])
MOUSECELLS_SOURCE = SEP.join([source_path, "MouseCells", "OTs.csv"])
BLESS_SOURCE = SEP.join([source_path, "BLESS", "nature14299-f3.csv"])
BLESS_SLAYMAKER_SOURCE = SEP.join([source_path, "BLESS", "slaymaker_data.csv"])
EPienTan_SOURCE = SEP.join([source_path, "HumanIPS_MouseES", "supptable#.csv"])
NEGATIVES_DIR = SEP.join([source_path, "final_potentials", ""])
NEGATIVES_W_FEATURES_DIR = SEP.join([OUTPUT_DIR, "negatives", ""])
POTENTIAL_LIBRARY = SEP.join([DATA_DIR, "potential_sites_library", ""])
BLESS_VALIDATION_DIR = SEP.join([DATA_DIR, "BLESS", "validation", ""])
HSU_SOURCE_DIR = SEP.join([DATA_DIR, "hsu", ""])
HAEUSSLER_SOURCE_FILE = SEP.join([DATA_DIR, "haeussler", "13059_2016_1012_MOESM14_ESM.tsv"])

#parsed db files
PARSED_GUIDESEQ_FILE = SEP.join([source_path, "guideSeq", "parsed_nbt.3117-S2.csv"])
PARSED_GUIDESEQ_HF1_FILE = SEP.join([source_path, "guideSeq", "parsed_nature16526-s5.csv"])
#PARSED_IDLV_FILE = SEP.join([source_path, "IDLV", "parsed_suppTable5-CR-4.csv"])
PARSED_IDLV_FILE = SEP.join([source_path, "IDLV", "parsed_suppTable2.3_CR4.5.csv"])
PARSED_CHO_SOURCE = SEP.join([source_path, "Cho", "parsed_data.csv"])
PARSED_HTGTS_FILE = SEP.join([source_path, "HTGTS", "parsed_htgts_ots.csv"])
TRUE_HTGTS_FILE = SEP.join([source_path, "HTGTS", "true_htgts_sites.csv"])
PARSED_MULTIPLEX_DIGENOMESEQ_FILE = SEP.join([source_path, "multiplex_digenomeseq", "parsed_muldiseq_data.csv"])
PARSED_DIGENOMESEQ_FILE = SEP.join([source_path, "digenomeSeq", "parsed_digenomeSeq_ots2.csv"])
PARSED_DIGENOMESEQ_VALIDATED_FILE = SEP.join([source_path, "digenomeSeq", "parsed_digenomeSeq_ots2_validated.csv"])
PARSED_KUSCO_FILE = SEP.join([source_path, "Kusco", "parsed_OTs.csv"])
PARSED_MOUSECELLS_FILE = SEP.join([source_path, "MouseCells", "parsed_OTs.csv"])
PARSED_BLESS_FILE = SEP.join([source_path, "BLESS", "parsed_OTs.csv"])
PARSED_SLAYMAKER_FILE = SEP.join([source_path, "BLESS", "parsed_slaymaker_OTs.csv"])
PARSED_EPienTan_ALL_FILE = SEP.join([source_path, "HumanIPS_MouseES", "parsed_ots_all.csv"])
PARSED_EPienTan_1Human_FILE = SEP.join([source_path, "HumanIPS_MouseES", "parsed_ots_assembly_unique.csv"])
PARSED_COSMID_DIR = SEP.join([source_path, "COSMID"])
PARSED_BLESS_VALIDATION = SEP.join([DATA_DIR, "BLESS", "parsed_bless_validation.csv"])
PARSED_HSU_SOURCE = SEP.join([DATA_DIR, "hsu", "parsed_data.csv"])
PARSED_HAEUSSLER_SOURCE_FILE = SEP.join([DATA_DIR, "haeussler", "parsed_haeussler.csv"])

PARSED_ALL_SITES_HG38 = SEP.join([source_path, "alldata", "all_4methods_sites_wo_filtering_hg38.csv"]) #I did liftOver online

ORIG_ALIGNED_DIFF_FILE = SEP.join([source_path, "alldata", "orig_aligned_diff.csv"])

SGRNAS_LIST_FILENAME = SEP.join([DATA_DIR, "sgRNAs_list.txt"])
SGRNAs_FEATURES_DICT = None

#Validation definitions
VALIDATION_DIR = SEP.join([DATA_DIR, "validationSgs1", ""])
CDS_COL = 108
EXON_COL = 94
START_POS_COL = 6
END_POS_COL = 7
CHROM_POS_COL = 8


LOCAL_DB_DIRS = {db: "/groups/itay_mayrose/shiranabad/CRISPR/reference_dbs/{0}/".format(db) for db in
				 ["hg19", "mm9", "ce6", "danRer10", "danRer7", "dm3", "rn5", "ci2", "mm10", "hg38"]}
CHOP_LENGTH = 10000000

PA_SCORE_THRES = 14.75 # as 95% of the data
MATCH_SCORE = 1
MISMATCH_PENALTY = 0
GAP_PENALTY = -1.25
NTREES = 1000

#MAX_OBSERVATION_FREQ = 1.0
#MIN_OBSERVATION_FREQ = float("inf")

DNASHAPE_TYPES = ["MGW", "HelT", "ProT", "Roll"]
DNASHAPE_DICT_FILE = DATA_DIR + "dnaShape.pkl"
DNASHAPE_DICT = None

#sources
GUIDESEQ_SN = "guideSeq"
GUIDESEQ_HEK_SN = "guideSeq_hek"
GUIDESEQ_HF1_SN = "guideSeq_HF1"
IDLV_SN = "IDLV"
CHO_SN = "CHO"
HTGTS_SN = "HTGTS"
DIGENOMESEQ_SN = "digenomeSeq"
MULTIPLEX_DIGENOMESEQ_SN = "multiplex_digenomeSeq"
KUSCO_SN = "KUSCO"
MOUSECELLS_SN = "MOUSECELLS"
BLESS_SN = "BLESS"
BLESSVALID_SN = "BLESS_validation"
EPienTan_SN = "ePienTan"
NEGATIVES = "negatives"
HSU_SN = "Hsu2013"
HAEUSSLER_SN = "haeussler"


GUIDESEQ_HEK_sgRNAs = ['HEK293_sgRNA1', 'HEK293_sgRNA2', 'HEK293_sgRNA3', 'HEK293_sgRNA4']
GUIDESEQ_U2OS_sgRNAs = ['EMX1', 'FANCF', 'RNF2', 'VEGFA_site1', 'VEGFA_site2', 'VEGFA_site3'] #'tru_EMX1', 'tru_VEGFA_site1', 'tru_VEGFA_site3',
GUIDESEQ_sgRNAs = GUIDESEQ_U2OS_sgRNAs + GUIDESEQ_HEK_sgRNAs
GUIDESEQ_hf1_sgRNAs = ["EMX1-2", "FANCF-2", "FANCF-3", "FANCF-4", "RUNX1-1", "ZSCAN2"]
GUIDESEQ_hf1_sgRNAs_redundant = ["EMX1_1", "FANCF-1", "VEGFA_2", "VEGFA_3"]
IDLV_sgRNAs = ["WAS_CR-5", "WAS_CR-3", "TAT_CR-4", "TAT_CR-1", "TAT_CR-6", "WAS_CR-4", "VEGFA_site_3"]
CHO_sgRNAs = ["CCR5", "CCR5-1", "CCR5-2", "CCR5-3", "CCR5-4", "CCR5-5", "CCR5-6", "CCR5-7", "CCR5-8", "C4BPB"]
HTGTS_sgRNAs = ["Chr12 47.0Mb OT2", "Chr19 DAZAP1 OT1", "Chr7 103.6Mb OT17", "RAG1A", "RAG1B", "RAG1C", "RAG1D"]
HTGTS_sgRNAs_redundant = ["EMX1","VEGFA"]
DIGENOMESEQ_sgRNAs = ['HBB']
DIGENOMESEQ_sgRNAs_redundant = ['VEGFA']
MULTIPLEX_DIGENOMESEQ_sgRNAs_redundant = ["MDS_" + x for x in ["VEGFA1", "VEGFA2", "VEGFA3", "EMX1", "FANCF", "RNF2",
                                              "HEK293-1", "HEK293-2", "HEK293-3", "HEK293-4"]]
MULTIPLEX_DIGENOMESEQ_sgRNAs_unique = ["MDS_HBB"]
KUSCO_sgRNAs = ["sgRNA1", "sgRNA2", "sgRNA3", "sgRNA4", "sgRNA5", "sgRNA6", "sgRNA7", "sgRNA8", "sgRNA9", "sgRNA10", "sgRNA11", "sgRNA12"]
MOUSECELLS_sgRNAs = ["nanog-sg2", "nanog-sg3", "phc1-sg1", "phc1-sg2"]
BLESS_sgRNAs = ["Sp_sg1", "Sp_sg2"]
SLAYMAKER_sgRNAs = ["VEGFA(1)", "EMX1(1)"]
BLESS_VALIDATION_sgRNAs = ["DNMT1-4", "EMX1-4", "GRIN2b-7", "VEGFA8"]
EPienTan_Human_iPS_sgRNAs = [HEK293_CELLINE + "_" + x for x in ["AAVS1", "AKT2", "ATPA", "ATPB", "CDKA", "CDKB", "SLCA", "SLCB"]]
EPienTan_Human_HEK_sgRNAs = [hIPS_CELLINE + "_" + x for x in ["AAVS1", "AKT2", "ATPA", "ATPB", "CDKA", "CDKB", "SLCA", "SLCB"]]
EPienTan_Mouse_sgRNAs = [mESCs_CELLINE + "_" + x for x in ["mAtpA", "mCdA", "mCdB", "mSlcA"]]

RAG1C_SEQ = "GCACCTAACATGATATATTAAGG"
RAG1D_SEQ = "GACCTTAAGGTTTTTGTGGAAGG"

SOURCES_GUIDES = {GUIDESEQ_SN: GUIDESEQ_U2OS_sgRNAs,
                  GUIDESEQ_HEK_SN: GUIDESEQ_HEK_sgRNAs,
                  GUIDESEQ_HF1_SN: GUIDESEQ_hf1_sgRNAs,
				  IDLV_SN: IDLV_sgRNAs,
				  HTGTS_SN: HTGTS_sgRNAs,
				  MULTIPLEX_DIGENOMESEQ_SN: MULTIPLEX_DIGENOMESEQ_sgRNAs_unique,
				  KUSCO_SN: KUSCO_sgRNAs,
				  MOUSECELLS_SN: MOUSECELLS_sgRNAs,
				  BLESS_SN: BLESS_sgRNAs,
				  BLESSVALID_SN: BLESS_VALIDATION_sgRNAs,
				  EPienTan_SN: EPienTan_Human_iPS_sgRNAs + EPienTan_Human_HEK_sgRNAs}

SOURCES_GUIDES_w_REDUNDANCY = {GUIDESEQ_SN: GUIDESEQ_U2OS_sgRNAs,
                               GUIDESEQ_HF1_SN: GUIDESEQ_hf1_sgRNAs + GUIDESEQ_hf1_sgRNAs_redundant,
                               GUIDESEQ_HEK_SN: GUIDESEQ_HEK_sgRNAs,
                               HTGTS_SN: HTGTS_sgRNAs + HTGTS_sgRNAs_redundant,
                               DIGENOMESEQ_SN: DIGENOMESEQ_sgRNAs + DIGENOMESEQ_sgRNAs_redundant,
                               MULTIPLEX_DIGENOMESEQ_SN: MULTIPLEX_DIGENOMESEQ_sgRNAs_unique + MULTIPLEX_DIGENOMESEQ_sgRNAs_redundant,
                               BLESS_SN: BLESS_sgRNAs + SLAYMAKER_sgRNAs,
                               BLESSVALID_SN: BLESS_VALIDATION_sgRNAs
                               }

CORRECTION_FACTORS = None
"""{HTGTS_SN: 10.8,
                      DIGENOMESEQ_SN: 259.58,
                      GUIDESEQ_SN: 1.0,
                      GUIDESEQ_HEK_SN: 4590/2136.0,
                      BLESS_SN: 220.0,
                      BLESSVALID_SN: 220.0
                      }
"""


# validation on-targets:
CHARI_DIRPATH = SEP.join([DATA_DIR, "Chari", ""])
CHARI_GUIDES_FILEPATH = CHARI_DIRPATH + "suppdata6_guides.csv"
CHARI_HEK293t_DATA_FILEPATH = CHARI_DIRPATH + "suppfigure11data_endogenous293t_SpCas.csv"
PARSED_CHARI_SOURCE = CHARI_DIRPATH + "parsed_suppfigure11data_endogenous293t_SpCas.csv"

DOENCH_DIRPATH = SEP.join([DATA_DIR, "doench", ""])
DOENCH14_SUPPT7_SOURCE = SEP.join([DOENCH_DIRPATH, "nbt.3026-S5.csv"])
PARSED_DOENCH14_SUPPT7_SOURCE_ALL = DOENCH_DIRPATH + "parsed_2014_suppt7.csv"
PARSED_DOENCH14_SUPPT7_SOURCE_WOMOUSE = DOENCH_DIRPATH + "parsed_womouse_2014_suppt7.csv"

DOENCH16_SUPPT20_SOURCE = DOENCH_DIRPATH + "doench2016_stable20_H2D_H2K.csv"
PARSED_DOENCH16_SUPPT20_SOURCE = DOENCH_DIRPATH + "parsed_doench2016_stable20_H2D_H2K.csv"

PARSED_DOENCH16_SUPPT18_SOURCE = DOENCH_DIRPATH + "parsed_2016_suppt18.csv"
DOENCH16_SUPPT18_SOURCE = DOENCH_DIRPATH + "doench2016_stable18_cd33offtargets.csv"
DOENCH16_SUPPT18_ANNOTATIONS = DOENCH_DIRPATH + "doench2016_stable18_annotations.csv"



DATASET_HEADER_SOURCES = ["source", "target_name", "sample_kind", "sgRNA_seq", "offtarget_seq",
                          "strand", "startpos", "endpos"]
DATASET_HEADER_NONGENOMIC_MAT = ["pa_score", "total_bulges", "rna_bulges", "dna_bulges", "PAM_2_first",
					"PAM_N_id", "last_nucleotide", "mms_cnt",
					"continuous_inconsistencies_cnt", "continuous_inconsistencies_length"] +\
					["mismatches_1_4", "mismatches_5_8", "mismatches_9_12", "mismatches_13_16", "mismatches_17_end"] +\
					["wobble_total", "YY_total", "RR_total", "Tv_total"] +\
					["PamProximalCouple1", "PamProximalCouple2", "PamProximalCouple3", "PamProximalCouple4"] +\
					["id_pos" + str(i) for i in range(5, 0, -1)] +\
					["extended_genomic_GC_content", #"upstream_50_extension_gc", "downstream_50_extension_gc",
					 "dna_enthalpy", "extended_dna_enthalpy",
					 'nA', 'nC', 'nT', 'nG'] +\
					 ["nuc_id_afterPAM1", "nuc_id_afterPAM2", "nuc_id_afterPAM3", "nuc_id_afterPAM4", "nuc_id_afterPAM5"] +\
					 ["min_MGW_N", "avg_HelT_NG", "avg_Roll_NG", "avg_ProT_N"] +\
					 ["MGW_NxxGxxG", "HelT_NxxGGxx", "Roll_NxxGGxx", "ProT_NxxGxxG"]
DATASET_HEADER_WGENOMIC_MAT = ["chromosome"] + DATASET_HEADER_NONGENOMIC_MAT + ["Distance_from_centromere", "Distance_from_telomere",
					  "DHS_signalValue", "distance_from_DHS", "nucleosome_mapping", "length_from_nucleosome",
					  "genomic_subcompartment",
					  "NGG_exon_signalValue", "opposite_exon_signalValue",
					  "in_NGG_exon", "in_opposite_exon", "in_tx", "in_cds"]

DATASET_HEADER_MAT = DATASET_HEADER_WGENOMIC_MAT

CCTOP_SCORE_HEADER = "CCTop_score"
OPTCD_SCORE_HEADER = "OptCD_score"
CFD_SCORE_HEADER = "CFD_score"
AZIMUTH_SCORE_HEADER = "azimuth_score"
# [model + "_" + str(i) for model in ["AAWEDGE", "TRIFONOV", "DESANTIS", "CALLADINE", "TRINUCLEOTIDE", "NUCLEOSOME", "REVERSED"] for i in range(16)] +\

DATASET_HEADER = DATASET_HEADER_SOURCES + DATASET_HEADER_MAT + ["intensity"]
DATASET_HEADER_WO_GENOMIC = DATASET_HEADER_SOURCES + DATASET_HEADER_NONGENOMIC_MAT + ["intensity"]

DATASET_HEADER_FOR_VALIDATION = DATASET_HEADER[:-1] + ["prediction", "NGG", "NAG", "CFD_score", "Azimuth_score", "cctop", "CRISPR_MIT"]
RF_PICKLE = "/groups/itay_mayrose/shiranabad/CRISPR/tests/final_model/test_20170108-153323_r100_f5_n1_lro_common_unique_woHtgts_woMds/0/RFpredictor.pkl" #"/groups/itay_mayrose/shiranabad/CRISPR/tests/final_model/test_20160817-180711_r100_f10_lro_common_unique_woHtgts/0/RFpredictor.pkl"
RF_MODEL_DIR = "/bioseq/crista/rf/"
RF_MODEL_DIR_WOGENOMIC = "/bioseq/crista/rf_nogenomic/"
RF_MODEL_DIR_WOFLANKING = "/bioseq/crista/rf_noflanking/"

FEATURES_NAMES_MAPPING_FILE = DATA_DIR + "features_display_name_mapping.csv"

#"histone_h3k4me_signalValue", "histone_h3k9me_signalValue", "histone_h3k27me_signalValue", "histone_h3k36me_signalValue",

#VALIDATION
TESTED_GENES_FOR_VALIDATION = ["Chr19_DAZAP1_OT1_OT6.csv", "HEK293_sgRNA4_OT12.csv", "Sp_sg2_OT3.csv", "VEGFA_site1_OT1.csv",
					"RAG1A_OT5.csv", "Chr12_47.0Mb_OT2_OT5.csv", "FANCF-2_OT1.csv", "HBB_OT1.csv", "FANCF-4_OT1.csv", "FANCF_OT1.csv", "FANCF-3_OT1.csv", "ZSCAN2_OT1.csv","Sp_sg2_OT4.csv",
			        "Chr12_47.0Mb_OT2_OT1.csv", "Chr7_103.6Mb_OT17_OT5.csv", "VEGFA_site2_OT13.csv", "RNF2_OT1.csv", "VEGFA_site1_OT2.csv"]

#REPEATS
NOREPEATS_PHANTOM_FILES_DIR = "/groups/itay_mayrose/shiranabad/CRISPR/data/validationSgs1/phantom/"
BLAT_EXE = "/share/apps/blat/blat"
repclass_dict = {'repClass': 1, 'Unknown': 2, 'srpRNA': 3, 'DNA': 4, 'Other': 5, 'scRNA': 6, 'RC': 11,
                 'LINE': 8, 'Simple_repeat': 12, 'SINE': 10, 'snRNA': 9, 'rRNA': 7, 'tRNA': 14,
                 'Low_complexity': 15, 'LTR': 16, 'Satellite': 13, 'RNA': 17, 'Retroposon': 18, "Unmasked": 0}

repfamily_dict = {'DNA/DNA': 1, 'DNA/hAT': 2, 'DNA/hAT-Ac': 3, 'DNA/hAT-Blackjack': 4, 'DNA/hAT-Charlie': 5,
                  'DNA/hAT-Tag1': 6, 'DNA/hAT-Tip100': 7, 'DNA/Merlin': 8, 'DNA/MuDR': 9, 'DNA/MULE-MuDR': 10,
                  'DNA/PIF-Harbinger': 11, 'DNA/PiggyBac': 12, 'DNA/TcMar': 13, 'DNA/TcMar-Mariner': 14,
                  'DNA/TcMar-Pogo': 15, 'DNA/TcMar-Tc2': 16, 'DNA/TcMar-Tigger': 17, 'LINE/CR1': 18, 'LINE/Dong-R4': 19,
                  'LINE/L1': 20, 'LINE/L2': 21, 'LINE/Penelope': 22, 'LINE/RTE': 23, 'LINE/RTE-BovB': 24,
                  'LINE/RTE-X': 25, 'Low_complexity/Low_complexity': 26, 'LTR/LTR': 27, 'LTR/ERV': 28, 'LTR/ERV1': 29,
                  'LTR/ERVK': 30, 'LTR/ERVL': 31, 'LTR/ERVL-MaLR': 32, 'LTR/Gypsy': 33, 'Other/Other': 34,
                  'RC/Helitron': 35, 'Retroposon/SVA': 36, 'RNA/RNA': 37, 'rRNA/rRNA': 38, 'Satellite/Satellite': 39,
                  'Satellite/acro': 40, 'Satellite/centr': 41, 'Satellite/telo': 42, 'scRNA/scRNA': 43,
                  'Simple_repeat/Simple_repeat': 44, 'SINE/SINE': 45, 'SINE/5S-Deu-L2': 46, 'SINE/Alu': 47,
                  'SINE/Deu': 48, 'SINE/MIR': 49, 'SINE/tRNA': 50, 'SINE/tRNA-Deu': 51, 'SINE/tRNA-RTE': 52,
                  'snRNA/snRNA': 53, 'srpRNA/srpRNA': 54, 'tRNA/tRNA': 55, 'Unknown/Unknown': 66, 'Unmasked/Unmasked': 0}

inv_repclass_dict = {v: k for k, v in repclass_dict.items()}
inv_repfamily_dict = {v: k for k, v in repfamily_dict.items()}

CRISTA_DISPLAY_NAME = r"$CRISTA\/$"
CRISTA_plus_DISPLAY_NAME = r"$CRISTA^+\/$"
CCTOP_DISPLAY_NAME = r"$CCTop\/$"
CrDESIGN_DISPLAY_NAME = r"$OptCD\/$"
CFD_DISPLAY_NAME = r"$CFD\/$"


def get_features_names_mapping(as_dict=False):
	"""
	:return: mapping (pd or dict) of features codes to display names
	"""
	import pandas as pd
	features_map_df = pd.read_csv(FEATURES_NAMES_MAPPING_FILE, index_col=0, sep=",")
	if as_dict:
		return features_map_df.to_dict
	return features_map_df


########## filenames getters ##########
def get_sg_negatives_with_features(sgseq, source_name):
	return NEGATIVES_W_FEATURES_DIR + source_name + "_" + sgseq[:-3] + "_negatives.csv"

def get_sg_positives_with_features(source_name, target_name):
	return SEP.join([OUTPUT_DIR, source_name, target_name + "_positives.csv"])


def get_sg_potential_targetsites(sgseq, genome_db):
	return NEGATIVES_DIR + genome_db + "_" + sgseq[:-3] + "_final.csv"


def get_sg_potential_targetsites_for_validation(sgseq, genome_db):
	return VALIDATION_DIR + "/offtargets/" + genome_db + "_" + sgseq[:-3] + "_final.csv"


def get_sg_potential_that_are_positives_log(sgseq, source_name):
	return NEGATIVES_DIR + source_name + "_" + sgseq[:-3] + "_log.txt"


def get_sg_all_potential_targetsites_dir(sgseq, genome_db):
	return POTENTIAL_LIBRARY + genome_db + "_" + sgseq[:-3] + SEP


def get_ucsc_filename(tablename, cell_line, db_assembly):
	return str.format(UCSC_FILES_DIR + "UCSC_{0}_{1}_{2}.txt", tablename, cell_line, db_assembly)


def get_refseqgene_filename(tablename, position):
	return str.format(SEP.join([UCSC_FILES_DIR, "Refseq", ""]) + "UCSC_{0}_{1}.txt", tablename, position)


def get_ucsc_picle_filename(table_type, cell_line, db_assembly):
	return str.format(UCSC_FILES_DIR + "UCSC_{0}_{1}_{2}.pl", table_type, cell_line, db_assembly)


def get_random_validation_targets_filename(chop_file_name):
	filename = chop_file_name.split(SEP)[-1]
	return VALIDATION_DIR + 'genomic_targets' + SEP + filename + "_targets.csv"


def get_ml_random_dataset(set_num):
	return "/groups/itay_mayrose/shiranabad/CRISPR/data/machineLearningRandomSets/set-" + str(set_num) + ".pkl"


def get_repeats_summary_filename(sgseq):
	"""
	:param sgseq: 20bp + .GG
	:return:
	"""
	return SEP.join([VALIDATION_DIR, "perfect_matches_per_guide", sgseq[:-3]])


def get_repeats_phantom_filename(sgseq):
	"""
	:param sgseq: 20bp + .GG
	:return: the filename of a rna with no repeats
	"""
	return NOREPEATS_PHANTOM_FILES_DIR + sgseq[:-3]



def get_matrix_header(include_genomic_features=True):
	if include_genomic_features:
		return DATASET_HEADER
	else:
		return DATASET_HEADER_WO_GENOMIC


RESULTS_DIR = r"D:\Dropbox\lab\CRISPR\data\results_files\\"
FEATURE_IMP_FILE = RESULTS_DIR + "feature_imp_avg.csv"
FWDSEL_FEATURE_IMP_FILE = RESULTS_DIR + "fwdsel_filtered_feature_imp_avg.csv"
FWDSEL_FEATURES_ACCURACY_SCORES_FILE = RESULTS_DIR + "selected_features.csv"
BEF_AFT_PA_FILE = RESULTS_DIR + "orig_aligned_diff.csv"
PREDICTION_RESCALING_FACTOR = 8.22
OBSERVED_RESCALING_FACTOR = 10.0


CHROMOSOME_HEADING = "chromosome"
STARTPOS_HEADING = "start position"
ENDPOS_HEADING = "end position"
STRAND_HEADING = "strand"
SGRNA_HEADING = "aligned sgRNA"
OFFTARGET_HEADING = "aligned site"

global ERROR_MESSAGE_ON_EXCEPTION
ERROR_MESSAGE_ON_EXCEPTION = ""

global WARNING_MESSAGES
WARNING_MESSAGES = []


def add_warning_message(msg):
	global WARNING_MESSAGES
	WARNING_MESSAGES.append(msg)


def get_warning_message():
	global WARNING_MESSAGES
	return "- " + "<br>- ".join(WARNING_MESSAGES)


def set_error_message(msg):
	global ERROR_MESSAGE_ON_EXCEPTION
	if ERROR_MESSAGE_ON_EXCEPTION == "":
		ERROR_MESSAGE_ON_EXCEPTION=msg
		return True
	else:
		return False

def get_error_message():
	global ERROR_MESSAGE_ON_EXCEPTION
	return ERROR_MESSAGE_ON_EXCEPTION