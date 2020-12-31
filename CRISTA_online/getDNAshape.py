import random
from definitions import *
from utils import *


def get_DNAshape_features(dna_seq, replace_N_randomly=False):

	global DNASHAPE_DICT
	if DNASHAPE_DICT is None:
		if os.path.exists(DNASHAPE_DICT_FILE):
			DNASHAPE_DICT = pickle.load(open(DNASHAPE_DICT_FILE, "rb"))
		#else:
			#return run_DNAshape(dna_seq, trim=True)

	mgw = [None]
	roll = [None]
	prot = [None]
	helt = [None]

	for i in range(2, len(dna_seq)-2):
		current_heptamer = dna_seq[i-2 : i+3]
		if replace_N_randomly:
			current_heptamer = re.sub("[^ACGT]", random.choice(["A", "C", "G", "T"]), current_heptamer)
		current_nucleotide = DNASHAPE_DICT[current_heptamer]
		mgw += current_nucleotide["MGW"]
		roll += current_nucleotide["Roll"]
		prot += current_nucleotide["ProT"]
		helt += current_nucleotide["HelT"]

	helt_modified = [helt[1]]
	for i in range(2, len(helt), 2):
		helt_modified.append(get_avg(helt[i:i+2]))
	roll_modified = [roll[1]]
	for i in range(2, len(roll), 2):
		roll_modified.append(get_avg(roll[i:i+2]))

	return {"MGW":mgw[1:], "ProT": prot[1:], "Roll": roll_modified, "HelT": helt_modified}

