import os, re
from definitions import *
from download_and_index_assembly import write_current_timestamp_for_assembly

if __name__ == '__main__':
	res_dirs = os.listdir("/bioseq/data/results/crista/")
	lst = []
	for dir in res_dirs:
		if os.path.exists("/bioseq/data/results/crista/" + dir + "/qsub.sh"):
			with open("/bioseq/data/results/crista/" + dir + "/qsub.sh") as fpr:
				qsub_cmd = fpr.readlines()[-1]

			genome_name = re.search("(?<=\-g ).*? ", qsub_cmd)
			if genome_name:
				lst.append(genome_name.group().strip())
	for genome in set(lst):
		write_current_timestamp_for_assembly(GENOMES_PATH + genome + "/")
