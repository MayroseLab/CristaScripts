import os, shutil
from definitions import *
from download_and_index_assembly import is_time_older_in_months


def delete_unused_genome_assemblies():
	for genome_dir in os.listdir(GENOMES_PATH):
		assembly_dirpath = GENOMES_PATH + genome_dir + SEP
		if (not os.path.exists(assembly_dirpath + "/protected")) and (not os.path.exists(assembly_dirpath + "/currently_processing"))\
				and (os.path.exists(assembly_dirpath + "/timestamp.txt") and is_time_older_in_months(assembly_dirpath, months=6)
					 or not os.path.exists(assembly_dirpath + "/timestamp.txt")):
			try:
				shutil.rmtree(assembly_dirpath)
				print("Deleted", assembly_dirpath, "succesfully")
			except:
				pass
		else:
			print(assembly_dirpath, "was not deleted.")
	print("Assemblies that were not deleted comply to one of the following:\n1) they are less than 6 months old\n2) they are protected from deletion\n3)there was an exeption while attempting to delete. Possibly, they were just downloaded.")


if __name__ == '__main__':
	delete_unused_genome_assemblies()