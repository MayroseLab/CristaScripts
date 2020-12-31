
import pickle
import sklearn
import shutil, os

location = "/bioseq/crista"

if __name__ == '__main__':
	for outer_dir in ["rf", "rf_noflanking", "rf_nogenomic"]:
		for ix in range(100):
			source = "/".join([location, outer_dir, str(ix)])
			dest = "/".join([location, outer_dir + "_0.19.1", str(ix)])
			os.makedirs(dest)

			with open(source + "/RFpredictor.pkl", "rb") as pklr:
				x = pickle.load(pklr)
			with open(dest + "/RFpredictor.pkl", "wb") as pklw:
				pickle.dump(x, pklw)
			shutil.copyfile(source + "/RF training.txt", dest + "/RF training.txt")

