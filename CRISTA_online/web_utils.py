import csv
import string
import subprocess, os
import pandas as pd
import regex as re


def get_email_recipient(results_path):
	if os.path.exists(results_path + "/user_email.txt"):

		with open(results_path + "/user_email.txt") as usr_fpr:
			recipient = usr_fpr.read().strip()
		return recipient


def send_email(recipient, run_id, script, logger=None):
	"""
	:param recipient:
	:param first: if True> sendFirstEmail.pl, else sendLastEmail.pl
	:return:
	"""
	if recipient is None:
		logger.info("email address not defined")
		return

	sendmail_cmd = "perl {0} -toEmail {1} -id {2}".format(script, recipient, run_id)
	proc = subprocess.Popen([sendmail_cmd], stdout=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	#os.system(sendmail_cmd)
	if logger is not None:
		logger.info("email sent")


def send_first_email(recipient, run_id, logger=None):
	script = "/bioseq/crista/CRISTA_online/sendFirstEmail.pl"
	send_email(recipient, run_id, script, logger=logger)


def send_last_email(recipient, run_id, logger=None):
	script = "/bioseq/crista/CRISTA_online/sendLastEmail.pl"
	send_email(recipient, run_id, script, logger=logger)


def send_error_email(recipient, run_id, logger=None):
	script = "/bioseq/crista/CRISTA_online/sendErrorEmail.pl"
	send_email(recipient, run_id, script, logger=logger)


def replace_line_in_html(results_path, old_line, new_line=""):
	html_path = results_path + "/output.php"

	with open(html_path) as fpr:
		html_content = fpr.read()
	html_content = re.sub(old_line, new_line, html_content)

	with open(html_path, "w") as fpw:
		fpw.write(html_content)


def add_to_results_html(results_path, add_text, header_message="Completed successfully"):
	html_path = results_path + "/output.php"

	with open(html_path) as fpr:
		html_content = fpr.read()

	html_content += "\n<p><font face=Verdana size=4>\n" \
					"<br>-------------------------------------------------------------------------------------------------\n" \
					"<br><br>\n" \
					"Results:<br>\n" \
					"</font></p>\n\n" \
					"<p><font face=Verdana size=2>\n"

	html_content += add_text


	if header_message != "Error occurred":
		html_content += "\n\n<b>Thank you for using CRISTA!</b><br><br><br>"
	html_content = re.sub("(?<=Your job status is.*?)Running", header_message, html_content)

	with open(html_path, "w") as fpw:
		fpw.write(html_content)
		fpw.flush()


def make_new_results_page(results_path, run_number, running_mode, running_params):
	with open('/bioseq/crista/CRISTA_online/results_html_template.html', "r") as fpr:
		content = fpr.read()
	content = re.sub("\{0\}",str(30), content)
	content = re.sub("\{1\}",run_number, content)
	content = re.sub("\{2\}",running_mode, content)
	content = re.sub("\{3\}",running_params, content)

	with open(results_path + "output.php", "w") as fpw:
		fpw.write(content)
		fpw.flush()

	return results_path + "output.php"


def stop_refreshing(results_path):
	html_path = results_path + "/output.php"

	with open(html_path) as fpr:
		html_content = fpr.read()
	html_content = re.sub("<HEAD> <META HTTP-EQUIV=\"REFRESH\" CONTENT=.*?> </HEAD>", "\n", html_content)

	with open(html_path, "w") as fpw:
		fpw.write(html_content)
		fpw.flush()


def load_error_page(results_path, error_message=""):
	error_message = "Error has occurred during processing of your job. <br><font color=\"red\">"+ \
					error_message +\
					"</font><br>Please contact us for assistance."
	add_to_results_html(results_path, error_message, "Error occurred")


def output_results_to_html(results_path, results_df):
	assert isinstance(results_df, pd.DataFrame)
	html_table_outpage = results_path + "/outtable.php"
	results_df.to_html(html_table_outpage, justify="left", float_format=lambda x: '%10.2f' % x)
	return html_table_outpage


def add_table_to_html_results_page(results_path, results_df, create_offtargets_link=False):

	if create_offtargets_link:
		results_df["off-targets search"] = "find off-targets!"
		results_df.rename(columns={"aligned sgRNA": "sgRNA", "aligned site": "DNA site"}, inplace=True)
	table_fpr = open(output_results_to_html(results_path, results_df))
	link_str = "<input id=\"{0}\" type=\"checkbox\" name=\"{0}\" value=\"{0}\" onchange='javascript:addOrRemoveFromList(this)'>"

	html_path = results_path + "/output.php"
	html_fpa = open(html_path, "a")

	for line in table_fpr:
		if create_offtargets_link and "NGG" in line and "<td>" in line:
			sgseq = re.search("(?<=\<td\>)[AGCT]{20}", line).group()
			content = link_str.format(sgseq)
			#line = re.sub(sgseq + "NGG", content, line)
		if create_offtargets_link and "find off-targets!" in line and "<td>" in line:
			line = re.sub("find off-targets!", content, line)
		html_fpa.write(line)
	html_fpa.write("\n</blockquote>\n</html>")

	html_fpa.close()
	table_fpr.close()

	return html_path


def is_csv(infile):
	try:
		with open(infile, newline='') as csvfile:
			csv.Sniffer().sniff(csvfile.read(4096))
			return True
	except:
		# Could not get a csv dialect -> probably not a csv.
		return False


if __name__ == '__main__':
	make_new_results_page("", "", "", "")
	add_table_to_html_results_page("/bioseq/data/results/crista/1487777465/", pd.read_csv("/bioseq/data/results/crista/1487777465/GAGTCCTAGCAGAAGAAGAANGG_CRISTA_offtargets_scores.csv"))