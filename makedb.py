#/usr/bin/env python3
import mysql.connector
db_host = "localhost"
with open(".dbpw", "r") as f:
	content = f.read().strip().split(",")
	db_user, db_pwd = content[0], content[1]
db_name = "boloa"

try:
	mydb = mysql.connector.connect(
		host=db_host,
		user=db_user,
		password=db_pwd,
		database=db_name
	)
	cur = mydb.cursor()
except:
	mydb = mysql.connector.connect(
                host=db_host,
                user=db_user,
                password=db_pwd
	)
	cur = mydb.cursor()
	cur.execute("CREATE DATABASE boloa")

def make_samples():
	cur.execute("DROP TABLE IF EXISTS samples;")
	cur.execute("CREATE TABLE samples ("
			"sample_hash VARCHAR(100) PRIMARY KEY, "
			"file_path VARCHAR(400), "
			"upload_date DATETIME, "
			"sample_name VARCHAR(100), "
			"chromatography_type TINYINT, "
			"original_file_name VARCHAR(400)"
		");"
	)


def make_jobs():
	cur.execute("DROP TABLE IF EXISTS jobs;")
	cur.execute("CREATE TABLE jobs ("
		"job_id INT AUTO_INCREMENT PRIMARY KEY, "
		"sample_hash VARCHAR(100), "
		"param_id INT, "
		"job_status VARCHAR(10), "
		"FOREIGN KEY (sample_hash)"
		"	REFERENCES samples(sample_hash),"
		"FOREIGN KEY (param_id)"
                "       REFERENCES parameters(param_id)"

		");"
	)

def make_parameters():
	cur.execute("DROP TABLE IF EXISTS parameters;")
	cur.execute("CREATE TABLE parameters ("
		"param_id INT AUTO_INCREMENT PRIMARY KEY, "
		"sample_type FLOAT, "
		"min_peakwidth FLOAT, "
		"max_peakwidth FLOAT, "
		"mzdiff FLOAT, "
		"snthresh FLOAT, "
		"bw FLOAT, "
		"peak_method FLOAT, "
		"ppm FLOAT, "
		"noise FLOAT, "
		"prefilter FLOAT, "
		"value_of_prefilter FLOAT, "
		"minfraction FLOAT, "
		"minsamples FLOAT, "
		"maxfeatures FLOAT, "
		"fitgauss TINYINT, "
		"mzcenterfun FLOAT, "
		"integrate FLOAT, "
		"extra FLOAT, "
		"span FLOAT, "
		"smooth FLOAT, "
		"family FLOAT, "
		"verbose_cols TINYINT, "
		"polarity FLOAT, "
		"perc_fwhm FLOAT, "
		"mz_abs_iso FLOAT, "
		"max_charge FLOAT, "
		"max_iso FLOAT, "
		"corr_eic_th FLOAT, "
		"mz_abs_add FLOAT, "
		"rmconts TINYINT"
		");"
	)
cur.execute("SET GLOBAL local_infile=1;")
cur.execute("SET FOREIGN_KEY_CHECKS=0;")
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
bash_remove = f"rm {dir_path}/massascans/*.mzXML"
os.system(bash_remove)
make_samples()
make_parameters()
make_jobs()
cur.execute("SET FOREIGN_KEY_CHECKS=1;")
