import mysql.connector

#EMPTY IN GIT
db_host = ""
db_user = ""
db_pwd = ""
db_name = ""

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
	cur.execute("DROP TABLE IF EXISTS samples")
	cur.execute("CREATE TABLE samples ("
			"sample_hash VARCHAR(100) PRIMARY KEY, "
			"file_path VARCHAR(100), "
			"upload_date DATETIME, "
			"sample_name VARCHAR(100)"
		");"
	)


def make_jobs():
	cur.execute("DROP TABLE IF EXISTS jobs")
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
	cur.execute("DROP TABLE IF EXISTS parameters; ")
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
		"verbose_cols_lcms TINYINT, "
		"polarity_lcms FLOAT, "
		"perc_fwhm_lcms FLOAT, "
		"mz_abs_iso_lcms TINYINT, "
		"max_charge_lcms FLOAT, "
		"max_iso_lcms FLOAT, "
		"corr_eic_th_lcms FLOAT, "
		"mz_abs_add_lcms FLOAT, "
		"rmconts_lcms TINYINT"
	");"
	)
make_samples()
make_parameters()
make_jobs()
