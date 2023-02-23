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

def make_sample():
	cur.execute("DROP TABLE IF EXISTS sample;")
	cur.execute("CREATE TABLE sample ("
			"sample_hash VARCHAR(100) PRIMARY KEY, "
			"file_path VARCHAR(400), "
			"upload_date DATETIME, "
			"sample_description VARCHAR(100), "
			"metadata VARCHAR(100),"
			"chromatography_type TINYINT, "
			"original_file_name VARCHAR(400), "
			"original_XCMSnExp_path VARCHAR(400)"
		");"
	)
	print("Sample table created...")


def make_job():
	cur.execute("DROP TABLE IF EXISTS job;")
	cur.execute("CREATE TABLE job ("
		"job_id INT AUTO_INCREMENT PRIMARY KEY, "
		"job_status VARCHAR(40), " #CHANGE 40
		"start_time DATETIME, "
		"end_time DATETIME, " #NEW
		"job_name VARCHAR(100) "
		");"
	)
	print("Job table created...")

def make_parameter():
	cur.execute("DROP TABLE IF EXISTS parameter;")
	cur.execute("CREATE TABLE parameter ("
		"job_id INT PRIMARY KEY, "
		"sample_type DOUBLE, "
		"Peak_method TINYINT, "
    "Ref_method TINYINT, "
    "Align_method TINYINT, "
    "Group_method TINYINT, "
    "absMz DOUBLE, "
    "absRt DOUBLE, "
    "baseValue DOUBLE, "
    "binSize DOUBLE, "
    "bw DOUBLE, "
    "centerSample DOUBLE, "
    "checkBack DOUBLE, "
    "consecMissedLimit DOUBLE, "
    "criticalValue DOUBLE, "
    "distance DOUBLE, "
    "distFun DOUBLE, "
    "expandMz DOUBLE, "
    "expandRt DOUBLE, "
    "extendLengthMSW DOUBLE, "
    "extraPeaks DOUBLE, "
    "factorDiag DOUBLE, "
    "factorGap DOUBLE, "
    "family TINYINT, "
    "firstBaselineCheck DOUBLE, "
    "fitgauss BOOL, "
    "fixedMz DOUBLE, "
    "fixedRt DOUBLE, "
    "fwhm DOUBLE, "
    "gapExtend DOUBLE, "
    "gapInit DOUBLE, "
    "impute TINYINT, "
    "`index` BOOL, "
    "initPenalty DOUBLE, "
    "integrate DOUBLE, "
    "kNN DOUBLE, "
    "localAlignment BOOL, " 
    "max DOUBLE, "
    "maxFeatures DOUBLE, "
    "minFraction DOUBLE, "
    "minProp DOUBLE, "
    "minSamples DOUBLE, "
    "mzCenterFun TINYINT, "
    "mzdiff DOUBLE, "
    "mzVsRtBalance DOUBLE, "
    "ncol DOUBLE, "
    "noise DOUBLE, "
    "nrow DOUBLE, "
    "nValues DOUBLE, "
    "peakGroupsMatrix DOUBLE, "
    "max_peakwidth DOUBLE, "
    "min_peakwidth DOUBLE, "
    "ppm DOUBLE, "
    "prefilter DOUBLE, "
    "value_of_prefilter DOUBLE, " 
    "response DOUBLE, "
    "roiList DOUBLE, "
    "roiScales DOUBLE, "
    "sampleGroups DOUBLE, "
    "sigma DOUBLE, "
    "smooth TINYINT, "
    "snthresh DOUBLE, "
    "span DOUBLE, "
    "steps DOUBLE, "
    "subset DOUBLE, "
    "subsetAdjust TINYINT, "
    "threshold DOUBLE, "
    "unions DOUBLE, "
    "value DOUBLE, "
    "verboseColumns BOOL, "
    "withWave BOOL, "
    "rtrmin DOUBLE, "
    "rtrmax DOUBLE, "
		"FOREIGN KEY (job_id) REFERENCES job(job_id) ON DELETE CASCADE"
		");"
	)
	print("Parameter table created...")
	
def make_processed_sample():
	cur.execute("DROP TABLE IF EXISTS processed_sample;")
	cur.execute("CREATE TABLE processed_sample ("
      "job_id INT PRIMARY KEY, "
      "file_path_rda VARCHAR(100), "
      "file_path_peaks VARCHAR(100), "
      "FOREIGN KEY (job_id) REFERENCES job(job_id) ON DELETE CASCADE"
    ");"
	)
	print("Processed sample table created...")

def make_sample_job():
  cur.execute("DROP TABLE IF EXISTS sample_job;")
  cur.execute("CREATE TABLE sample_job ("
      "job_id INT, "
      "sample_hash VARCHAR(100), "
      "sample_number INT, "
      "PRIMARY KEY(job_id, sample_hash), "
      "FOREIGN KEY(job_id) REFERENCES job(job_id) ON DELETE CASCADE, "
      "FOREIGN KEY(sample_hash) REFERENCES sample(sample_hash) ON DELETE CASCADE"
    ");"
  )
  print("Sample job table created...")

def make_peak():
  cur.execute("DROP TABLE IF EXISTS sample_job;")
  cur.exeute("CREATE TABLE peak"
    
  )


cur.execute("SET GLOBAL local_infile=1;")
cur.execute("SET FOREIGN_KEY_CHECKS=0;")
# import os
# dir_path = os.path.dirname(os.path.realpath(__file__))
# bash_remove = f"rm {dir_path}/massascans/*.mzXML"
# os.system(bash_remove)
# make_sample()
make_parameter()
make_job()
make_processed_sample()
make_sample_job()
cur.execute("SET FOREIGN_KEY_CHECKS=1;")
