#/usr/bin/env python3
'''
This script is responsible for the creation of the application database.
It is necessary to run this script's functions at least once in order to make
the application functional. Comment out tables which do not require to be reset.
Each time this script is run, the contents of tables are deleted.
'''

import mysql.connector
import json
import numpy as np
import subprocess
import os
db_host = "localhost"
with open(".dbpw", "r") as f:
	content = f.read().strip().split(",")
	db_user, db_pwd = content[0], content[1]
db_name = "boloa"

try:
  '''
  Try to connect to the database.
  '''
	mydb = mysql.connector.connect(
		host=db_host,
		user=db_user,
		password=db_pwd,
		database=db_name
	)
	cur = mydb.cursor()
except:
  '''
  If the database does not exist; create the database.
  '''
	mydb = mysql.connector.connect(
                host=db_host,
                user=db_user,
                password=db_pwd
	)
	cur = mydb.cursor()
	cur.execute("CREATE DATABASE boloa")

def make_sample():
  '''
  This function is responsible to create the table which will contain
  information about the samples which have been uploaded.
  '''
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
  '''
  This function is responsible to create the table which will contain
  information about the jobs which have been created by various users.
  '''
	cur.execute("DROP TABLE IF EXISTS job;")
	cur.execute("CREATE TABLE job ("
		"job_id INT AUTO_INCREMENT PRIMARY KEY, "
		"job_status VARCHAR(40), "
		"start_time DATETIME, "
		"end_time DATETIME, "
		"job_name VARCHAR(100) "
		");"
	)
	print("Job table created...")

def make_parameter():
  '''
  This function is responsible store parameters in relation to the job table.
  '''
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
    "fitgauss BOOLEAN, "
    "fixedMz DOUBLE, "
    "fixedRt DOUBLE, "
    "fwhm DOUBLE, "
    "gapExtend DOUBLE, "
    "gapInit DOUBLE, "
    "impute TINYINT, "
    "`index` BOOLEAN, "
    "initPenalty DOUBLE, "
    "integrate DOUBLE, "
    "kNN DOUBLE, "
    "localAlignment BOOLEAN, " 
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
    "verboseColumns BOOLEAN, "
    "withWave BOOLEAN, "
    "rtrmin DOUBLE, "
    "rtrmax DOUBLE, "
    "rsd_threshold DOUBLE, "
    "simthresh DOUBLE, "
		"FOREIGN KEY (job_id) REFERENCES job(job_id) ON DELETE CASCADE"
		");"
	)
	print("Parameter table created...")
	
def make_processed_sample():
  '''
  This function is responsible to store the paths to the processed samples.
  '''
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
  '''
  This function is responsible to create a junction table between the job and
  sample table.
  '''
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
  '''
  This creates the table that will contain the information about peaks which
  have been detected and processed by the application.
  '''
  cur.execute("DROP TABLE IF EXISTS peak;")
  cur.execute(
    "CREATE TABLE peak ("
    "job_id INT, "
    "peak_id VARCHAR(100), "
    "sample_hash  VARCHAR(100), "
    "peak_area DOUBLE, "
    "mol_id_coeffsim INT, "
    "mol_id_modcosinesim INT, "
    "modcosinesim DOUBLE, "
    "coeff_diff DOUBLE, "
    "coeff DOUBLE, "
    "mz DOUBLE, "
    "mzmin DOUBLE, "
    "mzmax DOUBLE, "
    "rt DOUBLE, "
    "rtmin DOUBLE, "
    "rtmax DOUBLE, "
    "apex_tic DOUBLE, "
    "splash VARCHAR(100), "
    "spectrum MEDIUMTEXT, "
    "PRIMARY KEY(job_id, peak_id, sample_hash), "
    "FOREIGN KEY (sample_hash) REFERENCES sample(sample_hash) ON DELETE CASCADE, "
    "FOREIGN KEY(job_id) REFERENCES job(job_id) ON DELETE CASCADE, "
    "FOREIGN KEY(mol_id_coeffsim) REFERENCES mol(mol_id) ON DELETE CASCADE, "
    "FOREIGN KEY(mol_id_modcosinesim) REFERENCES mol(mol_id) ON DELETE CASCADE"
    ");"
  )
  print("Peak table created...")

def make_mol():
  '''
  This function creates the molecule table. This table will be filled when
  this function called. In order to fill this table, the path to GCMS and LCMS
  MoNA files is needed. These files are subsequently parsed and relevant
  information is stored inside the molecule table.
  '''
  cur.execute("DROP TABLE IF EXISTS mol;")
  cur.execute(
    "CREATE TABLE mol ("
    "mol_id INT AUTO_INCREMENT PRIMARY KEY, "
    "splash VARCHAR(100), "
    "mol_name VARCHAR(1000), "
    "mass DOUBLE, "
    "pubid VARCHAR(1000), "
    "type VARCHAR(4), "
    "coeff DOUBLE, "
    "spectrum MEDIUMTEXT"
    ");"  
  )
  print("Molecule table created... Now filling mol, please be patient...")
  for chrom_method in ["GC", "LC"]: # Needed for both GC and LCMS files
    f = open(f'annot_data/MoNA-export-{chrom_method}-MS_Spectra.json')
    data = json.load(f)
    f.close()
    # For every compound inside the file: retrieve information
    for i in range(len(data)):
      splash = data[i]["splash"]["splash"]
      try:
          mol_name = data[i]["compound"][0]["names"][0]["name"]
          mol_name = mol_name.replace("'", r"\'")
      except:
          mol_name = "Unknown Compound"
      for md in data[i]["compound"][0]["metaData"]:
          if md["name"] == "total exact mass":
              mass = md["value"]
          if "pubchem" in md["name"]:
            pubid = md["value"]
      try:
        pubid
      except:
        pubid = "N/A"
      spectrum = data[i]["spectrum"]
      # Parse the spectrum in order to calculate a coefficient between
      # m/z and intensity. Might be used in a new annotation method.
      s = spectrum.split(" ")
      s = [bn.split(':') for bn in s]
      x = np.array([float(ms[0]) for ms in s])
      y = np.array([float(intens[1]) for intens in s])
      coeff = np.dot(x, y)/(np.linalg.norm(x) * np.linalg.norm(y))
      if coeff == "NA":
        coeff = 0
      # Store all the relevant variables inside the database.
      cur.execute(f"INSERT INTO mol (splash, mol_name, mass, pubid, type, coeff, spectrum) VALUES "
          f"('{splash}', '{mol_name}', '{mass}', '{pubid}', '{chrom_method}', '{coeff}', '{spectrum}');")
      del pubid
  print("Molecule table filled...")

cur.execute("SET GLOBAL local_infile=1;")
cur.execute("SET FOREIGN_KEY_CHECKS=0;")
# dir_path = os.path.dirname(os.path.realpath(__file__))
# bash_remove = f"rm {dir_path}/massascans/*.mzXML"
# os.system(bash_remove)

### Comment out tables where drops are unwanted.
# make_sample() # Warning: Will remove samples from database, requires manual reuploading
make_parameter()
make_job()
make_processed_sample()
make_sample_job()
# make_mol() # Warning: The creation and subsequent filling of this table takes a long time
make_peak()
cur.execute("SET FOREIGN_KEY_CHECKS=1;")
mydb.commit()
