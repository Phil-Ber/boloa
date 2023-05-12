#!/usr/bin/env Rscript

# This loads the variables containing the database username and password.
load('.hidden.RData')

# Below all libraries which are required to run the script are imported.
library(shiny)
library(shinyBS)
library(shinyjs)
library(bslib)
library(digest)
library(DBI)
library(RMySQL)
#library(tidyverse)
library(dplyr)
library(DT)
library(pacman)
library(devtools)
library(plot3D)
library(magrittr)
library(mzR)
library(rgl)
library(reticulate)
library(promises)
library(future)
library(waiter)
library(Spectra)
#library(deldir)
#library(xcms)
library(heatmaply)
pacman::p_load(c("Rserve", "ellipse", "scatterplot3d", "Cairo", "randomForest", "caTools", "e1071", "som", "impute", "pcaMethods", "RJSONIO", "ROCR", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "pheatmap", "sva", "Rcpp", "pROC", "data.table", "limma", "car", "fitdistrplus", "lars", "Hmisc", "magrittr", "methods", "xtable", "pls", "caret", "lattice", "igraph", "gplots", "KEGGgraph", "reshape", "RColorBrewer", "tibble", "siggenes", "plotly", "xcms", "CAMERA", "fgsea", "MSnbase", "BiocParallel", "multtest", "RBGL", "edgeR", "fgsea", "crmn", "progress", "qs", "glasso"), character.only = TRUE)
# install_github("berlinguyinca/spectra-hash", subdir="splashR")
library(splashR)

# This function simply returns columns from a selected table in the database.
get_query <- function(query){
	sqlconn <- dbConnect(
		drv = RMySQL::MySQL(),
		dbname='boloa',
		host="127.0.0.1",
		port=3306,
		user=db_usr,
		password=db_pwd
	)
	on.exit(dbDisconnect(sqlconn))
	dbGetQuery(sqlconn, query)
}

# This function allows R to send a query to the database. Used for more complex
# goals like joins.
send_query <- function(query){
  sqlconn <- dbConnect(
    drv = RMySQL::MySQL(),
    dbname='boloa',
    host="127.0.0.1",
    port=3306,
    user=db_usr,
    password=db_pwd
  )
  on.exit(dbDisconnect(sqlconn))
  dbSendQuery(sqlconn, query)
}

# This function allows for the simple insertion of data into the database.
# The function allows for the input of a data.frame object. The object is then
# inserted into the database under the table which has been inputted.
# Warning: The colnames() in the data.frame have to match 100% with those in
# the desired table of the database.
insert_query <- function(tb_name, data){
	sqlconn <- dbConnect(
		drv = RMySQL::MySQL(),
		dbname='boloa',
		host="127.0.0.1",
		port=3306,
		user=db_usr,
		password=db_pwd
	)
	on.exit(dbDisconnect(sqlconn))
	dbWriteTable(sqlconn, tb_name, data, append = TRUE, row.names = FALSE)
}

# Multisession allows for parallel processing. Necessary for the future package.
plan(multisession)

# Allows for larger file uploads. Necessary for mass spectrometry data.
options(shiny.maxRequestSize=10000*1024^2)

# Shiny configuration.
hostip <- "145.97.18.149"
portnr <- 7123
# portnr <- 9705


# This variable contains all the parameters which can possibly be used in peak
# detection. Adding new methods require adding the parameters to this array.
# Warning: Check first in this array if a parameter is already in use.
# If an overlap occurs you don't need to add a new parameter to the respective
# table in the "makedb.py" script.
all_params <- c("Peak_method", "Ref_method", "Align_method", "Group_method", "absMz", "absRt", "baseValue", "binSize", "bw", "centerSample", "checkBack", "consecMissedLimit", "criticalValue", "distance", "distFun", "expandMz", "expandRt", "extendLengthMSW", "extraPeaks", "factorDiag", "factorGap", "family", "firstBaselineCheck", "fitgauss", "fixedMz", "fixedRt", "fwhm", "gapExtend", "gapInit", "impute", "index", "initPenalty", "integrate", "kNN", "localAlignment", "max", "maxFeatures", "minFraction", "minProp", "minSamples", "mzCenterFun", "mzdiff", "mzVsRtBalance", "ncol", "noise", "nrow", "nValues", "peakGroupsMatrix", "max_peakwidth", "min_peakwidth", "ppm", "prefilter", "value_of_prefilter", "response", "roiList", "roiScales", "sampleGroups", "sigma", "smooth", "snthresh", "span", "steps", "subset", "subsetAdjust", "threshold", "unions", "value", "verboseColumns", "withWave", "rtrange", "rsd_threshold", "simthresh")
lcms_only <- tail(all_params, n=9)


# If an error occurs regarding the database, send this query.
# send_query(stringr::str_glue("SET GLOBAL local_infile=1;"))

ui <- fillPage(
  # This object contains all the hard coded front end elements.
  padding = NULL,
  autoWaiter(color = "white", html = spin_3()),
  
	useShinyjs(),
	theme = bs_theme(version = 4, bootswatch = "spacelab"),
	mainPanel(style = "padding: 0px;",
		titlePanel("Boloa"),
		tabsetPanel(id = "tabSwitch",
			tabPanel("Select data",  icon = icon("arrow-right"), style='width:100vw;height:90vh;overflow-y: scroll;', value = "p1",
				fluidRow(style='max-width:100%;padding:10px;',
					column(3,
						br(),
						textInput("sample_description", label = 'Sample description (e.g. "Rhino, feces" or "human, cancer, blood"):'),
						h4("Mass Spectrometry data:"),
						fileInput("msdata", label = "MS-data", multiple = TRUE, accept = c(".mzXML", ".raw")),#, ".CDF")),
						#DT::dataTableOutput("fileOverview"),
						uiOutput("fileOverview", style='margin-top:20px; border-top:1px solid #dfd7ca; margin-bottom:20px; border-bottom:1px solid #dfd7ca;overflow-y: scroll;'),
						actionButton("dataupl", label = "Submit", width = 180),
						verbatimTextOutput("upl_completed"),
						style='margin-bottom:30px;border-right:1px solid #dfd7ca;; padding: 10px;'
					),
					column(9,
						DT::dataTableOutput("uploaded_samples")
					)
				)
			),
			tabPanel("Process data", icon = icon("arrow-right"), style='width:100vw;height:90vh;overflow-y: scroll;', value = "p2",
				fluidRow(style='max-width:100%;padding:10px;',
					column(7,
  		      radioButtons("datatype", label = "Select correct chromatography type in order to properly display data:",
  		        choices = list("LC-MS" = 1, "GC-MS" = 2)),
						DT::dataTableOutput("selected_samples")
					),
					column(5,
						#plotOutput("ticplot"),
						plotlyOutput("ticplot")
					)
				),
				fluidRow(style='max-width:100%;padding:10px;',
					column(7,
						actionButton("preview", label = "\n Preview spectra")
						#actionButton("prevtic", label = "\n Preview TIC")
					),
					column(5,
						#uiOutput("mzrange"),
						uiOutput("previewplot")
					)
				),
				fluidRow(style='max-width:100%;padding:20px;margin-top:20px; border-top:1px solid #dfd7ca;',
				         h4("Parameters:")),
				fluidRow(style='max-width:100%;padding:10px;',
					column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
					  uiOutput("rtRangeSlider"),
						selectInput("preset", label = "Parameter presets", choices = list("None" = 0, "LC-MS" = 1, "GC-MS" = 2, "Automatic (centWave only)" = 3)),
						textInput("job_name", label = "Job name"),
						actionButton("submitJob", "Submit job", width = 180, icon=icon("play")),
						tableOutput("jobval")
					),
					column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
						h5("1. Detection"),
						selectInput("Peak_method", label = "Peak detection method", choices = list("centWave" = 0, "Massifquant" = 1, "MatchedFilter" = 2, "relsqdiff" = 3), selected = 3),
						uiOutput("peak_parameters")
					),
					column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
						h5("2. Refinement"),
						selectInput("Ref_method", label = "Peak refinement method", choices = list("MergeNeighboringPeaks" = 0, "FilterIntensity" = 1)),
						uiOutput("refinement_parameters")
					),
					column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
						h5("3. Alignment"),
						selectInput("Align_method", label = "Peak alignment method", choices = list("Obiwarp" = 0, "PeakGroups" = 1)),
						uiOutput("alignment_parameters"),
					),
					column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
					       h5("4. Grouping"),
					       selectInput("Group_method", label = "Peak grouping method", choices = list("PeakDensity" = 0, "MzClust" = 1, "NearestPeaks" = 2)),
					       uiOutput("grouping_parameters")
					),
					column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
					       h5("5. Filling"),
					       numericInput("fixedMz", label = "fixedMz", 0),
					       numericInput("fixedRt", label = "fixedRt", 0)
					       )
				)
			),
			tabPanel("Jobs", icon = icon("arrow-right"), style='width:100vw;height:90vh;overflow-y: scroll;', value = "p3",
			  fluidRow(style='max-width:100%;padding:10px;',
			    column(12,
  			    h4("Select a job to analyse contents:"),
  			    DT::dataTableOutput("jobs"),
  			    actionButton("updateJobsTable", label = "\n Update table", icon = icon("arrows-rotate")),
  			    actionButton("analyse", label = "\n Analyse job", icon=icon("play"))
			   )
			  )
			),
			tabPanel("Analysis",  icon = icon("arrow-right"), style='width:100vw;height:90vh;overflow-y: scroll;', value = "p4",
			         uiOutput("job_analysis")
		)
	)
)
)


# This function contains the back end for the web application.
server <- function(input, output, session) {
  nr_files <<- 0
  
  # The objects below contain peak processing methods. If the changing of
  # default parameters is desired: change them here.
  default_relsqdiff <- tagList(
    numericInput("noise", label = "noise", 10000),
    numericInput("rsd_threshold", "threshold", 5),
    numericInput("simthresh", "simthresh", 900)
  )
  
  default_cent <- tagList(
    numericInput("ppm", label = "ppm", 40.24726),
    numericInput("min_peakwidth", label = "min_peakwidth", 7.179200),
    numericInput("max_peakwidth", label = "max_peakwidth", 56.20481),
    numericInput("snthresh", label = "snthresh", 49.03500),
    numericInput("prefilter", label = "prefilter", 2),
    numericInput("value_of_prefilter", label = "value_of_prefilter", 113.6439),
    selectInput("mzCenterFun", label = "mzCenterFun", choices = list("wMean" = 0, "mean" = 1, "apex" = 2, "wMeanApex3" = 3, "meanApex3" = 4), selected = 3),
    numericInput("integrate", label = "integrate", 1),
    numericInput("mzdiff", label = "mzdiff", 0.033627882),
    checkboxInput("fitgauss", label = "fitgauss", 0),
    numericInput("noise", label = "noise", 5046.248),
    checkboxInput("verboseColumns", label = "verboseColumns", 0)
  )
  default_massif <- tagList(
    numericInput("ppm", label = "ppm", 74.77602),
    numericInput("min_peakwidth", label = "min_peakwidth", 3.727111),
    numericInput("max_peakwidth", label = "max_peakwidth", 29.69969),
    numericInput("snthresh", label = "snthresh", 175.7889),
    numericInput("prefilter", label = "prefilter", 15.51382),
    numericInput("value_of_prefilter", label = "value_of_prefilter", 3523.923),
    selectInput("mzCenterFun", label = "mzCenterFun", choices = list("wMean" = 0, "mean" = 1, "apex" = 2, "wMeanApex3" = 3, "meanApex3" = 4), selected = 0),
    numericInput("integrate", label = "integrate", 1),
    numericInput("mzdiff", label = "mzdiff", -1.305187),
    checkboxInput("fitgauss", label = "fitgauss", 1),
    numericInput("noise", label = "noise", 102536.0),
    checkboxInput("verboseColumns", label = "verboseColumns", 1),
    numericInput("criticalValue", label = "criticalValue", 1.125),
    numericInput("consecMissedLimit", label = "consecMissedLimit", 2),
    numericInput("unions", label = "unions", 1),
    numericInput("checkBack", label = "checkBack", 0),
    checkboxInput("withWave", label = "withWave", 0)
  )
  
  default_matchedfilter <- tagList(
    numericInput("binSize", label = "binSize", 0.1),
    selectInput("impute", label = "impute", choices = list("none" = 0, "lin" = 1, "linbase" = 2, "intlin" = 3, "imputeLinInterpol" = 4)),
    # numericInput("baseValue", label = "baseValue", 5),
    # numericInput("distance", label = "distance", 5),
    numericInput("fwhm", label = "fwhm", 30),
    numericInput("sigma", label = "sigma", 2.3548),
    numericInput("max", label = "max", 5),
    numericInput("snthresh", label = "snthresh", 10),
    numericInput("steps", label = "steps", 2),
    numericInput("mzdiff", label = "mzdiff", 0.8),
    checkboxInput("index", label = "index", 0)
  )
  
  default_mnp <- tagList(
    numericInput("expandRt", label = "expandRt", 0),
    numericInput("expandMz", label = "expandMz", 2),
    numericInput("minProp", label = "minProp", 0.75)
  )
  
  default_fi <- tagList(
    numericInput("threshold", label = "threshold", 0),
    # numericInput("nValues", label = "nValues", 5),
    # numericInput("value", label = "value", 5)
  )
  
  default_obiwarp <- tagList(
    #binSize = 1,
    numericInput("binSize", label = "binSize", 1),
    # numericInput("centerSample", label = "centerSample", 5),
    # numericInput("response", label = "response", 1),
    selectInput("distFun", label = "distFun", choices = list("cor" = 0, "cor_opt" = 1, "cov" = 2, "prd" = 3, "euc" = 4)),
    # numericInput("gapInit", label = "gapInit", 5),
    # numericInput("gapExtend", label = "gapExtend", 5),
    numericInput("factorDiag", label = "factorDiag", 2),
    numericInput("factorGap", label = "factorGap", 1),
    checkboxInput("localAlignment", label = "localAlignment", 0),
    numericInput("initPenalty", label = "initPenalty", 0),
    # numericInput("subset", label = "subset", 5),
    selectInput("subsetAdjust", label = "subsetAdjust", choices = list("average" = 0, "previous" = 1))
  )
  
  default_peakgroups <- tagList(
    numericInput("minFraction", label = "minFraction", 0.9),
    numericInput("extraPeaks", label = "extraPeaks", 1),
    selectInput("smooth", label = "smooth", choices = list("loess" = 0, "linear" = 1)),
    numericInput("span", label = "span", 0.2),
    selectInput("family", label = "family", choices = list("gaussian" = 0, "symmetric" = 1)),
    # numericInput("peakGroupsMatrix", label = "peakGroupsMatrix", 5),
    # numericInput("subset", label = "subset", 5),
    selectInput("subsetAdjust", label = "subsetAdjust", choices = list("average" = 0, "previous" = 1))
  )
  
  default_pd <- tagList(
    numericInput("bw", label = "bw", 30),
    numericInput("minFraction", label = "minFraction", 0.5),
    numericInput("minSamples", label = "minSamples", 1),
    numericInput("binSize", label = "binSize", 0.25),
    numericInput("maxFeatures", label = "maxFeatures", 50)
  )
  
  default_mzc <- tagList(
    numericInput("absMz", label = "absMz", 0),
    numericInput("minFraction", label = "minFraction", 0.5),
    numericInput("minSamples", label = "minSamples", 1)
  )
  
  default_np <- tagList(
    numericInput("mzVsRtBalance", label = "mzVsRtBalance", 10),
    numericInput("absMz", label = "absMz", 0.2),
    numericInput("absRt", label = "absRt", 15),
    numericInput("kNN", label = "kNN", 10)
  )
  
  # This is displayed in case the user proceeds to the analysis tab before
  # properly selecting a job to inspect.
  output$job_analysis <- renderUI(tagList(
    h5('Please select a finished job in the "Jobs" panel to analyse data'),
  ))
  
  # This function checks if the user changes the selection of a peak detection
  # method. This way alteration displays correct parameters for another method.
  observeEvent(input$Peak_method, {
    if (input$Peak_method == 0) {
      peak_options <- default_cent
    }
    else if (input$Peak_method == 1) {
      peak_options <- default_massif
    }
    else if (input$Peak_method == 2) {
      peak_options <- default_matchedfilter
    }
    else if (input$Peak_method == 3) {
      peak_options <- default_relsqdiff
    }
    output$peak_parameters <- renderUI(peak_options)
  })
  
  # This function checks if the user changes the selection of a peak refinement
  # method. This way alteration displays correct parameters for another method.
  observeEvent(input$Ref_method, {
    if (input$Ref_method == 0) {
      ref_options <- default_mnp
    }
    else if (input$Ref_method == 1) {
      ref_options <- default_fi
    }
    output$refinement_parameters <- renderUI(ref_options)
  })
  
  # This function checks if the user changes the selection of a peak alingment
  # method. This way alteration displays correct parameters for another method.
  observeEvent(input$Align_method, {
    if (input$Align_method == 0) {
      alig_options <- default_obiwarp
    }
    else if (input$Align_method == 1) {
      alig_options <- default_peakgroups
    }
    output$alignment_parameters <- renderUI(alig_options)
  })
  
  # This function checks if the user changes the selection of a peak grouping
  # method. This way alteration displays correct parameters for another method.
  observeEvent(input$Group_method, {
    if (input$Group_method == 0) {
      group_options <- default_pd
    }
    else if (input$Group_method == 1) {
      group_options <- default_mzc
    }
    else if (input$Group_method == 2) {
      group_options <- default_np
    }
    output$grouping_parameters <- renderUI(group_options)
  })
  
	# This allows for inputting metadata or labels for samples which the user
  # would like to upload.
	observeEvent(input$msdata, {
		output$fileOverview <- renderUI({
			file <- input$msdata
			ext <- tools::file_ext(file$datapath)
			req(file)
			# Checks if the filetype is correct.
			validate(need(ext %in% c("raw", "mzXML", "CDF"), "Please upload a .raw, .CDF or .mzXML file!"))
			nr_files <<- length(file[,3])
			# Returns text input boxes with a quantity equal to the number of samples
			# which have been uploaded. IDs of the boxes are generated with the format
			# "meta(n)" where n is the index of the sample.
			lapply(1:nr_files, function(i) {
			  if (i == 1) {
			    p('Please edit file metadata (REQUIRED!)')
			  }
			  return(textInput(inputId = paste("meta", i, sep = ""), label = paste(i, ": ", file[i, 1], sep = ""), placeholder = paste("Group sample", i)))
			})
		#}#, server = FALSE)
  	})
  })
  
	# Reactively updates content when a switch of a tab or the upload of a sample
	# is detected. Necessary for datatables. Especially the sample table.
	updateEvent <- reactive({
		list(input$tabSwitch, input$dataupl)
	})
	
	# Responsible for the creation of the sample table in the sample selection
	# screen.
	observeEvent(updateEvent, {
		query <- stringr::str_glue("SELECT * FROM sample ORDER BY upload_date DESC;")
		sample_table_content <<- get_query(query)
		sample_table_content$chromatography_type[sample_table_content$chromatography_type == 1] <- "Liquid"
		sample_table_content$chromatography_type[sample_table_content$chromatography_type == 2] <- "Gas"
		output$uploaded_samples <- DT::renderDataTable({
			DT::datatable(sample_table_content[, c(7, 3, 5, 4, 6)])
		}, server = FALSE)
		query <- stringr::str_glue("SELECT * FROM job ORDER BY start_time DESC;")
		job_table_content <<- get_query(query)
		output$jobs <- DT::renderDataTable({
		  DT::datatable(job_table_content[, c(1, 5, 3, 4, 2)], selection = 'single',
		                rownames= FALSE)
		}, server = FALSE)
	})
	
	# Responsible for the creation of the running and/or finished jobs table.
	observeEvent(input$updateJobsTable, {
	  query <- stringr::str_glue("SELECT * FROM job ORDER BY start_time DESC;")
	  job_table_content <<- get_query(query)
	  # if (length(job_table_content[, 1]) != 0) {
	  #   job_table_content$end_time[is.null(sample_table_content$end_time)] <- "-"
	  # }
	  output$jobs <- DT::renderDataTable({
	    DT::datatable(job_table_content[, c(1, 5, 3, 4, 2)], selection = 'single',
	                  rownames= FALSE)
	  }, server = FALSE)
	})
	
	# When the user selects a sample in the "process data" tab, and presses the
	# "preview" button. This event creates a 3-dimensional plot of the mass
	# spectra in the sample. This is the only function that uses metaboAnalystR
	# and optiLCMS. Might remove since it's fun but not really functional.
	observeEvent(input$preview, {
	  # Creation of a progress bar
	  progress <- shiny::Progress$new()
	  on.exit(progress$close())
	  progress$set(message = "Preview", value = 0)
		preview_file <- selsamples[input$selected_samples_rows_selected, ]$sample_hash
		query <- stringr::str_glue("SELECT file_path FROM sample WHERE sample_hash = '", preview_file, "';")
		progress$inc(1/4, detail = "Opening MS-file")
		aa <- openMSfile(get_query(query)$file_path)
		on.exit(close(aa))
		aahead <- header(aa)
		ms1 <- which(aahead$msLevel == 1)
		rtsel <- aahead$retentionTime[ms1] >= min(aahead$retentionTime) & aahead$retentionTime[ms1] <= max(aahead$retentionTime)
		progress$inc(1/4, detail = "Determining MZ-ranges")
		mzlow <- min(aahead$lowMZ)
		mzhigh <- max(aahead$highMZ)
		mzres <- (mzhigh - mzlow) / 50
		progress$inc(1/4, detail = "Making MSmap object")
		M <- MSmap(aa, ms1[rtsel], lowMz=mzlow, highMz=mzhigh, resMz=mzres, aahead, zeroIsNA = TRUE)
		close(aa)
		# Render 3d plot
		output$previewplot <- renderUI({
			open3d()
			plot3D(M, rgl = TRUE)#xlim = input$chromrange, ylim = input$mzrange, xlab = "Retention time (seconds)", ylab = "M/Z", zlab = "Intensity")
			rglwidget()
		  
		})
		progress$close()
	})
	
	# Clear plots on tab switch DEPRECATED
	observeEvent(input$tabSwitch, {
		# output$previewplot <- renderUI({
		# 	return(NULL)
		# })
		# output$ticplot <- renderUI({
		# 	return(NULL)
		# })
		# output$mzrange <- renderUI({
		# 	return(NULL)
		# })
		# output$chromrange <- renderUI({
		# 	return(NULL)
		# })
	})
	
	# This event is responsible for handling data processing. The code below runs
	# whenever the user submits newly uploaded samples.
	observeEvent(input$dataupl, {
	  # Checks if the user submitted a correct overarching name to the samples.
	  if (nchar(input$sample_description) == 0){
	    output$upl_completed <- renderText({
	      "Please input a valid Sample description"
	    })
	    return()
	  }
	  # Checks if files are present for uploading.
		if (nr_files == 0) {
			output$upl_completed <- renderText({
				'Please upload a file.'
			})
			return()
		}
	  # Checks if all samples correctly have labels assigned to them.
	  invalid_metadata <- FALSE
	  for (i in 1:nr_files) {
	    # Removes all characters that are not spaces and or alphanumeric.
	    metastring <- gsub("[^[:alnum:][:space:]]", "", input[[paste("meta", i, sep = "")]])
	    if (metastring == "") {
	      output$upl_completed <- renderText({
	        paste('Please assign a label to all files consisting of only letters [Aa-Zz] and numbers [0-9]. \nError: ', input[[paste("meta", i, sep = "")]], sep = "")
	      }) 
	      invalid_metadata <- TRUE
	      return()
	    }
	  }
	  output$upl_completed <- renderText({
	    "Files uploading, please be patient...."
	  })
	  # Collects the labels or "file_tags"
	  file_tags <- c()
	  for (i in 1:nr_files){
	    file_tags <- c(file_tags, gsub("[^[:alnum:][:space:]]", "", input[[paste("meta", i, sep = "")]]))
	  }
		dir <- getwd()
		dir <- paste(dir, "/massascans", sep = "")
		# Creation of metadata object
		file <- input$msdata
		time <- format(Sys.time() + 60*60, "%Y-%m-%d %X")
		# This function is responsible for uploading files to the server and
		# appending information to the database.
		run_file_upload <- function(rfuFile, rfuTime, rfuDir, rfuSample_description, rfuDb_usr, rfuD_pwd, rfuSample_table_content, rfuFile_tags){
		  count <- 0
		  # Because of a nested future call, multisession and database functions
		  # are redifined.
		  plan(multisession)
		  # Iterates through all uploaded files.
		  for (fileloc in rfuFile$datapath) {
		    upload_file <- function(ufCount, ufFileloc, ufDir, ufTime, ufSample_description, ufFile, ufFile_tag, ufDb_usr, ufDb_pwd, ufSample_table_content){
		      
		      # Only insert and send queries are required.
		      insert_query <- function(tb_name, data){
		        sqlconn <- dbConnect(
		          drv = RMySQL::MySQL(),
		          dbname='boloa',
		          host="127.0.0.1",
		          port=3306,
		          user=ufDb_usr,
		          password=ufDb_pwd
		        )

		        on.exit(dbDisconnect(sqlconn))
		        dbWriteTable(sqlconn, tb_name, data, append = TRUE, row.names = FALSE)
		      }
		      send_query <- function(query){
		        sqlconn <- dbConnect(
		          drv = RMySQL::MySQL(),
		          dbname='boloa',
		          host="127.0.0.1",
		          port=3306,
		          user=ufDb_usr,
		          password=ufDb_pwd
		        )
		        on.exit(dbDisconnect(sqlconn))
		        dbSendQuery(sqlconn, query)
		      }
		      
		      # Copies the file from a temporary to a persistent directory.
		      file.copy(ufFileloc, ufDir)
		      # If the user uploads a .raw file, msconvert docker container is used
		      # to convert it to .mzXML.
		      if (grepl('.raw', ufFileloc, fixed=TRUE)) {
		        system(paste("docker run --rm -v ", ufDir, ":/massascans chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /massascans/", ufCount, ".raw", " --mzXML -o /massascans", sep = "")) #change .raw files to .mzXML
		        system(paste("rm ", ufDir, "/", ufCount, ".raw", sep = ""))
		        filetype <- ".mzXML"
		      }
		      else if (grepl('.mzXML', ufFileloc, fixed=TRUE)) {
		        filetype <- ".mzXML"
		      }
		      # .CDF is not supported; this is here in case it changes.
		      else if (grepl('.CDF', ufFileloc, fixed=TRUE)) {
		        filetype <- ".CDF"
		      }
		      # If somehow an unsupported file extension was uploaded: removal.
		      else {
		        system(paste("rm ", ufFileloc, sep = ""))
		      }
		      # Converts the contents of the processed file into a sha224 hash.
		      # Necessary to prevent duplicate uploads. Linux overwrites the
		      # equally named files. Also used as a primary key in the sample table
		      # of the database.
		      hash <- system(paste("sha224sum ", paste(ufDir, "/", ufCount, filetype, sep = ""), " | awk '{print $1}'", sep = ""), intern=TRUE)
		      if (hash %in% ufSample_table_content) {
		        system(paste("rm -f ", ufDir, "/", ufCount, filetype, sep = ""))
		      }
		      else {
  	        file.rename(paste(ufDir, "/", ufCount, filetype, sep = ""), paste(ufDir, "/", hash, filetype, sep = ""))
  	        filepath <- toString(paste(ufDir, "/", hash, filetype, sep = ""))
  	        # Code below is wrapped in a tryCatch in case a hash is already
  	        # present as a primary key in the sample table.
  	        tryCatch(
  	          {
  	            # Retrieve information about the specifics of a file.
  	            aa <- openMSfile(filepath)
  	            aai <- instrumentInfo(aa)
  	            test.empty <- header(aa)
  	            # This removes files with a correct layout but no measured scans
  	            test.empty <- test.empty$lowMZ
  	            on.exit(close(aa))
  	            close(aa)
  	            if (length(test.empty) == 0) {
  	              file.remove(filepath)
  	            }
  	            # This assigns the correct chromatographic method to the sample.
  	            if (length(test.empty) > 1) {
  	              if (aai$ionisation == "electrospray ionization") {
  	                chromtype <- 1 # LCMS
  	              } else {
  	                chromtype <- 2 # GCMS
  	              }
  	              # Saves the sample to a separate ogXCMSset, in order to speed
  	              # up TIC creation in later steps.
  	              ogXCMSset <- readMSData(filepath, mode = "onDisk")
  	              saveRDS(ogXCMSset, paste(ufDir, "/", hash, ".rds", sep = ""))
  	              # Creation of a cataframe with corresponding column values
  	              # of the sample table in the sample table of the database.
  	              todf <- data.frame(
  	                sample_hash = toString(hash),
  	                file_path = toString(filepath),
  	                upload_date = toString(ufTime),
  	                sample_description = toString(ufSample_description),
  	                metadata = toString(ufFile_tag),
  	                chromatography_type = toString(chromtype),
  	                original_file_name = toString(tools::file_path_sans_ext(ufFile[ufCount + 1])),
  	                original_XCMSnExp_path = toString(paste(ufDir, "/", hash, ".rds", sep = ""))
  	              )
  	              # return(todf)
  	              send_query(stringr::str_glue("SET GLOBAL local_infile=1;"))
  	              insert_query("sample", todf)
  	            }
  	          },
  	          error = function(cnd){
  	            # Removes the sample if an error has occurred.
  	            file.remove(filepath)
  	            file.remove(paste(ufDir, "/", hash, ".rds", sep = ""))
  	            system("echo ERROR")
  	            print("REMOVED FILE")
  	            return(NA)
  	          }
  	        )
		      }
		    }
		    #for Nested asynchonous function.
		    future({upload_file(count, fileloc, rfuDir, rfuTime, rfuSample_description, rfuFile$name, rfuFile_tags[count + 1], rfuDb_usr, rfuD_pwd, rfuSample_table_content)}, seed = NULL)
		    #upload_file(count, fileloc, rfuDir, rfuTime, rfuSample_description, rfuFile$name, rfuFile_tags[count + 1], rfuDb_usr, rfuD_pwd, rfuSample_table_content)
		    count <- count + 1
		    # sink()
		  }
		}
		sample_descr <- toString(input$sample_description)
		sample_hashes_existing <- sample_table_content[, 1]
		# The "run_file_upload" function is called within furute, which allows for
		# parallel processing. Necessary to maintain full application functionality.
		future({run_file_upload(file, time, dir, sample_descr, db_usr, db_pwd, sample_hashes_existing, file_tags)}, seed = NULL)
		updateTextInput(session, "sample_description", value = "")
		shinyjs::alert('Data upload in progress... Check back later!\nWarning: Do not close the current tab!')
	})
  
	# Update table to selected chomatographic method
	observeEvent(input$datatype, {
	  output$selected_samples <- DT::renderDataTable({
	    selsamples <<- sample_table_content[input$uploaded_samples_rows_selected,]
	    selsamples <<- selsamples[selsamples$chromatography_type == input$datatype,]
	    selsamples$chromatography_type[selsamples$chromatography_type == 1] <- "Liquid"
	    selsamples$chromatography_type[selsamples$chromatography_type == 2] <- "Gas"
	    if (input$tabSwitch == "p2") {
	      if (length(selsamples$sample_hash) >= 1){
	        preview_files <- selsamples$sample_hash
	        filequery <- paste("sample_hash = '", preview_files[1], "'", sep = "")
	        if (length(selsamples >= 2)){
	          for (file in preview_files[2:length(preview_files)]) {
	            filequery <- paste(filequery, " OR sample_hash = '", file, "'", sep = "")
	          }
	        }
	        query <- stringr::str_glue("SELECT file_path, original_file_name, original_XCMSnExp_path FROM sample WHERE ", filequery, ";")
	        selection <- get_query(query)
	        # Render a TIC per sample in the same plotly plot.
	        minrts <- c()
	        maxrts <- c()
	        output$ticplot <- renderPlotly({
	          chrom <- plot_ly(x = NULL, y = NULL, type = 'scatter', mode = 'lines', fill = 'tozeroy') %>% layout(title = "Total Ion Current", xaxis = list(title = "Retention time (seconds)"), yaxis = list(title = "TotIonCurrent"))
	          for (i in 1:length(preview_files)){
	            aa <- readRDS(selection$original_XCMSnExp_path[i], refhook = NULL)
	            minrts <- c(min(rtime(aa)), minrts)
	            maxrts <- c(max(rtime(aa)), maxrts)
	            y <- aa@featureData@data$totIonCurrent
	            x <- aa@featureData@data$retentionTime
	            chrom <- chrom %>% add_lines(x = x, y = y, name = selection$original_file_name[i])
	          }
	          output$rtRangeSlider <- renderUI(sliderInput("rtrange", label = "Select a minimum and maximum retention time.", value = c(min(minrts), max(maxrts)), min = min(minrts), max = max(maxrts)))
	          chrom
	        })
	      }
	      # Removes the plots when not in use.
	      else{
	        output$ticplot <- renderPlotly({
	          return(NULL)
	        })
	      }
	      output$previewplot <- renderUI({
	        return(NULL)
	      })
	    }
	    else{
	      output$ticplot <- renderPlotly({
	        return(NULL)
	      })
	    }
	    output$previewplot <- renderUI({
	      return(NULL)
	    })
	    DT::datatable(selsamples[, c(7, 3, 5, 4, 6)], selection = 'single')
	  }, server = FALSE)
	})
	
	# Update parameter usage to datatype and preset selection
	observeEvent(input$preset, {
# 	#Disable all parameter inputs for automatic mode
  	for(app_element in all_params) {
  		if(input$preset == 3) {
  			shinyjs::disable(app_element)
  		}
  	  else {
      	shinyjs::enable(app_element)
  	  }
  	}
  	# 	else {
  	# 		if(input$datatype == 2) { #PAS AAN
  	# 			if(!(app_element %in% lcms_only)) {
  	# 				shinyjs::enable(app_element)
  	# 			}
  	# 		}
  	# 		else {
  	# 			shinyjs::enable(app_element)
  	# 		}
  	# 	}
  	# }
  	# })
  	# observeEvent(input$datatype, { #PAS AAN
  	# #Disable/enable buttons for LC-MS
  	# for(lc_element in lcms_only) {
  	# 	if(input$datatype == 2) { #PAS AAN
  	# 		shinyjs::disable(lc_element)
  	# 	}
  	# 	else {
  	# 		if(input$preset != 3) {
  	# 			shinyjs::enable(lc_element)
  	# 		}
  	# 	}
  	# }
	})
	
	# Conditionally disable the "Preview spectrum" button
	observeEvent(input$selected_samples_rows_selected, ignoreNULL = FALSE, {
	  output$previewplot <- renderUI({
	    return(NULL)
	  })
		if (is.null(input$selected_samples_rows_selected)) {
			shinyjs::disable("preview")
		}
		else {
			shinyjs::enable("preview")
		}
	})
	
	observeEvent(input$jobs_rows_selected, ignoreNULL = FALSE, {
	  if (is.null(input$jobs_rows_selected)) {
	    shinyjs::disable("analyse")
	  }
	  else if (job_table_content[input$jobs_rows_selected, ]$job_status != "Finished") {
	    shinyjs::disable("analyse")
	  }
	  else {
	    shinyjs::enable("analyse")
	  }
	  
	})
	
	# Change of parameter values to presets. Currently not in use.
	observeEvent(input$preset, {
		# if(input$preset == 1) {
		# 	updateCheckboxInput(session, "p1", value = 0)
		# 	updateCheckboxInput(session, "p2", value = 1)
		# 	#updateRadioButtons(session, "datatype", selected = 1)
		# 	updateCheckboxInput(session, "p3", value = 0) }
		# else if(input$preset == 2) {
		# 	updateCheckboxInput(session, "p1", value = 1)
		# 	updateCheckboxInput(session, "p2", value = 0)
		# 	#updateRadioButtons(session, "datatype", selected = 2)
		# 	updateCheckboxInput(session, "p3", value = 1) }
		# else if(input$preset == 0) {
		# 	updateCheckboxInput(session, "p1", value = 0)
		# 	updateCheckboxInput(session, "p2", value = 0)
		# 	updateCheckboxInput(session, "p3", value = 0) }
	})
	
	# Actions when job is submitted
	observeEvent(input$submitJob, {
	  output$jobval <- renderTable({
	    validate(need(nchar(input$job_name) != 0, "No job title avaiable"))
	    validate(need(length(selsamples$sample_hash) >= 3, "Please select at least 3 samples"))
	  })
	  if (nchar(input$job_name) != 0) {
	    if (length(selsamples$sample_hash) >= 1) {
	      params <- c()
	      used_param_names <- c()
	      for (p in all_params) {
	        if(is.null(input[[p]])) {
	        }
	        else {
	          if (p == "rtrange") {
	            params <- c(params, input[[p]][1], input[[p]][2])
	            used_param_names <- c(used_param_names, "rtrmin", "rtrmax")
	          }
	          else {
	            params <- c(params, input[[p]])
	            used_param_names <- c(used_param_names, p) 
	          }
	        }
	      }
	      names(params) <- used_param_names
	      
	      #Insert job with status
	      params <- c(params, sample_type=input$datatype)
	      start_time <- format(Sys.time() + 60*60, "%Y-%m-%d %X")
	      job_status <- "Initializing..."
	      preset <- input$preset
	      todf <- data.frame(
	        job_status = toString(job_status),
	        start_time = toString(start_time),
	        job_name = toString(input$job_name)
	      )
	      insert_query("job", todf)
	      
	      query <- stringr::str_glue("SELECT job_id FROM job ORDER BY job_id DESC LIMIT 1;")
	      job_id <- get_query(query)
	      params <- c(params, job_id=job_id)
	      # params$Peak_method <- c("centWave", "massifquant", "matchedFilter")[strtoi(input$Peak_method) + 1]
	      # params$RT_method <- c("loess", "obiwarp")[strtoi(input$RT_method) + 1]
	      # params$mzCenterFun <- c("wMean", "mean", "apex", "wMeanApex3", "meanApex3")[strtoi(input$mzCenterFun) + 1]
	      # params$smooth <- c("loess", "linear")[strtoi(input$smooth) + 1]
	      # params$family <- c("gaussian", "symmetric")[strtoi(input$family) + 1]
	      # params$polarity <- c("negative", "positive")[strtoi(input$polarity) + 1]
	      
	      params_insert <- do.call(rbind, params) %>% as.data.frame() %>% t() %>% as.data.frame()
	      names(params_insert)[names(params_insert) == 'job_id.job_id'] <- 'job_id'
	      insert_query("parameter", params_insert)
	      
	      jobfiles <- NULL
	      sample_number <- 1
	      for (hash in selsamples$sample_hash) {
	        query <- stringr::str_glue("SELECT * FROM sample WHERE sample_hash = '", hash, "';")
	        jobfiles <- rbind(jobfiles, get_query(query))
	        todf <- data.frame(
	          job_id = toString(job_id),
	          sample_hash = toString(hash),
	          sample_number = toString(sample_number)
	        )
	        insert_query("sample_job", todf)
	        sample_number <- sample_number + 1
	      }
	      shinyjs::alert(paste('Your job "', input$job_name, '" is running. Check progress in the "Jobs" tab.', sep = ""))
	      future({xcms_data_processing(jobfiles, params, preset, job_id, db_usr, db_pwd, 7)}, seed = NULL)
	      #xcms_data_processing(jobfiles, params, preset, job_id, db_usr, db_pwd, 7)
	      # if (input$debug == 0) {
	      #   #future({metaboanalyst_data_processing(jobfiles, params, preset, job_id, db_usr, db_pwd)}, seed = NULL)#, packages = .packages(TRUE))
	      #   future({xcms_data_processing(jobfiles, params, preset, job_id, db_usr, db_pwd, 7)}, seed = NULL)
	      # }
	      # else {
	      #   #metaboanalyst_data_processing(jobfiles, params, preset, job_id, db_usr, db_pwd) #debug
	      #   xcms_data_processing(jobfiles, params, preset, job_id, db_usr, db_pwd, 7)
	      # }
	      session$reload()
	      updateTabsetPanel(session, "tabSwitch",
	                        selected = "p3")
	    }
	  }
	})
	
	# Function to tweak parameters of a previously processed job. Ran if a resub
	# button was pressed.
	observeEvent(input$resub1 | input$resub2 | input$resub3 | input$resub4 | input$resub5, ignoreInit =  T, {
	  job_plan <- NA
	  # Defining the plan of the job that needs to be run. A job runs all plans
	  # greater then the job_plan variable.
	  if (input$resub1 == 1) {
	    job_plan <- 2
	  }
	  else if (input$resub2 == 1) {
	    job_plan <- 3
	  }
	  else if (input$resub3 == 1) {
	    job_plan <- 4
	  }
	  else if (input$resub4 == 1) {
	    job_plan <- 5
	  }
	  else if (input$resub5 == 1) {
	    job_plan <- 6
	  }
	  # If a valid job_plan has been assigned.
	  if (!is.na(job_plan)) {
	    job_id <- job_table_content[input$jobs_rows_selected, ]$job_id
	    params <- c()
	    used_param_names <- c()
	    for (p in all_params) {
	      if(is.null(input[[paste(p, "A", sep = "")]]) & p != "rtrange") {
	        params <- c(params, "NULL")
	        used_param_names <- c(used_param_names, p) 
	      }
	      else {
	        if (p == "rtrange") {
	          #params <- c(params, input[[paste(p, "A", sep = "")]][1], input[[paste(p, "A", sep = "")]][2])
	          #used_param_names <- c(used_param_names, "rtrmin", "rtrmax")
	        }
	        else {
	          params <- c(params, input[[paste(p, "A", sep = "")]])
	          used_param_names <- c(used_param_names, p) 
	        }
	      }
	    }
	    updatestring <- ""
	    for (i in 1:length(params)) {
	      if (is.na(params[i])) {
	        updatestring <- paste(updatestring, used_param_names[i], "=", "NULL", sep = "")
	      } else {
	        if (params[i] == "TRUE"){
	          params[i] <- 1
	        }
	        else if (params[i] == "FALSE") {
	          params[i] <- 0
	        }
	        if (used_param_names[i] == "index"){
	          updatestring <- paste(updatestring, "`index`", "=", params[i], sep = "")
	        }
	        else{
	          updatestring <- paste(updatestring, used_param_names[i], "=", params[i], sep = "")
	        }
	      }
	      if (i != length(params)) {
	        updatestring <- paste(updatestring, ",", sep = "")
	      }
	    }
	    send_query(stringr::str_glue(paste("UPDATE job SET job_status = 'Re-initializing...' WHERE job_id = ", job_id, ";", sep = "")))
	    send_query(stringr::str_glue(paste("UPDATE parameter SET ", updatestring, " WHERE job_id = ", job_id, ";", sep = "")))
	    query <- stringr::str_glue(paste("SELECT * FROM sample WHERE sample_hash IN (SELECT sample_hash FROM sample_job WHERE job_id = ", job_id, ");", sep = ""))
	    jobfiles <<- get_query(query)
	    names(params) <- used_param_names
	    params <- c(params, job_id=job_id)
	    params <- as.data.frame(t(params))
	    query <- stringr::str_glue(paste("SELECT * FROM parameter WHERE job_id = ", job_id, ";", sep = ""))
	    params <- get_query(query)
	    future({xcms_data_processing(jobfiles, params, 0, job_id, db_usr, db_pwd, job_plan)}, seed = NULL)
	    #xcms_data_processing(jobfiles, params, 0, job_id, db_usr, db_pwd, job_plan)
	    session$reload()
	    updateTabsetPanel(session, "tabSwitch",
	                      selected = "p3")
	  }
	})
	
	#Actions involving job analysis
	observeEvent(input$analyse, {
	  job_id <- job_table_content[input$jobs_rows_selected, ]$job_id
	  updateTabsetPanel(session, "tabSwitch",
	                    selected = "p4")
	  query1 <- stringr::str_glue(paste("SELECT * FROM processed_sample WHERE job_id = ", job_id, ";", sep = ""))
	  rda_path <- get_query(query1)
	  
	  query2 <- stringr::str_glue(paste("SELECT * FROM parameter WHERE job_id = ", job_id, ";", sep = ""))
	  used_parameters <- get_query(query2)
	  
	  query3 <- stringr::str_glue(paste("SELECT original_file_name, metadata FROM sample WHERE sample_hash IN (SELECT sample_hash FROM sample_job WHERE job_id = ", job_id, ");", sep = ""))
	  file_info <- get_query(query3)
	  
	  load(file=toString(rda_path$file_path_rda))
	  query4 <- stringr::str_glue(paste("SELECT * FROM sample WHERE sample_hash IN (SELECT sample_hash FROM sample_job WHERE job_id = ", job_id, ");", sep = ""))
	  selection <- get_query(query4)
	  
	  query5 <- stringr::str_glue(paste("SELECT peak.peak_id, peak.modcosinesim, mol.mol_name, mol.pubid, peak.spectrum FROM peak LEFT JOIN mol ON peak.mol_id_modcosinesim = mol.mol_id WHERE job_id = ",
	                                   job_id, ";", sep = ""))
	  detected_mols_cos <- get_query(query5)
	  query6 <- stringr::str_glue(paste("SELECT peak.peak_id, peak.coeff_diff, mol.mol_name, mol.pubid, peak.spectrum FROM peak LEFT JOIN mol ON peak.mol_id_coeffsim = mol.mol_id WHERE job_id = ",
	                                    job_id, ";", sep = ""))
	  detected_mols_coeff <- get_query(query6)
	  
	  detected_mols_cos <-detected_mols_cos[order(detected_mols_cos$peak_id),]
	  detected_mols_coeff <-detected_mols_coeff[order(detected_mols_coeff$peak_id),]
	  dfpeaks <- as.data.frame(chromPeaks(xset))
	  dfpeaks <- dfpeaks[order(rownames(dfpeaks)),]
	  dfpeaks[rownames(dfpeaks) %in% detected_mols_cos$peak_id,"mol_name_cosd"] <- detected_mols_cos$mol_name
	  dfpeaks[rownames(dfpeaks) %in% detected_mols_cos$peak_id,"pubid_cosd"] <- detected_mols_cos$pubid
	  dfpeaks[rownames(dfpeaks) %in% detected_mols_cos$peak_id,"cossim"] <- detected_mols_cos$modcosinesim
	  # dfpeaks[rownames(dfpeaks) %in% detected_mols_coeff$peak_id,"mol_name_coeff"] <- detected_mols_coeff$mol_name
	  # dfpeaks[rownames(dfpeaks) %in% detected_mols_coeff$peak_id,"pubid_coeff"] <- detected_mols_coeff$pubid
	  # dfpeaks[rownames(dfpeaks) %in% detected_mols_coeff$peak_id,"coeff_diff"] <- detected_mols_coeff$coeff_diff
	  # Add sample name to dataframe
	  dfpeaks$sample_name <- NULL
	  # Order the dfpeaks according to original xmcs ordering
	  dfpeaks <- dfpeaks[rownames(chromPeaks(xset)),]
	  feature_mols <- data.frame(mol = NULL)
	  # Add a compound to each feature
	  for (feature in rownames(featureDefinitions(xset))) {
	    # Retrieve peak ids for a feature
	    dfpeaks <- dfpeaks[rownames(chromPeaks(xset)),]
	    peakids_for_mols <- featureDefinitions(xset)[feature, "peakidx"][[1]]
	    # Retrieve most frequent compound
	    most_frequent <- names(which.max(table(na.omit(dfpeaks[peakids_for_mols,"mol_name_cosd"]))))
	    # Replace NULL values by Unknown Compound
	    if(is.null(most_frequent)) {
	      most_frequent <- "Unknown Compound"
	    }
	    # Create a new dataframe from the most frequent compound
	    add_frame <- data.frame(mol = most_frequent)
	    rownames(add_frame) <- feature
	    # Add the frame to the annotation frame
	    feature_mols <- rbind(feature_mols, add_frame)
	  }
	  # Responsible for everything in the model when the user clicks a peak
	  observeEvent(input$peaktable_rows_selected, {
	    pspec <- detected_mols_cos[detected_mols_cos$peak_id == rownames(dfpeaks)[input$peaktable_rows_selected],"spectrum"]
	    # Split the spectrum ([mz]:[intensity]) string into mz and intensity values
	    # Dissect the string into a vector of mzs and intensities
	    bins <- strsplit(pspec, " ")
	    if (!is.null(bins)) {
	      # Split into two vectors to allow plotting
	      ints <- c()
	      mzs <- c()
	      if (length(bins) != 0) {
	        for (section in bins[[1]]){
	          mzs <- c(mzs, as.double(strsplit(section, ":")[[1]][1]))
	          ints <- c(ints, as.double(strsplit(section, ":")[[1]][2]))
	        }
	        pspec <- data.frame(mzs, ints)
	        colnames(pspec) <- c("Mass (m/z)", "Intensity (relative to highest peak)")
	        # Allows for downloadign the spectrum to a .txt file.
	        output$download_spectr <- downloadHandler(
	          filename = function(){"spectrum.txt"}, 
	          content = function(file){
	            writeLines(paste(detected_mols_cos[detected_mols_cos$peak_id == rownames(dfpeaks)[input$peaktable_rows_selected],"spectrum"], collapse = ", "), file)
	          }
	        )
	        # Rendering the modal
	        showModal(modalDialog(
	          tagList(
	            renderPlot({
	              plot(pspec, type = 'h')
	              text(x = pspec[pspec[,2] > 10,1], y = pspec[pspec[,2] > 10,2] + 1, round(pspec[pspec[,2] > 10,1], 3), col = 'blue')
	            }),
	            # Add button to download spectrum
	            downloadButton('download_spectr', "Download spectrum [mz]:[intensity]")
	          ), easyClose = TRUE
	        ))
	      }
	    }
	  })
	  # Add file name to datatable
	  for (sample_index in 1:length(selection$original_file_name)) {
	    dfpeaks[dfpeaks[,"sample"] == sample_index,"sample_name"] <- selection$original_file_name[sample_index]
	  }
	  # Render peak table
	  output$peaktable <- DT::renderDataTable({
	    DT::datatable(dfpeaks, extensions = "Buttons", selection = 'single',
	                  rownames= TRUE, options = list(scrollX = TRUE, autoWidth = TRUE, dom = 'Bfrtip',
	                                                 buttons = c('csv')))
	  }, server = FALSE)
	  
	  # Needed in order to render parameters at bottom of page
	  observeEvent(input$Ref_methodA, ignoreInit =  T, {
	    ref_options <- c(default_mnp, default_fi)[is.na(input$Ref_methodA) + 1]
	  })
	  # Peak refinement
	  if (is.null(input$Ref_methodA)){
	    if (used_parameters$Ref_method == 0) {
	      ref_options <- tagList(
	        numericInput("expandRtA", label = "expandRt", value = used_parameters$expandRt),
	        numericInput("expandMzA", label = "expandMz", value = used_parameters$expandMz),
	        numericInput("minPropA", label = "minProp", value = used_parameters$minProp)
	      )
	    }
	    else if (used_parameters$Ref_method == 1){
	      ref_options <- tagList(
	        numericInput("thresholdA", label = "threshold", value = used_parameters$threshold)
	      )
	    }
	  } else {
	    if (used_parameters$Ref_method == 0 & input$Ref_methodA == 0) {
	      ref_options <- tagList(
	        numericInput("expandRtA", label = "expandRt", value = used_parameters$expandRt),
	        numericInput("expandMzA", label = "expandMz", value = used_parameters$expandMz),
	        numericInput("minPropA", label = "minProp", value = used_parameters$minProp)
	      )
	    }
	    else if (used_parameters$Ref_method == 1 & input$Ref_methodA == 1){
	      ref_options <- tagList(
	        numericInput("thresholdA", label = "threshold", value = used_parameters$threshold)
	      )
	    }
	  }
	  
	  # Peak alignment
	  observeEvent(input$Align_methodA, ignoreInit =  T, {
	    aln_options <- c(default_obiwarp, default_peakgroups)[is.na(input$Align_methodA) + 1]
	  })
	  if (is.null(input$Align_methodA)) {
	    if (used_parameters$Align_method == 0) {
	      aln_options <- tagList(
	        #binSize = 1,
	        numericInput("binSizeA", label = "binSize", used_parameters$binSize),
	        # numericInput("centerSample", label = "centerSample", 5),
	        # numericInput("response", label = "response", 1),
	        selectInput("distFunA", label = "distFun", choices = list("cor" = 0, "cor_opt" = 1, "cov" = 2, "prd" = 3, "euc" = 4), selected = used_parameters$distFun),
	        # numericInput("gapInit", label = "gapInit", 5),
	        # numericInput("gapExtend", label = "gapExtend", 5),
	        numericInput("factorDiagA", label = "factorDiag", used_parameters$factorDiag),
	        numericInput("factorGapA", label = "factorGap", used_parameters$factorGap),
	        checkboxInput("localAlignmentA", label = "localAlignment", used_parameters$localAlignment),
	        numericInput("initPenaltyA", label = "initPenalty", used_parameters$initPenalty),
	        # numericInput("subset", label = "subset", 5),
	        selectInput("subsetAdjustA", label = "subsetAdjust", choices = list("average" = 0, "previous" = 1), selected = used_parameters$subsetAdjust)
	      )
	    }
	    else if (used_parameters$Align_method == 1) {
	      aln_options <- tagList(
	        numericInput("minFractionA", label = "minFraction", used_parameters$minFraction),
	        numericInput("extraPeaksA", label = "extraPeaks", used_parameters$extraPeaks),
	        selectInput("smoothA", label = "smooth", choices = list("loess" = 0, "linear" = 1), selected = used_parameters$smooth),
	        numericInput("spanA", label = "span", used_parameters$span),
	        selectInput("familyA", label = "family", choices = list("gaussian" = 0, "symmetric" = 1), selected = used_parameters$family),
	        # numericInput("peakGroupsMatrix", label = "peakGroupsMatrix", 5),
	        # numericInput("subset", label = "subset", 5),
	        selectInput("subsetAdjustA", label = "subsetAdjust", choices = list("average" = 0, "previous" = 1), selected = used_parameters$subsetAdjust)
	      )
	    }
	  } else {
	    if (used_parameters$Align_method == 0 & input$Align_methodA == 0) {
	      aln_options <- tagList(
	        #binSize = 1,
	        numericInput("binSizeA", label = "binSize", used_parameters$binSize),
	        # numericInput("centerSample", label = "centerSample", 5),
	        # numericInput("response", label = "response", 1),
	        selectInput("distFunA", label = "distFun", choices = list("cor" = 0, "cor_opt" = 1, "cov" = 2, "prd" = 3, "euc" = 4), selected = used_parameters$distFun),
	        # numericInput("gapInit", label = "gapInit", 5),
	        # numericInput("gapExtend", label = "gapExtend", 5),
	        numericInput("factorDiagA", label = "factorDiag", used_parameters$factorDiag),
	        numericInput("factorGapA", label = "factorGap", used_parameters$factorGap),
	        checkboxInput("localAlignmentA", label = "localAlignment", used_parameters$localAlignment),
	        numericInput("initPenaltyA", label = "initPenalty", used_parameters$initPenalty),
	        # numericInput("subset", label = "subset", 5),
	        selectInput("subsetAdjustA", label = "subsetAdjust", choices = list("average" = 0, "previous" = 1), selected = used_parameters$subsetAdjust)
	      )
	    }
	    else if (used_parameters$Align_method == 1 & input$Align_methodA == 1) {
	      aln_options <- tagList(
	        numericInput("minFractionA", label = "minFraction", used_parameters$minFraction),
	        numericInput("extraPeaksA", label = "extraPeaks", used_parameters$extraPeaks),
	        selectInput("smoothA", label = "smooth", choices = list("loess" = 0, "linear" = 1), selected = used_parameters$smooth),
	        numericInput("spanA", label = "span", used_parameters$span),
	        selectInput("familyA", label = "family", choices = list("gaussian" = 0, "symmetric" = 1), selected = used_parameters$family),
	        # numericInput("peakGroupsMatrix", label = "peakGroupsMatrix", 5),
	        # numericInput("subset", label = "subset", 5),
	        selectInput("subsetAdjustA", label = "subsetAdjust", choices = list("average" = 0, "previous" = 1), selected = used_parameters$subsetAdjust)
	      )
	    }
	  }
	  # Peak grouing
	  observeEvent(input$Group_methodA, ignoreInit =  T, {
	    grp_options <- c(default_pd, default_mzc, default_np)[is.na(input$Group_methodA) + 1]
	  })
	  if (is.null(input$Group_methodA)) {
	    if (used_parameters$Group_method == 0) {
	      grp_options <- tagList(
	        numericInput("bwA", label = "bw", used_parameters$bw),
	        numericInput("minFractionA", label = "minFraction", used_parameters$minFraction),
	        numericInput("minSamplesA", label = "minSamples", used_parameters$minSamples),
	        numericInput("binSizeA", label = "binSize", used_parameters$binSize),
	        numericInput("maxFeaturesA", label = "maxFeatures", used_parameters$maxFeatures)
	      )
	    }
	    else if (used_parameters$Group_method == 1) {
	      grp_options <- tagList(
	        numericInput("absMzA", label = "absMz", used_parameters$absMz),
	        numericInput("minFractionA", label = "minFraction", used_parameters$minFraction),
	        numericInput("minSamplesA", label = "minSamples", used_parameters$minSamples)
	      )
	    }
	    else if (used_parameters$Group_method == 2) {
	      grp_options <- tagList(
	        numericInput("mzVsRtBalanceA", label = "mzVsRtBalance", used_parameters$mzVsRtBalance),
	        numericInput("absMzA", label = "absMz", used_parameters$absMz),
	        numericInput("absRtA", label = "absRt", used_parameters$absRt),
	        numericInput("kNNA", label = "kNN", used_parameters$kNN)
	      )
	    }
	  } else {
	    if (used_parameters$Group_method == 0 & input$Group_methodA == 0) {
	      grp_options <- tagList(
	        numericInput("bwA", label = "bw", used_parameters$bw),
	        numericInput("minFractionA", label = "minFraction", used_parameters$minFraction),
	        numericInput("minSamplesA", label = "minSamples", used_parameters$minSamples),
	        numericInput("binSizeA", label = "binSize", used_parameters$binSize),
	        numericInput("maxFeaturesA", label = "maxFeatures", used_parameters$maxFeatures)
	      )
	    }
	    else if (used_parameters$Group_method == 1 & input$Group_methodA == 1) {
	      grp_options <- tagList(
	        numericInput("absMzA", label = "absMz", used_parameters$absMz),
	        numericInput("minFractionA", label = "minFraction", used_parameters$minFraction),
	        numericInput("minSamplesA", label = "minSamples", used_parameters$minSamples)
	      )
	    }
	    else if (used_parameters$Group_method == 2 & input$Group_methodA == 2) {
	      grp_options <- tagList(
	        numericInput("mzVsRtBalanceA", label = "mzVsRtBalance", used_parameters$mzVsRtBalance),
	        numericInput("absMzA", label = "absMz", used_parameters$absMz),
	        numericInput("absRtA", label = "absRt", used_parameters$absRt),
	        numericInput("kNNA", label = "kNN", used_parameters$kNN)
	      )
	    }
	  }
	  ### Code necessary in order to create a stacked bar plot of all detected compounds
	  molmatches <- data.frame(compound = character(0), sample = numeric(0), area = numeric(0))
	  for (sample in 1:length(unique(chromPeaks(xset)[,"sample"]))) {
	    comp <- dfpeaks[dfpeaks[,"sample"] == sample,]
	    for (row in 1:nrow(comp)){
	      mol_name <- comp[row,"mol_name_cosd"]
	      if (mol_name == '' | is.na(mol_name) == TRUE){
	        mol_name <- "Unknown Compound"
	      }
	      newrow <- data.frame(mol_name, sample, comp[row, "into"])
	      names(newrow) <- c("compound", "sample", "area")
	      molmatches <- rbind(molmatches, newrow)
	    }
	  }
	  molmatches$area <- log10(molmatches$area)
	  ## Splits the frame for the stacked plot into separate sub-frame. Necessary for plotly.
	  stacked_samps <- split(molmatches, molmatches$sample)
	  analysisTL <- renderUI(tagList(
	    h4(paste("Output of job", job_id, ":\n", job_table_content[input$jobs_rows_selected, ]$job_name)),
        fluidRow(style='max-width:100%;padding:10px;',
  			         column(12,
  			                uiOutput("jobInformation")
  			                )
  			         ),
  			fluidRow(style='max-width:100%;padding:10px;',
  			         column(6,
  			                # Render stacked bar plot of detected compounds
  			                # renderPlotly({
  			                #   stacked_mols <- plot_ly(x = NULL, y = NULL, type = 'bar', name = 'SF Zoo')
  			                #   for (sample in 1:length(unique(chromPeaks(xset)[,"sample"]))) {
  			                #     stacked_mols <- stacked_mols %>% add_trace(stacked_samps[[sample]], x = stacked_samps[[sample]]$compound, y = stacked_samps[[sample]]$area, name =  selection$original_file_name[sample])
  			                #   }
  			                #   stacked_mols <- stacked_mols %>% layout(yaxis = list(title = 'Log10 of Detected Compounds AUC'), barmode = 'stack')
  			                #   stacked_mols
  			                # })
  			         ),
  			         column(6,
  			                div(DT::renderDataTable({
  			                  file_info["Number_of_detected_peaks"] <- unlist(lapply(1:nrow(file_info), function(i){nrow(chromPeaks(xset, bySample = TRUE)[[i]])}))
  			                  DT::datatable(file_info, selection = 'none',
  			                                rownames= TRUE, options = list(scrollY = TRUE))
  			                }, server = FALSE))#, style = "font-size: 75%; width: 50%")

  			         )
  			),
  			fluidRow(style='max-width:100%;padding:10px;',
  			         column(6,
  			                renderPlotly({
  			                  # Take the log2 of the feature values
  			                  #ft_fv <- na.omit(log2(featureValues(xset, filled = TRUE)))
  			                  # Normalize the dataset using scale
  			                  ft_fv <- scale(featureValues(xset, filled = TRUE), scale = TRUE)
  			                  pc <- prcomp(t(na.omit(ft_fv)), center = TRUE)
  			                  pca <- plot_ly(x = NULL, y = NULL, type = "scatter") %>% layout(title  = "PCA", scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2")))
  			                  #hashes <- substr(names(pc$x[,1]), 0, 56)
  			                  meta_selection <- selection$metadata
  			                  for (metda in unique(meta_selection)){
  			                    pca <- pca %>% add_trace(y = pc$x[,2][meta_selection == metda], x = pc$x[,1][meta_selection == metda], name = paste(metda), mode = 'markers', hovertext = selection$original_file_name[meta_selection == metda], hoverinfo = "text", marker = list(size = 10))
  			                  }
  			                  pca
  			                })
  			         ),
  			         column(6,
  			                # TIC plot and detected peaks
  			                numericInput("sample_nr_peaks", label = "Select a sample to view its detected peaks.", min = 1, max = nrow(selection), value = 1),
  			                renderPlotly({
  			                  chrom <- plot_ly(x = NULL, y = NULL, type = 'scatter', mode = 'markers') %>% layout(title = "Peaks over TIC", xaxis = list(title = "Retention time (seconds)"), yaxis = list(title = "TotIonCurrent"))
  			                  #featspec <- featureSpectra(xset, msLevel = 1, return.type="Spectra")
  			                  load(file = rda_path$file_path_peaks) #loads variable annot_spectra
  			                  tt <- spectraData(annot_spectra)
  			                  #aa <- readRDS(selection$original_XCMSnExp_path[i], refhook = NULL)
  			                  for (i in round(input$sample_nr_peaks)){
  			                    aa <- filterFile(xset, i)
  			                    # Omit filled peaks
  			                    is_filled <-chromPeakData(aa)[,2]
  			                    # Get ids of peaks in sample
  			                    tp <- rownames(chromPeaks(aa))[!is_filled]
  			                    for (peak in 1:length(tp)) {
  			                      p <- unique(tt[tt[, "peak_id"] == tp[peak], c("peak_id", "scanIndex", "rtime", "totIonCurrent")])
  			                      xp <- p[,"rtime"]
  			                      yp <- p[,"totIonCurrent"]
  			                      label1 <- dfpeaks[tp[peak], "mol_name_cosd"]
  			                      label2 <- dfpeaks[tp[peak], "mol_name_coeff"]
  			                      chrom <- chrom %>% add_lines(x = xp, y = yp, name = paste(tp[peak]), fill = 'tozeroy', hovertext = paste(peak, ": ", label1, "|", label2), hoverinfo = "text")
  			                    }
  			                    y <- tic(aa)
  			                    x <- rtime(aa)
  			                    chrom <- chrom %>% add_lines(x = x, y = y, name = paste(selection$original_file_name[i], "|", selection$metadata[i]), hoverinfo = 'skip')
  			                  }
  			                  chrom
  			                })
  			         )
  			),
  			fluidRow(style='max-width:100%;padding:10px;',
  			         column(6,
  			                # Potential future addition
  			                uiOutput("analysisCompounds")
  			         ),
  			         column(6,
  			                # Potetial future addition
  			                uiOutput("analysisDifferential")
  			         )
  			),
  			fluidRow(style='max-width:100%;padding:10px;',
  			         column(12,
  			                # Group correlation selection
  			                renderUI({
  			                    checkboxGroupInput(inputId = "matrixselect", label = "Select groups to alter view of heatmap", choices <- unique(selection$metadata), selected = unique(selection$metadata))
  			                }),
  			                ### HEATMAP
  			                renderUI({
  			                  heatframe <- featureValues(xset, filled = TRUE)
  			                  # heatframe <- data.frame(matrix(nrow = length(unique(molmatches$compound)), ncol = length(selection$original_file_name)))
  			                  # rownames(heatframe) <- unique(molmatches$compound)
  			                  colnames(heatframe) <- selection$original_file_name
  			                  heatframe <- heatframe[,selection$metadata %in% input$matrixselect]
  			                  # for (df in stacked_samps){
  			                  #   agg <- aggregate(. ~ compound, data=df[,c(1,3)], FUN=sum)
  			                  #   heatframe[agg$compound, selection$original_file_name[df$sample[1]]] <- agg$area
  			                  # }
  			                  heatframe[is.na(heatframe)] <- 0
  			                  heatframe[is.na(heatframe)] <- 0
  			                  heatframe <- heatframe[!is.infinite(rowSums(heatframe)),]
  			                  # heatframe <- heatframe[rownames(heatframe) != "Unknown Compound",]
  			                  heatframe <<- heatframe[,sort(colnames(heatframe))]
  			                  # Allows for dataframe export to .csv
  			                  output$download_featvals <- downloadHandler(
  			                    filename = paste("featureValues", job_id, ".csv", sep = ""),
  			                    content = function(file) {
  			                      write.csv(cbind(heatframe, feature_mols), file, row.names = TRUE)
  			                    }
  			                  )
  			                  heatmaply(heatframe, xlab = "Samples", ylab = "Features", main = "Clustered heatmap of features between samples")
  			                }),
  			                # Rendering of download button
  			                downloadButton("download_featvals", "Download feature values")
  			         )
  			),
	      # Creation of correlation matrix and heatmap between features
	      fluidRow(style='max-width:100%;padding:10px;font-size:75%;', column(12, 
	                                                                          renderUI({
	                                                                            selection$metadata %in% input$matrixselect
	                                                                            heatmaply(cor(na.omit(heatframe)), main = "Pearson correlation between sample features")
	                                                                          })
	                                                                          )),
  	    fluidRow(style='max-width:100%;padding:10px;font-size:75%;',
  	             column(6,
  	                    # Rendering of peak table
  	                    p("Detected peak specifics:"),
  	                    DTOutput('peaktable')
  	             ),
  	             column(6,
  	                    # Rendering of feature table
  	                    p("Detected features:"),
  	                    div(DT::renderDataTable({
  	                      DT::datatable(as.data.frame(cbind(featureDefinitions(xset), feature_mols)), extensions = "Buttons", selection = 'none',
  	                                    rownames= TRUE, options = list(scrollX = TRUE, autoWidth = TRUE, dom = 'Bfrtip',
  	                                                                   buttons = c('csv')))
  	                    }, server = FALSE))
  	             )
  	    ),
	      # Row with the utilized parameters and rerun buttons
  	    fluidRow(style='max-width:100%;padding:10px;',
  	             column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  	             ),
  	             column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  	                    # Peak detection
  	                    h5("1. Detection"),
  	                    selectInput("Peak_methodA", label = "Peak detection method", choices = list("centWave" = 0, "Massifquant" = 1, "MatchedFilter" = 2, "relsqdiff" = 3), selected = used_parameters$Peak_method),
                      renderUI({
                        # Choose utilized method to render correct parameters
                        observeEvent(input$Peak_methodA, ignoreInit =  F, {
                          peak_options <- c(default_cent, default_massif, default_matchedfilter, default_relsqdiff)[is.na(input$Peak_methodA) + 1]
                        })
                        # Render parameters with utilized values
                        if (is.null(input$Peak_methodA)) {
                          if (used_parameters$Peak_method == 0) {
                            peak_options <- tagList(
                              numericInput("ppmA", label = "ppm", value = used_parameters$ppm),
                              numericInput("min_peakwidthA", label = "min_peakwidth", value = used_parameters$min_peakwidth),
                              numericInput("max_peakwidthA", label = "max_peakwidth", value = used_parameters$max_peakwidth),
                              numericInput("snthreshA", label = "snthresh", value = used_parameters$snthresh),
                              numericInput("prefilterA", label = "prefilter", value = used_parameters$prefilter),
                              numericInput("value_of_prefilterA", label = "value_of_prefilter", value = used_parameters$value_of_prefilter),
                              selectInput("mzCenterFunA", label = "mzCenterFun", choices = list("wMean" = 0, "mean" = 1, "apex" = 2, "wMeanApex3" = 3, "meanApex3" = 4), selected = used_parameters$mzCenterFun),
                              numericInput("integrateA", label = "integrate", value = used_parameters$integrate),
                              numericInput("mzdiffA", label = "mzdiff", value = used_parameters$mzdiff),
                              checkboxInput("fitgaussA", label = "fitgauss", value = used_parameters$fitgauss),
                              numericInput("noiseA", label = "noise", value = used_parameters$noise),
                              checkboxInput("verboseColumnsA", label = "verboseColumns", value = used_parameters$verboseColumns)
                            )
                          }
                          else if (used_parameters$Peak_method == 1) {
                            peak_options <- tagList(
                              numericInput("ppmA", label = "ppm", value = used_parameters$ppm),
                              numericInput("min_peakwidthA", label = "min_peakwidth", value = used_parameters$min_peakwidth),
                              numericInput("max_peakwidthA", label = "max_peakwidth", value = used_parameters$max_peakwidth),
                              numericInput("snthreshA", label = "snthresh", value = used_parameters$snthresh),
                              numericInput("prefilterA", label = "prefilter", value = used_parameters$prefilter),
                              numericInput("value_of_prefilterA", label = "value_of_prefilter", value = used_parameters$value_of_prefilter),
                              selectInput("mzCenterFunA", label = "mzCenterFun", choices = list("wMean" = 0, "mean" = 1, "apex" = 2, "wMeanApex3" = 3, "meanApex3" = 4), selected = used_parameters$mzCenterFun),
                              numericInput("integrateA", label = "integrate", value = used_parameters$integrate),
                              numericInput("mzdiffA", label = "mzdiff", value = used_parameters$mzdiff),
                              checkboxInput("fitgaussA", label = "fitgauss", value = used_parameters$fitgauss),
                              numericInput("noiseA", label = "noise", value = used_parameters$noise),
                              checkboxInput("verboseColumnsA", label = "verboseColumns", value = used_parameters$verboseColumns),
                              numericInput("criticalValueA", label = "criticalValue", value = used_parameters$criticalValue),
                              numericInput("consecMissedLimitA", label = "consecMissedLimit", value = used_parameters$consecMissedLimit),
                              numericInput("unionsA", label = "unions", value = used_parameters$unions),
                              numericInput("checkBackA", label = "checkBack", value = used_parameters$checkBack),
                              checkboxInput("withWaveA", label = "withWave", value = used_parameters$withWave)
                            )
                          }
                          else if (used_parameters$Peak_method == 2) {
                            peak_options <- tagList(
                              numericInput("binSizeA", label = "binSize", value = used_parameters$binSize),
                              selectInput("imputeA", label = "impute", choices = list("none" = 0, "lin" = 1, "linbase" = 2, "intlin" = 3, "imputeLinInterpol" = 4), selected = used_parameters$impute),
                              # numericInput("baseValue", label = "baseValue", 5),
                              # numericInput("distance", label = "distance", 5),
                              numericInput("fwhmA", label = "fwhm", value = used_parameters$fwhm),
                              numericInput("sigmaA", label = "sigma", value = used_parameters$sigma),
                              numericInput("maxA", label = "max", value = used_parameters$max),
                              numericInput("snthreshA", label = "snthresh", value = used_parameters$snthresh),
                              numericInput("stepsA", label = "steps", value = used_parameters$steps),
                              numericInput("mzdiffA", label = "mzdiff", value = used_parameters$mzdiff),
                              checkboxInput("indexA", label = "index", value = used_parameters$index)
                            )
                          }
                        } else {
                          if (used_parameters$Peak_method == 0 & input$Peak_methodA == 0) {
                            peak_options <- tagList(
                              numericInput("ppmA", label = "ppm", value = used_parameters$ppm),
                              numericInput("min_peakwidthA", label = "min_peakwidth", value = used_parameters$min_peakwidth),
                              numericInput("max_peakwidthA", label = "max_peakwidth", value = used_parameters$max_peakwidth),
                              numericInput("snthreshA", label = "snthresh", value = used_parameters$snthresh),
                              numericInput("prefilterA", label = "prefilter", value = used_parameters$prefilter),
                              numericInput("value_of_prefilterA", label = "value_of_prefilter", value = used_parameters$value_of_prefilter),
                              selectInput("mzCenterFunA", label = "mzCenterFun", choices = list("wMean" = 0, "mean" = 1, "apex" = 2, "wMeanApex3" = 3, "meanApex3" = 4), selected = used_parameters$mzCenterFun),
                              numericInput("integrateA", label = "integrate", value = used_parameters$integrate),
                              numericInput("mzdiffA", label = "mzdiff", value = used_parameters$mzdiff),
                              checkboxInput("fitgaussA", label = "fitgauss", value = used_parameters$fitgauss),
                              numericInput("noiseA", label = "noise", value = used_parameters$noise),
                              checkboxInput("verboseColumnsA", label = "verboseColumns", value = used_parameters$verboseColumns)
                            )
                          }
                          else if (used_parameters$Peak_method == 1 & input$Peak_methodA == 1) {
                            peak_options <- tagList(
                              numericInput("ppmA", label = "ppm", value = used_parameters$ppm),
                              numericInput("min_peakwidthA", label = "min_peakwidth", value = used_parameters$min_peakwidth),
                              numericInput("max_peakwidthA", label = "max_peakwidth", value = used_parameters$max_peakwidth),
                              numericInput("snthreshA", label = "snthresh", value = used_parameters$snthresh),
                              numericInput("prefilterA", label = "prefilter", value = used_parameters$prefilter),
                              numericInput("value_of_prefilterA", label = "value_of_prefilter", value = used_parameters$value_of_prefilter),
                              selectInput("mzCenterFunA", label = "mzCenterFun", choices = list("wMean" = 0, "mean" = 1, "apex" = 2, "wMeanApex3" = 3, "meanApex3" = 4), selected = used_parameters$mzCenterFun),
                              numericInput("integrateA", label = "integrate", value = used_parameters$integrate),
                              numericInput("mzdiffA", label = "mzdiff", value = used_parameters$mzdiff),
                              checkboxInput("fitgaussA", label = "fitgauss", value = used_parameters$fitgauss),
                              numericInput("noiseA", label = "noise", value = used_parameters$noise),
                              checkboxInput("verboseColumnsA", label = "verboseColumns", value = used_parameters$verboseColumns),
                              numericInput("criticalValueA", label = "criticalValue", value = used_parameters$criticalValue),
                              numericInput("consecMissedLimitA", label = "consecMissedLimit", value = used_parameters$consecMissedLimit),
                              numericInput("unionsA", label = "unions", value = used_parameters$unions),
                              numericInput("checkBackA", label = "checkBack", value = used_parameters$checkBack),
                              checkboxInput("withWaveA", label = "withWave", value = used_parameters$withWave)
                            )
                          }
                          else if (used_parameters$Peak_method == 2 & input$Peak_methodA == 2) {
                            peak_options <- tagList(
                              numericInput("binSizeA", label = "binSize", value = used_parameters$binSize),
                              selectInput("imputeA", label = "impute", choices = list("none" = 0, "lin" = 1, "linbase" = 2, "intlin" = 3, "imputeLinInterpol" = 4), selected = used_parameters$impute),
                              # numericInput("baseValue", label = "baseValue", 5),
                              # numericInput("distance", label = "distance", 5),
                              numericInput("fwhmA", label = "fwhm", value = used_parameters$fwhm),
                              numericInput("sigmaA", label = "sigma", value = used_parameters$sigma),
                              numericInput("maxA", label = "max", value = used_parameters$max),
                              numericInput("snthreshA", label = "snthresh", value = used_parameters$snthresh),
                              numericInput("stepsA", label = "steps", value = used_parameters$steps),
                              numericInput("mzdiffA", label = "mzdiff", value = used_parameters$mzdiff),
                              checkboxInput("indexA", label = "index", value = used_parameters$index)
                            )
                          }
                          else if (used_parameters$Peak_method == 3 & input$Peak_methodA == 3) {
                            peak_options <- tagList(
                              numericInput("noise", label = "noise", used_parameters$noise),
                              numericInput("rsd_threshold", "threshold", used_parameters$rsd_threshold),
                              numericInput("simthresh", "simthresh", used_parameters$simthresh)
                            )
                          }
                        }
                        peak_options
                      }),
                      actionButton("resub1", label = "\n Rerun")
  	             ),
  	             column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  	                    h5("2. Refinement"),
  	                    # Render utilized refinement options
  	                    selectInput("Ref_methodA", label = "Peak refinement method", choices = list("MergeNeighboringPeaks" = 0, "FilterIntensity" = 1), selected = used_parameters$Ref_method),
  	                    renderUI(ref_options),
  	                    actionButton("resub2", label = "\n Rerun")
  	             ),
  	             column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  	                    h5("3. Alignment"),
  	                    # Render utilized alignment options
  	                    selectInput("Align_methodA", label = "Peak alignment method", choices = list("Obiwarp" = 0, "PeakGroups" = 1), selected = used_parameters$Align_method),
  	                    renderUI(aln_options),
  	                    actionButton("resub3", label = "\n Rerun")
  	             ),
  	             column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  	                    h5("4. Grouping"),
  	                    # Render utilized grouping options
  	                    selectInput("Group_methodA", label = "Peak grouping method", choices = list("PeakDensity" = 0, "MzClust" = 1, "NearestPeaks" = 2), selected = used_parameters$Group_method),
  	                    renderUI(grp_options),
  	                    actionButton("resub4", label = "\n Rerun")
  	             ),
  	             column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  	                    h5("5. Filling"),
  	                    # Render utilized filling options
  	                    # numericInput("expandMzA", label = "expandMz", 0, value = used_parameters$expandMz),
  	                    # numericInput("expandRtA", label = "expandRt", 0, value = used_parameters$expandRt),
  	                    numericInput("fixedMzA", label = "fixedMz", 0, value = used_parameters$fixedMz),
  	                    numericInput("fixedRtA", label = "fixedRt", 0, value = used_parameters$fixedRt),
  	                    actionButton("resub5", label = "\n Rerun")
  	             )
  	    )
	  ))
	  # Insert analysis taglist to the front end
	  output$job_analysis <- analysisTL
	})
	
	#################
	###########       ASYNC FUNCTIONS
	# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
	#Function to process the selected MS-files
  xcms_data_processing <- function(massfiles, parameters, preset, job_id, db_usr, db_pw, job_plan){
    library(xcms)
    register(bpstart(MulticoreParam(8)))
    ##################### Import RMySQL query functions
    send_query <- function(query){
      sqlconn <- dbConnect(
        drv = RMySQL::MySQL(),
        dbname='boloa',
        host="127.0.0.1",
        port=3306,
        user=db_usr,
        password=db_pwd
      )
      on.exit(dbDisconnect(sqlconn))
      dbSendQuery(sqlconn, query)
    }
    insert_query <- function(tb_name, data){
      sqlconn <- dbConnect(
        drv = RMySQL::MySQL(),
        dbname='boloa',
        host="127.0.0.1",
        port=3306,
        user=db_usr,
        password=db_pwd
      )
      on.exit(dbDisconnect(sqlconn))
      dbWriteTable(sqlconn, tb_name, data, append = TRUE, row.names = FALSE)
    }
    get_query <- function(query){
      sqlconn <- dbConnect(
        drv = RMySQL::MySQL(),
        dbname='boloa',
        host="127.0.0.1",
        port=3306,
        user=db_usr,
        password=db_pwd
      )
      on.exit(dbDisconnect(sqlconn))
      dbGetQuery(sqlconn, query)
    }
    ################
    tryCatch(
      {
        def_params <- parameters
        files <- massfiles$file_path
        if (preset == 3) {
          # optiLCMS Parameter optimization (depricated)
          param_initial <- SetPeakParam(platform = "general")
          jstep <- 1
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '1/9 Performing ROI extraction...' WHERE job_id = ", job_id, ";", sep = "")))
          raw_train <- PerformROIExtraction(files, rt.idx = 0.2, rmConts = FALSE, plot = FALSE)
          jstep <- 2
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '2/9 Performing parameter optimization...' WHERE job_id = ", job_id, ";", sep = "")))
          def_params <- PerformParamsOptimization(raw_train, param = param_initial, ncore = 5)
          send_query(stringr::str_glue(paste("UPDATE parameter SET min_peakwidth = ", toString(def_params$min_peakwidth), ", max_peakwidth = ", toString(def_params$max_peakwidth), ", mzdiff = ", toString(def_params$mzdiff), ", snthresh = ", toString(def_params$snthresh), ", bw = ", toString(def_params$bw), ", Peak_method = '", toString(def_params$Peak_method), "', ppm = ", toString(def_params$ppm), ", noise = ", toString(def_params$noise), ", prefilter = ", toString(def_params$prefilter), ", value_of_prefilter = ", toString(def_params$value_of_prefilter), ", minFraction = ", toString(def_params$minFraction), ", minSamples = ", toString(def_params$minSamples), ", maxFeatures = ", toString(def_params$maxFeatures), ", fitgauss = ", toString(def_params$fitgauss), ", mzCenterFun = '", toString(def_params$mzCenterFun), "', integrate = ", toString(def_params$integrate), ", extra = ", toString(def_params$extra), ", span = ", toString(def_params$span), ", smooth = '", toString(def_params$smooth), "', family = '", toString(def_params$family), "', polarity = '", toString(def_params$polarity), "', perc_fwhm = ", toString(def_params$perc_fwhm), ", max_charge = ", toString(def_params$max_charge), ", max_iso = ", toString(def_params$max_iso), ", corr_eic_th = ", toString(def_params$corr_eic_th), ", mz_abs_add = ", toString(def_params$mz_abs_add), ", rmConts = ", toString(def_params$rmConts), ", RT_method = '", toString(def_params$RT_method), "' WHERE job_id = ", job_id, ";", sep = "")))
        }
        else {
          def_params[is.na(def_params)] <- 0
        }
        # Processing pipeline
        dir <- getwd()
        dir <- paste(dir, "/processed_data", sep = "")
        if (job_plan == 7){
          jstep <- 3
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '3/9 Reading MS data...' WHERE job_id = ", job_id, ";", sep = "")))
          # Load raw files into XCMSnExp-object
          pd <- data.frame(sample_name = massfiles$original_file_name,
                           sample_group = massfiles$metadata,
                           stringsAsFactors = FALSE)
          xset <- readMSData(files = files, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")
          # Select user-specified RT-range
          xset <- filterRt(xset, c(as.double(def_params$rtrmin), as.double(def_params$rtrmax)))
        } else {
          # Whenever the job is rerun, choose previously initialized xset
          rdaquer <- stringr::str_glue(paste("SELECT * FROM processed_sample WHERE job_id = ", job_id, ";", sep = ""))
          rda_path <- get_query(rdaquer)
          load(file=toString(rda_path$file_path_rda))
        }
        if (job_plan == 2 | job_plan == 7){
          # Functions utilized for peak detection
          jstep <- 4
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '4/9 Finding peaks...' WHERE job_id = ", job_id, ";", sep = "")))
          
          # centWave
          if (as.integer(def_params$Peak_method) == 0){
            peak_params <- CentWaveParam(
              ppm = as.double(def_params$ppm),
              peakwidth = c(as.double(def_params$min_peakwidth), as.double(def_params$max_peakwidth)),
              snthresh = as.double(def_params$snthresh),
              prefilter = c(as.double(def_params$prefilter), as.double(def_params$value_of_prefilter)),
              mzCenterFun = c("wMean", "mean", "apex", "wMeanApex3", "meanApex3")[as.integer(def_params$mzCenterFun) + 1],
              integrate = as.double(def_params$integrate),
              mzdiff = as.double(def_params$mzdiff),
              fitgauss = as.logical(def_params$fitgauss)  %>% replace(is.na(.), FALSE),
              noise = as.double(def_params$noise),
              verboseColumns = as.logical(def_params$verboseColumns) %>% replace(is.na(.), FALSE)
              #roiList = list(),
              #firstBaselineCheck = TRUE,
              #roiScales = numeric(),
              #extendLengthMSW = FALSE
            )
          }
          # Massifquant
          else if (as.integer(def_params$Peak_method) == 1){
            peak_params <- MassifquantParam(
              ppm = as.double(def_params$ppm),
              peakwidth = c(as.double(def_params$min_peakwidth), as.double(def_params$max_peakwidth)),
              snthresh = as.double(def_params$snthresh),
              prefilter = c(as.double(def_params$prefilter), as.double(def_params$value_of_prefilter)),
              mzCenterFun = c("wMean", "mean", "apex", "wMeanApex3", "meanApex3")[as.integer(def_params$mzCenterFun) + 1],
              integrate = as.double(def_params$integrate),
              mzdiff = as.double(def_params$mzdiff),
              fitgauss = as.logical(def_params$fitgauss),
              noise = as.double(def_params$noise),
              verboseColumns = as.logical(def_params$verboseColumns),
              criticalValue = as.double(def_params$criticalValue),
              consecMissedLimit = as.double(def_params$consecMissedLimit),
              unions = as.double(def_params$unions),
              checkBack = as.double(def_params$checkBack),
              withWave = as.logical(def_params$withWave)
            )
          }
          # MatchedFilter
          else if (as.integer(def_params$Peak_method) == 2){
            peak_params <- MatchedFilterParam(
              binSize = as.double(def_params$binSize),
              impute = c("none", "lin", "linbase", "intlin", "imputeLinInterpol")[as.integer(def_params$impute) + 1],
              baseValue = as.double(def_params$baseValue),
              distance = as.double(def_params$distance),
              fwhm = as.double(def_params$fwhm),
              sigma = as.double(def_params$fwhm)/as.double(def_params$sigma),
              max = as.double(def_params$max),
              snthresh = as.double(def_params$snthresh),
              steps = as.double(def_params$steps),
              mzdiff = as.double(def_params$mzdiff) - as.double(def_params$binSize) * as.double(def_params$steps),
              index = as.logical(def_params$index)
            )
          }
          ## New method for peak detection, square difference
          else if (as.integer(def_params$Peak_method) == 3){
            def_params$ppm <- 40 # Necessary for other steps! Does not utilize ppm
            xset <- as(xset, "XCMSnExp")
            for (sampleN in 1:length(files)) {
              raw_data <- readMSData(files = files[sampleN], pdata = new("NAnnotatedDataFrame"), mode = "onDisk")
              raw_data <- filterRt(raw_data, c(as.double(def_params$rtrmin), as.double(def_params$rtrmax)))
              # Subtract squared tic from original
              altered_tic <- tic(raw_data) - sqrt(tic(raw_data))
              # Calcualte of difference between subtracted tic and original tic
              diffalt <- (altered_tic - tic(raw_data))/tic(raw_data)*100
              # Create a local linear regression using supsmu
              trend <- supsmu(rtime(raw_data), diffalt, bass = -100)
              # Calculation of threshold based on a local regression model
              trend$y <- trend$y - mean(diffalt) / as.double(def_params$rsd_threshold)
              newvals <- tic(raw_data)
              # Remove scans below threshold
              newvals[diffalt < trend$y] <- NA
              # Initialisation of vectors for 'manual' peak assignment
              minrts <- c()
              maxrts <- c()
              minmzs <- c()
              maxmzs <- c()
              simthresh <- as.double(def_params$simthresh) # Cosine similarty threshold for new peak splitting
              peakmzs <- c()
              sim <- 0
              sampleSpec <- spectra(raw_data)
              # Loop through the peak TIC in order to retieve information and separate overlapping peaks.
              for (i in 2:length(newvals)) {
                # If the current is above and the previous is below threshold:
                # Initialise new peak
                if (!is.na(newvals[i]) & is.na(newvals[i - 1])) {
                  minrts <- c(minrts, rtime(raw_data)[i])
                  scan_info <- intensity(sampleSpec[[i]])
                  names(scan_info) <- mz(sampleSpec[[i]])
                  scan_info <- scan_info[scan_info > as.double(def_params$noise)]
                  peakmzs <- append(peakmzs, as.double(names(scan_info)))
                }
                # If the current is below and the previous is above threshold or if both scans are above but similarity is below
                # End of peak. Wrap up values.
                if ((is.na(newvals[i]) & !is.na(newvals[i - 1])) | (sim < simthresh & !is.na(newvals[i]) & !is.na(newvals[i - 1]))) {
                  maxrts <- c(maxrts, rtime(raw_data)[i - 1])
                  minmzs <- c(minmzs, min(peakmzs))
                  maxmzs <- c(maxmzs, max(peakmzs))
                  peakmzs <- c()
                }
                # If the current and previous scans are above threshold
                # Append values to peak or initialise new peak with low cosim.
                if (!is.na(newvals[i]) & !is.na(newvals[i - 1])) {
                  # If the previously calculated similarity was below simtrhesh: Initialise new peak
                  if (sim < simthresh) {
                    minrts <- c(minrts, rtime(raw_data)[i])
                    scan_info <- intensity(sampleSpec[[i]])
                    names(scan_info) <- mz(sampleSpec[[i]])
                    scan_info <- scan_info[scan_info > as.double(def_params$noise)]
                    peakmzs <- append(peakmzs, as.double(names(scan_info)))
                  }
                  # Retrieve information about the current and previous scan for comparison
                  scan_info <- intensity(sampleSpec[[i]])
                  names(scan_info) <- mz(sampleSpec[[i]])
                  scan_info <- scan_info[scan_info > as.double(def_params$noise)]
                  peakmzs <- append(peakmzs, as.double(names(scan_info)))
                  scan_info_prev <- intensity(sampleSpec[[i - 1]])
                  names(scan_info_prev) <- mz(sampleSpec[[i - 1]])
                  scan_info_prev <- scan_info_prev[scan_info_prev > as.double(def_params$noise)]
                  
                  # Pad x and y to make lengths equal
                  x <- as.double(names(scan_info))
                  y <- as.double(names(scan_info_prev))
                  x <- c(x, rep(0, abs(length(y) - length(x))))
                  y <- c(y, rep(0, abs(length(x) - length(y))))
                  
                  # Calculate cosine similarity (composite score)
                  sim <- 999*(sum(sqrt(x) * sqrt(y))^2)/(sum(x)*sum(y))
                }
              }
              # Add peaks to xset object
              new_ranges <- matrix(c(minmzs, maxmzs, minrts, maxrts), ncol = 4, nrow = length(maxrts), byrow = FALSE)
              colnames(new_ranges) <- c('mzmin', 'mzmax', 'rtmin', 'rtmax')
              xset <- manualChromPeaks(xset, new_ranges, samples = sampleN)
            }
          }
          else {
            return()
          }
          if (as.integer(def_params$Peak_method) != 3){
            # Perform peak detection with xcms-methods
            print(peak_params)
            xset <-findChromPeaks(xset, param = peak_params)
          }
        }
        if (job_plan %in% c(2, 3) | job_plan == 7){
          # Functions utilized for peak refinement
          jstep <- 5
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '5/9 Refining peaks...' WHERE job_id = ", job_id, ";", sep = "")))
          # When optimization was chosen, utilize Merge Neighboring Peaks
          if (preset == 3) {
            refmeth <- MergeNeighboringPeaksParam()
          }
            else {
              # MergeNeighboringPeaks
            if (as.integer(def_params$Ref_method) == 0){
                refmeth <- MergeNeighboringPeaksParam(
                                                  expandRt = as.double(def_params$expandRt),
                                                  expandMz = as.double(def_params$expandMz),
                                                  ppm = as.double(def_params$ppm),
                                                  minProp = as.double(def_params$minProp)
                                                  )
              }
            else if (as.integer(def_params$Ref_method) == 1) {
              # FilterIntensity
              refmeth <-  FilterIntensityParam(
                  threshold = as.double(def_params$threshold),
                  nValues = 1L,
                  value = "maxo"
                  )
            }
            }
          print(refmeth)
          # Perform peak refinement
          xset <- refineChromPeaks(xset, refmeth)
        }
        if (job_plan %in% c(2,3,4) | job_plan == 7){
          # Functions for peak alignment
          jstep <- 6
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '6/9 Aligning peaks...' WHERE job_id = ", job_id, ";", sep = "")))
          if (preset == 3) {
            align_param <- ObiwarpParam()
          }
          else {
            # Obiwarp
            if (as.integer(def_params$Align_method) == 0){
              align_param <- ObiwarpParam(
                binSize = as.double(def_params$binSize),
                # centerSample = numeric(),
                response = 1L,
                distFun = c("cor", "cor_opt", "cov", "prd", "euc")[as.integer(def_params$distFun) + 1],
                # gapInit = numeric(),
                # gapExtend = numeric(),
                factorDiag = as.double(def_params$factorDiag),
                factorGap = as.double(def_params$factorGap),
                localAlignment = as.logical(def_params$localAlignment),
                initPenalty = as.double(def_params$initPenalty),
                # subset = integer(),
                subsetAdjust = c("average", "previous")[as.integer(def_params$subsetAdjust) + 1]
              )
            }
            # PeakGroups
            else if (as.integer(def_params$Align_method) == 1){
              align_param <- PeakGroupsParam(
                minFraction = as.double(def_params$minFraction),
                extraPeaks = as.double(def_params$extraPeaks),
                smooth = c("loess", "linear")[as.integer(def_params$smooth) + 1],
                span = as.double(def_params$span),
                family = c("gaussian", "symmetric")[as.integer(def_params$family) + 1],
                # peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
                # subset = integer(),
                subsetAdjust = c("average", "previous")[as.integer(def_params$subsetAdjust) + 1]
              )
            }
          }
          print(align_param)
          # Perform peak alignment
          xset <- adjustRtime(xset, param = align_param)
        }
        if (job_plan %in% c(2,3,4,5) | job_plan == 7){
          #Perform peak grouping
          jstep <- 7
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '7/9 Grouping peaks...' WHERE job_id = ", job_id, ";", sep = "")))
          if (preset == 3) {
            peak_group_param <- PeakDensityParam()
          }
          else {
            # PeakDensity
            if (as.integer(def_params$Group_method) == 0){
              peak_group_param <- PeakDensityParam(
                sampleGroups = massfiles$metadata,
                bw = as.double(def_params$bw),
                minFraction = as.double(def_params$minFraction),
                minSamples = as.double(def_params$minSamples),
                binSize = as.double(def_params$binSize),
                maxFeatures = as.double(def_params$maxFeatures)
              )
            }
            # MzClust
            else if (as.integer(def_params$Group_method) == 1){
              peak_group_param <- MzClustParam(
                sampleGroups = massfiles$metadata,
                ppm = as.double(def_params$ppm),
                absMz = as.double(def_params$absMz),
                minFraction = as.double(def_params$minFraction),
                minSamples = as.double(def_params$minSamples)
              )
            }
            # NearestPeaks
            else if (as.integer(def_params$Group_method) == 2){
              peak_group_param <- NearestPeaksParam(
                sampleGroups = massfiles$metadata,
                mzVsRtBalance = as.double(def_params$mzVsRtBalance),
                absMz = as.double(def_params$absMz),
                absRt = as.double(def_params$absRt),
                kNN = as.double(def_params$kNN)
              )
            }
          }
          print(peak_group_param)
          # Perform peak grouping
          xset <- groupChromPeaks(xset, param = peak_group_param)
        }
        if (job_plan %in% c(2,3,4,5,6) | job_plan == 7){
          # Peak filling functions
          jstep <- 8
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '8/9 Filling peaks...' WHERE job_id = ", job_id, ";", sep = "")))
          # fillparam <- ChromPeakAreaParam(
          #     # mzmin = function(z) quantile(z, probs = 0.25),
          #     # mzmax = function(z) quantile(z, probs = 0.75),
          #     # rtmin = function(z) quantile(z, probs = 0.25),
          #     # rtmax = function(z) quantile(z, probs = 0.75)
          #   )
          if (preset == 3) {
            fillparam <- FillChromPeaksParam()
          }
          # FillChromPeaks
          else {
            fillparam <- FillChromPeaksParam(
              expandMz = as.double(def_params$expandMz),
              expandRt = as.double(def_params$expandRt),
              ppm = 10,
              fixedMz = as.double(def_params$fixedMz),
              fixedRt = as.double(def_params$fixedRt)
            )
          }
          print(fillparam)
          # Perform peak filling
          xset <- fillChromPeaks(xset, param = fillparam)
          # Save xset to .rda for further alterations and actions
          save(xset, file = paste(dir, "/", job_id, ".rda", sep = ""))
        }
        if (job_plan %in% c(2,3,4,5,6) | job_plan == 7){
          # Functions for peak annotations
          jstep <- 9
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '9/9 Annotating peaks...' WHERE job_id = ", job_id, ";", sep = "")))
          # Remove previous annotation for same job, necessary for reruns.
          send_query(stringr::str_glue(paste("DELETE FROM peak WHERE job_id = ", job_id, ";", sep = "")))
          annot_spectra <- featureSpectra(xset, msLevel = 1, return.type = "Spectra", skipFilled = TRUE)
          save(annot_spectra, file = paste(dir, "/", job_id, "_annot.rda", sep = ""))
          chromtype <- c("LC", "GC")[as.integer(def_params$sample_type)]
          # Store information about the detected peaks
          tt <- spectraData(annot_spectra)
          # Loop through all the indexes of filenames used in the job.
          for (i in 1:length(files)) {
            # Split the xset into a single file.
            aa <- filterFile(xset, i)
            is_filled <-chromPeakData(aa)[,2]
            tp <- rownames(chromPeaks(aa))[!is_filled]
            # Take various columns from the tt variable.
            p <- unique(tt[tt[, "peak_id"] %in% tp, c("peak_id", "scanIndex", "rtime", "totIonCurrent")])
            # Loop through each peak detected in a single file in order to retrieve information
            print(length(tp))
            for (peak in 1:length(tp)) {
              # Retrieve information about a single peak
              p <- unique(tt[tt[, "peak_id"] == tp[peak], c("peak_id", "scanIndex", "rtime", "totIonCurrent")])
              if (nrow(p) >= 1) {
                apex <- p[order(p[,"totIonCurrent"], decreasing = TRUE ),][1,] #Scan number of peak apex
                if (apex[,"scanIndex"] > length(annot_spectra)) {
                  next
                }
                apexspectra <- peaksData(annot_spectra[apex[,"scanIndex"]])[[1]]
                # Subtract noise from intensities
                apexspectra[,2] <- apexspectra[,2] - as.double(def_params$noise)
                # Keep only positive intensities
                apexspectra <- apexspectra[apexspectra[,2] > 0,]
                # Calculate relative abundance for intensity:
                #   Divide all intensities by maximum intensity.
                apexspectra[,2] <- (apexspectra[,2]/max(apexspectra[,2])) * 100
                splash <- getSplash(apexspectra) #Splash generated from peak apex
                search_splash <- paste(strsplit(splash , split = "-")[[1]][1:3], collapse='-')
                # calculate the pearson correlation between mz and intensity
                
                spectrum_formatted <- c()
                pspec <- apexspectra
                colnames(pspec) <- NULL
                # Format data to [mz]:[intensity]
                rownames(pspec) <- NULL
                for (specrow in 1:nrow(pspec)) {
                  spectrum_formatted <- c(spectrum_formatted, paste(pspec[specrow,], collapse = ':'))
                }
                spectrum_formatted <- paste(spectrum_formatted, collapse = " ")
                # Molecule annotation !DEPRECATED! NOT MY RESPONSIBILITY, FRAMEWORK IS HERE
                # Splash filtering and modified cosine similarity classification: METHOD 1
                mol_annot_splash <- stringr::str_glue(paste("SELECT * FROM mol WHERE type = '", chromtype, "' AND (splash LIKE '%", search_splash, "%');", sep = ""))
                mol_annot_splash <- get_query(mol_annot_splash)
                if (nrow(mol_annot_splash) == 0) {
                  highest_cosim_mol <- "NULL"
                  highest_cosim <- "NULL"
                } else {
                  # Iterate through all matched compounds
                  scan1 <- apexspectra
                  for (mspectr in 1:length(mol_annot_splash)) {
                    scan2 <- mol_annot_splash$spectrum[mspectr]
                    # Split the scan2 string into mz and intensity values
                    bins <- strsplit(scan2, " ")
                    ints <- c()
                    mzs <- c()
                    for (section in bins[[1]]){
                      mzs <- c(mzs, as.double(strsplit(section, ":")[[1]][1]))
                      ints <- c(ints, as.double(strsplit(section, ":")[[1]][2]))
                    }
                    scan2 <- data.frame(mzs, ints)

                    # minlength <- min(nrow(scan1), nrow(scan2))
                    x <- ((scan1[,2] / max(scan1[,2])) * 100)
                    y <- ((scan2[,2] / max(scan2[,2])) * 100)
                    # scan1in2 <- round(scan1[,1]) %in% scan2[,1]
                    # scan2in1 <- scan2[,1] %in% round(scan1[,1])
                    # print(scan1)
                    # print(scan2)
                    # print(scan1in2)
                    # print(scan2in1)
                    # print(length(scan1in2))
                    # print(length(scan2in1))
                    # This will align the relative abundances of m/z bins by
                    # maximum value, which is always 100.
                    # Example:
                    # x = m m m M m m
                    # y = m m M m
                    # --
                    # x = m m m M m m
                    # y = 0 m m M m 0
                    if (!any(is.na(x)) & !any(is.na(y))){
                      while (which.max(x) != which.max(y)) {
                        max_indexes <- c(which.max(x), which.max(y))
                        min_max <- which.min(max_indexes)
                        if (min_max == 1) {
                          x <- c(0, x)
                        }
                        else if (min_max == 2) {
                          y <- c(0, y)
                        }
                      }
                      x <- c(x, rep(0, abs(length(y) - length(x))))
                      y <- c(y, rep(0, abs(length(x) - length(y))))
                      # + 1 to remove 0 values and keep ratio's similar
                      x <- x + 1
                      y <- y + 1
                      # Keep only similar mz values between spectra
                      # x <- x[scan1in2]
                      # y <- x[scan2in1]
                      # Calculate the cosine similarity between peak and compound
                      sim <- 999*(sum(sqrt(x) * sqrt(y))^2)/(sum(x)*sum(y))
                    }
                    else {
                      sim <- 0
                    }
                    mol_annot_splash[mspectr,"sim"] <- sim
                  }
                  highest_cosim <- mol_annot_splash[order(mol_annot_splash$sim, decreasing = TRUE),][1,]
                  highest_cosim_mol <- highest_cosim$mol_id
                  highest_cosim <- highest_cosim$sim
                  if (highest_cosim < 900) { # If highest is below threshold
                    highest_cosim_mol <- "NULL"
                    highest_cosim <- "NULL"
                  }
                }

                # Coeff similarity search METHOD: 2
                # v1 <- apexspectra[,1]
                # v2 <- apexspectra[,2]
                # coeff <- sum(v1 * v2)/sqrt(sum(v1^2)*sum(v2^2))
                # # Calculate range limits for confidence intervals
                # corr_UL_range <- 0.05 #Edit this value for a preferred limit
                # max_diff <- coeff * corr_UL_range
                #
                #
                # if (coeff >= 0) {
                #   mol_annot <- stringr::str_glue(paste("SELECT * FROM mol WHERE type = '", chromtype, "' AND coeff < ", sub("−", "-", paste(coeff + max_diff)), " AND coeff > ", sub("−", "-", paste(coeff - max_diff)), ";", sep = ""))
                #   mol_annot <- get_query(mol_annot)
                # }
                # else {
                #   mol_annot <- stringr::str_glue(paste("SELECT * FROM mol WHERE type = '", chromtype, "' AND coeff < ", sub("−", "-", paste(coeff - max_diff)), " AND coeff > ", sub("−", "-", paste(coeff + max_diff)), ";", sep = ""))
                #   mol_annot <- get_query(mol_annot)
                # }
                # # If no matching compounds are available in the database, annotation is not possible
                # if (nrow(mol_annot) == 0) {
                #   mol_id_coeffsim <- "NULL"
                #   coeff_diff <- "NULL"
                # } else {
                #   # Calculate the difference between correlation of peak and matched compounds in database
                #   mol_annot$diff <- abs((coeff - mol_annot$coeff)/mol_annot$coeff)
                #   best_match <- mol_annot[order(mol_annot$diff),][1,]
                #   if (best_match$diff > corr_UL_range) { # Outdated
                #     mol_id_coeffsim <- "NULL"
                #     coeff_diff <- "NULL"
                #   }
                #   else {
                #     mol_id_coeffsim <- best_match$mol_id
                #     coeff_diff <- best_match$diff
                #   }
                # }
                ##### TESTING
                # highest_cosim_mol <- "NULL"
                # highest_cosim <- "NULL"
                # 
                # 
                mol_id_coeffsim <- "NULL"
                coeff_diff <- "NULL"
                coeff <- "NULL"
                chromaa <- as.data.frame(chromPeaks(aa))
                send_query(stringr::str_glue(paste("INSERT INTO peak VALUES (",
                                                   job_id, ", '", #job_id
                                                   tp[peak], "', '", #peak_id
                                                   massfiles$sample_hash[i], "', ", #sample_hash
                                                   sum(p[,"totIonCurrent"]), ", ", #peak_area
                                                   mol_id_coeffsim, ", ", #mol_id_coeffsim
                                                   highest_cosim_mol, ", ", #mol_id_modcosinesim
                                                   highest_cosim, ", ", #modcosinesim
                                                   coeff_diff, ", ", #coeff_diff
                                                   coeff, ", ", #coeff
                                                   chromaa[peak,"mz"], ", ", #mz
                                                   chromaa[peak,"mzmin"], ", ", #mzmin
                                                   chromaa[peak,"mzmax"], ", ", #mzmax
                                                   chromaa[peak,"rt"], ", ", #rt
                                                   chromaa[peak,"rtmin"], ", ", #rtmin
                                                   chromaa[peak,"rtmax"], ", ", #rtmax
                                                   apex[,"totIonCurrent"], ", '", #apex_tic
                                                   splash, "', '", #splash key
                                                   spectrum_formatted, "'", #spectrum
                                                   ");", sep = "")))
              }
            }
          }
        }
        if (job_plan == 7) {
          # Job finished for the first time
          todf <- data.frame(
            job_id = toString(job_id),
            file_path_rda = toString(paste(dir, "/", job_id, ".rda", sep = "")),
            file_path_peaks = toString(paste(dir, "/", job_id, "_annot.rda", sep = ""))
          )
          # Insert processed job to DB
          insert_query("processed_sample", todf)
        }
        # Retrieve end time of job.
        end_time <- format(Sys.time() + 60*60, "%Y-%m-%d %X")
        send_query(stringr::str_glue(paste("UPDATE job SET job_status = 'Finished' WHERE job_id = ", job_id, ";", sep = "")))
        send_query(stringr::str_glue(paste("UPDATE job SET end_time = '", end_time, "' WHERE job_id = ", job_id, ";", sep = "")))
      },
      error = function(cnd){
        # Catch errors and set status to "CRASHED"
        print(cnd)
        send_query(stringr::str_glue(paste("UPDATE job SET job_status = 'CRASHED AT ", jstep, "/9', end_time = '", format(Sys.time() + 60*60, "%Y-%m-%d %X"), "' WHERE job_id = ", job_id, ";", sep = "")))
        return(NA)
      }
    )
  }
}

runApp(shinyApp(ui = ui, server = server), host = hostip, port = portnr)
