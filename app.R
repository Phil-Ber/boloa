#!/bin/Rscript

load('.hidden.RData')

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

#if (!require("BiocManager", quietly = TRUE))
#BiocManager::install("OptiLCMS")
tryCatch(
	expr = {
		library("MetaboAnalystR")
	  library("OptiLCMS")
	  #library(shinyvalidate)
	},
	error = {
	#   install.packages("BiocManager")
	# 	#pacman::p_load(c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR"), character.only = TRUE)
		pacman::p_load(c("Rserve", "ellipse", "scatterplot3d", "Cairo", "randomForest", "caTools", "e1071", "som", "impute", "pcaMethods", "RJSONIO", "ROCR", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "pheatmap", "sva", "Rcpp", "pROC", "data.table", "limma", "car", "fitdistrplus", "lars", "Hmisc", "magrittr", "methods", "xtable", "pls", "caret", "lattice", "igraph", "gplots", "KEGGgraph", "reshape", "RColorBrewer", "tibble", "siggenes", "plotly", "xcms", "CAMERA", "fgsea", "MSnbase", "BiocParallel", "multtest", "RBGL", "edgeR", "fgsea", "crmn", "progress", "qs", "glasso"), character.only = TRUE)
	# 	BiocManager::install("mtbls2")
	  # install_github("berlinguyinca/spectra-hash", subdir="splashR")
	#   #install.packages("shinyvalidate", repos="https://mirrors.evoluso.com/CRAN/")
	# 	devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE)
	# 	devtools::install_github("xia-lab/OptiLCMS", build = TRUE, build_vignettes = FALSE, build_manual =TRUE)
	},
	finally = {
		library("MetaboAnalystR")
	  library("OptiLCMS")
	  #library(shinyvalidate)
	}
)
# install_github("berlinguyinca/spectra-hash", subdir="splashR")
library(splashR)
#print(warnings())
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

plan(multisession)

options(shiny.maxRequestSize=10000*1024^2)
hostip <- "145.97.18.149"
portnr <- 7123


# all_params <- c("Peak_method", "RT_method", "mzdiff", "snthresh", "bw", "ppm", "min_peakwidth", "max_peakwidth", "noise", "prefilter", "value_of_prefilter", "minFraction", "minSamples", "maxFeatures", "mzCenterFun", "integrate", "extra", "span", "smooth", "family", "fitgauss", "polarity", "perc_fwhm", "mz_abs_iso", "max_charge", "max_iso", "corr_eic_th", "mz_abs_add", "rmConts", "verboseColumns")
all_params <- c("Peak_method", "Ref_method", "Align_method", "Group_method", "absMz", "absRt", "baseValue", "binSize", "bw", "centerSample", "checkBack", "consecMissedLimit", "criticalValue", "distance", "distFun", "expandMz", "expandRt", "extendLengthMSW", "extraPeaks", "factorDiag", "factorGap", "family", "firstBaselineCheck", "fitgauss", "fixedMz", "fixedRt", "fwhm", "gapExtend", "gapInit", "impute", "index", "initPenalty", "integrate", "kNN", "localAlignment", "max", "maxFeatures", "minFraction", "minProp", "minSamples", "mzCenterFun", "mzdiff", "mzVsRtBalance", "ncol", "noise", "nrow", "nValues", "peakGroupsMatrix", "max_peakwidth", "min_peakwidth", "ppm", "prefilter", "value_of_prefilter", "response", "roiList", "roiScales", "sampleGroups", "sigma", "smooth", "snthresh", "span", "steps", "subset", "subsetAdjust", "threshold", "unions", "value", "verboseColumns", "withWave", "rtrange")
lcms_only <- tail(all_params, n=9)

# send_query(stringr::str_glue("SET GLOBAL local_infile=1;"))

ui <- fillPage(
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
						checkboxInput("debug", label = "Debug mode", 0),
						tableOutput("jobval")
					),
					column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
						h5("1. Detection"),
						selectInput("Peak_method", label = "Peak detection method", choices = list("centWave" = 0, "Massifquant" = 1, "MatchedFilter" = 2)),
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
  			# fluidRow(style='max-width:100%;padding:10px;',
  			#          column(12,
  			#                 uiOutput("jobInformation")
  			#                 )
  			#          ),
  			# fluidRow(style='max-width:100%;padding:10px;',
  			#          column(6,
  			#                 div(DT::dataTableOutput("analysisParams"))#, style = "font-size: 75%; width: 25%")
  			#          ),
  			#          column(6,
  			#                 div(DT::dataTableOutput("analysisSamples"))#, style = "font-size: 75%; width: 50%")
  			#                 
  			#          )
  			# ),
  			# fluidRow(style='max-width:100%;padding:10px;',
  			#          column(6,
  			#                 plotOutput("analysisPCA")
  			#          ),
  			#          column(6,
  			#                 plotlyOutput("analysisPeaksTIC")
  			#          )
  			# ),
  			# fluidRow(style='max-width:100%;padding:10px;',
  			#          column(6,
  			#                 uiOutput("analysisCompounds")
  			#          ),
  			#          column(6,
  			#                 uiOutput("analysisDifferential")
  			#          )
  			# ),
  			# fluidRow(style='max-width:100%;padding:10px;',
  			#          column(6,
  			#                 uiOutput("analysisVenn")
  			#          ),
  			#          column(6,
  			#                 uiOutput("analysisHeatmap")
  			#          )
  			# )#,
  			# fluidRow(style='max-width:100%;padding:10px;',
  			#   column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  			#          h5("1. Detection"),
  			#          uiOutput("used_peak_parameters")
  			#   ),
  			#   column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  			#          h5("2. Refinement"),
  			#          uiOutput("used_refinement_parameters")
  			#   ),
  			#   column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  			#          h5("3. Alignment"),
  			#          uiOutput("used_alignment_parameters")
  			#   ),
  			#   column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  			#          h5("4. Grouping"),
  			#          uiOutput("used_grouping_parameters")
  			#   ),
  			#   column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  			#          h5("5. Filling"),
  			#          uiOutput("used_grouping_parameters"),
  			#   ),
  			#   column(2)
  			# )
		)
	)
)
)


server <- function(input, output, session) {
  nr_files <<- 0
  
  default_cent <- tagList(
    numericInput("ppm", label = "ppm", 32.9),
    numericInput("min_peakwidth", label = "min_peakwidth", 3.783993),
    numericInput("max_peakwidth", label = "max_peakwidth", 52.19042),
    numericInput("snthresh", label = "snthresh", 72.97847),
    numericInput("prefilter", label = "prefilter", 2),
    numericInput("value_of_prefilter", label = "value_of_prefilter", 247.2712),
    selectInput("mzCenterFun", label = "mzCenterFun", choices = list("wMean" = 0, "mean" = 1, "apex" = 2, "wMeanApex3" = 3, "meanApex3" = 4), selected = 4),
    numericInput("integrate", label = "integrate", 1),
    numericInput("mzdiff", label = "mzdiff", 0.6590461),
    checkboxInput("fitgauss", label = "fitgauss", 0),
    numericInput("noise", label = "noise", 8778.373),
    checkboxInput("verboseColumns", label = "verboseColumns", 0)
  )
  default_massif <- tagList(
    numericInput("ppm", label = "ppm", 93.07),
    numericInput("min_peakwidth", label = "min_peakwidth", 5),
    numericInput("max_peakwidth", label = "max_peakwidth", 30),
    numericInput("snthresh", label = "snthresh", 10),
    numericInput("prefilter", label = "prefilter", 2),
    numericInput("value_of_prefilter", label = "value_of_prefilter", 328.63),
    selectInput("mzCenterFun", label = "mzCenterFun", choices = list("wMean" = 0, "mean" = 1, "apex" = 2, "wMeanApex3" = 3, "meanApex3" = 4)),
    numericInput("integrate", label = "integrate", 1),
    numericInput("mzdiff", label = "mzdiff", 0.01),
    checkboxInput("fitgauss", label = "fitgauss", 0),
    numericInput("noise", label = "noise", 0),
    checkboxInput("verboseColumns", label = "verboseColumns", 0),
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
    numericInput("expandRt", label = "expandRt", 10),
    numericInput("expandMz", label = "expandMz", 10),
    numericInput("minProp", label = "minProp", 1)
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
  
  output$job_analysis <- renderUI(tagList(
    h5('Please select a finished job in the "Jobs" panel to analyse data')
  ))
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
    output$peak_parameters <- renderUI(peak_options)
  })
  
  observeEvent(input$Ref_method, {
    if (input$Ref_method == 0) {
      ref_options <- default_mnp
    }
    else if (input$Ref_method == 1) {
      ref_options <- default_fi
    }
    output$refinement_parameters <- renderUI(ref_options)
  })
   
  observeEvent(input$Align_method, {
    if (input$Align_method == 0) {
      alig_options <- default_obiwarp
    }
    else if (input$Align_method == 1) {
      alig_options <- default_peakgroups
    }
    output$alignment_parameters <- renderUI(alig_options)
  })
  
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
  
	#Weergave van geuploade bestand(en)
	observeEvent(input$msdata, {
		output$fileOverview <- renderUI({
			file <- input$msdata
			ext <- tools::file_ext(file$datapath)
			req(file)
			validate(need(ext %in% c("raw", "mzXML", "CDF"), "Please upload a .raw, .CDF or .mzXML file!"))
			nr_files <<- length(file[,3])
			lapply(1:nr_files, function(i) {
			  if (i == 1) {
			    p('Please edit file metadata (REQUIRED!)')
			  }
			  return(textInput(inputId = paste("meta", i, sep = ""), label = paste(i, ": ", file[i, 1], sep = ""), placeholder = paste("Group sample", i)))
			})
		#}#, server = FALSE)
  	})
  })
  
	#Update van tabel op "Upload data" tab
	updateEvent <- reactive({
		list(input$tabSwitch, input$dataupl)
	})
	
	#Sample selection table
	observeEvent(updateEvent, {
		query <- stringr::str_glue("SELECT * FROM sample ORDER BY upload_date DESC;")
		sample_table_content <<- get_query(query)
		sample_table_content$chromatography_type[sample_table_content$chromatography_type == 1] <- "Liquid"
		sample_table_content$chromatography_type[sample_table_content$chromatography_type == 2] <- "Gas"
		output$uploaded_samples <- DT::renderDataTable({
			DT::datatable(sample_table_content[, c(1, 7, 3, 5, 4, 6)])
		}, server = FALSE)
		query <- stringr::str_glue("SELECT * FROM job ORDER BY start_time DESC;")
		job_table_content <<- get_query(query)
		output$jobs <- DT::renderDataTable({
		  DT::datatable(job_table_content[, c(1, 5, 3, 4, 2)], selection = 'single',
		                rownames= FALSE)
		}, server = FALSE)
	})
	
	#Job table
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
	
	#On click preview button
	observeEvent(input$preview, {
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
		#Render 3d plot
		output$previewplot <- renderUI({
			open3d()
			plot3D(M, rgl = TRUE)#xlim = input$chromrange, ylim = input$mzrange, xlab = "Retention time (seconds)", ylab = "M/Z", zlab = "Intensity")
			rglwidget()
		  
		})
		progress$close()
	})
	
	#Clear plots on tab switch
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
	
	#Handelingen als data ingevoerd is ---
	observeEvent(input$dataupl, {
	  if (nchar(input$sample_description) == 0){
	    output$upl_completed <- renderText({
	      "Please input your name."
	    })
	    return()
	  }
		if (nr_files == 0) {
			output$upl_completed <- renderText({
				'Please upload a file.'
			})
			return()
		}
	  invalid_metadata <- FALSE
	  for (i in 1:nr_files) {
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
	  file_tags <- c()
	  for (i in 1:nr_files){
	    file_tags <- c(file_tags, gsub("[^[:alnum:][:space:]]", "", input[[paste("meta", i, sep = "")]]))
	  }
		dir <- getwd()
		dir <- paste(dir, "/massascans", sep = "")
		#aanmaak metadata object
		file <- input$msdata
		time <- format(Sys.time() + 60*60, "%Y-%m-%d %X")
		#progress <- shiny::Progress$new()
		#on.exit(progress$close())
		#progress$set(message = "Uploading", value = 0)s
		run_file_upload <- function(rfuFile, rfuTime, rfuDir, rfuSample_description, rfuDb_usr, rfuD_pwd, rfuSample_table_content, rfuFile_tags){
		  count <- 0
		  plan(multisession)
		  for (fileloc in rfuFile$datapath) {
		    #progress$inc(1/length(file$datapath), detail = paste("File:", file$name[count + 1]))
		    upload_file <- function(ufCount, ufFileloc, ufDir, ufTime, ufSample_description, ufFile, ufFile_tag, ufDb_usr, ufDb_pwd, ufSample_table_content){
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
		      file.copy(ufFileloc, ufDir)
		      if (grepl('.raw', ufFileloc, fixed=TRUE)) {
		        system(paste("docker run --rm -v ", ufDir, ":/massascans chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /massascans/", ufCount, ".raw", " --mzXML -o /massascans", sep = "")) #change .raw files to .mzXML
		        system(paste("rm ", ufDir, "/", ufCount, ".raw", sep = ""))
		        filetype <- ".mzXML"
		      }
		      else if (grepl('.mzXML', ufFileloc, fixed=TRUE)) {
		        filetype <- ".mzXML"
		      }
		      else if (grepl('.CDF', ufFileloc, fixed=TRUE)) {
		        filetype <- ".CDF"
		      }
		      else {
		        system(paste("rm ", ufFileloc, sep = ""))
		      }
		      hash <- system(paste("sha224sum ", paste(ufDir, "/", ufCount, filetype, sep = ""), " | awk '{print $1}'", sep = ""), intern=TRUE)
		      if (hash %in% ufSample_table_content) {
		        system(paste("rm -f ", ufDir, "/", ufCount, filetype, sep = ""))
		      }
		      else {
  	        file.rename(paste(ufDir, "/", ufCount, filetype, sep = ""), paste(ufDir, "/", hash, filetype, sep = ""))
  	        filepath <- toString(paste(ufDir, "/", hash, filetype, sep = ""))
  	        tryCatch(
  	          {
  	            aa <- openMSfile(filepath)
  	            aai <- instrumentInfo(aa)
  	            test.empty <- header(aa)
  	            test.empty <- test.empty$lowMZ
  	            on.exit(close(aa))
  	            close(aa)
  	            if (length(test.empty) == 0) {
  	              file.remove(filepath)
  	            }
  	            if (length(test.empty) > 1) {
  	              if (aai$ionisation == "electrospray ionization") {
  	                chromtype <- 1
  	              } else {
  	                chromtype <- 2
  	              }
  	              ogXCMSset <- readMSData(filepath, mode = "onDisk")
  	              saveRDS(ogXCMSset, paste(ufDir, "/", hash, ".rds", sep = ""))
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
  	            file.remove(filepath)
  	            file.remove(paste(ufDir, "/", hash, ".rds", sep = ""))
  	            system("echo ERROR")
  	            print("REMOVED FILE")
  	            return(NA)
  	          }
  	        )
		      }
		    }
		    #for potential async useage
		    future({upload_file(count, fileloc, rfuDir, rfuTime, rfuSample_description, rfuFile$name, rfuFile_tags[count + 1], rfuDb_usr, rfuD_pwd, rfuSample_table_content)}, seed = NULL)
		    #upload_file(count, fileloc, rfuDir, rfuTime, rfuSample_description, rfuFile$name, rfuFile_tags[count + 1], rfuDb_usr, rfuD_pwd, rfuSample_table_content)
		    count <- count + 1
		    # sink()
		  }
		}
		sample_descr <- toString(input$sample_description)
		sample_hashes_existing <- sample_table_content[, 1]
		future({run_file_upload(file, time, dir, sample_descr, db_usr, db_pwd, sample_hashes_existing, file_tags)}, seed = NULL)
		# run_file_upload(file, time, dir, input$sample_description, db_usr, db_pwd, sample_table_content, file_tags)
		updateTextInput(session, "sample_description", value = "")
		shinyjs::alert('Data upload in progress... Check back later!\nWarning: Do not close the current tab!')
		# session$reload()
	})
  
	#Update table to selected datatype
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
	        #Render TIC(s)
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
	
	#Updata parameter usage to datatype and preset selection
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
	
	#Conditionally disable the "Preview spectrum" button
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
	
	#Verandering van parameters naar preset
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
	
	#Show relevant parameters for peak method
	observeEvent(input$peak_method, {
	  
	})
	
	#Actions when job gets submitted
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
	      if (input$debug == 0) {
	        #future({metaboanalyst_data_processing(jobfiles, params, preset, job_id, db_usr, db_pwd)}, seed = NULL)#, packages = .packages(TRUE))
	        future({xcms_data_processing(jobfiles, params, preset, job_id, db_usr, db_pwd, 7)}, seed = NULL)
	      }
	      else {
	        #metaboanalyst_data_processing(jobfiles, params, preset, job_id, db_usr, db_pwd) #debug
	        xcms_data_processing(jobfiles, params, preset, job_id, db_usr, db_pwd, 7)
	      }
	      session$reload()
	      updateTabsetPanel(session, "tabSwitch",
	                        selected = "p3")
	    }
	  }
	})
	
	#Function to tweak parameters of a previously processed job
	observeEvent(input$resub1 | input$resub2 | input$resub3 | input$resub4 | input$resub5, ignoreInit =  T, {
	  job_plan <- NA
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
	  if (is.na(job_plan) == FALSE) {
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
	    send_query(stringr::str_glue(paste("UPDATE parameter SET ", updatestring, " WHERE job_id = ", job_id, ";", sep = "")))
	    query <- stringr::str_glue(paste("SELECT * FROM sample WHERE sample_hash IN (SELECT sample_hash FROM sample_job WHERE job_id = ", job_id, ");", sep = ""))
	    jobfiles <<- get_query(query)
	    names(params) <- used_param_names
	    params <- c(params, job_id=job_id)
	    params <- as.data.frame(t(params))
	    #future({xcms_data_processing(jobfiles, params, 0, job_id, db_usr, db_pwd, job_plan)}, seed = NULL)
	    xcms_data_processing(jobfiles, params, 0, job_id, db_usr, db_pwd, job_plan)
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
	  
	  output$job_analysis <- renderUI(tagList(
	    h4(paste("Output of job", job_id, ":\n", job_table_content[input$jobs_rows_selected, ]$job_name)),
        fluidRow(style='max-width:100%;padding:10px;',
  			         column(12,
  			                uiOutput("jobInformation")
  			                )
  			         ),
  			fluidRow(style='max-width:100%;padding:10px;',
  			         column(6,
  			                div(DT::renderDataTable({
  			                  DT::datatable(used_parameters, selection = 'none',
  			                                rownames= FALSE, options = list(scrollX = TRUE, autoWidth = TRUE, dom = 't'))
  			                }, server = FALSE))#, style = "font-size: 75%; width: 25%")
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
  			                  ft_fv <- log2(featureValues(xset, filled = TRUE))
  			                  ft_fv <- scale(featureValues(xset), scale = TRUE)
  			                  pc <- prcomp(t(na.omit(ft_fv)), center = TRUE)
  			                  pca <- plot_ly(x = NULL, y = NULL, z = NULL, type = "scatter3d") %>% layout(title  = "PCA", scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2"), zaxis = list(title = "PC3")))
  			                  for (i in 1:nrow(selection)){
  			                    pca <- pca %>% add_trace(y = pc$x[,2][i], x = pc$x[,1][i], z = pc$x[,3][i], name = paste(selection$original_file_name[i], "|", selection$metadata[i]), mode = 'markers')
  			                  }
  			                  pca
  			                })
  			         ),
  			         column(6,
  			                numericInput("sample_nr_peaks", label = "Select a sample to view its detected peaks.", min = 1, max = nrow(selection), value = 1),
  			                renderPlotly({
  			                  chrom <- plot_ly(x = NULL, y = NULL, type = 'scatter', mode = 'markers') %>% layout(title = "Peaks over TIC", xaxis = list(title = "Retention time (seconds)"), yaxis = list(title = "TotIonCurrent"))
  			                  #featspec <- featureSpectra(xset, msLevel = 1, return.type="Spectra")
  			                  load(file = rda_path$file_path_peaks) #loads variable annot_spectra
  			                  tt <- spectraData(annot_spectra)
  			                  #aa <- readRDS(selection$original_XCMSnExp_path[i], refhook = NULL)
  			                  for (i in round(input$sample_nr_peaks)){
  			                    aa <- filterFile(xset, i)
  			                    #Get ids of peaks in sample
  			                    tp <- rownames(chromPeaks(aa))
  			                    
  			                    #Get scanids for peaks in sample
  			                    # p <- unique(tt[tt[, "peak_id"] %in% tp, c("scanIndex", "rtime", "totIonCurrent")])
  			                    # y <- tic(aa)
  			                    # x <- rtime(aa)
  			                    # xp <- p[,"rtime"]
  			                    # yp <- p[,"totIonCurrent"]
  			                    # peak <- x %in% xp
  			                    for (peak in 1:length(tp)) {
  			                      p <- unique(tt[tt[, "peak_id"] == tp[peak], c("peak_id", "scanIndex", "rtime", "totIonCurrent")])
  			                      xp <- p[,"rtime"]
  			                      yp <- p[,"totIonCurrent"]
  			                      chrom <- chrom %>% add_lines(x = xp, y = yp, name = paste(selection$original_file_name[i], "| peak", peak), fill = 'tozeroy')
  			                    }
  			                    y <- tic(aa)
  			                    x <- rtime(aa)
  			                    chrom <- chrom %>% add_lines(x = x, y = y, name = paste(selection$original_file_name[i], "|", selection$metadata[i]))
  			                  }
  			                  # y <- tic(featspec)
  			                  # x <- rtime(featspec)
  			                  # x <- featureDefinitions(xset)[,"rtmed"]
  			                  # y <- featureDefinitions(xset)[,"mzmed"]
  			                  
  			                  
  			                  # for (i in 1:length(xset)) {
  			                  #   print(i)
  			                  # #   print(length(xset))
  			                  #   y <- tic(xset)
  			                  #   x <- rtime(xset)
  			                  #   chrom <- chrom %>% add_trace(y = y, x = x, name = 'peak', mode = 'markers')
  			                  # # }
  			                  chrom
  			                })
  			         )
  			),
  			fluidRow(style='max-width:100%;padding:10px;',
  			         column(6,
  			                uiOutput("analysisCompounds")
  			         ),
  			         column(6,
  			                uiOutput("analysisDifferential")
  			         )
  			),
  			fluidRow(style='max-width:100%;padding:10px;',
  			         column(6,
  			                uiOutput("analysisVenn")
  			         ),
  			         column(6,
  			                uiOutput("analysisHeatmap")
  			         )
  			),
  	    fluidRow(style='max-width:100%;padding:10px;',
  	             column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  	             ),
  	             column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  	                    h5("1. Detection"),
  	                    selectInput("Peak_methodA", label = "Peak detection method", choices = list("centWave" = 0, "Massifquant" = 1, "MatchedFilter" = 2), selected = used_parameters$Peak_method),
                      renderUI(peak_options),
                      actionButton("resub1", label = "\n Rerun")
  	             ),
  	             column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  	                    h5("2. Refinement"),
  	                    selectInput("Ref_methodA", label = "Peak refinement method", choices = list("MergeNeighboringPeaks" = 0, "FilterIntensity" = 1), selected = used_parameters$Ref_method),
  	                    renderUI(ref_options),
  	                    actionButton("resub2", label = "\n Rerun")
  	             ),
  	             column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  	                    h5("3. Alignment"),
  	                    selectInput("Align_methodA", label = "Peak alignment method", choices = list("Obiwarp" = 0, "PeakGroups" = 1), selected = used_parameters$Align_method),
  	                    renderUI(aln_options),
  	                    actionButton("resub3", label = "\n Rerun")
  	             ),
  	             column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  	                    h5("4. Grouping"),
  	                    selectInput("Group_methodA", label = "Peak grouping method", choices = list("PeakDensity" = 0, "MzClust" = 1, "NearestPeaks" = 2), selected = used_parameters$Group_method),
  	                    renderUI(grp_options),
  	                    actionButton("resub4", label = "\n Rerun")
  	             ),
  	             column(2, style='margin-bottom:30px;border-left:1px solid #dfd7ca;; padding: 10px;',
  	                    h5("5. Filling"),
  	                    numericInput("expandMzA", label = "expandMz", 0, value = used_parameters$expandMz),
  	                    numericInput("expandRtA", label = "expandRt", 0, value = used_parameters$expandRt),
  	                    numericInput("fixedMzA", label = "fixedMz", 0, value = used_parameters$fixedMz),
  	                    numericInput("fixedRtA", label = "fixedRt", 0, value = used_parameters$fixedRt),
  	                    actionButton("resub5", label = "\n Rerun")
  	             )
  	    )
	  ))
	})
	
  mass_data_annotation <- function(){
    #SPLASH FUTURE
  }
  
  mass_data_visualisation <- function(){
    
  }
  
	#################
	###########       ASYNC FUNCTIONS
	# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
	#Function to process the selected MS-files
	metaboanalyst_data_processing <- function(massfiles, parameters, preset, job_id, db_usr, db_pwd){ #https://cran.r-project.org/web/packages/future.batchtools/future.batchtools.pdf
	  #####################
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
    	  files <- massfiles$file_path
    	  # preset <- 3 # !!! CUSTOM PARAMETERS ARE NONFUNCTIONAL, AUTOMATIC HARDCODED. REMOVE WHEN FUNCTIONAL !!!
    	  if (preset == 3) {
    	    param_initial <- SetPeakParam(platform = "general")
    	    send_query(stringr::str_glue(paste("UPDATE job SET job_status = '1/5 Performing ROI extraction...' WHERE job_id = ", job_id, ";", sep = "")))
    	    raw_train <- PerformROIExtraction(files, rt.idx = 0.2, rmConts = FALSE, plot = FALSE)
    	    send_query(stringr::str_glue(paste("UPDATE job SET job_status = '2/5 Performing parameter optimization...' WHERE job_id = ", job_id, ";", sep = "")))
    	    def_params <- PerformParamsOptimization(raw_train, param = param_initial, ncore = 5)
    	    send_query(stringr::str_glue(paste("UPDATE parameter SET min_peakwidth = ", toString(def_params$min_peakwidth), ", max_peakwidth = ", toString(def_params$max_peakwidth), ", mzdiff = ", toString(def_params$mzdiff), ", snthresh = ", toString(def_params$snthresh), ", bw = ", toString(def_params$bw), ", Peak_method = '", toString(def_params$Peak_method), "', ppm = ", toString(def_params$ppm), ", noise = ", toString(def_params$noise), ", prefilter = ", toString(def_params$prefilter), ", value_of_prefilter = ", toString(def_params$value_of_prefilter), ", minFraction = ", toString(def_params$minFraction), ", minSamples = ", toString(def_params$minSamples), ", maxFeatures = ", toString(def_params$maxFeatures), ", fitgauss = ", toString(def_params$fitgauss), ", mzCenterFun = '", toString(def_params$mzCenterFun), "', integrate = ", toString(def_params$integrate), ", extra = ", toString(def_params$extra), ", span = ", toString(def_params$span), ", smooth = '", toString(def_params$smooth), "', family = '", toString(def_params$family), "', polarity = '", toString(def_params$polarity), "', perc_fwhm = ", toString(def_params$perc_fwhm), ", max_charge = ", toString(def_params$max_charge), ", max_iso = ", toString(def_params$max_iso), ", corr_eic_th = ", toString(def_params$corr_eic_th), ", mz_abs_add = ", toString(def_params$mz_abs_add), ", rmConts = ", toString(def_params$rmConts), ", RT_method = '", toString(def_params$RT_method), "' WHERE job_id = ", job_id, ";", sep = "")))
    	  }
    	  else {
    	    #fwhm = parameters$fwhm, steps= parameters$steps, peakBinSize = parameters$peakbinsize, criticalValue = parameters$criticalvalue, consecMissedLimit = parameters$consecmissedlimit, unions = parameters$unions, checkBack = parameters$checkback, withWave = parameters$withwave, profStep = parameters$profstep,
    	    #print(t(parameters))
    	    parameters[is.na(parameters)] <- 0
    	    def_params <- SetPeakParam(Peak_method = parameters$Peak_method, RT_method = parameters$RT_method, mzdiff = as.double(parameters$mzdiff), snthresh = as.double(parameters$snthresh), bw = as.double(parameters$bw), ppm = as.double(parameters$ppm), min_peakwidth = as.double(parameters$min_peakwidth), max_peakwidth = as.double(parameters$max_peakwidth), noise = as.double(parameters$noise), prefilter = as.double(parameters$prefilter), value_of_prefilter = as.double(parameters$value_of_prefilter), minFraction = as.double(parameters$minFraction), minSamples = as.double(parameters$minSamples), maxFeatures = as.double(parameters$maxFeatures), mzCenterFun = parameters$mzCenterFun, integrate = as.double(parameters$integrate), extra = as.double(parameters$extra), span = as.double(parameters$span), smooth = parameters$smooth, family = parameters$family, fitgauss = as.logical(parameters$fitgauss), polarity = parameters$polarity, perc_fwhm = as.double(parameters$perc_fwhm), mz_abs_iso = as.double(parameters$mz_abs_iso), max_charge = as.double(parameters$max_charge), max_iso = as.double(parameters$max_iso), corr_eic_th = as.double(parameters$corr_eic_th), mz_abs_add = as.double(parameters$mz_abs_add), rmConts = parameters$rmConts) #verboseColumns
    	  }
    	  send_query(stringr::str_glue(paste("UPDATE job SET job_status = '3/5 Importing raw spectra...' WHERE job_id = ", job_id, ";", sep = "")))
    	  rawData <- ImportRawMSData(path = files, plotSettings = SetPlotParam(Plot = FALSE)) #ontbreekt ppm, min_peakwidth, max_peakwidth, mzdiff, snthresh, noise, prefilter, value_of_prefilter
    	  send_query(stringr::str_glue(paste("UPDATE job SET job_status = '4/5 Performing peak profiling...' WHERE job_id = ", job_id, ";", sep = "")))
    	  mSet <- PerformPeakProfiling(rawData, def_params, plotSettings = SetPlotParam(Plot = FALSE))
    	  send_query(stringr::str_glue(paste("UPDATE job SET job_status = '5/5 Performing peak annotation...' WHERE job_id = ", job_id, ";", sep = "")))
    	  annParams <- SetAnnotationParam(polarity = def_params$polarity, mz_abs_add = def_params$mz_abs_add)
    	  annotPeaks <- PerformPeakAnnotation(mSet, annParams)
    	  maPeaks <- FormatPeakList(annotPeaks, annParams, filtIso =F, filtAdducts = FALSE, missPercent = 1)
    	  end_time <- format(Sys.time() + 60*60, "%Y-%m-%d %X")
    	  send_query(stringr::str_glue(paste("UPDATE job SET job_status = 'Finished' WHERE job_id = ", job_id, ";", sep = "")))
    	  send_query(stringr::str_glue(paste("UPDATE job SET end_time = '", end_time, "' WHERE job_id = ", job_id, ";", sep = "")))
    	  dir <- getwd()
    	  dir <- paste(dir, "/msets/", toString(job_id), sep = "")
    	  dir.create(dir)
    	  Export.PeakTable(mSet = maPeaks, path = dir)
    	  save(mSet, file = paste(dir, "/mSet.rda", sep = ""))
    	  todf <- data.frame(
    	    job_id = toString(job_id),
    	    file_path_rda = toString(paste(dir, "/mSet.rda", sep = "")),
    	    file_path_peaks = toString(paste(dir, "/", "metaboanalyst_input.csv", sep = ""))
    	  )
    	  insert_query("processed_sample", todf)
    	  
	    },
	    error = function(cnd){
	        print(cnd)
	        send_query(stringr::str_glue(paste("UPDATE job SET job_status = 'CRASHED', end_time = '", format(Sys.time() + 60*60, "%Y-%m-%d %X"), "' WHERE job_id = ", job_id, ";", sep = "")))
	      return(NA)
	    }
	  )
	}

  xcms_data_processing <- function(massfiles, parameters, preset, job_id, db_usr, db_pw, job_plan){ #https://cran.r-project.org/web/packages/future.batchtools/future.batchtools.pdf
    register(bpstart(MulticoreParam(8)))
    #parameters[] <- lapply(parameters, function(x) as.double(as.numeric(x)))
    # query <- stringr::str_glue("SELECT * FROM parameter WHERE job_id = '", job_id, "';")
    # def_params <- get_query(query)
    #####################
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
          print(def_params)
          def_params[is.na(def_params)] <- 0
          # def_params <- SetPeakParam(Peak_method = parameters$Peak_method, RT_method = parameters$RT_method, mzdiff = as.double(parameters$mzdiff), snthresh = as.double(parameters$snthresh), bw = as.double(parameters$bw), ppm = as.double(parameters$ppm), min_peakwidth = as.double(parameters$min_peakwidth), max_peakwidth = as.double(parameters$max_peakwidth), noise = as.double(parameters$noise), prefilter = as.double(parameters$prefilter), value_of_prefilter = as.double(parameters$value_of_prefilter), minFraction = as.double(parameters$minFraction), minSamples = as.double(parameters$minSamples), maxFeatures = as.double(parameters$maxFeatures), mzCenterFun = parameters$mzCenterFun, integrate = as.double(parameters$integrate), extra = as.double(parameters$extra), span = as.double(parameters$span), smooth = parameters$smooth, family = parameters$family, fitgauss = as.logical(parameters$fitgauss), polarity = parameters$polarity, perc_fwhm = as.double(parameters$perc_fwhm), mz_abs_iso = as.double(parameters$mz_abs_iso), max_charge = as.double(parameters$max_charge), max_iso = as.double(parameters$max_iso), corr_eic_th = as.double(parameters$corr_eic_th), mz_abs_add = as.double(parameters$mz_abs_add), rmConts = parameters$rmConts) #verboseColumns
          #def_params <- parameters
        }
        #processing
        dir <- getwd()
        dir <- paste(dir, "/processed_data", sep = "")
        if (job_plan == 7){
          jstep <- 3
          #range meegeven
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '3/9 Reading MS data...' WHERE job_id = ", job_id, ";", sep = "")))
          #Load raw files into XCMSnExp-object
          pd <- data.frame(sample_name = massfiles$original_file_name,
                           sample_group = massfiles$metadata,
                           stringsAsFactors = FALSE)
          xset <- readMSData(files = files, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")
          xset <- filterRt(xset, c(as.double(def_params$rtrmin), as.double(def_params$rtrmax)))
        } else {
          rdaquer <- stringr::str_glue(paste("SELECT * FROM processed_sample WHERE job_id = ", job_id, ";", sep = ""))
          rda_path <- get_query(rdaquer)
          load(file=toString(rda_path$file_path_rda))
        }
        if (job_plan >= 2 | job_plan == 7){
          jstep <- 4
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '4/9 Finding peaks...' WHERE job_id = ", job_id, ";", sep = "")))
          
          #Find peaks
          if (as.integer(def_params$Peak_method) == 0){
            peak_params <- CentWaveParam(
              ppm = as.double(def_params$ppm),
              peakwidth = c(as.double(def_params$min_peakwidth), as.double(def_params$max_peakwidth)),
              snthresh = as.double(def_params$snthresh),
              prefilter = c(as.double(def_params$prefilter), as.double(def_params$value_of_prefilter)),
              mzCenterFun = c("wMean", "mean", "apex", "wmMeanApex3", "meanApex3")[as.integer(def_params$mzCenterFun) + 1],
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
          else {
            return()
          }
          print(peak_params)
          xset <-findChromPeaks(xset, param = peak_params)
        }
        if (job_plan >= 3 | job_plan == 7){
          jstep <- 5
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '5/9 Refining peaks...' WHERE job_id = ", job_id, ";", sep = "")))
          #Peak refinement
          if (preset == 3) {
            refmeth <- MergeNeighboringPeaksParam()
          }
            else {
            if (as.integer(def_params$Ref_method) == 0){
                refmeth <- MergeNeighboringPeaksParam(
                                                  expandRt = as.double(def_params$expandRt),
                                                  expandMz = as.double(def_params$expandMz),
                                                  ppm = as.double(def_params$ppm),
                                                  minProp = as.double(def_params$minProp)
                                                  )
              }
            else if (as.integer(def_params$Ref_method) == 1) {
            # cpp <- CleanPeaksParam()
              refmeth <-  FilterIntensityParam(
                  threshold = as.double(def_params$threshold),
                  nValues = 1L,
                  value = "maxo"
                  )
            }
          }
          xset <- refineChromPeaks(xset, refmeth)
        }
        if (job_plan >= 4 | job_plan == 7){
          jstep <- 6
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '6/9 Aligning peaks...' WHERE job_id = ", job_id, ";", sep = "")))
          #Peak alignment
          if (preset == 3) {
            align_param <- ObiwarpParam()
          }
          else {
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
          xset <- adjustRtime(xset, param = align_param)
        }
        if (job_plan >= 5 | job_plan == 7){
          jstep <- 7
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '7/9 Grouping peaks...' WHERE job_id = ", job_id, ";", sep = "")))
          if (preset == 3) {
            peak_group_param <- PeakDensityParam()
          }
          else {
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
            else if (as.integer(def_params$Group_method) == 1){
              peak_group_param <- MzClustParam(
                sampleGroups = massfiles$metadata,
                ppm = as.double(def_params$ppm),
                absMz = as.double(def_params$absMz),
                minFraction = as.double(def_params$minFraction),
                minSamples = as.double(def_params$minSamples)
              )
            }
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
          xset <- groupChromPeaks(xset, param = peak_group_param)
        }
        if (job_plan >= 6 | job_plan == 7){
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
          else {
            fillparam <- FillChromPeaksParam(
              expandMz = as.double(def_params$expandMz),
              expandRt = as.double(def_params$expandRt),
              ppm = as.double(def_params$ppm),
              fixedMz = as.double(def_params$fixedMz),
              fixedRt = as.double(def_params$fixedRt)
            )
          }
          xset <- fillChromPeaks(xset, param = fillparam)
          save(xset, file = paste(dir, "/", job_id, ".rda", sep = ""))
        }
        if (job_plan == 7){
          jstep <- 9
          send_query(stringr::str_glue(paste("UPDATE job SET job_status = '9/9 Annotating peaks...' WHERE job_id = ", job_id, ";", sep = "")))
          # intmzs <- intensity(xset)
          # intmzs <- split(intmzs, f = fromFile(xset))
          # mzs <- mz(xset)
          # mzs <-split(mzs, f = fromFile(xset))
          # #intmzs[[1]]$F1.S0106
          # peaknr <- 1
          # chrom_peakids <- featureDefinitions(xset)[[peaknr,"peakidx"]]
          # chrompeaks_values <- chromPeaks(xset)[chrom_peakids,]
          
          
          annot_spectra <- featureSpectra(xset, msLevel = 1, return.type = "Spectra", skipFilled = FALSE)
          save(annot_spectra, file = paste(dir, "/", job_id, "_annot.rda", sep = ""))
        }
        todf <- data.frame(
          job_id = toString(job_id),
          file_path_rda = toString(paste(dir, "/", job_id, ".rda", sep = "")),
          file_path_peaks = toString(paste(dir, "/", job_id, "_annot.rda", sep = ""))
        )
        # length(chromPeaks(xset)[,chromPeaks(xset)[,"sample"] == 1])
        insert_query("processed_sample", todf)
        end_time <- format(Sys.time() + 60*60, "%Y-%m-%d %X")
        send_query(stringr::str_glue(paste("UPDATE job SET job_status = 'Finished' WHERE job_id = ", job_id, ";", sep = "")))
        send_query(stringr::str_glue(paste("UPDATE job SET end_time = '", end_time, "' WHERE job_id = ", job_id, ";", sep = "")))
      },
      error = function(cnd){
        print(cnd)
        send_query(stringr::str_glue(paste("UPDATE job SET job_status = 'CRASHED AT ", jstep, "/9', end_time = '", format(Sys.time() + 60*60, "%Y-%m-%d %X"), "' WHERE job_id = ", job_id, ";", sep = "")))
        return(NA)
      }
    )
  }
}

runApp(shinyApp(ui = ui, server = server), host = hostip, port = portnr)
