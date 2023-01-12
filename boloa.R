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

all_params <- c("Peak_method", "RT_method", "mzdiff", "snthresh", "bw", "ppm", "min_peakwidth", "max_peakwidth", "noise", "prefilter", "value_of_prefilter", "minFraction", "minSamples", "maxFeatures", "mzCenterFun", "integrate", "extra", "span", "smooth", "family", "fitgauss", "polarity", "perc_fwhm", "mz_abs_iso", "max_charge", "max_iso", "corr_eic_th", "mz_abs_add", "rmConts", "verboseColumns")
lcms_only <- tail(all_params, n=9)


ui <- fluidPage(
  autoWaiter(color = "white", html = spin_3()),
  
	useShinyjs(),
	theme = bs_theme(version = 4, bootswatch = "spacelab"),
	mainPanel(
		titlePanel("Boloa"),
		tabsetPanel(id = "tabSwitch",
			tabPanel("Select data",  icon = icon("arrow-right"),
				fluidRow(
					column(3,
						br(),
						textInput("title", label = "Sample description(s)"),
						h4("Mass Spectrometry data:"),
						fileInput("msdata", label = "MS-data", multiple = TRUE, accept = c(".mzXML", ".raw", ".CDF")),
						tableOutput("bestanden"),
						actionButton("dataupl", label = "Submit", width = 180),
						br(),
						verbatimTextOutput("upl_completed"),
						style='margin-bottom:30px;border-right:1px solid #dfd7ca;; padding: 10px;'
					),
					column(9,
						DT::dataTableOutput("uploaded_samples")
					)
				)
			),
			tabPanel("Process data", icon = icon("arrow-right"),
				fluidRow(
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
				fluidRow(
					column(7,
						actionButton("preview", label = "\n Preview spectra")
						#actionButton("prevtic", label = "\n Preview TIC")
					),
					column(5,
						#uiOutput("mzrange"),
						uiOutput("previewplot")
					)
				),
				fluidRow(
					style='margin-top:20px; border-top:1px solid #dfd7ca;',
					column(2,
						selectInput("preset", label = "Parameter presets", choices = list("None" = 0, "LC-MS" = 1, "GC-MS" = 2, "Automatic (centWave only)" = 3), width = 160),
						textInput("job_name", label = "Job name"),
						actionButton("submitJob", "Submit job", width = 180, icon=icon("play")),
						checkboxInput("debug", label = "Debug mode", 0),
						tableOutput("jobval")
					),
					column(2, #zet deze in een nieuwe tab, de 1e wordt voor alleen geuploade samples.
						h4("Parameters:"),
						h5("1. Peak detection"),
						numericInput("min_peakwidth", label = "Minimum peak width", 5),
						numericInput("max_peakwidth", label = "Maximum peak width", 30),
						numericInput("mzdiff", label = "mzdiff", 0.01),
						numericInput("snthresh", label = "snthresh", 10),
						numericInput("bw", label = "bw", 10),
						selectInput("Peak_method", label = "Peak_method", choices = list("centWave" = 0, "Massifquant" = 1, "MatchedFilter" = 2)),
						selectInput("RT_method", label = "RT_method", choices = list("loess" = 0, "obiwarp" = 1)),
						numericInput("ppm", label = "ppm", 93.07),
						numericInput("noise", label = "noise", 0),
						numericInput("prefilter", label = "prefilter", 2),
						numericInput("value_of_prefilter", label = "value_of_prefilter", 328.63),
					),
					column(2,
						br(),br(),
						h5("2. Alignment"),
						numericInput("minFraction", label = "minFraction", 0.8),
						numericInput("minSamples", label = "minSamples", 1),
						numericInput("maxFeatures", label = "maxFeatures", 100),
						checkboxInput("fitgauss", label = "fitgauss", 0),
						selectInput("mzCenterFun", label = "mzCenterFun", choices = list("wMean" = 0, "mean" = 1, "apex" = 2, "wMeanApex3" = 3, "meanApex3" = 4)),
						numericInput("integrate", label = "integrate", 1),
						numericInput("extra", label = "extra", 1),
						numericInput("span", label = "span", 0.4),
						selectInput("smooth", label = "smooth", choices = list("Loess" = 0, "Linear" = 1)),
						selectInput("family", label = "family", choices = list("gaussian" = 0, "symmetric" = 1))
					),
					column(2,
						br(),br(),
						h5("3. Annotation"),
						checkboxInput("verboseColumns", label = "verboseColumns", 0),
						selectInput("polarity", label = "polarity", choices = list("negative" = 0, "positive" = 1)),
						numericInput("perc_fwhm", label = "perc_fwhm", 0.6),
						numericInput("mz_abs_iso", label = "mz_abs_iso", 0.005),
						numericInput("max_charge", label = "max_charge", 2),
						numericInput("max_iso", label = "max_iso", 2),
						numericInput("corr_eic_th", label = "corr_eic_th", 0.85),
						numericInput("mz_abs_add", label = "mz_abs_add", 0.001),
						checkboxInput("rmConts", label = "rmConts", 1),
						br()
					),
					column(4)
				)
			),
			tabPanel("Jobs", icon = icon("arrow-right"),
			  fluidRow(
			    h4("Select a job to analyse contents:"),
			    DT::dataTableOutput("jobs"),
			    actionButton("updateJobsTable", label = "\n Update table", icon = icon("arrows-rotate")),
			    actionButton("analyse", label = "\n Analyse job", icon=icon("play"))
			  )
			),
			tabPanel("Analysis",  icon = icon("arrow-right"))
		)
	)
)



server <- function(input, output, session) {
  
	#Weergave van geuploade bestand(en)
	observeEvent(input$msdata, {
		output$bestanden <- renderTable({
			file <- input$msdata
			ext <- tools::file_ext(file$datapath)
			req(file)
			validate(need(ext %in% c("raw", "mzXML", "CDF"), "Please upload a .raw, .CDF or .mzXML file!"))
		})
	})
  
	#Update van tabel op "Upload data" tab
	updateEvent <- reactive({
		list(input$tabSwitch, input$dataupl)
	})
	
	#Sample selection table
	observeEvent(updateEvent, {
		query <- stringr::str_glue("SELECT upload_date, original_file_name, sample_hash, sample_name, chromatography_type FROM sample ORDER BY upload_date DESC;")
		sample_table_content <<- get_query(query)
		sample_table_content$chromatography_type[sample_table_content$chromatography_type == 1] <- "Liquid"
		sample_table_content$chromatography_type[sample_table_content$chromatography_type == 2] <- "Gas"
		output$uploaded_samples <- DT::renderDataTable({
			DT::datatable(sample_table_content[, c(1, 2, 4, 5)])
		}, server = FALSE)
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
		if (is.null(input$msdata)) {
			output$upl_completed <- renderText({
				'Please upload a file'
			})
		}
		if (nchar(input$title) != 0 && is.null(input$msdata) != TRUE) {
			dir <- getwd()
			dir <- paste(dir, "/massascans", sep = "")
			#aanmaak metadata object
			file <- input$msdata
			time <- format(Sys.time() + 60*60, "%Y-%m-%d %X")
			count <- 0
			progress <- shiny::Progress$new()
			on.exit(progress$close())
			progress$set(message = "Uploading", value = 0)
			for (fileloc in file$datapath) {
			  progress$inc(1/length(file$datapath), detail = paste("File:", file$name[count + 1]))
				file.copy(fileloc, dir)
				upload_file <- function(countUF, filelocUF, dirUF, timeUF, titleUF, filenamesUF){
  				if (grepl('.raw', filelocUF, fixed=TRUE)) {
  					system(paste("docker run -it --rm -v ", dirUF, ":/massascans chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /massascans/", countUF, ".raw", " --mzXML -o /massascans", sep = "")) #change .raw files to .mzXML
  					system(paste("rm ", dirUF, "/", countUF, ".raw", sep = ""))
  				}
				  if (grepl('.mzXML', filelocUF, fixed=TRUE)) {
				    filetype <- ".mzXML"
				  }
				  if (grepl('.CDF', filelocUF, fixed=TRUE)) {
				    filetype <- ".CDF"
				  }
  				hash <- system(paste("sha224sum ", paste(dirUF, "/", countUF, filetype, sep = ""), " | awk '{print $1}'", sep = ""), intern=TRUE)
  				file.rename(paste(dirUF, "/", countUF, filetype, sep = ""), paste(dirUF, "/", hash, filetype, sep = ""))
  				filepath <- toString(paste(dirUF, "/", hash, filetype, sep = ""))
  				tryCatch(
  					{
  						aa <- openMSfile(filepath)
  						aai <- instrumentInfo(aa)
  						test.empty <- header(aa)
  						test.empty <- test.empty$lowMZ
  						on.exit(close(aa))
  						close(aa)
  						if (length(test.empty) == 0) {
  						  shinyjs::alert(paste(toString(file$name[countUF + 1]), " is empty! Please refrain from uploading empty MS-files!"))
  						  file.remove(filepath)
  						}
  						if (length(test.empty) > 1) {
  						  if (aai$ionisation == "electrospray ionization") {
  						    chromtype <- 1
  						  } else {
  						    chromtype <- 2
  						  }
  						  todf <- data.frame(
  						    sample_hash = toString(hash),
  						    file_path = toString(filepath),
  						    upload_date = toString(timeUF),
  						    sample_name = toString(titleUF),
  						    chromatography_type = toString(chromtype),
  						    original_file_name = toString(filenamesUF[countUF + 1])
  						  )
  						  insert_query("sample", todf)
  						}
  					},
  					error = function(cnd){
  			  		file.remove(filepath)
  					  system("echo ERROR")
  					  print("REMOVED FILE")
  			  		return(NA)
  					}
  				)
				}
				upload_file(count, fileloc, dir, time, input$title, file$name) #for potential async useage
				count <- count + 1
			}
			shinyjs::alert('Data upload completed.\nSelect samples from the table in order to continue analysis.')
			session$reload()
  	}
	})
  
	#Update table to selected datatype
	observeEvent(input$datatype, {
	  output$selected_samples <- DT::renderDataTable({
	    selsamples <<- sample_table_content[input$uploaded_samples_rows_selected,]
	    selsamples <<- selsamples[selsamples$chromatography_type == input$datatype,]
	    selsamples$chromatography_type[selsamples$chromatography_type == 1] <- "Liquid"
	    selsamples$chromatography_type[selsamples$chromatography_type == 2] <- "Gas"
	    DT::datatable(selsamples[, c(1, 2, 4, 5)], selection = 'single')
	  }, server = FALSE)
	  output$previewplot <- renderUI({
	    return(NULL)
	  })
	  output$ticplot <- renderPlotly({
	    return(NULL)
	  })
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
	  output$ticplot <- renderPlotly({
	    return(NULL)
	  })
		if (is.null(input$selected_samples_rows_selected)) {
			shinyjs::disable("preview")
		}
		else {
		  preview_file <- selsamples[input$selected_samples_rows_selected, ]$sample_hash
		  query <- stringr::str_glue("SELECT file_path FROM sample WHERE sample_hash = '", preview_file, "';")
		  aa <- openMSfile(get_query(query)$file_path)
		  on.exit(close(aa))
		  aahead <- header(aa)
		  y <- aahead$totIonCurrent
		  x <- aahead$retentionTime
		  close(aa)
		  output$chromrange <- renderUI({
		    sliderInput("chromrange", "Retention time range:", min = min(x), max = max(x), value = c(min(x), max(x)))
		  })
		  #Render TIC
		  output$ticplot <- renderPlotly({
		    plot_ly(data = data.frame(x,y), x = ~x, y = ~y, type = 'scatter', mode = 'lines', fill = 'tozeroy') %>% layout(title = "Total Ion Current", xaxis = list(title = "Retention time (seconds)"), yaxis = list(title = "TotIonCurrent"))
		  })
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
	    if (length(selsamples$sample_hash) >= 3) {
	      params <- c()
	      for (p in all_params) {
	        params <- c(params, input[[p]])
	      }
	      names(params) <- all_params
	      
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
	      params$Peak_method <- c("centWave", "massifquant", "matchedFilter")[strtoi(input$Peak_method) + 1]
	      params$RT_method <- c("loess", "obiwarp")[strtoi(input$RT_method) + 1]
	      params$mzCenterFun <- c("wMean", "mean", "apex", "wMeanApex3", "meanApex3")[strtoi(input$mzCenterFun) + 1]
	      params$smooth <- c("loess", "linear")[strtoi(input$smooth) + 1]
	      params$family <- c("gaussian", "symmetric")[strtoi(input$family) + 1]
	      params$polarity <- c("negative", "positive")[strtoi(input$polarity) + 1]
	      
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
	        future({metaboanalyst_data_processing(jobfiles, params, preset, job_id, db_usr, db_pwd)}, seed = NULL)#, packages = .packages(TRUE))
	      }
	      else {
	        metaboanalyst_data_processing(jobfiles, params, preset, job_id, db_usr, db_pwd) #debug
	      }
	      session$reload()
	    }
	  }
	})
	
	#Actions involving job analysis
	observeEvent(input$analyse, {
	  
	})
	
  mass_data_annotation <- function(){
    
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
    	    send_query(stringr::str_glue(paste("UPDATE job SET job_status = '1/4 Performing ROI extraction...' WHERE job_id = ", job_id, ";", sep = "")))
    	    raw_train <- PerformROIExtraction(files, rt.idx = 0.2, rmConts = FALSE, plot = FALSE)
    	    send_query(stringr::str_glue(paste("UPDATE job SET job_status = '2/4 Performing parameter optimization...' WHERE job_id = ", job_id, ";", sep = "")))
    	    def_params <- PerformParamsOptimization(raw_train, param = param_initial, ncore = 5)
    	    send_query(stringr::str_glue(paste("UPDATE parameter SET min_peakwidth = ", toString(def_params$min_peakwidth), ", max_peakwidth = ", toString(def_params$max_peakwidth), ", mzdiff = ", toString(def_params$mzdiff), ", snthresh = ", toString(def_params$snthresh), ", bw = ", toString(def_params$bw), ", Peak_method = '", toString(def_params$Peak_method), "', ppm = ", toString(def_params$ppm), ", noise = ", toString(def_params$noise), ", prefilter = ", toString(def_params$prefilter), ", value_of_prefilter = ", toString(def_params$value_of_prefilter), ", minFraction = ", toString(def_params$minFraction), ", minSamples = ", toString(def_params$minSamples), ", maxFeatures = ", toString(def_params$maxFeatures), ", fitgauss = ", toString(def_params$fitgauss), ", mzCenterFun = '", toString(def_params$mzCenterFun), "', integrate = ", toString(def_params$integrate), ", extra = ", toString(def_params$extra), ", span = ", toString(def_params$span), ", smooth = '", toString(def_params$smooth), "', family = '", toString(def_params$family), "', polarity = '", toString(def_params$polarity), "', perc_fwhm = ", toString(def_params$perc_fwhm), ", max_charge = ", toString(def_params$max_charge), ", max_iso = ", toString(def_params$max_iso), ", corr_eic_th = ", toString(def_params$corr_eic_th), ", mz_abs_add = ", toString(def_params$mz_abs_add), ", rmConts = ", toString(def_params$rmConts), ", RT_method = '", toString(def_params$RT_method), "' WHERE job_id = ", job_id, ";", sep = "")))
    	  }
    	  else {
    	    #fwhm = parameters$fwhm, steps= parameters$steps, peakBinSize = parameters$peakbinsize, criticalValue = parameters$criticalvalue, consecMissedLimit = parameters$consecmissedlimit, unions = parameters$unions, checkBack = parameters$checkback, withWave = parameters$withwave, profStep = parameters$profstep,
    	    #print(t(parameters))
    	    parameters[is.na(parameters)] <- 0
    	    def_params <- SetPeakParam(Peak_method = parameters$Peak_method, RT_method = parameters$RT_method, mzdiff = as.double(parameters$mzdiff), snthresh = as.double(parameters$snthresh), bw = as.double(parameters$bw), ppm = as.double(parameters$ppm), min_peakwidth = as.double(parameters$min_peakwidth), max_peakwidth = as.double(parameters$max_peakwidth), noise = as.double(parameters$noise), prefilter = as.double(parameters$prefilter), value_of_prefilter = as.double(parameters$value_of_prefilter), minFraction = as.double(parameters$minFraction), minSamples = as.double(parameters$minSamples), maxFeatures = as.double(parameters$maxFeatures), mzCenterFun = parameters$mzCenterFun, integrate = as.double(parameters$integrate), extra = as.double(parameters$extra), span = as.double(parameters$span), smooth = parameters$smooth, family = parameters$family, fitgauss = as.logical(parameters$fitgauss), polarity = parameters$polarity, perc_fwhm = as.double(parameters$perc_fwhm), mz_abs_iso = as.double(parameters$mz_abs_iso), max_charge = as.double(parameters$max_charge), max_iso = as.double(parameters$max_iso), corr_eic_th = as.double(parameters$corr_eic_th), mz_abs_add = as.double(parameters$mz_abs_add), rmConts = parameters$rmConts) #verboseColumns
    	  }
    	  send_query(stringr::str_glue(paste("UPDATE job SET job_status = '3/4 Importing raw spectra...' WHERE job_id = ", job_id, ";", sep = "")))
    	  rawData <- ImportRawMSData(path = files, plotSettings = SetPlotParam(Plot = FALSE)) #ontbreekt ppm, min_peakwidth, max_peakwidth, mzdiff, snthresh, noise, prefilter, value_of_prefilter
    	  send_query(stringr::str_glue(paste("UPDATE job SET job_status = '4/4 Performing peak profiling...' WHERE job_id = ", job_id, ";", sep = "")))
    	  mSet <- PerformPeakProfiling(rawData, def_params, plotSettings = SetPlotParam(Plot = FALSE))
    	  end_time <- format(Sys.time() + 60*60, "%Y-%m-%d %X")
    	  send_query(stringr::str_glue(paste("UPDATE job SET job_status = 'Finished' WHERE job_id = ", job_id, ";", sep = "")))
    	  send_query(stringr::str_glue(paste("UPDATE job SET end_time = '", end_time, "' WHERE job_id = ", job_id, ";", sep = "")))
    	  dir <- getwd()
    	  dir <- paste(dir, "/msets/", toString(job_id), sep = "")
    	  save(mSet, file = paste(dir, ".rda", sep = ""))
    	  todf <- data.frame(
    	    job_id = toString(job_id),
    	    file_path = toString(paste(dir, ".rda", sep = ""))
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

}

runApp(shinyApp(ui = ui, server = server), host = hostip, port = portnr)
