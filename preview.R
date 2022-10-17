#!/bin/Rscript
library(shiny)
library(shinyBS)
library(shinyjs)
library(bslib)
library(digest)
library(RMySQL)

options(shiny.maxRequestSize=100*4096^2)
hostip <- "145.97.18.149"
portnr <- 7000

all_params <- c("min_peakwidth", "max_peakwidth", "mzdiff", "snthresh", "bw", "Peak_method", "RT_method", "ppm", "noise", "prefilter", "value_of_prefilter", "minFraction", "minSamples", "maxFeatures", "fitgauss", "mzCenterFun", "integrate", "pextrapm", "span", "smooth", "gaussian", "verbose.columns", "polarity", "perc_fwhm", "mz_abs_iso", "max_charge", "max_iso", "corr_eic_th", "mz_abs_add", "rmConts")
lcms_only <- tail(all_params, n=9)

ui <- fluidPage(
	useShinyjs(),
	theme = bs_theme(version = 4, bootswatch = "sandstone"),
	mainPanel(
		titlePanel("Boloa"),
		tabsetPanel(id = "tabSwitch",
			tabPanel("Upload data",
				fluidRow(
					column(3,
						br(),
						textInput("title", label = "Job title"),
						h4("Mass Spectrometry data:"),
						fileInput("msdata", label = "MS-data", multiple = TRUE, accept = ".raw"),
						tableOutput("bestanden"),
						radioButtons("datatype", label = "Method",
							choices = list("LC-MS" = 1, "GC-MS" = 2)),
						selectInput("preset", label = "Parameter presets", choices = list("None" = 0, "LC-MS" = 1, "GC-MS" = 2, "Automatic" = 3)),
						actionButton("parsubm", label = "Submit", width = 180),
						br(),
						verbatimTextOutput("upl_completed"),
						style='margin-bottom:30px;border-right:1px solid #dfd7ca;; padding: 10px;'
					),
					column(3, #zet deze in een nieuwe tab, de 1e wordt voor allen geuploade scans.
						h4("Parameters:"),
						h5("1. Peak detection"),
						numericInput("min_peakwidth", label = "Minimum peak width", 6),
						numericInput("max_peakwidth", label = "Maximum peak width", 27),
						numericInput("mzdiff", label = "mzdiff", 27),
						numericInput("snthresh", label = "snthresh", 8),
						numericInput("bw", label = "bw", 2),
						selectInput("Peak_method", label = "Peak_method", choices = list("centWave" = 0, "Massifquant" = 1, "MatchedFilter" = 2)),
						selectInput("RT_method", label = "RT_method", choices = list("loess" = 0, "obiwarp" = 1)),
						numericInput("ppm", label = "ppm", 3.6),
						numericInput("noise", label = "noise", 6000.1),
						numericInput("prefilter", label = "prefilter", 2.0),
						numericInput("value_of_prefilter", label = "value_of_prefilter", 13043.25),
					),
					column(3,
						br(),br(),
						h5("2. Alignment"),
						numericInput("minFraction", label = "minFraction", 0.6),
						numericInput("minSamples", label = "minSamples", 1),
						numericInput("maxFeatures", label = "maxFeatures", 100),
						checkboxInput("fitgauss", label = "fitgauss", 0),
						selectInput("mzCenterFun", label = "mzCenterFun", choices = list("nognietkloppend" = 0, "obiwarp" = 1)),
						numericInput("integrate", label = "integrate", 1),
						numericInput("pextrapm", label = "extra", 1),
						numericInput("span", label = "span", 0.4),
						selectInput("smooth", label = "smooth", choices = list("nognietkloppend" = 0, "obiwarp" = 1)),
						selectInput("gaussian", label = "gaussian", choices = list("nognietkloppend" = 0, "obiwarp" = 1))
					),
					column(3,
						br(),br(),
						h5("3. LC-MS settings"),
						checkboxInput("verbose.columns", label = "verbose.columns", 0),
						selectInput("polarity", label = "polarity", choices = list("negative" = 0, "positive" = 1)),
						numericInput("perc_fwhm", label = "perc_fwhm", 100),
						checkboxInput("mz_abs_iso", label = "mz_abs_iso", 0),
						numericInput("max_charge", label = "max_charge", 1),
						numericInput("max_iso", label = "max_iso", 1),
						numericInput("corr_eic_th", label = "corr_eic_th", 1),
						numericInput("mz_abs_add", label = "mz_abs_add", 0.4),
						checkboxInput("rmConts", label = "rmConts", 1),
						br()
					)
				)
			),
			tabPanel("Mass spectra",
				tableOutput("finished_jobs")
			),
			tabPanel("Analysis")
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
			validate(need(ext == "raw", "Please upload a .raw file"))
		})
	})
	#Update van tabel op "Mass spectra" tab
	observeEvent(input$tabSwitch, {
		output$finished_jobs <- renderTable({
			fj <- data.frame(
				date = c(Sys.time()),
				JobID = c('Test'),
				Spectrum = list.files("massascans")
			)
			data.frame(JobNr = rownames(fj), fj[rownames(fj),])
		})
	})
	#Handelingen als data ingevoerd is ---
	observeEvent(input$parsubm, {
		if (is.null(input$msdata)) {
			output$upl_completed <- renderText({
				'Please upload a file'
			})
		}
		if (nchar(input$title) != 0 && is.null(input$msdata) != TRUE) {
			genhash <- sapply(paste(toString(format(Sys.time(), "%S")), toString(input$datatype), toString(input$title), toString(runif(1))), digest, algo="md5")
			dir <- getwd()
			dir <- paste(dir, "/massascans", sep = "")
			parameters <- c()
			for(app_param in all_params) {
				parameters <- c(parameters, input[[app_param]])
			}
			#aanmaak metadata object
			print(input$msdata[1])
			metadata <- c("title"=toString(input$title), "date"=c(Sys.time()), "params"=parameters, "hash"=toString(genhash))
			file <- input$msdata
			file.copy(file$datapath, dir)
			file.rename(paste(dir, "/0.raw", sep = ""), paste(dir, "/", toString(genhash), ".raw", sep = "")) 
			print(metadata)
			#output$upl_completed <- renderText({
			#	print(metadata)
			#	'Data upload completed.\nCheck the "Mass spectra" tab to track progress.'
			#})
			shinyjs::alert('Data upload completed.\nCheck the "Mass spectra" tab to track progress.')
		}
	})
	observeEvent(input$preset, {
		#Disable all parameter inputs for automatic mode
		for(app_element in all_params) {
			if(input$preset == 3) {
				shinyjs::disable(app_element)
			}
			else {
				if(input$datatype == 2) {
					if(!(app_element %in% lcms_only)) {
						shinyjs::enable(app_element)
					}
				}
				else {
					shinyjs::enable(app_element)
				}
			}
		}
	})
	observeEvent(input$datatype, {
	#Disable/enable buttons for LC-MS
		for(lc_element in lcms_only) {
			if(input$datatype == 2) {
				shinyjs::disable(lc_element)
			}
			else {
				if(input$preset != 3) {
					shinyjs::enable(lc_element)
				}
			}
		}
	})
	#Verandering van parameters naar preset
	observeEvent(input$preset, {
		if(input$preset == 1) {
			updateCheckboxInput(session, "p1", value = 0)
			updateCheckboxInput(session, "p2", value = 1)
			updateRadioButtons(session, "datatype", selected = 1)
			updateCheckboxInput(session, "p3", value = 0) }
		else if(input$preset == 2) {
			updateCheckboxInput(session, "p1", value = 1)
			updateCheckboxInput(session, "p2", value = 0)
			updateRadioButtons(session, "datatype", selected = 2)
			updateCheckboxInput(session, "p3", value = 1) }
		else if(input$preset == 0) {
			updateCheckboxInput(session, "p1", value = 0)
			updateCheckboxInput(session, "p2", value = 0)
			updateCheckboxInput(session, "p3", value = 0) }
	})
}


runApp(shinyApp(ui = ui, server = server), host = hostip, port = portnr)
