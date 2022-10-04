library(shiny)
library(bslib)
library(digest)

options(shiny.maxRequestSize=100*4096^2)
hostip <- "145.97.18.149"
portnr <- 6999

ui <- fluidPage(
	theme = bs_theme(version = 4, bootswatch = "minty"),
	navbarPage(
		titlePanel("MS-analyser"),
		tabPanel("Input",
			mainPanel(
			fluidRow(
				column(4),
				column(6,
					h5("Mass Spectrometry data:"),
					fileInput("msdata", label = "MS-data", multiple = TRUE, accept = ".raw"),
					tableOutput("bestanden"),
					radioButtons("datatype", label = "Method",
						choices = list("GC-MS" = 1, "LC-MS" = 2)),
					h5("Parameters:"),
					selectInput("preset", label = "Parameter presets", choices = list("None" = 0, "preset1" = 1, "preset2" = 2)),
					checkboxInput("p1", label = "p1"),
					checkboxInput("p2", label = "p2"),
					checkboxInput("p3", label = "p3"),
					textInput("email", label = "E-mail"),
					actionButton("parsubm", label = "Submit"),
					br(),
					verbatimTextOutput("hash")
				),
				column(2,
					textInput("jobid", label = "Existing JobID", width = "200"),
					actionButton("idsubm", label = "Submit")
				)
			)
		)
	),
		tabPanel("Job status",
			h2("Job status line")
		),
		tabPanel("Output")
	)
)

server <- function(input, output, session) {
	output$bestanden <- renderTable({
		file <- input$msdata
   		ext <- tools::file_ext(file$datapath)

   		req(file)
    		validate(need(ext == "raw", "Please upload a .raw file"))
	})
	observeEvent(input$parsubm, {
		genhash <- sapply(paste(toString(format(Sys.time(), "%S")), toString(input$datatype), toString(input$email), toString(runif(1))), digest, algo="md5")
		output$hash <- renderText({
			genhash 
		})
	})
	observeEvent(input$preset, {
		if(input$preset == 1) {
			updateCheckboxInput(session, "p1", value = 0)
			updateCheckboxInput(session, "p2", value = 1)
			updateCheckboxInput(session, "p3", value = 0) }
		else if(input$preset == 2) {
			updateCheckboxInput(session, "p1", value = 1)
			updateCheckboxInput(session, "p2", value = 0)
			updateCheckboxInput(session, "p3", value = 1) }
		else if(input$preset == 0) {
			updateCheckboxInput(session, "p1", value = 0)
			updateCheckboxInput(session, "p2", value = 0)
			updateCheckboxInput(session, "p3", value = 0) }
	})
}


runApp(shinyApp(ui = ui, server = server), host = hostip, port = portnr)
