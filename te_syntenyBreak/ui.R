library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(plotly)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
        titlePanel("TE synteny app"),
        sidebarLayout(
            sidebarPanel(
                width=3,
                tabsetPanel(
                    tabPanel(
                        "Step 1: Basic msynt",
                        # inputs
                        "Select data set and load the basic mysnt before proceeding to TE plot.",
                        shinyWidgets::pickerInput("speciesA",label = "Species A", choices = species_name, options = list(title = "Select species")), # species A
                        shinyWidgets::pickerInput("speciesB", label = "Species B", choices = species_name, options = list(title = "Select species")), # species B
                        numericInput("lim",label = "min. #paralog per species", value = 4),
                        numericInput("nmin",label = "min. #msynt per scaffold", value = 10),
                        shinyWidgets::prettySwitch(inputId = "shared", label = "show shared OG only",value = TRUE, fill=TRUE),
                        shinyWidgets::prettySwitch(inputId = "showlines", label = "plot scaffold dot lines",value = TRUE, fill=TRUE),
                        actionButton("basic_cal", label = "1. Calculate and Plot"),
                        hr(),
                        h5("PNG size:"),
                        numericInput("plot_width",label = "width (mm)", value = 250),
                        numericInput("plot_height",label = "height (mm)", value = 250)
                    ),
                    tabPanel(
                        "Step 2: Add TE density",
                        "Calculate the TE content in flanking size for the TE plot. This may take long time depending on the size of genome and #msynt.
                        Once the calculataion is done, subsetting TE for visualization is quick.",
                        numericInput("flank_size",label = "Flanking size", value = 5000),
                        shinyWidgets::radioGroupButtons("upstream",label = "Flanking (+ strand):",choices = c("5'","3'"), justified = TRUE, checkIcon = list(yes=icon("ok",lib = "glyphicon"))),
                        actionButton("update", label = "2. Update count (Take long time)"), # calculate TE in window
                        hr(),
                        shinyWidgets::prettyRadioButtons(inputId = "te_species", label = "Plot te density of:", choices = c("x-axis","y-axis"), inline=TRUE, fill=TRUE),
                        shinyWidgets::prettyRadioButtons(inputId = "counttype", label = "Count type:", choices = c("bp","n"), inline=TRUE, fill=TRUE, selected="bp"),
                        uiOutput("teclass_picker"),
                        actionButton("update_plot", label = "3. Update TE plot (quick)"),
                        actionButton("update_plotly", label = "Refresh interactive plot")
                    )
                )
                
            ),
            
            mainPanel(
                # output
                width = 9,
                tabsetPanel(
                    tabPanel(
                        "basic msynt plot",
                        plotOutput("basic_plot",width = "100%", height = "900px") %>% withSpinner(),
                        downloadButton('basic_plot_dl', 'Download plot (.png)')
                    ),
                    tabPanel(
                        "static TE plot",
                        plotOutput("static_te_plot", width="100%", height = "900px") %>%  withSpinner(),
                        downloadButton('ggplot_dl', 'Download plot (.png)')
                    ),
                    tabPanel(
                        "interactive TE plot",
                        plotlyOutput("plotly_o", width="100%", height = "900px") %>%  withSpinner()
                    )
                ) # close tabsetPanel
            ) # close mainPanel
        )
    )
    
)
