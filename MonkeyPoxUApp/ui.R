#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("University Monkeypox Model"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            h3("University inputs"),
            sliderInput("Yale_undegrad_pop",
                        "Undergraduate student population:",
                        min = 1000,
                        max = 20000,
                        value = 6500),
            
            numericInput("isocapacity",
                         "Isolation capacity, daily", 0),
            
            numericInput("quarcapacity",
                         "Quarantine capacity, daily", 0),
            
            sliderInput("studentcontacts",
                        "Quarantined contacts per diagnosed student",
                        min = 0,
                        max = 20,
                        value = 0),
            
            sliderInput("diagrate",
                        "Daily probability of diagnosis of symptomatic infection",
                        min = 0,
                        max = 1,
                        value = 0),
            
            h3("Disease inputs"),  
            sliderInput("initial_inf",
                        "Infected students at start of semester",
                        min = 0,
                        max = 100,
                        value = 0),
            
            sliderInput("exograte",
                        "Number of infections coming in from outside of university community per month",
                        min = 0,
                        max = 30,
                        value = 0),

            sliderInput("R0_h",
                        "R0",
                        min = 0.8,
                        max = 1.5,
                        value = 1.1),




            
            submitButton(text = "Apply Changes", icon = NULL, width = NULL)),
        
        
        


        # Show a plot of the generated distribution
        mainPanel(
            #tabsetPanel(type="tabs",
                       # tabPanel("Summary results",
                                    fluidRow(valueBoxOutput("isocaplikelihood"),
                                             valueBoxOutput("isocaptime")),
                                    fluidRow(valueBoxOutput("quarcaplikelihood"),
                                            valueBoxOutput("quarcaptime")),
                                    fluidRow(valueBoxOutput("medianinfections")),
                                    fluidRow(splitLayout(cellWidths = c("33%", "33%","33%"), plotOutput("DPlot"),plotOutput("QPlot"), plotOutput("IPlot"), width=8)),
                        #tabPanel("Stochastic results",fluidRow(splitLayout(cellWidths = c("50%", "50%"), plotOutput("D1Plot"),plotOutput("Q1Plot"))),
                                    #fluidRow(splitLayout(cellWidths = c("50%", "50%"), plotOutput("I1Plot"),plotOutput("R1Plot"))),

                                    plotOutput("I_plot")
        ))
  #  )
))


