#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("MonkeyPox U"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("Yale_undegrad_pop",
                        "Undergraduate student population:",
                        min = 1000,
                        max = 20000,
                        value = 6500),

            sliderInput("initial_inf",
                        "Infected students at start of semester",
                        min = 0,
                        max = 100,
                        value = 0),

            sliderInput("R0_h",
                        "R0",
                        min = 0.8,
                        max = 1.5,
                        value = 1.1),

            sliderInput("studentcontacts",
                        "Quarantined contacts per diagnosed student",
                        min = 0,
                        max = 20,
                        value = 6),

            sliderInput("quarduration",
                        "Quarantine duration (days)",
                        min = 0,
                        max = 20,
                        value = 14),

            sliderInput("diagrate",
                        "Daily probability of diagnosis of symptomatic infection",
                        min = 0,
                        max = 1,
                        value = 0.1),

            sliderInput("exograte",
                        "Number of infections coming in from outside of university community per week",
                        min = 0,
                        max = 10,
                        value = 0),
            
            submitButton(text = "Apply Changes", icon = NULL, width = NULL)),
        
        
        


        # Show a plot of the generated distribution
        mainPanel(
            # plotOutput("QPlot"),
            # plotOutput("DPlot"),
            # plotOutput("IPlot"),
            # plotOutput("RPlot")
            plotOutput("I_plot")
        )
    )
))
