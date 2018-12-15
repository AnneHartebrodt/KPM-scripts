#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
require(plotly)
SAVE_PLOT <- TRUE
# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Old Faithful Geyser Data"),
 
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      fileInput("jaccard_table", "Evaluation", multiple = T),
      fileInput("annotation_table", "Annotation", multiple = T),
      fileInput("fdr", "Thresholds", multiple = T),
      actionButton("action", label = "Action"),
      uiOutput("check"),
      actionButton("update", label = "Update Plot")
      
    ),
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
       #textOutput("value"),
      tabPanel("Jaccards", plotOutput("eval")),
      tabPanel("Best", plotOutput("best")),
      tabPanel("Thresholds", plotOutput("thresholds")),
      tabPanel("Dots", plotOutput("weirdDots"))
      #plotOutput("plo"),
      #textOutput("pr")
    )
    )
  )
))
