#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# ----------------------------------------------------- #
# Load library
# ----------------------------------------------------- #
library(shiny)
library(tidyverse)


# ----------------------------------------------------- #
# Preparation
# ----------------------------------------------------- #
data_set_2sp_bi <- function(t_length = 200,
                            inits = c(1, 1),
                            parms = c(3.8, -3.8, 3.5, -3.5, -0.02, -0.1),
                            time_series_plot = F)
{
  t <- t_length			# length of time series
  Y <- as.data.frame(matrix(rep(NA,3*t),ncol=3))
  colnames(Y) <- c("t","X","Y")
  
  Y[,1] <- 1:t
  Y[1,2:3] <- inits
  
  r1 <- parms[1]
  s1 <- parms[2]
  r2 <- parms[3]
  s2 <- parms[4]
  beta1 <- parms[5]
  beta2 <- parms[6]
  
  for(i in 1:(t-1)){
    Y[i+1,2] <- Y[i,2]*(r1 +   s1*Y[i,2] + beta1*Y[i,3])
    Y[i+1,3] <- Y[i,3]*(r2 + beta2*Y[i,2] + s2*Y[i,3])
    if (Y[i+1,2] <= 0) Y[i+1,2] <- 0
    if (Y[i+1,3] <= 0) Y[i+1,3] <- 0
  }
  
  if(time_series_plot){
    plot(Y[,1], Y[,2], type="l", col="royalblue", lwd=2, ylim=range(0,1.5),
         ylab="Abundance", xlab="Time", main="Dynamics of X, Y")
    lines(Y[,1], Y[,3], type="l", col="red3", lwd=1)
  }
  
  return(Y)
}


# ----------------------------------------------------- #
# UI for application
# ----------------------------------------------------- #
ui <- fluidPage(

    # Application title
    titlePanel("OCES3160 L16: Discrete time competition model"),

    fluidRow(column(12, p(HTML("2023-11-10 Masayuki Ushio"),
                          style = "font-size:14px; color:gray;"))),
    fluidRow(column(12, withMathJax("$$X_{t+1} = ax X_t - bx X_t^2 - cx X_t Y_t$$"))),
    fluidRow(column(12, withMathJax("$$Y_{t+1} = ay Y_t - by Y_t^2 - cy Y_t X_t$$"))),
    
    # ----------------------------------------------- #
    # Show time series
    # ----------------------------------------------- #
    fluidRow(column(12,  plotOutput('ts_plot'))),
    
    # ----------------------------------------------- #
    # Growth rate
    # ----------------------------------------------- #
    fluidRow(column(12, h3("1. Specify growth rates"))),
    fluidRow(
      column(6,
             wellPanel(
               # Input: Simple integer interval for the time series plot
               sliderInput("ax",
                           h5(withMathJax("$$ax : \\text{Growth rate of } X$$")),
                           min = 0, max = 5, value = 3.8, step = 0.1)
            )       
        ),
      
      column(6,
             wellPanel(
               # Input: Simple integer interval for the time series plot
               sliderInput("ay",
                           h5(withMathJax("$$ay : \\text{Growth rate of } Y$$")),
                           min = 0, max = 5, value = 3.5, step = 0.1)
             )       
      )
    ),

    # ----------------------------------------------- #
    # Density dependence
    # ----------------------------------------------- #
    fluidRow(column(12, h3("2. Specify density dependence effects"))),
    fluidRow(
      column(6,
             wellPanel(
               # Input: Simple integer interval for the time series plot
               sliderInput("bx",
                           h5(withMathJax("$$bx : \\text{Negative density dependence of } X$$")),
                           min = 0, max = 5, value = 3.8, step = 0.1)
             )       
      ),
      
      column(6,
             wellPanel(
               # Input: Simple integer interval for the time series plot
               sliderInput("by",
                           h5(withMathJax("$$by : \\text{Negative density dependence of } Y$$")),
                           min = 0, max = 5, value = 3.5, step = 0.1)
             )       
      )
    ),

    # ----------------------------------------------- #
    # Competition strength
    # ----------------------------------------------- #
    fluidRow(column(12, h3("3. Specify the strength of competition"))),
    fluidRow(
      column(6,
             wellPanel(
               # Input: Simple integer interval for the time series plot
               sliderInput("cx",
                           h5(withMathJax("$$cx : \\text{Competition from } X \\text{ to } Y$$")),
                           min = 0, max = 5, value = 0.02, step = 0.1)
             )       
      ),
      
      column(6,
             wellPanel(
               # Input: Simple integer interval for the time series plot
               sliderInput("cy",
                           h5(withMathJax("$$cy : \\text{Competition from } Y \\text{ to } X$$")),
                           min = 0, max = 5, value = 0.1, step = 0.1)
             )       
      )
    ),
    
)



# ----------------------------------------------------- #
# Define server logic
# ----------------------------------------------------- #
server <- function(input, output) {

    output$ts_plot <- renderPlot({
        # Collect parameters
        ## Growth
        ax <- input$ax
        ay <- input$ay
        ## Density dependence
        bx <- -input$bx
        by <- -input$by
        ## Competition
        cx <- input$cx
        cy <- input$cy
        
        # Generate population dynamics
        d <- data_set_2sp_bi(t_length = 200,
                             inits = c(0.1, 0.1),
                             #parms = c(3.8, -3.8, 3.5, -3.5, -0.02, -0.1), # => Chaotic
                             parms = c(ax, bx, ay, by, cx, cy),
                             time_series_plot = F)
        
        # Visualize dynamics
        d_long <- tidyr::pivot_longer(d, cols = c(-t))
        ggplot(d_long, aes(x = t, y = value, color = name)) +
          geom_point() + geom_line() +
          scale_color_manual(values = c("red3", "royalblue"), name = "Species") +
          xlab("Time") + ylab("Population size") +
          theme(legend.position = "top", text = element_text(size = 14)) + coord_cartesian(ylim = c(0, 3)) +
          #ggtitle("Population size of X and Y") +
          geom_hline(yintercept = 0, linetype = 2)
        
        # plot(d[,1], d[,2], type="l", col="royalblue", lwd=2, ylim=range(0,3),
        #      ylab="Abundance", xlab="Time", main="Dynamics of X, Y")
        # lines(d[,1], d[,3], type="l", col="red3", lwd=1)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)


