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
                            parms = c(4, 2, 3, 2, 0.1, 0.1),
                            time_series_plot = F)
{
  t <- t_length			# length of time series
  Y <- as.data.frame(matrix(rep(NA,3*t),ncol=3))
  colnames(Y) <- c("t","N","P")

  Y[,1] <- 1:t
  Y[1,2:3] <- inits

  rn <- parms[1]
  bn <- parms[2]
  rp <- parms[3]
  bp <- parms[4]
  s1 <- parms[5]
  e1 <- parms[6]

  for(i in 1:(t-1)){
    Y[i+1,2] <- rn*Y[i,2] - bn*Y[i,2]^2 - s1*Y[i,2]*Y[i,3]
    Y[i+1,3] <- rp*Y[i,3] - bp*Y[i,3]^2 + e1*s1*Y[i,2]*Y[i,3]
    if (Y[i+1,2] <= 0) Y[i+1,2] <- 0
    if (Y[i+1,3] <= 0) Y[i+1,3] <- 0
  }

  if(time_series_plot){
    plot(Y[,1], Y[,2], type="l", col="royalblue", lwd=2, ylim=range(0,1.5),
         ylab="Abundance", xlab="Time", main="Dynamics of N, P")
    lines(Y[,1], Y[,3], type="l", col="red3", lwd=1)
  }

  return(Y)
}


# ----------------------------------------------------- #
# UI for application
# ----------------------------------------------------- #
ui <- fluidPage(

  # Application title
  titlePanel("OCES3160 L19: Discrete time prey-predator model"),

  fluidRow(column(12, p(HTML("2023-11-15 Masayuki Ushio"),
                        style = "font-size:14px; color:gray;"))),
  fluidRow(column(12, withMathJax("$$N_{t+1} = r_n N_t - b_n N_t^2 - s P_t N_t$$"))),
  fluidRow(column(12, withMathJax("$$P_{t+1} = r_p P_t - b_p P_t^2 + es P_t N_t$$"))),

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
             sliderInput("rn",
                         h5(withMathJax("$$r_n : \\text{Growth rate of } N$$")),
                         min = 0, max = 5, value = 4, step = 0.1)
           )
    ),

    column(6,
           wellPanel(
             # Input: Simple integer interval for the time series plot
             sliderInput("rp",
                         h5(withMathJax("$$r_p : \\text{Growth rate of } P$$")),
                         min = 0, max = 5, value = 3, step = 0.1)
           )
    )
  ),

  # ----------------------------------------------- #
  # Density dependence
  # ----------------------------------------------- #
  fluidRow(column(12, h3("2. Specify the density-dependent effect"))),
  fluidRow(
    column(6,
           wellPanel(
             # Input: Simple integer interval for the time series plot
             sliderInput("bn",
                         h5(withMathJax("$$b_n : \\text{Negative density dependence of } N$$")),
                         min = 0, max = 5, value = 2, step = 0.1)
           )
    ),

    column(6,
           wellPanel(
             # Input: Simple integer interval for the time series plot
             sliderInput("bp",
                         h5(withMathJax("$$b_p : \\text{Negative density dependence of } P$$")),
                         min = 0, max = 5, value = 2, step = 0.1)
           )
    )
  ),

  # ----------------------------------------------- #
  # Competition strength
  # ----------------------------------------------- #
  fluidRow(column(12, h3("3. Specify the prey-predator relationship"))),
  fluidRow(
    column(6,
           wellPanel(
             # Input: Simple integer interval for the time series plot
             sliderInput("s1",
                         h5(withMathJax("$$s : \\text{Consumption rate of } N \\text{ by } P$$")),
                         min = 0, max = 5, value = 0.1, step = 0.1)
           )
    ),

    column(6,
           wellPanel(
             # Input: Simple integer interval for the time series plot
             sliderInput("e1",
                         h5(withMathJax("$$e : \\text{Conversion efficiency of } P$$")),
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
    rn <- input$rn
    rp <- input$rp
    ## Density dependence
    bn <- input$bn
    bp <- input$bp
    ## Competition
    s1 <- input$s1
    e1 <- input$e1

    # Generate population dynamics
    d <- data_set_2sp_bi(t_length = 200,
                         inits = c(1, 1),
                         parms = c(rn, bn, rp, bp, s1, e1),
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


