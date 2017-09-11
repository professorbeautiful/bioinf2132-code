library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Adaptive Rejection Sampling"), 
  br(), 
  actionButton("next_point", label = "sample point(s)"),
  numericInput("nPointsToAdd", "# Points to sample", 1),
  fluidRow( column(width = 6,
      h2("g = p^(a-1) q^(b-1), proportional to the beta distribution"),
      plotOutput("samplingPlot")
    ),
    column(width = 6,
           plotOutput("qqPlot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  a = 2
  b = 4
  x = 5
  n = 20
  p.seq <- seq (0.01, 0.99, 0.01)
  # g = function(p, a.=a, b.=b) dbeta(p,a.,b.) * dbinom(x=x, size=n, p) 
  g = function(p, a.=a, b.=b) dbeta(p,a.,b.)  
  #### distribution we want to sample from
  mode_g = (a-1)/(a+b-2)
  max_g = g(mode_g)
  glog = function(...) log(g(...))   
  #### log of the  distribution we want to sample from
  glogdot = function(p, a.=a, b.=b) 
    ((a.-1)/p - (b.-1)/(1-p) ) 
  ####Derivative of log of g
  h = function(p) dbeta(p, 1, 1)  ## uniform
  r_h = function() rbeta(1, 1, 1) 
  max_h = 1
  M = max_g/max_h
  
  rValues = reactiveValues(Tset=c(0.25, 0.75),
                           acceptance = numeric(0),
                           proposalSet  = numeric(0),
                           h.u.Set = numeric(0)
  )
  observe({
    if(input$next_point > 0) isolate({  ## reaction
      for(i in 1:input$nPointsToAdd) {
      proposal = rbeta(1,1,1)
      rValues$proposalSet = c(rValues$proposalSet, proposal)
      acceptProb = g(proposal)/h(proposal)/M
      accept.u = runif(1)
      h.u = h(proposal) * accept.u
      accepted = (accept.u <= acceptProb)
      rValues$h.u.Set = c(rValues$h.u.Set, h.u)
      rValues$acceptance = c(rValues$acceptance, accepted)
      #cat('proposal ', proposal, '\n')
      }
    })
  })
  output$logSamplingPlot <- renderPlot({
    plot(p.seq, glog(p.seq), type="l", ylim = c(-3, 2))
    points(Tset, glog(Tset))
    title ("log of g, the likelihood (* prior) of interest")
  })
  output$samplingPlot <- renderPlot({
    plot(p.seq, h(p.seq), type="l", ylim = c(0, 1), lwd=2, xlim=c(0,1))
    lines(p.seq, g(p.seq)/M*max_h, type="l", col="green")
    #points(Tset, g(Tset))
    rug(rValues$proposalSet[rValues$acceptance==TRUE], col="blue")
    points(rValues$proposalSet, rValues$h.u.Set,
           col=c("red","blue")
           [rValues$acceptance+1]) 
  })
  output$qqPlot = renderPlot({
    if(length(rValues$proposalSet) > 0) {
      plot(col="blue", main="green=Truth, blue=sample smooth",
           density(rValues$proposalSet[rValues$acceptance==TRUE]))
      rug(rValues$proposalSet[rValues$acceptance==TRUE], col="blue")
      lines(p.seq, dbeta(p.seq, a, b), col="green")
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

