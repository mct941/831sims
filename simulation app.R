library(dplyr)
library(ggplot2)
library(shiny)


read.nonmem <- function(file, n=-1) {
  ## auxiliary function to split text lines at blanks
  my.split <- function(line, numeric=FALSE) {
    pieces <- unlist(strsplit(line, split=" +"))[-1]
    if( numeric )
      return(as.numeric(pieces))
    else
      return(pieces)
  }
  
  cat(sprintf("Reading NONMEM data from '%s'\n", file))
  lines <- readLines(file, n) # read file as text
  cat(sprintf("- %d lines\n", length(lines)))
  
  idx <- substring(lines,1,1)!="T"
  cat(sprintf("- %d tables\n", sum(!idx)))
  lines <- lines[idx] # strip lines starting with T (TABLE NO ...)
  
  ## do we have header lines??
  if( length(grep("^ +[A-Za-z]", lines[1]))>0 ) { # yes!
    data <- sapply(lines[lines!= lines[1]], my.split, numeric=TRUE)
    header <-  my.split(lines[1])
    cat(sprintf("- file has column names (%s)\n", paste(header,collapse=", ")))
  } else {                                        # no
    data <- sapply(lines, my.split, numeric=TRUE)
    header <- sprintf("Column%02d", 1:nrow(data)) # make fake column names
    cat("- file has NO header names - creating names Column01, Column02, ...\n")
  }
  cat(sprintf("- %d columns\n", nrow(data)))
  
  ## transpose data and make a data.frame
  df <- data.frame(data[1,])
  for( i in 2:nrow(data))
    df <- cbind(df, data[i,])
  
  ## set column and row names
  rownames(df) <- NULL
  colnames(df) <- header
  cat("ok.\n")
  return(df)
}

sim.out <- read.nonmem(file.path("U:","Projects","TAK-831","1001","tak831_2-cmpt_dose-v2-cl-ka_v2-sim.tbl")) %>% 
  filter(EVID==0) %>% unite(SUBJ,REP,ID) %>% group_by(SUBJ,DOSE) 

## app.R ##
server <- function(input, output) {
  simdta <- reactive({
    if(input$max==TRUE) dta <- sim.out %>% mutate(RO=IPRED/(IPRED+input$ec50)*100)  %>%
                                           summarise(CMAX=max(IPRED),RMAX=max(RO),CMIN=min(IPRED),RMIN=min(RO))
    else  dta <- sim.out %>% mutate(RO=IPRED/(IPRED+input$ec50)*100) %>% filter(TIME==input$filter) 
    return(dta)
  })
  
  output$distPlot <- renderPlot({
    if(input$pwd=="takeda") {
      if(input$max==TRUE) p <- ggplot(simdta(), aes(x=factor(DOSE),y=RMAX,group=factor(DOSE)))  else
        p <- ggplot(simdta(), aes(x=factor(DOSE),y=RO,group=factor(DOSE))) 
      p <- p + geom_boxplot() + slide_theme() + xlab("Dose (mg)") + ylab("Enzyme Occupancy (%)")
      if(input$ref==TRUE) p <- p + geom_hline(yintercept=input$cutoff, color="red", linetype="dashed", size=rel(2))
      return(p)
    }
  })
}

ui <- fluidPage(titlePanel("TAK-831 Simulations"),
  textInput("pwd","Enter password:"),              
  fluidRow(
  column(2, wellPanel(
    sliderInput("ec50", "EC50:", min = 0, max = 100, value = 12.7),
    checkboxInput("max","Use Cmax",value=TRUE),
    conditionalPanel("input.max==false",sliderInput("filter", "Use concentration @ time=",min=0,max=24,value=2)),
    checkboxInput("ref","Add reference line",value=FALSE),
    conditionalPanel("input.ref==true",sliderInput("cutoff", "Set threshold to ",min=0,max=100,value=50))
    ),
    tags$small(paste0(
      "Note: PK model was developed, based on SRD/MRD study data and concentration-occupancy relationship was established from PET study data.",
      " Model-based simulations were performed to evaluate steady state brain levels of DAO enzyme occupancy at various doses (n=100).",
      " Sensitivity analysis can be evaluated within this application to evaluate the effects of EC50 and timepoint selection."
    ))
  ),
  column(10,
    wellPanel(plotOutput("distPlot"))
  )
))

shinyApp(ui = ui, server = server)

