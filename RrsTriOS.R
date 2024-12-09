# Specify required packages
required_packages <- c("shiny", "shinydashboard", "data.table", "DT")

# Install missing packages
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Load required packages
lapply(required_packages, library, character.only = TRUE)

# Prevent scientific notation
options(scipen=999)

ui <- fluidPage(
       titlePanel("Process Rrs above water TriOS Data"),
       sidebarLayout(
         sidebarPanel(
           textInput("Ed_device","Ed_device",""),
           textInput("Lw_device","Lw_device",""),
           fileInput("Rrs.file","Upload TRIOS Data"),
           textInput("outdir", "Output Directory", ""),
           textInput("Experiment","Experiment","Missouri_Reservoirs_RSWQ"),
           textInput("Cruise","Cruise","MO_RSWQ_2024"),
           textInput("Station","Station","NA"),
           textInput("Latitude","Latitude","-99"),
           textInput("Longitude","Longitude","-99"),
           textInput("Secchi","Secchi Depth","-99"),
           textInput("Depth","Water Depth","-99"),
           textInput("Comments","Comments","-99"),
           actionButton("SaveFile","Save Data",icon = icon("play-circle"))),
         mainPanel(uiOutput("tb"))
       )
     )


server <- function(input, output, session) {

################################################################################
#  Functions called by code

# read.TRIOS reads in raw TRIOS file and returns a list of
# Spectra, TimeStamp, Pressure (if present), and Inclination (if present)

read.TRIOS <- function(TRIOS.file, Sensor){
  
  con = file(TRIOS.file, "r")
  
  dat  = readLines(con)
  
  # Find indices of Calibrated Spectrum
  spec = which(dat=="IDDataTypeSub1     = CALIBRATED", arr.ind=T)
  for (i in 1:length(spec)){ 
    foo = spec[i] - 1
    bar = spec[i] - 2
    if (dat[foo] != "IDDataType         = SPECTRUM") {spec[i] = NA}
    IDDevice  = gsub("IDDevice           = ", "", dat[bar])
    if (IDDevice != Sensor) {spec[i] = NA}
  }
  
  spec = spec[is.finite(spec)]
  
  # Declare variables
  n           = length(spec)
  wv          = c(350:950)
  TimeStamp   = rep(NA, n)
  data        = array(NA, c(n, length(wv)))
  
  # Retrieve Device, Timestamp, and Spectral Data
  for (i in 1:n){
    foo          = spec[i] + 3
    TimeStamp[i] = gsub("DateTime           = ", "", dat[foo])
    
    # Loop through file until beginning of spectral data
    while (dat[foo]!= "[DATA]"){
      foo = foo + 1
    }
    
    foo          = matrix(as.numeric(unlist(strsplit(dat[(foo+2):(foo+256)], split=" "))),
                          nrow=255, byrow = TRUE)
    data[i,]    = approx(x=foo[,2], y=foo[,3], xout=wv)$y/10 # unit conversion (uW/cm2/nm)
    
  }
  
  close(con)
  
  return(list(Spectra = data, TimeStamp = TimeStamp))
}

calc.RRs <- function(){
  
  if(is.null(input$Rrs.file)){return()}
  
  Ed = read.TRIOS(input$Rrs.file$datapath, input$Ed_device)
  Lw = read.TRIOS(input$Rrs.file$datapath, input$Lw_device)

  Ed.spectra = Ed$Spectra
  Lw.spectra = Lw$Spectra
  
  # Find TimeStamps where both sensors are measuring
  Ed.ts      = Ed$TimeStamp
  Lw.ts      = Lw$TimeStamp
  timestamps = c(Ed.ts, Lw.ts)
  foo        = duplicated(timestamps)
  timestamps = timestamps[foo]
  
  #Populate Data Frame [Sensors, Inclination, Spectra(350:950)]
  RRs = array(NA, c(length(timestamps),601))
  
  for (i in 1:length(timestamps)){
    foo          = which(Lw.ts==timestamps[i], arr.ind=T)
    bar          = which(Ed.ts==timestamps[i], arr.ind=T)
    RRs[i,1:601] = round(Lw.spectra[foo,] / Ed.spectra[bar,],5)
  }
  
  df = data.frame(RRs, row.names=timestamps)
  names(df) = c(paste(c(350:950)))
  return(df) 
  
}

################################################################################
## Below code to display the structure of the input file object

# output$plot_RRs finds timestamps where both Ed and Lu are present,
# then computes and plots RRs

output$plot_Ed <- renderPlot({
  
  if(is.null(input$Rrs.file)){return()}
  Ed = read.TRIOS(input$Rrs.file$datapath, input$Ed_device)$Spectra
  wv = c(350:950)
  
  s1 = input$Rrs_data_rows_selected
  s2 = input$Rrs_data_rows_all 
  
  if (length(s1)>0){
  
    ind = s1[1]
    maxEd = max(Ed[s1,],na.rm=TRUE)
    
    plot(wv, Ed[ind,1:601], type="l", axes=T,ylim=c(0,maxEd),
         xlab="Wavelength (nm)", ylab=expression(E[D]))
    
    if (length(s1)>1) {
      for (i in 1:length(s1)){
        ind = s1[i]
        lines(wv, Ed[ind,])
      }
    }
  }else{
    plot(wv, Ed[1,], type="l", axes=T,
         xlab="Wavelength (nm)", ylab=expression(E[D]))
    
    if (dim(Ed)[1]>1) {
      for (i in 1:dim(Ed)[1]){
        lines(wv, Ed[i,])
      }
    }
  }  
})

output$plot_Lw <- renderPlot({
  
  if(is.null(input$Rrs.file)){return()}
  
  Lw = read.TRIOS(input$Rrs.file$datapath, input$Lw_device)$Spectra
  Lw = data.frame(Lw)
  wv = c(350:950)
  
  selected_rows = input$Rrs_data_rows_selected
  if (length(selected_rows)==0){
    selected_rows = input$Rrs_data_rows_all 
  }

  if (length(selected_rows)==0){
    Lw_mean = 0
    Lw_sd = 0
  } else {  
    Lw_mean = apply(Lw[selected_rows, 401:451], 2, mean,na.rm=TRUE) 
    Lw_sd = apply(Lw[selected_rows, 401:451], 2, sd,na.rm=TRUE)
  }  
  maxLw = ceiling(max(Lw[selected_rows,],na.rm=TRUE))
  
  # Empty plot to start  
  plot(wv, rep(NA, length(wv)), type="n", axes=F, ylim=c(0, maxLw),
       xlab="Wavelength (nm)", ylab=expression(L[W]))
  
  for (i in selected_rows){
    
    if (any(Lw[i,401:451] > Lw_mean + (3*Lw_sd)) |
        any(Lw[i,401:451] < Lw_mean - (3*Lw_sd))) { 
      lines(wv, Lw[i,], col='red', lwd=2)
    } else { # Good data
      lines(wv, Lw[i,], col='black')
    }
    if (any(Lw[i,251:351] <= 0.0005)){ # Highlighting Lw < 0.001
      lines(wv, Lw[i,], col='orange', lwd=2)
    }
  }
  
  axis(side=1, at=seq(350, 950, by=50), pos=0)
  axis(side=2, at=seq(0, maxLw, by=1), pos=350)  
  
})

# highlight selected rows in the scatterplot
output$plot_Rrs = renderPlot({
  
  if(is.null(input$Rrs.file)){return()}
  
  wv  = c(350:950)
  Rrs = calc.RRs()
  
  selected_rows = input$Rrs_data_rows_selected
  if (length(selected_rows)==0){
    selected_rows = input$Rrs_data_rows_all 
  }
  
  if (length(selected_rows)==0){
    Rrs_mean = 0
    Rrs_sd = 0
  } else {  
    Rrs_mean = apply(Rrs[selected_rows, 401:601], 2, mean, na.rm=TRUE) 
    Rrs_sd = apply(Rrs[selected_rows, 401:601], 2, sd, na.rm=TRUE)
  }  
  
  # For nice axis in increments of 0.0025
  maxRrs = ceiling(max(Rrs[selected_rows,], na.rm=TRUE)/0.0025)  * 0.0025
  
  # Empty plot to start  
  plot(wv, rep(NA, length(wv)), type="n", axes=F, ylim=c(0, maxRrs),
       xlab="Wavelength (nm)", ylab=expression(R[Rs]~(sr^-1)))
  
  for (i in selected_rows){
    
    if (any(Rrs[i,401:451] > Rrs_mean + 3*(Rrs_sd)) |
        any(Rrs[i,401:451] < Rrs_mean - 3*(Rrs_sd))) { # Highlighting bad data
      lines(wv, Rrs[i,], col='red', lwd=2)
    } else { # Good data
      lines(wv, Rrs[i,], col='black')
    }
  }
  
 axis(side=1, at=seq(350,950, by=50), pos=0)
 axis(side=2, at=seq(0, maxRrs, by=0.0025), pos=350)    
 
  
})


output$Rrs_data <- renderDataTable(datatable(calc.RRs(), extensions = 'Scroller', 
                          options=list(columnDefs = list(list(width = '200px', targets = c(1, 201))),
                                       dom = 't',
                                       server = TRUE,
                                       deferRender = TRUE,
                                       scrollX = TRUE,
                                       scrollY = 500,
                                       scroller = TRUE)) %>%
                            formatRound(c(1:601),5))


observeEvent(input$SaveFile, {
  
  if(is.null(input$Rrs.file)){return()}
  
  if(is.null(input$outdir)){
    outdir = getwd('')
  } else {
    outdir = input$outdir
  }
  
  wv  = c(350:950)
  RRs = calc.RRs()
  Ed = read.TRIOS(input$Rrs.file$datapath, input$Ed_device)$Spectra 
  Lw = read.TRIOS(input$Rrs.file$datapath, input$Lw_device)$Spectra
  
  s1 = input$Rrs_data_rows_selected
  s2 = input$Rrs_data_rows_all 

  if (length(s1)>0){
    
    Rrs = RRs[s1,]
    Rrs_mean = apply(Rrs, 2, mean, na.rm=T)
    Rrs_sdev = apply(Rrs, 2, sd, na.rm=T)
    Ed = Ed[s1,]
    Ed_mean = apply(Ed, 2, mean, na.rm=T)
    Ed_sdev = apply(Ed, 2, sd, na.rm=T)
    Lw = Lw[s1,]
    Lw_mean = apply(Lw, 2, mean, na.rm=T)
    Lw_sdev = apply(Lw, 2, sd, na.rm=T)
    out = cbind(wv, Rrs_mean, Rrs_sdev, Ed_mean, Ed_sdev, Lw_mean, Lw_sdev)

    timestamp0 = min(row.names(Rrs))
    timestamp1 = max(row.names(Rrs))

    filename=paste0(input$Experiment,"_",input$Cruise,"_",input$Station,"_Rrs_above_water_",gsub("-","",substr(timestamp0,1,10)),"_",gsub(":","",substr(timestamp0,12,19)),'_R1.csv')
    
    fileConn<-file(paste0(outdir,filename))
    writeLines(c("/begin_header",
                 "/investigators=Greg_Silsbe,Lorena_Silva,Rebecca_North",
                 "/affiliations=University_of_Maryland_Center_for_Environmental_Science_and_University_of_Missouri_Columbia",
                 "/contact=gsilsbe@umces.edu",
                 paste0("/experiment=",input$Experiment),
                 paste0("/cruise=",input$Cruise),
                 paste0("/station=",input$Station), 
                 paste0("/data_file_name=",filename),
                 "/documents=Rrs_Seabass_methodology.pdf,Calibration_certificate_850C.pdf,Calibration_certificate_8760.pdf,Missouri_Reservoirs_RSWQ_2023_checklist.txt",
                 "/data_type=above_water",
                 "/calibration_files=Back_SAM_850C.dat,Cal_SAM_850C.dat,SAM_850C.ini,Back_SAM_8760.dat,Cal_SAM_8760.dat,SAM_8760.ini",
                 "/calibration_date=20220214",
                 paste0("/start_date=",gsub("-","",substr(timestamp0,1,10))),
                 paste0("/end_date=",gsub("-","",substr(timestamp1,1,10))),
                 paste0("/start_time=",gsub(" ", "",substr(timestamp0,11,19)),"[GMT]"),
                 paste0("/end_time=",gsub(" ", "",substr(timestamp1,11,19)),"[GMT]"),
                 paste0("/north_latitude=",input$Latitude,"[DEG]"),
                 paste0("/south_latitude=",input$Latitude,"[DEG]"),
                 paste0("/east_longitude=",input$Longitude,"[DEG]"),
                 paste0("/west_longitude=",input$Longitude,"[DEG]"),
                 "/missing=-99",
                 "/delimiter=comma",
                 "/cloud_percent=-99",
                 "/data_status=final",
                 "/instrument_manufacturer=TRIOS",
                 "/instrument_model=RAMSES",
                 paste0("/secchi_depth=",input$Secchi),
                 paste0("/water_depth=",input$Depth),
                 "/measurement_depth=0",
                 "/wave_height=-99",
                 paste0("!comments=",input$Comments),
                 "/fields=wavelength,Rrs,Rrs_sd,Ed,Ed_sd,Lw,Lw_sd",
                 "/units=nm,1/sr,1/sr,uW/cm^2/nm,uW/cm^2/nm,uW/cm^2/nm/sr,uW/cm^2/nm/sr",
                 "/end_header"), 
                 fileConn)
    close(fileConn)
    write.table(out, paste0(input$outdir, filename), quote=F, sep=",", row.names = F,  append = TRUE,
                col.names = F)
  }
})

################################################################################
## MainPanel tabset renderUI code ##
# the following renderUI is used to dynamically generate the tabsets when the file is loaded. 
# Until the file is loaded, app will not show the tabset.


output$tb <- renderUI({
   fluidPage(
    fluidRow(column(12, DT::dataTableOutput('Rrs_data')), style='50px'),
    fluidRow(column(3, plotOutput('plot_Ed'), style='100px'),
             column(3, plotOutput('plot_Lw'), style='100px'),
             column(6, plotOutput('plot_Rrs'), style='100px'))
    )
})
}

shinyApp(ui = ui, server = server)

