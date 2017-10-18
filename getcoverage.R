library(ballgown)
library(ggplot2)
library(gtable)
library(grid)

### path_transcriptsID: transcriptsID is a file generated after parsing "*.tracking" output file 
### from gffcompare using the collowing bash command:
  # cat cuffcmp-C.tracking | cut -f1,2,4 > cuffcmp-C.ref.transcriptID.txt

### whole_tx_table and transcript_cov are dataframes obtained using ballgown, as below (copied from Ballgown tutorial) :
  # library(ballgown)
  # data_directory <- "stringTie_ballgown/stringTie-example" # path to coverage data obtained using stringtie
  # pattern <- "ballgown-BGI_ILU_"    #pattern of "stringtie -B" output folders, used as input for Ballgown. 
  # bg = ballgown(dataDir=data_directory, samplePattern=pattern, meas='all')
  # transcript_cov = texpr(bg, 'cov')
  # whole_tx_table = texpr(bg, 'all') 

### xmax and stp are plotting options, defining length of x axis (xmax) and setp (stp) 

getCoverage <- function(whole_tx_table, transcript_cov, path_transcriptsID, outname, complete = FALSE, xmax = 4000, stp = 10){
  gffcompareclass = read.table(path_transcriptsID)
  colnames(gffcompareclass) = c("trID", "geneID","class")
  gffcompareclass = gffcompareclass[order(gffcompareclass$trID),]
  #create a coverage table
  covTable = data.frame(transcript=whole_tx_table[,6], cov=rowSums(transcript_cov), class = as.character(""), 
                        stringsAsFactors = FALSE)
  covTable= covTable[covTable$cov > 0,]
  covTable=covTable[order(covTable$transcript),]
  covTable$class = gffcompareclass[gffcompareclass$trID %in% covTable$transcript,3]
  covTable$cov = as.numeric(covTable$cov)
  theme_set(theme_gray(base_size = 12))
  ## combining two graphs - histogram and cumulative
  for (i in c("=", "j", "u")) {
    if (i == "=") {
      a = "conserved"
    } 
    else if (i == "j") {
      a = "novel isoforms"
    }
    else if (i == "u") {
      a = "unknown"
    } else {a = i}
    plotcum= ggplot(data=covTable[covTable$class == i,][,c(3,2)], aes(cov)) +
      coord_cartesian(xlim =  c(0, xmax)) +
      stat_ecdf(geom = "step", pad = FALSE) + 
      labs(title="") +
      labs(x="Reads per bp", y="Cumulative transcr. density")
    plothist= ggplot(data=covTable[covTable$class == i,][,c(3,2)], aes(cov)) + 
      geom_histogram(breaks=seq(0, xmax, by = stp), aes(fill=..count..)) +
      labs(title=paste("Transfrag class:", a)) + labs(x="", y="Transcripts") 
    g1 <- ggplotGrob(plothist)
    g2 <- ggplotGrob(plotcum)
    g2 <- gtable_add_cols(g2, unit(0,"mm"))
    g2 <- gtable_add_cols(g2, unit(0,"mm")) # add a column for missing legend
    g <- rbind(g1, g2, size="first") # stack the two plots
    g$widths <- unit.pmax(g1$widths, g2$widths) # use the largest widths
    # center the legend vertically
    g$layout[grepl("guide", g$layout$name),c("t","b")] <- c(1,nrow(g))
    grid.newpage()
    grid.draw(g)
    ggsave(g, filename=paste("Cumulative_densities_histogram_", outname,"_", a, ".png", sep = ""), width = 7, height = 5, units = "in")
  }
  #theme_set(theme_gray(base_size = 14))
  plot_all=ggplot(data=covTable[covTable$class %in% c("=", "j", "u"),][,c(3,2)], aes(cov, colour = class))  +
    stat_ecdf(geom = "line", pad = FALSE) + coord_cartesian(xlim =  c(0, 10000)) + 
    geom_vline(xintercept = xmax, linetype = "dotted", colour = "gray50") +
    labs(x="Reads per bp", y="Cumulative transcript density")
  ggsave(plot_all, filename=paste("Cumulative_densities_", outname, "_select.png" ), width = 7, height = 5, units = "in")
  
  if (complete) {
    for (i in c("=", "j", "e", "i", "o", "p", "s", "u", "x", "c")) {
    print(paste("histogram_", i, ".png", sep = ""))
    if (i == "=") {a = "conserved"} else {a = i}
    histg = ggplot(data=covTable[covTable$class == i,][,c(3,2)], aes(cov)) + 
      geom_histogram(breaks=seq(0, xmax, by =stp), aes(fill=..count..)) +
      labs(title=paste("Transfrag class:", a)) + labs(x="Reads per bp", y="Transcripts") 
    ggsave(histg,
           filename=paste("histogram_", outname,"_", a, ".png", sep = ""), 
           width = 7, height = 5, units = "in")
    }
    
    plot_all= ggplot(data=covTable[,c(3,2)], aes(cov, colour = class)) + 
      coord_cartesian(xlim =  c(0, 10000)) +
      stat_ecdf(geom = "step", pad = FALSE) + 
      geom_vline(xintercept = xmax, linetype = "dotted", colour = "gray50") +
      labs(x="Reads per bp", y="Cumulative transcript density")
    ggsave(plot_all, filename= paste("Cumulative_densities_", outname,"_","All.png"), 
           width = 7, height = 5, units = "in")
    
    for (i in c("=", "j", "e", "i", "o", "p", "s", "u", "x", "c")) {
      print(paste("histogram_strt_", i, ".png", sep = ""))
      if (i == "=") {a = "conserved"} else {a = i}
      plot_i= ggplot(data=covTable[covTable$class == i,][,c(3,2)], aes(cov)) + 
        coord_cartesian(xlim =  c(0, 10000)) +
        geom_vline(xintercept = xmax, linetype = "dotted", colour = "gray50") +
        stat_ecdf(geom = "step", pad = FALSE) + 
        labs(title=paste("Transfrag class:", a)) +
        labs(x="Reads per bp", y="Cumulative transcript density")
      ggsave(plot_i,filename=paste("cumulative_density_", outname,"_", a, ".png", sep = ""), 
             width = 7, height = 5, units = "in")
    } 
  }
  print(paste("plots saved in ", getwd()))
  return(covTable)
}

# usage example
setwd("/path/for/workingDirectory/")

coverageTable = getCoverage(whole_tx_table, transcript_cov,
                 "/path/for/ref.transcriptID.txt",
                 "testName", complete = FALSE)
coverageTable = getCoverage(whole_tx_table, transcript_cov,
                 "/path/for/ref.transcriptID.txt",
                 "testName", complete = TRUE, xmax = 8000, stp = 20)
