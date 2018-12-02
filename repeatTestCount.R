#******************************************************************************************
#
#   Evaluate how many samples it would take to capture all clones and what the accuracy
#   of the samples is for Z09132 day 104 
#
#   1 uL blood 4,000 - 11,000 white blood cells
#
#  Z09132 - 6kg - 366mk
#      from excel ave = 7,600  WBC / uL
#      7,600 WBC/uL * 1000 ul/ml * 366ml -> 2,781,600,000 WBC in animal
#      29.7% gene marking                ->   826,140,000 GFP+ WBC in animal
#      blood sample ~ 3.5ml              ->    26,600,000 WBC in sample
#      29.7% marking                     ->     7,900,200 GFP+ in sample
#      want 500k genomes                 ->       500,000 / 26,600,000 = 0.018 fraction of 3.5 ml
#      PCR                               ->       500,000 WBC genomes
#      29.7% marking                     ->       148,000 GFP+ genomes (aT LATEST TIME POINT)
#      one sample                        ->       500,000 / 2,781,000    = 1.8e-4
#      one sample                        ->       148,000 / 826,140,000  = 1.8e-4
#   
#   M Enstrom 5-30-2018
#
#******************************************************************************************
#
#
#source("InitRisBC.R")
#isa9132@labels
#bc9132@labels

cumulative_9132 <- function(doPlot=TRUE)
{
  bc118gfp = bc9132@data[[3]]
  bc104    = bc9132@data[[2]]
  isa104   = isa9132@data[[4]]
  #
  # make up sample population. at day 104, Z09132 at about 11% gene marking
  #
  geneModPop = 0.113 * 826140000
  geneModPop   # ~ 90,000,000
  #
  # use 118 GFP to sim
  #
  fixedList = list(round(bc118gfp[,5] * geneModPop))
  head(fixedList[[1]])
  tail(fixedList[[1]])
  cellPop = generatePopulationFromDistribution(fixedList = fixedList)
  cellPop = round(cellPop)
  head(cellPop)
  tail(cellPop)
  length(cellPop)
  
  
  cumulativeClones = c()
  lengths = c()
  
  for (i in c(1:10)) {
    #
    # build 104 blood composition entirely from 118GFP
    #
    
    #
    # sim barcode   1
    #               3.5% GM = 17500
    #
    #simBarcode_ISA <- function(cellPop, 
    #                           popSample = 500000, 
    #                           captureEfficiency = 0.17, 
    #                           pcrEff = 1.0, 
    #                           pcrCycles = 18,
    #                           sequenceDepth=4000000) 
      
    res = simBarcode_ISA(cellPop,
                         popSample = round(500000 * 0.113),
                         captureEfficiency = 0.435,
                         pcrEff = 0.85,
                         pcrCycles = 10,
                         sequenceDepth = 3173781
                         )
    res
    #res = simBarcodeID(cellPop,popSample = round(500000 * 0.113),useSimple = T, BCeff = 0.435,pcrEff = 0.85,pcrCycles = 10, sequenceDepth = 3173781)
    #
    # resR$values are clones captured
    #
    length(res$values)
    length(cellPop)
    head(res$values)
    tail(res$values)
    for (clone in res$values) {
      old = which(cumulativeClones == clone)
      if (length(old) == 0) {
        cumulativeClones = c(cumulativeClones,clone)
        #print(paste0(clone," ",which(cumulativeClones == clone)))
      }
    }
    
    head(cumulativeClones)
    tail(cumulativeClones)
    a = length(cumulativeClones)
    a
    
    lengths = c(lengths,a)
    print(lengths)
  }
  bcLengths = lengths
  
  #-------------------------------------------------------------------------------------
  #
  # ISA - use same population
  #
  #
  #-------------------------------------------------------------------------------------
  
  cumulativeClones = c()
  lengths = c()
  
  for (i in c(1:10)) {
    #res = simISA_ID(cellPop,popSample = round(500000 * 0.113), ISAeff = 0.081,pcrEff = 0.85,sequenceDepth = 64964)
    #resT = table(res)
    #resS = sort(res)
    resR = simBarcode_ISA(cellPop,
                         popSample = round(500000 * 0.113),
                         captureEfficiency = 0.07,
                         pcrEff = 0.85,
                         pcrCycles = 10,
                         sequenceDepth = 64964
    )
    #
    # resR$values are clones captured
    #
    for (clone in resR$values) {
      old = which(cumulativeClones == clone)
      if (length(old) == 0) {
        cumulativeClones = c(cumulativeClones,clone)
        #print(paste0(clone," ",which(cumulativeClones == clone)))
      }
    }
    
    head(cumulativeClones)
    tail(cumulativeClones)
    a = length(cumulativeClones)
    
    lengths = c(lengths,a)
    print(lengths)
  }
  
  # if (doPlot) {
  #   pdf(paste0(plotDir,"cumulative_9132.pdf"))
  # }
  # 
  # par(mfrow=c(1,1))
  # par(mar=c(3,3,5,2))
  # 
  # plot(bcLengths,type='l',ylim = c(0,10000),axes = F,xlab = "",ylab = "")
  # axis(1,labels = seq(1,10,1),at=(as.numeric(seq(1,10,1)))   )
  # axis(2)
  # 
  # 
  # lines(lengths,lty=2)
  # 
  # 
  # title("Z09132 DBS and ISA Cumulative Clones Captured",
  #       outer = F,
  #       line=2.0,
  #       xlab = "Number of Samples Processed",
  #       ylab="Total Clones")
  # 
  # legend(7,
  #        9000,
  #        c("DBS","ISA"),
  #        cex=1.0,
  #        col=c('black','black'),
  #        #pch=rep(19,3),
  #        lty = c(1,2),
  #        bty = 'n')
  # 
  # abline(h=0.5 * nrow(bc118gfp),lty=3)
  # 
  # if (doPlot) {
  #   dev.off()
  # }
  # 
  goal = floor(0.5 * nrow(bc118gfp))
  goal = rep(goal,10)
  
  df <- data.frame(
    samples = rep(c(1:10),3),
    clones  = c(goal,bcLengths,lengths),
    subject = c(rep('goal',10),rep('DBS',10),rep('ISA',10)),
    ls      = c(rep(1,10),rep(1,10),rep(1,10))
    )
  
  #levels(df$subject) <- c('Goal','DBS','ISA')
  levels(df$subject)
  df$subject = factor(c(rep('goal',10),rep('DBS',10),rep('ISA',10)),levels=c('goal','DBS','ISA'))
  df$subject
  head(df)  
  
  g <- ggplot(df,aes(samples,clones,group=subject,linetype=subject,size=ls))
  
  base <- g + geom_line(aes(color=subject)) +
    labs(x = "Number of Samples Processed",y = "Total Clones", title = "Z09132") +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title = element_text(size=14),
          plot.title = element_text(size=16),
          legend.text = element_text(size=14),
          legend.title=element_blank(),
          legend.position =  c(0.2,0.85)
    ) +
    scale_size(range = c(1, 2), guide="none") +
    scale_color_manual(values=c('#e0e0e0',"#b0b0b0","#404040")) +
    scale_linetype_manual(values = c('dotted','solid','solid')) +
    guides(linetype = guide_legend(override.aes = list(size = 2))) +
    scale_x_continuous(breaks = c(1:10),labels = c("1","2","3","4","5","6","7","8","9","10"))

  base
    
  pdf(paste0(plotDir,"cumulative_9132.pdf"))
  base
  dev.off()
    

}
dev.off()
cumulative_9132()



#-------------------------------------------------------------------------------------
bc8103@labels
cumulative_8103 <- function(doPlot = TRUE)
{
  #
  # from highest one time clone detection
  #
  bc644    = bc8103@data[[30]]
  #
  # make up sample population. at day 104, Z09132 at about 11% gene marking
  #
  geneModPop = 0.0174 * 826140000
  geneModPop   # ~ 90,000,000
  #
  # use 118 GFP to sim
  #
  fixedList = list(round(bc644[,5] * geneModPop))
  head(fixedList[[1]])
  tail(fixedList[[1]])
  cellPop = generatePopulationFromDistribution(fixedList = fixedList)
  cellPop = round(cellPop)
  head(cellPop)
  tail(cellPop)
  length(cellPop)
  
  
  cumulativeClones = c()
  lengths = c()
  
  for (i in c(1:20)) {
    #
    # build 104 blood composition entirely from 118GFP
    #
    
    #
    # sim barcode   1
    #               3.5% GM = 17500
    #
    #simBarcode_ISA <- function(cellPop, 
    #                           popSample = 500000, 
    #                           captureEfficiency = 0.17, 
    #                           pcrEff = 1.0, 
    #                           pcrCycles = 18,
    #                           sequenceDepth=4000000) 
    
    res = simBarcode_ISA(cellPop,
                         popSample = round(500000 * 0.0174),
                         captureEfficiency = 0.435,
                         pcrEff = 0.85,
                         pcrCycles = 10,
                         sequenceDepth = 3173781
    )
    res
    #res = simBarcodeID(cellPop,popSample = round(500000 * 0.113),useSimple = T, BCeff = 0.435,pcrEff = 0.85,pcrCycles = 10, sequenceDepth = 3173781)
    #
    # resR$values are clones captured
    #
    length(res$values)
    length(cellPop)
    head(res$values)
    tail(res$values)
    for (clone in res$values) {
      old = which(cumulativeClones == clone)
      if (length(old) == 0) {
        cumulativeClones = c(cumulativeClones,clone)
        #print(paste0(clone," ",which(cumulativeClones == clone)))
      }
    }
    
    head(cumulativeClones)
    tail(cumulativeClones)
    a = length(cumulativeClones)
    a
    
    lengths = c(lengths,a)
    print(lengths)
  }
  bcLengths = lengths
  
  #-------------------------------------------------------------------------------------
  #
  # ISA - use same population
  #
  #
  #-------------------------------------------------------------------------------------
  
  cumulativeClones = c()
  lengths = c()
  
  for (i in c(1:20)) {
    #res = simISA_ID(cellPop,popSample = round(500000 * 0.113), ISAeff = 0.081,pcrEff = 0.85,sequenceDepth = 64964)
    #resT = table(res)
    #resS = sort(res)
    resR = simBarcode_ISA(cellPop,
                          popSample = round(500000 * 0.0174),
                          captureEfficiency = 0.07,
                          pcrEff = 0.85,
                          pcrCycles = 10,
                          sequenceDepth = 64964
    )
    #
    # resR$values are clones captured
    #
    for (clone in resR$values) {
      old = which(cumulativeClones == clone)
      if (length(old) == 0) {
        cumulativeClones = c(cumulativeClones,clone)
        #print(paste0(clone," ",which(cumulativeClones == clone)))
      }
    }
    
    head(cumulativeClones)
    tail(cumulativeClones)
    a = length(cumulativeClones)
    
    lengths = c(lengths,a)
    print(lengths)
  }
  
  goal = 0.5*nrow(bc644)
  goal = rep(goal,20)
  
  df <- data.frame(
    samples = rep(c(1:20),3),
    clones  = c(goal,bcLengths,lengths),
    subject = c(rep('goal',10),rep('DBS',10),rep('ISA',10)),
    ls      = c(rep(1,20),rep(1,20),rep(1,20))
  )
  
  levels(df$subject)
  df$subject = factor(c(rep('goal',20),rep('DBS',20),rep('ISA',20)),levels=c('goal','DBS','ISA'))
  df$subject
  head(df)  
  
  g <- ggplot(df,aes(samples,clones,group=subject,linetype=subject,size=ls))
  
  base <- g + geom_line(aes(color=subject)) +
    labs(x = "Number of Samples Processed",y = "Total Clones", title = "Z08103") +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title = element_text(size=14),
          plot.title = element_text(size=16),
          legend.text = element_text(size=14),
          legend.title=element_blank(),
          legend.position =  c(0.2,0.90)
    ) +
    scale_size(range = c(1, 2), guide="none") +
    scale_color_manual(values=c('#e0e0e0',"#b0b0b0","#404040")) +
    scale_linetype_manual(values = c('dotted','solid','solid')) +
    guides(linetype = guide_legend(override.aes = list(size = 2))) +
    scale_x_continuous(breaks = c(1:20),labels = c("1","2","3","4","5","6","7","8","9","10",
                                                   '11','12','13','14','15','16','17','18','19','20'))
  
  base
  
  pdf(paste0(plotDir,"cumulative_8103.pdf"))
  base
  dev.off()
  #
  # old plot code
  #
  # 
  # if (doPlot) {
  #   pdf(paste0(plotDir,"Cumulative_8103.pdf"))
  # }
  # 
  # par(mfrow=c(1,1))
  # par(mar=c(3,3,5,2))
  # 
  # plot(bcLengths,type='l',ylim = c(0,20000),axes = F,xlab = "",ylab = "")
  # axis(1,labels = seq(1,20,2),at=(as.numeric(seq(1,20,2)))   )
  # axis(2)
  # 
  # 
  # lines(lengths,lty=2)
  # 
  # 
  # title("Z08103 DBS and ISA Cumulative Clones Captured",
  #       outer = F,
  #       line=2.0,
  #       xlab = "Number of Samples Processed",
  #       ylab="Total Clones")
  # 
  # legend(14,
  #        20000,
  #        c("DBS","ISA"),
  #        cex=1.0,
  #        col=c('black','black'),
  #        #pch=rep(19,3),
  #        lty = c(1,2),
  #        bty = 'n')
  # 
  # abline(h=0.5*nrow(bc644),lty=3)
  # 
  # if (doPlot) {
  #   dev.off()
  # }
}

cumulative_9132()
cumulative_8103()
bc8103@labels
