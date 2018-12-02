#******************************************************************************************
#
#   Look at Z09132 118 DPT GFP barcode and compare counts to Z09132 104 DPT barcode and RIS
#
#   For a first approximation of barcode efficiency use 118 GFP as the model for actual
#   blood composition of clones...since sampling many more clones than 104_dpt. Sample the 
#   118dpt_gfp population then apply a barcode efficiency estimate and solve. Do the same
#   for RIS at 104dpt.
#   
#   M Enstrom 4-11-2018
#
#******************************************************************************************
isa9132@labels
bc9132@labels

#source("InitRisBC.R")
library(reshape2)

#******************************************************************************************
#
#   Comare DBS results for Z09132 104 and 118 gfp
#   
#   M Enstrom 4-11-2018
#
#******************************************************************************************
# 
# bc9103_104_118g <- function(plotName="") {
#   #
#   # load Z09132 104dpt and 118 dpt gfp for comparison
#   #
#   d104   = bc9132@data[[2]]
#   d118   = bc9132@data[[3]]
#   isa104 = isa9132@data[[4]]  
#   #
#   # to show that the blood composition is very similar from day 104 to day 118, show plot of 
#   # both days together
#   #
#   # plot row data for Z9132 104 v 118 gfp
#   #
#   if (plotName != "") {
#     pdf(paste0(plotDir,"Z09132_104_118gfp.pdf"))
#   }
# 
#   df1 <- data.frame(x = c(1:nrow(d118)),
#                     frac=d118[,6],
#                     lab=rep("DBS GFP Sort",times=nrow(d118)))
#   head(df1,1)
#   
#   df2 <- data.frame(x = c(1:nrow(d104)),
#                     frac=d104[,6],
#                     lab=rep("DBS",times=nrow(d104)))
#   
#   head(df2,1)
#   
#   dfG <- rbind(df1,df2)
# 
#   g <- ggplot(dfG,aes(x=x,y=frac,color=lab))
#   g + geom_line()
#   
#   if (plotName != "") {
#     dev.off()
#   }
#   
#   par(mfrow=c(2,1))
#   par(mar=c(2,5,1,1))
#   
#   plot(d118[1:1000,5],type='l',col='green',axes=F,ylim=c(0,0.010))
#   lines(d104[1:1000,5],col='black')
#   axis(1,labels = seq(0,1000,200),at=(as.numeric(seq(0,1000,200)))   )
#   axis(2)
#   
#   legend(600,
#          0.006,
#          c("barcode 118dpt GFP","barcode 104dpt"),
#          #cex=0.2,
#          col=c('green','black'),
#          #pch=c(20,2),
#          lty = c(1,1),
#          bty = 'n')
#   
#   
#   pr = c(1000:5000)
#   plot(d118[pr,5],type='l',col='green',axes = F,ylim=c(0,0.0001))
#   lines(d104[pr,5],col='black')
#   
#   axis(1,labels = seq(1000,5000,500),at=(as.numeric(seq(1000,5000,500))-1000)   )
#   axis(2)
#   
#   
#   
#   title("Z09132 Barcode 118_DPT_GFP -v- 104_DPT",outer = T,line=-2.0)
# 
#   
# }
#----------------------------------------------------------------------------------------------------
#
# build a simulation of 104 DPT using clone contribution of 118DPT GFP. Fin d the barcdeo capture 
# efficiency based on total number of clones captured
#
#
#----------------------------------------------------------------------------------------------------
sim_cap_9132_104 <- function(doPlot = FALSE) {
  doPlot=TRUE
  d104 = bc9132@data[[2]]
  nrow(d104)
  d118 = bc9132@data[[3]]
  nrow(d118)
  popFrac = d118[,5]
  isa104 = isa9132@data[[4]]
  head(popFrac)        # 0.008  0.008  0.006 ...
  sum(popFrac)         # 1.0
  l = length(popFrac) 
  l
  l104 = nrow(d104)
  l104
  l104_isa = nrow(isa104)
  l104_isa
  #
  # simulate sample of animal to match 104dpt : 500,000 genomes at 29.7% gene marking = 148,000 clones
  #
  gmCells = round(500000 * 0.113)
  gmCells #56,500
  

  par(mfrow=c(1,1))
  par(mar=c(2,2,5,2))
  #
  # horizontal line for 104DPT total count
  #
  dfG = data.frame(x=c(1:10),y=rep(l104,10),lab="104")
  #
  # sim different capture rates
  #
  eff=c()
  for (i in c(1:20)) {
  
    bloodSample = sample(length(popFrac),gmCells,replace = T, prob = popFrac)
    #
    # check
    #
    bloodTable = table(bloodSample)
    length(bloodTable) 
    #
    # actual value of 104DPT is 3205 so need to simulate less than 100% efficiency of capturing barcodes
    #
  
    ld = lapply(c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1), function(x) {
      bcSample = sample(bloodSample,length(bloodSample) * x,replace = F)
      t = table(bcSample)
      return(length(t))
    })
    unlist(ld)
    
    df = data.frame(x=c(1:10),y=unlist(ld),lab=paste0(i))
  

    dfG = rbind(dfG,df)
    
    #
    # estimate intercept
    #
    
    y = unlist(ld)
    j = 1
    while (y[j] > l104) {
      j = j + 1
    }
    dx = 1
    dy = y[j-1] - y[j]
    xInt = j-1 + (y[j-1]-l104) * dx / dy
    #
    # convert axis [1-10] -> [1.0 - 0.1]
    #
    xInt = 1.1 - (xInt/10)
    eff = c(eff,xInt)
    message(paste0(j," ",dx," ",dy," ",xInt))
       

    #lines(unlist(ld))
    #points(9,3205)
    #abline(h=3149,lty=2)
    #abline(h=3250,col='blue')
  }
  effMean = mean(eff)
  effMean = sprintf("mean = %0.2f",effMean)
  #axis(1,labels = seq(1.0,0.1,-0.1),at=(as.numeric(seq(1,10,1))))
  #axis(2)
  #title("Total Clones Identified -v- Capture Efficiency",outer = T,line=-2.0)
  if (doPlot) {
    pdf(paste0(plotDir,"Z09132_BarcodeCaptureEff.pdf"))
  }
  
  g <- ggplot(dfG,aes(x=x,y=y,color=lab,linetype=lab))
  g + geom_line(size=0.8) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
        axis.title = element_text(size=14),
        plot.title = element_text(size=16),
        legend.title = element_text(size=32)
        ) +
    
    scale_x_continuous(name="Capture Efficiency",breaks = c(1:10),labels = c("1.0","0.9","0.8","0.7","0.6","0.5","0.4","0.3","0.2","0.1")) +
    guides(color=FALSE) +
    scale_color_manual(values=rep(c("#000000"),21),guide=FALSE) + 
    scale_linetype_manual(values=c("dotted",rep(c("solid"),20)),guide=FALSE) +
    labs(y="Clones Captured",title="Z09132 Simulation of DBS Capture Efficiency") +
    annotate("text",x=8, y = 4000,label=effMean,size=6)
  if (doPlot) {
    dev.off()
  }
  #
  # Now ISA...horizontal line for 104DPT clone count
  #
  dfG = data.frame(x=c(1:10),y=rep(l104_isa,10),lab="104")
  eff = c()
  #
  # sim capture rates
  #
  for (i in c(1:20)) {
    
    bloodSample = sample(length(popFrac),gmCells,replace = T, prob = popFrac)
    #
    # check
    #
    bloodTable = table(bloodSample)
    length(bloodTable) 
    #
    # actual value of 104DPT is 3205 so need to simulate less than 100% efficiency of capturing barcodes
    #
    
    ld = lapply(c(0.2,0.18,0.16,0.14,0.12,0.1,0.08,0.06,0.04,0.02), function(x) {
      bcSample = sample(bloodSample,length(bloodSample) * x,replace = F)
      t = table(bcSample)
      return(length(t))
    })
    unlist(ld)
    
    df = data.frame(x=c(1:10),y=unlist(ld),lab=paste0(i))
    
    
    dfG = rbind(dfG,df)
    #
    # estimate intercept
    #
    y = unlist(ld)
    j = 1
    while (y[j] > l104_isa) {
      j = j + 1
    }
    dx = 1
    dy = y[j-1] - y[j]
    xInt = j-1 + (y[j-1]-l104_isa) * dx / dy
    #
    # convert axis
    #
    xInt = 0.22 - (xInt/50)
    eff = c(eff,xInt)
    message(paste0(j," ",dx," ",dy," ",xInt))
    
    #lines(unlist(ld))
    #points(9,3205)
    #abline(h=3149,lty=2)
    #abline(h=3250,col='blue')
  }
  effMean = mean(eff)
  effMean = sprintf("mean = %0.2f",effMean)
  #axis(1,labels = seq(1.0,0.1,-0.1),at=(as.numeric(seq(1,10,1))))
  #axis(2)
  #title("Total Clones Identified -v- Capture Efficiency",outer = T,line=-2.0)
  if (doPlot) {
    pdf(paste0(plotDir,"Z09132_ISACaptureEff.pdf"))
  }
  
  g <- ggplot(dfG,aes(x=x,y=y,color=lab,linetype=lab))
  g + geom_line(size=0.8) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title = element_text(size=14),
          plot.title = element_text(size=16),
          legend.title = element_text(size=16)
    ) +
    
    scale_x_continuous(name="Capture Efficiency",breaks = c(1:10),labels = c("0.2","0.18","0.16","0.14","0.12","0.10","0.08","0.06","0.04","0.02")) +
    guides(color=FALSE) +
    scale_color_manual(values=rep(c("#000000"),21),guide=FALSE) + 
    scale_linetype_manual(values=c("dotted",rep(c("solid"),20)),guide=FALSE) +
    labs(y="Clones Captured",title="Z09132 Simulation ISA Capture Efficiency") + 
    annotate("text",x=8, y = 2000,label=effMean,size=6)
  
  if (plotName != "") {
    dev.off()
  }
  
}
sim_cap_9132_104(doPlot = TRUE)

#******************************************************************************************
#
#   Look at Z09132 day 104 (blood composition from 118)
#
#   1 uL blood 4,000 - 11,000 white blood cells
#
#  Z09132 - 6kg - 366mk
#      from excel ave = 7,600  WBC / uL
#      7,600 WBC/uL * 1000 ul/ml * 366ml -> 2,781,600,000 WBC in animal
#      11.3% gene marking                ->   314,320,800 GFP+ WBC in animal (11.3 on day 104)
#      blood sample ~ 3.5ml              ->    26,600,000 WBC in sample
#      11.3% marking                     ->     3,005,800 GFP+ in sample
#      want 500k genomes                 ->       500,000 / 26,600,000 = 0.018 fraction of 3.5 ml
#      PCR                               ->       500,000 WBC genomes
#      11.3% marking                     ->       148,000 GFP+ genomes
#      one sample                        ->       500,000 / 2,781,000    = 1.8e-4
#      one sample                        ->        56,500 / 826,140,000  = 1.8e-4
#   
#   M Enstrom 4-11-2018
#
#******************************************************************************************
#***********************************************************************************************
# add in ris
#
#    Z09132 104 dpt new data
#
#
#
#***********************************************************************************************
#
# sample RIS efficiency as a fixed probability of shearing dna at proper length + binding
# of linker. Use same blood sample of 148000 cells from barcode
#
#
sim_verify_cap_ris_9132_104 <- function(plotName="") {
  d104 = bc9132@data[[2]]
  d118 = bc9132@data[[3]]
  isa104 = isa9132@data[[4]]
  #
  # blood composition for simulatio from 118gfp
  #
  popFrac = d118[,5]
  head(popFrac)        # 0.008  0.008  0.006 ...
  sum(popFrac)         # 1.0
  l = length(popFrac)  # 11272
  l104 = nrow(d104)    # 3205
  #
  #  build sample population
  #
  cellPop = sample(l,500000,replace = T,prob = popFrac)
  r = rle(sort(cellPop))
  r
  #
  # run DBS and ISA sim
  #
  simRle = simBarcode_ISA(r$lengths,
                          popSample = round(500000 * .113),
                          captureEfficiency = .41,
                          pcrEff = 0.5,
                          pcrCycles = 12, 
                          sequenceDepth = 9999791)
  
  dbs_sim_104 = sort(simRle$lengths, decreasing = TRUE)
  dbs_sim_104 = log(dbs_sim_104)
  #
  # ISA sim
  #
  simRle = simBarcode_ISA(r$lengths,
                               popSample = round(500000 * .113),
                               captureEfficiency = .07,
                               pcrEff = 0.5,
                               pcrCycles = 12, 
                               sequenceDepth = 99000)
  isa_sim_104 = sort(simRle$lengths, decreasing = TRUE)
  isa_sim_104 = log(isa_sim_104)
  #
  # plot results  
  #
  plotName = "Z09132_ISA_Sim.pdf"
  if (plotName != "") {
    name = paste0(plotDir,plotName)
    message(name)
    pdf(name)
  }
  
  par(mfrow=c(1,1))
  par(mar=c(2,5,1,1))
  
  
  df0 <- data.frame(x = c(1:nrow(d118)),
                    frac=log(d118[,4]),
                    ID=rep("DBS GFP",times=nrow(d118)))
  
  df1 <- data.frame(x = c(1:nrow(d104)),
                    frac=log(d104[,4]),
                    ID=rep("DBS",times=nrow(d104)))
  
  df2 <- data.frame(x = c(1:length(dbs_sim_104)),
                    frac=dbs_sim_104,
                    ID=rep("Sim DBS",times=length(dbs_sim_104)))

  df3 <- data.frame(x = c(1:nrow(isa104)),
                    frac=log(isa104[,4]),
                    ID=rep("ISA",times=nrow(isa104)))
  
  df4 <- data.frame(x = c(1:length(isa_sim_104)),
                    frac=isa_sim_104,
                    ID=rep("Sim ISA",times=length(isa_sim_104)))
  
  
  dfG <- rbind(df0,df1)
  dfG <- rbind(dfG,df2)
  dfG <- rbind(dfG,df3)
  dfG <- rbind(dfG,df4)
  
  
  g <- ggplot(dfG,aes(x=x,y=frac,color=ID,linetype=ID))
  g + geom_line(size=1)  + 
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title = element_text(size=14),
          plot.title = element_text(size=16),
          legend.title = element_blank(),
          legend.text = element_text(size=16)
    ) +
    theme(legend.position = c(0.8,0.75)) +
    scale_color_manual(values=c("#686868", "#505050", "#909090","#343434","#d0d0d0")) + 
    scale_linetype_manual(values=c("dotted","dotdash","solid","dashed","solid")) +
    labs(x="Clone Rank",y="Log of Normalized Sequence Reads",title="Simulation of DBS and ISA")# + coord_fixed(ratio=300)
  
  
  if (plotName != "") {
    dev.off()
  }
}

sim_verify_cap_ris_9132_104("Z09132_ISA_Sim.pdf")
