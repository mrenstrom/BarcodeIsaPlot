#******************************************************************************************
#
#   Evaluate how many samples it would take to capture all clones and what the accuracy
#   of the samples is for Z08103
#
#  Z08103 - 1555 kg - 975 ml blood
#       5,7600,000 infused, 53% GFP+      ->     3,050,000 maximum possible clones
#       from excel: ~5000 WBC/ul
#       5000 WBC/ul * 1000 ul/ml * 975 ml -> 4,875,000,000  WBC in animal
#       3.5% gene marking                 ->   170,630,000  GFP+ WBC in animal
#       blood sample ~ 3.5 ml             ->    17,500,000  WBC   in whole sample
#       3.5% gene marking                 ->       612,500  GFP+  in whole sample
#       want 500k geneomes                ->       500,000 / 17,500,000 -> 0.028 fraction of 3.5 ml
#       PCR                               ->       500,000 genomes
#       3.5% gene marking                 ->        17,500  GFP+ in PCR sample
#       one pcr sample                    ->   500,000 / 4,875,000,000 = 1.0e-4 of animal
#
#   
#   M Enstrom 5-30-2018
#
#******************************************************************************************
#
#
#source("InitRisBC.R")
#isa8103@labels
#bc8103@labels
#isa8103@labels



#******************************************************************************************
#
# Setup and simulation for barcode and ISA day 160 and 145GPF
#
#   
#   M Enstrom 5-30-2018
#
#******************************************************************************************
findVals = function(cellBase,p1,m1,d1,p2,m2,d2,plotName="") {
  # debug
  cellBase = bc8103@data[[30]]
  #f=2500
  p1 = 1
  m1 = 1
  d1 = 1
  p2 = 1
  m2 = 1
  d2 = 1
  f = nrow(cellBase)
  plotName = "8103_145_160_sim.pdf"
 
  bc145gfp = bc8103@data[[3]]
  bc160    = bc8103@data[[4]]
  bc223    = bc8103@data[[6]]
  isa160   = isa8103@data[[4]]
  isa223   = isa8103@data[[5]]
  geneModPop = 0.0174 * 4875000000
  par(mfrow=c(1,1))
  par(mar=c(2,2,5,2))
  #
  # set up cells from an actual sample
  #
  fixedList = list(cellBase[1:f,5] * geneModPop)
  #
  # make up rest of random population
  #
  attrList = list(list(p1,m1,d1),list(p2,m2,d2))
  cellPop = generatePopulationFromDistribution(fixedList = fixedList,attrList = attrList)
  #
  # barcode sim
  #
  simRle = simBarcode_ISA(cellPop,
                       popSample = round(500000 * .0174),
                       captureEfficiency = .44,
                       pcrEff = 0.5,
                       pcrCycles = 12, 
                       sequenceDepth = 4344719)
  
  sim160 = sort(simRle$lengths,decreasing = T)
  #
  # set up plot
  #
  #res = sort(res, decreasing = FALSE)
  #resS = sort(res)
  #resR = rle(resS)
  #resL = sort(resR$lengths, decreasing = TRUE)
  #print("tail resL  =")
  #print(tail(resL))
  #sim160 = resL
  #print("nrow(bc160) = ")
  #print(nrow(bc160))
  #print("nrow(sim160) = ")
  #print(length(resR$lengths))
  #
  # 145gfp based on day 160
  #
  simRle = simBarcode_ISA(cellPop,
                       popSample = round(225000),
                       captureEfficiency = 0.44,
                       pcrEff = 0.5,
                       pcrCycles = 12, 
                       sequenceDepth = 3697580)
  sim145gfp = sort(simRle$lengths,decreasing = TRUE)
  #
  # set up plot
  #
  #res145 = sort(res145, decreasing = FALSE)
  #resS_145 = sort(res145)
  #resR_145 = rle(resS_145)
  #resL_145 = sort(resR_145$lengths, decreasing = TRUE)
  #print("tail resL  =")
  #print(tail(resL_145))
  #sim145gfp = resL_145
  #print("nrow(bc145gfp) = ")
  #print(nrow(bc145gfp))
  #print("nrow(sim145) = ")
  #print(length(resR_145$lengths))
  #
  # ISA sim
  #
  simRle = simBarcode_ISA(cellPop,
                       popSample = round(500000 * .0174),
                       captureEfficiency = 0.04,
                       pcrEff = 0.5,
                       pcrCycles = 20, 
                       sequenceDepth = 227363)
  sim160_ISA = sort(simRle$lengths,decreasing = TRUE)
  #
  #
  #
  #resISA = sort(resISA, decreasing = FALSE)
  #resS_ISA = sort(resISA)
  #resR_ISA = rle(resS_ISA)
  #resL_ISA = sort(resR_ISA$lengths, decreasing = TRUE)
  #print("tail resL_ISA  =")
  #print(tail(resL_ISA))
  #sim160_ISA = resL_ISA
  
  df1 <- data.frame(x = c(1:nrow(bc145gfp)),
                    frac=log(bc145gfp[,3]),
                    lab=rep("DBS GFP",times=nrow(bc145gfp)))
  
  
  df2 <- data.frame(x = c(1:nrow(bc160)),
                    frac=log(bc160[,3]),
                    lab=rep("DBS",times=nrow(bc160)))

  df3 <- data.frame(x = c(1:length(sim160)),
                    frac=log(sim160),
                    lab=rep("Sim DBS",times=length(sim160)))
  

  df4 <- data.frame(x = c(1:nrow(isa160)),
                    frac=log(isa160[,3]),
                    lab=rep("ISA",times=nrow(isa160)))

  df5 <- data.frame(x = c(1:length(sim160_ISA)),
                    frac=log(sim160_ISA),
                    lab=rep("Sim ISA",times=length(sim160_ISA)))


  #
  # not used
  #
  # df6 <- data.frame(x = c(1:length(sim145gfp)),
  #                   frac=log(sim145gfp),
  #                   lab=rep("Sim DBS GFP",times=length(sim145gfp)))
  # 
  # 
  # df7 <- data.frame(x = c(1:nrow(bc223)),
  #                   frac=log(bc223[,3]),
  #                   lab=rep("DBS 223",times=nrow(bc223)))
  
  
  dfG <- rbind(df1,df2)
  dfG <- rbind(dfG,df3)
  dfG <- rbind(dfG,df4)
  dfG <- rbind(dfG,df5)
  #dfG <- rbind(dfG,df6)
  #dfG <- rbind(dfG,df7)
  
  if (plotName != "") {
    pdf(paste0(plotDir,plotName))
  }
  
  g <- ggplot(dfG,aes(x=x,y=frac,color=lab,linetype=lab))
  g + geom_line(size=1) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title = element_text(size=14),
          plot.title = element_text(size=16),
          legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.position = c(0.8,0.75)
    ) +
    scale_color_manual(values=c("#686868", "#505050", "#909090","#343434","#d0d0d0")) + 
    scale_linetype_manual(values=c("dotted","dotdash","solid","dashed","solid")) +
      
    #scale_color_manual(values=c("#686808","#500050","#009090","#343434","#f00000","#00ff00","#0000ff","#ffff00")) + 
    #scale_linetype_manual(values=c("dotted","dotdash","solid","dashed","solid","dotted","dashed","solid")) +
    labs(x="Clone Rank",y="Log of Normalized Sequence Reads",title="Simulation of Z08103 DBS and ISA")# + coord_fixed(ratio=300)
  
  
  if (plotName != "") {
    dev.off()
  }
  return(list(cellPop,sim160,sim160_ISA))
}
#
# run the tests
#
#
# match from 644 - all conditions
#
l = findVals(cellBase = bc8103@data[[30]],
             1,1,1,
             1,1,1,
             "8103_145_160_sim.pdf")

#----------------------------------------------------------------------------------------------------
#
# build a sim population for Z08103 and sample this population to compare barcode and
# ISA capture
#
#----------------------------------------------------------------------------------------------------
bc8103@labels
sim_cap_80103 <- function(doPlot = FALSE) {
  bc160 = bc8103@data[[4]]
  nrow(bc160)
  bc145gfp = bc8103@data[[3]]
  nrow(bc145gfp)
  isa160 = isa8103@data[[4]]
  nrow(isa160)
  #
  # from combined 644 dpt
  #
  popFrac = bc8103@data[[30]][,5]
  head(popFrac)        # 0.008  0.008  0.006 ...
  sum(popFrac)         # 1.0
  l = length(popFrac) 
  l
  l160 = nrow(bc160)
  l160
  l160_isa = nrow(isa160)
  l160_isa
  #
  # simulate sample of animal to match 104dpt : 500,000 genomes at 29.7% gene marking = 148,000 clones
  #
  gmCells = round(500000 * 0.0174)
  gmCells #8700
  
  
  par(mfrow=c(1,1))
  par(mar=c(2,2,5,2))
  
  dfG = data.frame(x=c(1:10),y=rep(l160,10),lab="160")
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
    while (y[j] > l160) {
      j = j + 1
    }
    dx = 1
    dy = y[j-1] - y[j]
    xInt = j-1 + (y[j-1]-l160) * dx / dy
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
    pdf(paste0(plotDir,"Z08103_BarcodeCaptureEff.pdf"))
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
    labs(y="Clones Captured",title="Simulation of Z08103 DBS Capture Efficiency") +
    annotate("text",x=8, y = 4000,label=effMean,size=6)
  if (doPlot) {
    dev.off()
  }
  #
  # Now ISA
  #
  dfG = data.frame(x=c(1:10),y=rep(l160_isa,10),lab="160")
  eff = c()
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
    while (y[j] > l160_isa) {
      j = j + 1
    }
    dx = 1
    dy = y[j-1] - y[j]
    xInt = j-1 + (y[j-1]-l160_isa) * dx / dy
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
    pdf(paste0(plotDir,"Z08103_ISACaptureEff.pdf"))
  }
  
  g <- ggplot(dfG,aes(x=x,y=y,color=lab,linetype=lab))
  g + geom_line(size=0.8) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title = element_text(size=14),
          plot.title = element_text(size=16),
          legend.title = element_text(size=32)
    ) +
    scale_x_continuous(name="Capture Efficiency",breaks = c(1:10),labels = c("0.2","0.18","0.16","0.14","0.12","0.10","0.08","0.06","0.04","0.02")) +
    guides(color=FALSE) +
    scale_color_manual(values=rep(c("#000000"),21),guide=FALSE) + 
    scale_linetype_manual(values=c("dotted",rep(c("solid"),20)),guide=FALSE) +
    labs(y="Clones Captured",title="Simulation of Z08103 ISA Capture Efficiency") + 
    annotate("text",x=8, y = 1100,label=effMean,size=6)
  
  if (plotName != "") {
    dev.off()
  }
}

sim_cap_80103(doPlot = TRUE)


#
# from bc9132@118 gfp
#
#l = findVals(cellBase = bc9132@data[[3]],
#             11000,
#             1,1,1,
#             1,5,1.5)
# 
# 
# cellPop = l[[1]]
# sim160 = l[[2]]
# sim160_ISA = l[[3]]
# l502 = nrow(bc8103@data[[20]])
# l160 = nrow(bc8103@data[[4]])
# 
# cellFreq = cellPop/sum(cellPop)
# bloodSample = sample(length(cellPop), 500000 * 0.0174, replace = T, prob = cellFreq)
# 
# ld = lapply(c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1), function(x) {
#   capture  = sample(bloodSample,length(bloodSample) * x,replace = F)
#   #bcSample = sample(capture,2498599,replace = T)
#   #t = table(bcSample)
#   t = table(capture)
#   return(length(t))
# })
# unlist(ld)
# 
# 
# pdf(paste0(plotDir,"8103_bc_eff.pdf"))
# 
# plot(1:10,type='n',ylim=c(1000,5000),axes=F)
# lines(unlist(ld))
# #points(9,3205)
# abline(h=l502,col='pink')
# abline(h=l160,col='blue')
# 
# axis(1,labels = seq(1.0,0.1,-0.1),at=(as.numeric(seq(1,10,1))))
# axis(2)
# dev.off()
# 
# #
# # ISA
# #
# l160_isa = nrow(isa8103@data[[4]])
# ld = lapply(c(0.2,0.18,0.16,0.14,0.12,0.1,0.08,0.06,0.04,0.02), function(x) {
#   capture  = sample(bloodSample,length(bloodSample) * x,replace = F)
#   t = table(capture)
#   return(length(t))
# })
# unlist(ld)
# pdf(paste0(plotDir,"8103_isa_eff.pdf"))
# 
# plot(1:10,type='n',ylim=c(0,1500),axes=F)
# lines(unlist(ld))
# abline(h=l160_isa,col='blue')
# 
# axis(1,labels = seq(0.2,0.02,-0.02),at=(as.numeric(seq(1,10,1))))
# axis(2)
# dev.off()

