#
#   Simulation functions
#     common ISA and barcode simulation functions
# 
#   Mark Enstrom 5-10-18
#
#
#



#******************************************************************************************
#
# Common sim function - given a clone population simulate blood sample, ISA/Barcode
# capture and sequencing
#
# cellPop = array of clone clounts  1000, 800, 400 ...
#
#
#  returns clone rle
#      r$values  = clone ID
#      r$lengths = clone counts
#
#   M Enstrom 5-30-2018
#
#******************************************************************************************
simBarcode_ISA <- function(cellPop, 
                           popSample = 500000, 
                           captureEfficiency = 0.17, 
                           pcrEff = 1.0, 
                           pcrCycles = 18,
                           sequenceDepth=4000000) 
{
  # debug
  #popSample = round(500000 * .113)
  #BCeff = 0.435
  #pcrEff = 0.85
  #pcrCycles = 14
  #sequenceDepth = 4344719
  #
  # draw blood sample  (popSample, 500,000 if 100% gene modified)
  #
  bloodSample = sample(length(cellPop),popSample,replace = TRUE, prob = cellPop/sum(cellPop))
  bt = table(bloodSample)
  #print(head(bt))
  #print(tail(bt))
  #print(paste0("l(bt) = ",length(bt)))
  #
  #
  #
  r = rle(sort(bloodSample))
  #print(r)
  length(cellPop)
  #
  # barcode or isa capture
  #
  capture = sample(bloodSample,length(bloodSample) * captureEfficiency,replace = FALSE)
  ct = table(capture)
  #print(head(ct))
  #print(tail(ct))
  #print(paste0("l(ct) = ",length(ct)))
  #
  # pcr
  #
  if (pcrCycles > 1) {
    for (i in c(1:pcrCycles)) {
      cp2 = sample(capture,length(capture)*pcrEff,replace = FALSE)
      capture = c(capture,cp2)
      #print(length(capture))
    }
    #print(length(capture))
  }
  #
  # seq to deptth, equal prob every pick
  #
  seq = sample(capture,sequenceDepth,replace = T)
  seq = sort(seq)
  st = rle(seq)
  #print(st)
  return(st)
}
#---------------------------------------------------------------------------------------------------------
#
# simBarcode
#
#
#
#---------------------------------------------------------------------------------------------------------
simBarcode <- function(cellPop, 
                       popSample = 500000, 
                       useSimple = T, 
                       BCeff = 0.17, 
                       pcrEff = 1.0, 
                       pcrCycles = 18,
                       sequenceDepth=4000000
) {
  #message(paste0("simBarcode ",length(cellPop)," bc eff = ",BCeff," pcreff = ",pcrEff," cycles ",pcrCycles," seq depth = ",sequenceDepth))
  cellFreq = cellPop/sum(cellPop)
  
  bloodSample = sample(length(cellFreq),popSample,replace = T, prob = cellFreq)
  sampleCount = table(bloodSample)
  sampleCount = as.numeric(sampleCount)
  sampleCount = sampleCount[order(-sampleCount)]
  
  #message(length(sampleCount))
  #print(head(sampleCount))
  #print(tail(sampleCount))
  #
  # a few rounds of amplification
  #
  pcrCap = c()
  
  if (useSimple) {
    #
    # 
    #
    odds = runif(length(bloodSample))
    oi = which(odds <= BCeff)
    bnew = bloodSample[oi]
    pcrCap = c(bnew)
    
    for (i in c(1:pcrCycles)) {
      odds = runif(length(pcrCap))
      oi = which(odds < pcrEff) 
      pcrnew = pcrCap[oi]
      pcrCap = c(pcrCap,pcrnew)
    }
  } else {
    for (i in c(1:pcrCycles)) {
      odds = runif(length(bloodSample))
      oi = which(odds < BCeff) #0.04
      bnew = bloodSample[oi]
      odds = runif(length(pcrCap))
      oi = which(odds < pcrEff) #80
      pcrnew = pcrCap[oi]
      pcrCap = c(pcrCap,bnew,pcrnew)
    }
  }
  bloodCap = pcrCap
  #
  # sequence
  #
  seq = sample(bloodCap,sequenceDepth,replace=F)
  
  simResult = table(seq)
  simResult = as.numeric(simResult)
  simResult = simResult[order(-simResult)]
  simFreq   = simResult / sum(simResult)
  return(list(simResult,simFreq,sampleCount))
}


#---------------------------------------------------------------------------------------------------------
#
# simISA
#
#
#
#---------------------------------------------------------------------------------------------------------
simISA <- function(cellPop, popSample = 500000, ISAeff = 0.17, pcrEff = 1.0, sequenceDepth=200000) {
  #message(paste0("simISA ",length(cellPop)," isa eff = ",ISAeff," pcreff = ",pcrEff," seq depth = ",sequenceDepth))
  cellFreq = cellPop/sum(cellPop)
  
  bloodSample = sample(length(cellFreq),popSample,replace = T, prob = cellFreq)
  sampleCount = table(bloodSample)
  #
  # a few rounds of amplification
  #
  pcrCap = c()
  #
  # 
  #
  odds = runif(length(bloodSample))
  oi = which(odds <= ISAeff)
  bnew = bloodSample[oi]
  pcrCap = c(bnew)
  
  for (i in c(1:14)) {
    odds = runif(length(pcrCap))
    oi = which(odds < pcrEff) 
    pcrnew = pcrCap[oi]
    pcrCap = c(pcrCap,pcrnew)
  }
  
  bloodCap = pcrCap
  #
  # sequence
  #
  seq = sample(bloodCap,sequenceDepth,replace=F)
  
  simResult = table(seq)
  simResult = as.numeric(simResult)
  simResult = simResult[order(-simResult)]
  simFreq = simResult / sum(simResult)
  return(list(simResult,simFreq,sampleCount))
}


#---------------------------------------------------------------------------------------------------------
#
# simBarcodeID   (return clode ID)
#
#
#
#---------------------------------------------------------------------------------------------------------
simBarcodeID <- function(cellPop, 
                         popSample = 500000, 
                         useSimple = T, 
                         BCeff = 0.17, 
                         pcrEff = 1.0, 
                         pcrCycles = 18,
                         sequenceDepth=4000000
) {
  #message(paste0("simBarcode ",length(cellPop)," bc eff = ",BCeff," pcreff = ",pcrEff," cycles ",pcrCycles," seq depth = ",sequenceDepth))
  cellFreq = cellPop/sum(cellPop)
  
  bloodSample = sample(length(cellFreq),popSample,replace = T, prob = cellFreq)
  #
  # a few rounds of amplification
  #
  pcrCap = c()
  
  if (useSimple) {
    #
    # 
    #
    odds = runif(length(bloodSample))
    oi = which(odds <= BCeff)
    bnew = bloodSample[oi]
    pcrCap = c(bnew)
    
    for (i in c(1:pcrCycles)) {
      odds = runif(length(pcrCap))
      oi = which(odds < pcrEff) 
      pcrnew = pcrCap[oi]
      pcrCap = c(pcrCap,pcrnew)
    }
  } else {
    for (i in c(1:pcrCycles)) {
      odds = runif(length(bloodSample))
      oi = which(odds < BCeff) #0.04
      bnew = bloodSample[oi]
      odds = runif(length(pcrCap))
      oi = which(odds < pcrEff) #80
      pcrnew = pcrCap[oi]
      pcrCap = c(pcrCap,bnew,pcrnew)
    }
  }
  bloodCap = pcrCap
  #
  # sequence
  #
  seq = sample(bloodCap,sequenceDepth,replace=F)
  
  return(seq)
}

#---------------------------------------------------------------------------------------------------------
#
# simISA
#
#
#
#---------------------------------------------------------------------------------------------------------
simISA_ID <- function(cellPop, popSample = 500000, ISAeff = 0.17, pcrEff = 1.0, sequenceDepth=200000) {
  #message(paste0("simISA ",length(cellPop)," isa eff = ",ISAeff," pcreff = ",pcrEff," seq depth = ",sequenceDepth))
  cellFreq = cellPop/sum(cellPop)
  
  bloodSample = sample(length(cellFreq),popSample,replace = T, prob = cellFreq)
  #
  # a few rounds of amplification
  #
  pcrCap = c()
  #
  # 
  #
  odds = runif(length(bloodSample))
  oi = which(odds <= ISAeff)
  bnew = bloodSample[oi]
  pcrCap = c(bnew)
  
  for (i in c(1:14)) {
    odds = runif(length(pcrCap))
    oi = which(odds < pcrEff) 
    pcrnew = pcrCap[oi]
    pcrCap = c(pcrCap,pcrnew)
  }
  
  bloodCap = pcrCap
  #
  # sequence
  #
  seq = sample(bloodCap,sequenceDepth,replace=F)
  
  return(seq)
}

br = seq(0,30,0.5)

#******************************************************************************************
# generatePopulationFromDistribution
#
#   given a list of population attributes [[cell#, mean, sd]...]
#
#   1) generate a blood population
#
#******************************************************************************************
generatePopulationFromDistribution <- function(fixedList=list(),attrList=list()) {
  cellPop = c()
  #
  # fixed cells
  #
  for (fix in fixedList) {
    cellPop = c(cellPop,fix)
  }
  #
  # random dist cells
  #
  for (attr in attrList) {
    cells = attr[[1]]
    mean = attr[[2]]
    sd   = attr[[3]]
    top = rnorm(cells,mean,sd)
    top = top[which(top > 0)]
    cCount = 2^top
    cCount = floor(cCount + 0.5)
    cellPop = c(cCount,cellPop)
  }
  cellPop = cellPop[order(-cellPop)]
  names(cellPop) = c(1:length(cellPop))
  popSize = sum(cellPop)
  #message(paste0("Total population = ",floor(popSize/1000000)," million"))
  return(round(cellPop))
}


#******************************************************************************************
# sampleForBarcode
#
#   2) sample bloodSampleSize cells from the population
#   3) PCR (with efficieant pcrEff)
#   4) sequence with sample size seqSampleSize
#
#******************************************************************************************
sampleForBarcode <- function(cellPop,bloodSampleSize,seqSampleSize) {
  #
  # sample from blood
  #
  probVec = cellPop/sum(cellPop)
  pop = sample(length(cellPop),size=bloodSampleSize,replace = T,prob = probVec)
  #pop = pop[order(pop)]
  #
  # pop now represents the full dna of single cells, with the value being
  # the clone id of the parent.
  #
  # To sim PCR, each round there is a chance of replicating the barcode from the
  # clone dna and also a chance of replicating the much smaller products of
  # previous PCR rounds
  #
  pcr_product_eff = 0.88
  full_dna_eff = 0.017
  pcrProduct = c()
  for (i in c(1:11)) {
    #message(i)
    bcResult = sample(pop, size = full_dna_eff * length(pop),replace=F)
    prResult = sample(pcrProduct, size = pcr_product_eff * length(pcrProduct),replace=F)
    pcrProduct = c(pcrProduct,bcResult,prResult)
  }
  #
  # new sequence, more rounds of PCR are run in reality but numbers hard to simulate so
  # pick based on equal probability with replacement
  #
  seqResult = sample(pcrProduct,size = seqSampleSize,replace = TRUE)
  #
  # rle sum...must be sorted because rle counts runs of same value
  #
  b <- rle(sort(seqResult))
  #
  # build data frame to match actual values
  #
  df = data.frame(b$values,row.names = b$values)
  df = cbind(df,c(1:length(b$lengths)))
  df = cbind(df,b$lengths)
  df = cbind(df,b$lengths)
  df = cbind(df,b$lengths/sum(b$lengths))
  df = cbind(df,log2(b$lengths))
  df = df[order(-df[,3]),]
  colnames(df) <- c('clodeID','index','count','nCount','freq','log2Count')
  #
  # return all sampled values
  #
  blood = rle(sort(pop))
  bl = blood$lengths
  names(bl) = blood$values
  pcr   = rle(sort(pcrProduct))
  p = pcr$lengths
  names(p) = pcr$values
  return( list(df,bl,p))
}
#******************************************************************************************
# runMultipleBarcodeSamples
#
#   given a list of test populations:
#     call genAndSample multiple times, save and plot results 
#
#  bcTest     - actual test results to use for comparison
#  popLists = list[   
#                     list[ list[cell#,mean,sd] ... ] 
#                 ]
# 
# sampleBlood - number of samples to draw from blood
# pcr         - pcr Efficiency
# sampleSeq   - # sequencing samples
# ranges      - ranges for multiple plots
# labels      - labels for plot legend
# title       - plot title
# file
# 
#
#******************************************************************************************
runMultipleBarcodeSamples <- function(bcTest,popLists,sampleBlood,pcr,sampleSeq,ranges,lables,titleStr="",file=""){
  #
  # run sims
  #
  plotList = lapply(popLists,function(l) {
    gRet = genAndSampleBarcode(l,sampleBlood,pcr,sampleSeq)
    gRet[[1]] = sort(gRet[[1]],decreasing = T)
    return(gRet[[1]])
  })
  
  if (file != "") {
    pdf(file)
  }
  
  par(mfrow=c(2,2))
  par(mar=c(2,2,5,2))
  
  plotchar <- seq(18,18+3,1)
  plotchar  <- rep(plotchar,10)
  colors = c(rep('red',4),rep('green',4),rep('blue',4),rep('purple',4),rep('orange',4))
  
  maxY = bcTest[1,3]
  #  for (j in c(1:length(plotList))) {
  #    if (plotList[[j]][1] > maxY) {
  #      maxY = plotList[[j]][1]
  #    }
  #  }
  
  plot(bcTest[ranges[[1]],3],type='l',lwd=2,ylim=c(0,maxY))
  
  for (j in c(1:length(plotList))) {
    lines(plotList[[j]][ranges[[1]]],type='l',col=colors[j],cex=0.5)
  }
  lines(bcTest[ranges[[1]],3],lwd=2) # put black back on top
  
  legend(30,
         floor(maxY * 0.98),
         labels,
         cex=0.8,
         col=c('red','green','blue','purple','orange'),
         pch=rep(20,5),
         lty = rep(1,5),
         bty = 'n')
  
  
  plot(bcTest[ranges[[2]],3],type='l',col='black',lwd=2,ylim = c(0,10000))
  for (j in c(1:length(plotList))) {
    lines(plotList[[j]][ranges[[2]]],type='l',col=colors[j],cex=0.5)
  }
  lines(bcTest[ranges[[2]],3],lwd=2) # put black back on top
  
  
  plot(bcTest[ranges[[3]],3],type='l',col='black',lwd=2,ylim = c(0,2000))
  for (j in c(1:length(plotList))) {
    lines(plotList[[j]][ranges[[3]]],type='l',col=colors[j],cex=0.5)
  }
  lines(bcTest[ranges[[3]],3],lwd=2) # put black back on top
  
  plot(bcTest[ranges[[4]],3],type='l',col='black',lwd=2,ylim=c(0,500))
  for (j in c(1:length(plotList))) {
    lines(plotList[[j]][ranges[[4]]],type='l',col=colors[j],cex=0.5)
  }
  lines(bcTest[ranges[[4]],3],lwd=2) # put black back on top
  
  if (titleStr != "") {
    title(titleStr, outer=TRUE,line=-4.5,cex=0.75)  
  }
  
  if (file != "") {
    dev.off()
  }
}
#******************************************************************************************
# sampleForRIS
#
#   1) cellPop - full blood population
#   2) sample bloodSampleSize cells from the population
#   3) RIS with efficiency
#   4) sequence with sample size seqSampleSize
#   5) align to genome
#
#******************************************************************************************
sampleForRIS <- function(cellPop,bloodSampleSize,seqSampleSize) {
  #
  # sample from blood
  #
  popSize = sum(cellPop)
  probVec = cellPop/popSize
  pop = sample(length(cellPop),size=bloodSampleSize,replace = T,prob = probVec)
  #
  # RIS : some percent of cells survive
  # + PCR : may need random function
  #
  risEff = 0.05
  risSize   = floor((length(pop) * risEff) + 0.5)
  risResult = sample(pop,size=risSize,replace=T) # replace  = T...sampling from PCR pool
  #message(paste0(length(risResult)," survive RIS ",risEff))
  #
  # sim random pcr amp
  #
  pcr_product_eff = 0.88
  pcrProduct = risResult
  for (i in c(1:11)) {
    #message(i)
    prResult = sample(pcrProduct, size = pcr_product_eff * length(pcrProduct),replace=F)
    pcrProduct = c(pcrProduct,prResult)
  }
  #
  # sequence : pcr really run amy more cycles so simulate with sample with replacement
  #
  seqResult = sample(pcrProduct,size=seqSampleSize,replace=TRUE)
  #message(paste0("seq sample size = ",length(seqResult)))
  #
  # combine clones
  #
  b <- rle(sort(seqResult))
  #
  # sample complete clones for alignment
  #
  #message(paste0("Clones before align ",length(b$lengths)))
  alignEff = 0.75
  #
  #
  #
  seqArray = b$lengths
  names(seqArray) = b$values
  alignResult = sample(seqArray,size=floor(alignEff * length(seqArray)),replace = FALSE)
  #
  # 
  #
  #message(paste0("Clones after align ",length(alignResult)))
  #print(head(alignResult))
  #
  # build data frame to match actual values
  #
  df = data.frame(b$values,row.names = b$values)
  df = cbind(df,c(1:length(b$lengths)))
  df = cbind(df,b$lengths)
  df = cbind(df,b$lengths)
  df = cbind(df,b$lengths/sum(b$lengths))
  df = cbind(df,log2(b$lengths))
  df = df[order(-df[,3]),]
  colnames(df) <- c('clodeID','index','count','nCount','freq','log2Count')
  
  blood = rle(sort(pop))
  bl = blood$lengths
  names(bl) = blood$values
  pcr   = rle(sort(pcrProduct))
  p = pcr$lengths
  names(p) = pcr$values
  
  return( list(df,bl,p))
}
#******************************************************************************************
# runMultipleRisSamples
#
#  
#
#  risTest    - actual test results to use for comparison
#  pop        - list[ list[cell#,mean,sd] ... ] - assume population already determined
#               
# 
# sampleBlood - [] number of samples to draw from blood these are raw cells
# risEff      - [] RIS Efficiency : how many blood cells are gene modified and complete fragmentation and linker ligation
# sampleSeq   - sequencing samples - samples from good RIS events
# alignEff    - [] alignment efficiency
# ranges      - ranges for multiple plots
# labels      - labels for plot legend
# title       - plot title
# file        - plot file
#
#******************************************************************************************
runMultipleRisSamples <- function(risTest,pop,sampleBlood,risEff,sampleSeq,alignEff,ranges,lables,titleStr="",file=""){
  #message("runMultipleRisSamples")
  #
  # run sims, risEff, alingEff could be list
  #
  if (length(risEff) > 1) {
    
    plotList = lapply(risEff,function(l) {
      message(paste0("ris ",l))
      gRet = genAndSampleRIS(pop,sampleBlood[[1]],l,sampleSeq,alignEff[[1]])
      gRet = sort(gRet,decreasing = T)
      return(gRet)
    })
  } else if (length(alignEff) > 1)  {
    plotList = lapply(alignEff,function(l) {
      message("align")
      gRet = genAndSampleRIS(pop,sampleBlood[[1]],risEff[[1]],sampleSeq,l)
      gRet = sort(gRet,decreasing = T)
      return(gRet)
    })
  } else if (length(sampleBlood) > 1) {
    plotList = lapply(sampleBlood,function(l) {
      message("bloodSample")
      gRet = genAndSampleRIS(pop,l,risEff[[1]],sampleSeq,alignEff[[1]])
      gRet = sort(gRet,decreasing = T)
      return(gRet)
    })
  }
  
  if (file != "") {
    pdf(file)
  }
  
  par(mfrow=c(2,2))
  par(mar=c(2,2,5,2))
  
  plotchar <- seq(18,18+3,1)
  plotchar  <- rep(plotchar,10)
  
  colors = c('red','green','blue','purple','orange','pink','grey','#ffa0a0','#a0ffa0','#a0a0ff')
  
  maxY = risTest[1,3]
  
  plot(risTest[ranges[[1]],3],type='l',lwd=2,ylim=c(0,maxY))
  
  for (j in c(1:length(plotList))) {
    lines(plotList[[j]][ranges[[1]]],type='l',col=colors[j],cex=0.5)
  }
  lines(risTest[ranges[[1]],3],lwd=2) # put black back on top
  
  legend(80,
         floor(maxY * 0.98),
         labels,
         cex=0.4,
         col=colors,
         pch=rep(20,5),
         lty = rep(1,5),
         bty = 'n')
  
  
  plot(risTest[ranges[[2]],3],type='l',col='black',lwd=2)
  for (j in c(1:length(plotList))) {
    lines(plotList[[j]][ranges[[2]]],type='l',col=colors[j],cex=0.5)
  }
  lines(risTest[ranges[[2]],3],lwd=2) # put black back on top
  
  
  plot(risTest[ranges[[3]],3],type='l',col='black',lwd=2)
  for (j in c(1:length(plotList))) {
    lines(plotList[[j]][ranges[[3]]],type='l',col=colors[j],cex=0.5)
  }
  lines(risTest[ranges[[3]],3],lwd=2) # put black back on top
  
  plot(risTest[ranges[[4]],3],type='l',col='black',lwd=2,ylim=c(0,100))
  for (j in c(1:length(plotList))) {
    lines(plotList[[j]][ranges[[4]]],type='l',col=colors[j],cex=0.5)
  }
  lines(risTest[ranges[[4]],3],lwd=2) # put black back on top
  
  if (titleStr != "") {
    title(titleStr, outer=TRUE,line=-4.5,cex=0.75)  
  }
  
  if (file != "") {
    dev.off()
  }
}
#-----------------------------------------------------------------------------------------------
#
# track clones for given bcObj and list index
#
#-----------------------------------------------------------------------------------------------
CloneTrack <- function(bcObj,dataIndex) {
  ml = 0
  divider = ""
  for (l in bcObj@labels) {
    s = nchar(l)
    if (s > ml) {
      ml = s
    }
    divider = paste0(divider,"------")
  }
  
  
  for (i in c(1:ml)) {
    sp = ""
    for (l in bcObj@labels) {
      ch = substr(l,i,i)
      if (nchar(ch) == 0) {
        ch = " "
      }
      sp = paste0(sp,"     ",ch)
    }
    print(sp)
  }
  
  print(divider)
  
  for (i in c(1:100)) {
    bci = c()
    sp = " "
    
    for (j in bcObj@data) {
      index = which(j[,1] == bcObj@data[[dataIndex]][i,1])
      if (length(index)==0) {
        index = -1
      }
      bci = c(bci,index)
      sp = paste0(sp,sprintf("%5d ",index))
    }
    print(sp)
  }
}


#-----------------------------------------------------------------------------------------------
#
# GetExpressionMatrix for a given clone in all samples
#
#-----------------------------------------------------------------------------------------------
getCloneRankExpression <- function(bcObj,dataIndex,cloneIndex) {
  
  target = bcObj@data[[dataIndex]][cloneIndex,1]
  
  bci = c()
  
  for (j in bcObj@data) {
    index = which(j[,1] == target)
    if (length(index)==0) {
      exp = 0.0
    } else {
      exp = j[index,5]
    }
    bci = c(bci,exp)
  }
  return(bci)
}