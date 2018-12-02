# 
# CloneBarplots
#  11-5-2018 Mark Enstrom
#
# quick generation of clone tracking barplosts
#
timesSeen <- function(bcObj,cl_id,numConditions) {
  ls = lapply(c(1:numConditions),function(y) {
    f = bcObj@data[[y]][cl_id,5]
    return (!(is.na(f)))
  })
  return(sum(unlist(ls)))
}
#----------------------------------------------------------------------
# 
# getDiffColor - generate a random color with at least a 
# minumum difference
#
#----------------------------------------------------------------------
getDiffColor <- function(r,g,b) {
  while (TRUE) {
    fr = runif(1)
    fg = runif(1)
    fb = runif(1)
    #
    # color saturation
    #
    cMax = max(fr,fg,fb)
    cMin = min(fr,fg,fb)
    if (cMax > 0.0) {
      sat = (cMax-cMin)/cMax
    } else {
      sat = 0
    }
    
    r1 = floor(fr*255)
    g1 = floor(fg*255)
    b1 = floor(fb*255)
    #
    # absolute difference from prev color
    #
    d = sqrt( (r1-r)^2 + (g1-g)^2 + (b1-b)^2)
    #
    # color intensity
    #
    intensity = r1 + g1 + g1
    
    #message(paste0(intensity," ",sat))
    #
    #
    #
    if ((d > 300.0) & (sat > 0.6) & (intensity > 300) & (intensity < 600)) {
      return(c(r1,g1,b1))
    }
  }
}
#----------------------------------------------------------------------
#  getNewClones - given a list of clones, a bcData object and a list
#  if indexes find unique new clones
#----------------------------------------------------------------------
getNewClones <- function(clones,dataObj,newCloneIndex) {
  lt = lapply(newCloneIndex,function(x){
    c = dataObj[x,1]
    if (c %in% clones) {
      #message("already in")
    } else {
      #message("new")
      return(c)
    }
  })
  return(unlist(lt))
}

#----------------------------------------------------------------------
#  buildCloneBarplot
#     given bcObj and dpt[index] ...
#         get clones from each time point >= minFraction
#         combine into one clone list
#         build dataframe compatible with ggplot
#         build bar plot
#
#----------------------------------------------------------------------



buildCloneBarplot<- function(bcObj,
                             dpt,
                             minF = 0.01,
                             maxF = 1.0,
                             plotType = "Bar",
                             bkColor = "#505050",
                             fixedClones = c()) {
  # debug start
  #bcObj = isa8103
  #dpt = c(1,3,4,8,10)
  #minF = 0.005
  #
  # sel limits
  #
  numTests = length(dpt)
  numConditions = length(bcObj@labels)
  #
  # get list of clones with expression >= 1% for each data object
  #
  print(length(fixedClones))
  if (length(fixedClones) == 0) {
    print("Build Clone List from Data")
    
    clones = c()
    for (i in dpt) {
      l1p = which(bcObj@data[[i]][,5] >= minF)
      nc  = getNewClones(clones,bcObj@data[[i]],l1p)
      clones = c(clones,nc)
    }
    # small is for the remainder (all clones < 1%)
    clones = c(clones,"small")
    iSmall = length(clones)
    #
    # now revers order so lowest starts on top
    #
    clones
    clones = clones[iSmall:1]
    clones
  } else {
    clones = fixedClones
    iSmall = length(clones)
    # reverse order
    clones = clones[iSmall:1]
    
    
    print("Init from fixed")
    print("length fixed = ")
    print(length(fixedClones))
    print(length(clones))
    print(paste0("iSmall = ",iSmall))
  }
  #
  # how many times is each clone seen?
  #
  lt = lapply(c(1:(iSmall)),function(x) {
    #
    # for clone x, how many times was it seen
    #
    cl_id = clones[x]
    ls = lapply(c(1:numConditions),function(y) {
      f = bcObj@data[[y]][cl_id,5]
      return (!(is.na(f)))
    })
    #
    # was clone seen only one time
    #
    return(sum(unlist(ls)))
  })
  seen = unlist(lt)
  seen
  which(seen == 1)
  #
  # build from fractions
  #
  fracs = c()
  for (i in dpt) {
    fr = bcObj@data[[i]][clones,5]
    fr[is.na(fr)] <- 0
    #
    # the remainder of clones below threshold
    #
    theRest = maxF - sum(fr)
    # all the remainder
    if (theRest < 0) {
      theRest = 0
    }
    fr[1] = theRest
    # accumulate
    fracs = c(fracs,fr)
  }
  #
  # get labels
  #
  dptLabels = bcObj@labels[dpt]
  message(paste0(" ",dptLabels))
  #
  # shorten to just the number will break with > 999 dpt
  #
  for (i in c(1:numTests)) {
    s = dptLabels[i]
    pos = regexpr('DPT',s)
    if (pos == -1) {
      pos = regexpr('_',s)
    }
    keep = substr(s,0,pos-1)
    keep = sprintf("%03d",as.numeric(keep))
    dup = 1
    keep_base = keep
    while (TRUE) {
      if (keep %in% dptLabels) {
        keep = sprintf("%s_%d",keep_base,dup)
        dup = dup + 1
      } else {
        break
      }
    }
    dptLabels[i] = keep
  }
  
  #message(paste0(" ",dpt))
  #message(paste0(" ",dptLabels))
  #
  # build barcode data
  #
  DPT   = rep(dptLabels,each = iSmall)
  #message(paste0(" ",DPT))
  #message(length(clones))
  #
  # fix factor order: order is reversed
  #
  cloneFactor = factor(clones,levels = clones[1:iSmall])
  seenOneTime = which(seen == 1)
  seenOneTime
  cloneFacRep = rep(cloneFactor,times=numTests)
  message(length(DPT))
  message(length(fracs))
  message(length(cloneFacRep))
  #
  # make plot database
  #
  df <- data.frame(labs = cloneFacRep, clones = fracs, dpt = DPT)
  #
  # generate random color map
  #
  r = 0
  g = 255
  b = 128
  lt = lapply(c(1:iSmall), function(x) {
    newColor = getDiffColor(r,g,b)
    r = newColor[1]
    g = newColor[2]
    b = newColor[3]
    sr = as.hexmode(floor(runif(1) * 255))
    sg = as.hexmode(floor(runif(1) * 255))
    sb = as.hexmode(floor(runif(1) * 255))
    sc = sprintf("#%02x%02x%02x",sr,sg,sb)
    return(sc)
  })
  colorScale = unlist(lt)
  # change constant colors
  colorScale[1] = bkColor
  colorScale[seenOneTime] = "#ffffff"
  #colorScale[368] = "#ff0000"
  
  #colorScale[294] = "#000000"
  #colorScale[295] = "#ff0000"
  #colorScale[296] = "#ffffff"
  #colorScale[297] = "#0000ff"
  #colorScale[298] = "#000000"
  
  print(paste("length of color scale = ",length(colorScale)))
  
  #
  # labels for top of bars = Total number of clones
  #
  lt = lapply(dpt,function(x) {
    nrow(bcObj@data[[x]])
  })
  sizeText = unlist(lt)
  
  dfp = data.frame(x=c(1:numTests),y=rep(c(maxF),times=numTests),labs=rep(c('small'),times=numTests),txt=sizeText)
  
  g = ""
  
  if (plotType == "Bar") {
    #
    # barplot 
    #
    g <- ggplot(df,aes(dpt,clones,fill=labs)) + geom_bar(stat = "identity") + 
      scale_fill_manual(values=colorScale) + 
      guides(fill=FALSE) +
      theme(text = element_text(size=10)) #+
      #position = position_dodge(width=0.9), 
      #geom_text(dfp,mapping=aes(x=x,y=y,label=txt), size=3,nudge_y = +0.02)
  } else if (plotType == "Area") {
    df$dpt = as.numeric(df$dpt)
    g <- ggplot(df, aes(x=dpt, y=clones, fill=labs)) + 
      geom_area(stat='identity') + 
      scale_fill_manual(values=colorScale) + 
      scale_x_continuous(breaks = c(1:length(dptLabels)),labels=dptLabels) + 
      guides(fill=FALSE) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+
      #geom_text(dfp,mapping=aes(x=x,y=y,label=txt), size=3,nudge_y = +0.02)
    
    #g
  }
  
  
  return(g)
}

#----------------------------------------------------------------------
#  Build master list of clones >= 1% for Z08103 for 
#  ISA 33,160,233,368,502 DPT
#
#
# Defined as a function so "Source" works on this file but
# meant to be run as script
#
#----------------------------------------------------------------------

runTests <- function() {
  
  isa8103@labels
  
  g <- buildCloneBarplot(isa8103,c(1,4,5,9,11),minF = 0.005)
  pdf(paste0(plotDir,"Z08103_ISA_cloneBar.pdf"))
  g + 
    labs(x= "Days post transplant",y = "Frequency of early repopulating clones") +
    scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
    ggtitle("Z08103 ISA") + 
    theme(plot.title=element_text(size='16',hjust = 0.5))
  dev.off()
  
  
  bc8103@labels
  g <- buildCloneBarplot(bc8103,c(1,4,6,18,20),minF = 0.005)
  pdf(paste0(plotDir,"Z08103_BC_cloneBar.pdf"))
  g + 
    labs(x= "Days post transplant",y = "Frequency of early repopulating clones") +
    scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
    ggtitle("Z08103 DBS") + 
    theme(plot.title=element_text(size='16',hjust = 0.5))
  
  dev.off()
  #
  # Z09132 ISA
  #
  isa9132@labels
  g <- buildCloneBarplot(isa9132,c(1,4,8,11,12),minF = 0.005)
  pdf(paste0(plotDir,"Z09132_ISA_cloneBar.pdf"))
  g + 
    labs(x= "Days post transplant",y = "Frequency of early repopulating clones") +
    scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
    ggtitle("Z09132 ISA") + 
    theme(plot.title=element_text(size='16',hjust = 0.5))
  
  dev.off()
  #
  # Z09132 DBS
  #
  bc9132@labels
  g <- buildCloneBarplot(bc9132,c(1,2,4,16,17),minF = 0.005,plotType = "Bar")
  pdf(paste0(plotDir,"Z09132_BC_cloneBar.pdf"))
  g + 
    labs(x= "Days post transplant",y = "Frequency of early repopulating clones") +
    scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
    ggtitle("Z09132 DBS") + 
    theme(plot.title=element_text(size='16',hjust = 0.5))
  
  dev.off()
  #
  # try area
  #
  
  #
  # Z08103 ISA
  #
  isa8103@labels
  g <- buildCloneBarplot(isa8103,c(1,4,5,9,11),minF = 0.005,plotType = "Area")
  pdf(paste0(plotDir,"Z08103_ISA_area.pdf"))
  g + 
    labs(x= "Days post transplant",y = "Frequency of early repopulating clones") +
    scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
    ggtitle("Z08103 ISA") + 
    theme(plot.title=element_text(size='16',hjust = 0.5),
          axis.text=element_text(size=14),
          axis.title = element_text(size=14),
          legend.title = element_text(size=16))
  dev.off()
  
  
  bc8103@labels
  g <- buildCloneBarplot(bc8103,c(1,4,6,18,20),minF = 0.005,plotType = "Area")
  pdf(paste0(plotDir,"Z08103_BC_area.pdf"))
  g + 
    labs(x= "Days post transplant",y = "Frequency of early repopulating clones") +
    scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
    ggtitle("Z08103 DBS") + 
    theme(plot.title=element_text(size='16',hjust = 0.5),
          axis.text=element_text(size=14),
          axis.title = element_text(size=14),
          legend.title = element_text(size=16))
  
  dev.off()
  
  #
  # Z09132 ISA
  #
  g <- buildCloneBarplot(isa9132,c(1,4,8,11,12),minF = 0.005,plotType = "Area")
  pdf(paste0(plotDir,"Z09132_ISA_area.pdf"))
  g + 
    labs(x= "Days post transplant",y = "Frequency of early repopulating clones") +
    scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
    ggtitle("Z09132 ISA") + 
    theme(plot.title=element_text(size='16',hjust = 0.5),
          axis.text=element_text(size=14),
          axis.title = element_text(size=14),
          legend.title = element_text(size=16))
  
  dev.off()
  #
  # Z09132 DBS
  #
  g <- buildCloneBarplot(bc9132,c(1,2,4,16,17),minF = 0.005,plotType = "Area")
  pdf(paste0(plotDir,"Z09132_BC_area.pdf"))
  g + 
    labs(x= "Days post transplant",y = "Frequency of early repopulating clones") +
    scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
    ggtitle("Z09132 DBS") + 
    theme(plot.title=element_text(size='16',hjust = 0.5),
          axis.text=element_text(size=14),
          axis.title = element_text(size=14),
          legend.title = element_text(size=16))
  dev.off()
  
  
  
  #
  # test
  #
#   g <- buildCloneBarplot(bc9132,c(1,2,3,4,5,17,18),minF = 0.005,maxF = 0.5,plotType = "Area")
#   g + 
#     labs(x= "Days post transplant",y = "Frequency of early repopulating clones") +
#     scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
#     ggtitle("Z09132 DBS") + 
#     theme(plot.title=element_text(size='16',hjust = 0.5))
#   
#   
#   
#   bc8103@labels
#   
#   g <- buildCloneBarplot(bc8103,c(18,19,20,21,22,23,24,25,26,27,28,29),minF = 0.001,maxF = 1.0,plotType = "Area")
#   g + 
#     labs(x= "Days post transplant",y = "Frequency of early repopulating clones") +
#     scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
#     ggtitle("Z08103 DBS") + 
#     theme(plot.title=element_text(size='16',hjust = 0.5))
#   
#   g <- buildCloneBarplot(bc8103,c(1,4,6,18,22,23,24,29),minF = 0.001,maxF = 1.0,plotType = "Bar")
#   g + 
#     labs(x= "Days post transplant",y = "Frequency of early repopulating clones") +
#     scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
#     ggtitle("Z08103 DBS") + 
#     theme(plot.title=element_text(size='16',hjust = 0.5))
#   
#   
#   isa8103@labels
#   
#   g <- buildCloneBarplot(isa8103,c(1,3,4,8,16,17,18,19),minF = 0.001,maxF = 1.0,plotType = "Area")
#   g + 
#     labs(x= "Days post transplant",y = "Frequency of early repopulating clones") +
#     scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
#     ggtitle("Z08103 ISA") + 
#     theme(plot.title=element_text(size='16',hjust = 0.5))
}


