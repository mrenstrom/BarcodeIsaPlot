#
# for a given list of barcodes, find the barcode in each condition
# and build a matrix. plot the matrix
# 
#
#
#----------------------------------------------------------------------------------------------------
# mapEntries
#
#  given list if dataframes of form
#     [id] [data1] ....  [dataN]
#     ----------------------------
#  
#  build a new matrix of (targetName) in ALL DFs
#
#
#  for all targetNames
#     make a new matrix row
#     search each df for id = targetName and add [dataColumn] to matrix if id exists
#
#----------------------------------------------------------------------------------------------------
mapEntries <- function(targetNames,dfList,dfNames,dataColumn) 
{
  #
  # mixBarcodes are the barcodes to find
  #
  bc_row = length(targetNames)
  bc_col = length(dfList)
  #
  # for log use 5
  # for raw counts use 2
  #
  fred <- matrix(nrow=0,ncol=bc_col)
  CN   <- matrix(nrow=0,ncol=bc_col)
  for (i in c(1:bc_row)) {
    newRow = matrix(nrow=1,ncol=bc_col) 
    newCN  = matrix(nrow=1,ncol=bc_col) 
    #
    # for each condition, find the index of given barcode
    #
    for (j in c(1:bc_col)) {
      cName = "NA"
      count = NA
      index = which(dfList[[j]][,1] == targetNames[i])
      if (length(index) != 0) {
        count = dfList[[j]][index,dataColumn]
        cName = paste0(index)
      } 
      newRow[j] = as.numeric(count)
      newCN[j]  = cName
    }
    fred <- rbind(fred,newRow)
    CN   <- rbind(CN,newCN)
  }
  rownames(fred) <- targetNames
  colnames(fred) <- dfNames
  
  return(list(fred,CN))
}
#----------------------------------------------------------------------------------------------------
# followBarcode
#
#
#
#----------------------------------------------------------------------------------------------------

followBarcode <- function(mainLab="",mixBarcodes,rowIndexNames,conds,names,column,colorMap,outName="",doPlot=TRUE,showKey=FALSE)
{
  print(paste0("FollowBarcode, showKey = ",showKey))
  
  fredList <- mapEntries(mixBarcodes,conds,names,column)
  
  fred <- fredList[[1]]
  cn  <- fredList[[2]]
  
  cn[which(cn != '1')] = ""
  cn[which(cn == '1')] = "*"
  
  
  
  bc_row = length(mixBarcodes)
  bc_col = length(conds)
  
  par(mar=c(2,2,12,2))
  #
  #
  #
  if (doPlot == TRUE)
  {
    par(mfrow=c(1,1))
    par(mar=c(5,5,2,2))
    if (outName != "") {
      message(paste0("send plot to ",outName))
      
      pdf(outName,
          width = 10, 
          height = 7
      )
    }
    
    message("Call heatmap")
    
    #
    # lmat stands for layout matrix, and in heatmap.2 the default is a 2x2 matrix used to organize four components 
    # (1=heatmap, 2=row dendogram, 3=col dendogram, 4= key ) 
    # 
    heatmap.2(fred,
              cellnote=cn,
              labRow = rowIndexNames,
              cexRow = 1.25,
              dendrogram = 'none',
              Rowv = 'NA',
              Colv = 'NA',
              notecol = 'black',
              key=showKey,
              #    uncomment to make a giant key   keysize = 10,
              density.info = 'none',
              main = mainLab,
              trace='none',
              margins =c(10,8), #12,6
              colsep = c(1:dim(fred)[2]),  # comment these two lines out to make color scale
              rowsep=c(1:dim(fred)[1]),    # when showing top 250 clones
              sepcolor = "#c0c0c0",
              #ylab="Barcode index",
              #scale = 'row',
              notecex = 2.0,
              na.color = "#000000",
              #lmat=rbind(c(4,3),c(2,1)),
              #lhei = c(.25,1),
              #lwid = c(.05,1),
              srtCol = 90,
              cexCol = 1.0,
              col=colorMap
    )
    #title("Title",cex=0.75,outer = T,line=-2.0)
    
    if (outName != "") {
      dev.off()  
    }
  }
  return(list(fred,cn))
}

clone_track <- function(bcobj,index=0,n=25,findLocal=FALSE,doPlot=FALSE,inputConds=c(),showKey=FALSE)
{
  #debug
  #bcobj = bc9132
  #index = 0
  #findLocal = FALSE
  #doPlot = FALSE
  #n = 25
  #inputConds=c()
  #
  # get color map !!! why isn't this global?
  #
  print(paste0("clone_track, showKey = ",showKey))
  iArray = buildIndexColorMap2()
  if (length(inputConds) == 0) {
    goodTests = bcobj@testsToRun
  } else {
    goodTests = inputConds
  }
  conds     = bcobj@data[goodTests]
  namesC    = bcobj@short[goodTests]
  #
  # !!! until this is added to main code:
  # add index column
  #
  for (j in c(1:length(conds))) {
    nr = nrow(conds[[j]])
    conds[[j]]$ic = log2(c(1:nr)+1)
  }
  #
  # use global max barcodes if index == 0
  #
  if (index == 0) {
    plotName = paste0(bcobj@objName," Global MAX")
    maxBarcodes = bcobj@global[1:n,1]
  } else if (findLocal == TRUE) {
    plotName = paste0(namesC[index]," by global")
    fred = conds[[index]][,1]
    bob = sort(fred)
    maxBarcodes = bob[1:n]
    
  } else {
    plotName = paste0(bcobj@objName," ",namesC[index])
    maxBarcodes = conds[[index]][1:n,1]
  }
  #
  # rank or full name
  #
  rowIndexNames = c(1:n)
  #rowIndexNames = maxBarcodes
  if (doPlot==TRUE) {
    name = paste0(plotDir,bcobj@objName,"_clonetrack",index,".pdf")
    message(name)
  } else {
    name=""
    message("plot to screen")
  }
  #
  # for followBarcde....doPlot == FALSE meanse just return matricies
  #    name == "" means plot to screen
  #    name != "" means plot to file
  # 
  r = followBarcode(plotName,maxBarcodes,rowIndexNames,conds,namesC,11,iArray,outName =  name,doPlot = TRUE, showKey=showKey)
}

#------------------------------------------------------------------------------------------
#
# track a fixed list of clones
#
#
#
#------------------------------------------------------------------------------------------

clone_track_fixed <- function(bcobj,maxBarcodes,doPlot=FALSE,inputConds=c(),showKey=FALSE)
{
  n = length(maxBarcodes)
  if (n < 1) {
    print("Length of clone vector is 0")
    return()
  }
  print(paste0("clone_track, showKey = ",showKey))
  iArray = buildIndexColorMap2()
  if (length(inputConds) == 0) {
    goodTests = bcobj@testsToRun
  } else {
    goodTests = inputConds
  }
  conds     = bcobj@data[goodTests]
  namesC    = bcobj@short[goodTests]
  #
  # !!! until this is added to main code:
  # add index column
  #
  for (j in c(1:length(conds))) {
    nr = nrow(conds[[j]])
    conds[[j]]$ic = log2(c(1:nr)+1)
  }
  #
  # rank or full name
  #
  rowIndexNames = c(1:n)
  #rowIndexNames = maxBarcodes
  if (doPlot==TRUE) {
    name = paste0(plotDir,bcobj@objName,"_fixed_clone_track.pdf")
    message(name)
  } else {
    name=""
    message("plot to screen")
  }
  plotName = "Fixed Clone Track"
  #
  # for followBarcde....doPlot == FALSE meanse just return matricies
  #    name == "" means plot to screen
  #    name != "" means plot to file
  # 
  r = followBarcode(plotName,maxBarcodes,rowIndexNames,conds,namesC,11,iArray,outName =  name,doPlot = TRUE, showKey=showKey)
}


knownBarcodes = c(
  "GCGTTGGGCGTTATGTCCTT", #1
  "GGGCGGTGTCTTGGCCCAGA", #2
  "ACAAAAGAGAACTCGGAGGC", #3
  "AATTGAACTGTTCAGTCAGC", #4
  "TGTCCCTCTAGGTCATACGT", #5
  "TATTTAACTAATGTCTGGCA", #6
  #"TGGTGGGTTACGGAGGAGAT", #7
  "GATGCGCTAATGATGTGACA", #8
  #"GTAGCTTCTACTATCAGCCA", #9
  #"GTCTTAGGACTGTCAGCTAT", #10
  #"GTTCTGTAGTGTGACGGGAC", #11
  "TGGTAGGATCCAGATGATTG", #12
  "TAAAATGGTATCCCAGGCGC", #13
  "AGTGTAGGGGATTAACAACC", #14
  #"TTAGTTTAGGTTCTGGCTTG", #15
  #"GGGACCTAGACTGGAAGCTA", #16
  #"TAATTCACGGTGCGATGCGT", #17
  "AGTAGAACGAAGCTGGACAC", #18
  #"GGAGAACCGAGGGCAACTCG", #19
  #"AGGAGTTTGCGTCAGGATGC", #20
  #"ACTGGAAGCGTGATAGAACA", #21
  #"CTCGGGTATAAGGTTCGGCC", #22
  #"CACATGAGGCGCCTCAGGGA", #23
  "TATTTCTCTAGCCGGTTGAA" #24
  #"CGCAGTTGGCTCTACTACAA" #25
)


clone_track_fixed(bc8103,knownBarcodes,doPlot = TRUE)

knownISA = c(
  "Z08103_000020", #1
  "Z08103_000049", #2
  "Z08103_000039", #3
  #"Z08103_000051", #4
  "Z08103_000041", #5
  "Z08103_000057", #6
  "Z08103_000056", #7
  #"Z08103_000054", #8
  #"Z08103_000066", #9
  "Z08103_000069", #10
  #"Z08103_000068", #11
  #"Z08103_000064", #12
  #"Z08103_000013", #13
  "Z08103_000072", #14
  #"Z08103_000062", #15
  #"Z08103_000069", #16
  #"Z08103_000070", #17
  "Z08103_000074", #18
  "Z08103_000082", #19
  "Z08103_000029", #20
  "Z08103_000094" #21
  #"Z08103_000084", #22
  #"Z08103_000065", #23
  #"Z08103_000001", #24
  #"Z08103_000025"  #25
)
clone_track_fixed(isa8103,knownISA,doPlot = TRUE)
