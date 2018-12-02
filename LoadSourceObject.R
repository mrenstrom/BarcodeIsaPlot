#----------------------------------------------------------------------------------------------------
# initBarcodeObject
#   creates a "BCOBJ" for barcde data
#
#
#----------------------------------------------------------------------------------------------------

initBarcodeObject <- function(bcDir,baseDir,objName)
{
  #debug
  #baseDir = "/Users/mark_enstrom/Barcode/"
  #bcDir   = "/Users/mark_enstrom/Barcode/RIS/Z09132/final_files/"
  #objName = "Z09132"
  
  #bcDir = "/Users/mark_enstrom/Barcode/Z08103/final_files/"
  #baseDir = "/Users/mark_enstrom/Barcode/"
  #objName = "Z08103"
  
  
  bcObj = new("BCObj")  
  bcObj@objName   = objName
  #
  # load files and labels from bc directory
  #
  df <- read.csv(paste0(bcDir,"RInputList.txt"),stringsAsFactors = FALSE)
  bcObj@fileNames = df$filename
  bcObj@labels = df$label
  bcObj@short  = df$short
  #
  # read data files
  #
  files <- vector("list",length(bcObj@fileNames))
  #
  #
  #
  for (i in c(1:length(bcObj@fileNames))) {
    name = paste0(bcDir,bcObj@fileNames[i])
    files[[i]] <- read.csv(name,stringsAsFactors = FALSE,comment.char = '#',header = T)
    rownames(files[[i]]) <- files[[i]][,1]
    bcObj@sizes[i] = dim(files[[i]])[1]
    bcObj@labels[i] <- paste0(bcObj@labels[i]," ",dim(files[[i]])[1])
  }
  bcObj@data = files
  #
  # load global insertion list
  #
  bcObj@global = read.table(paste0(bcDir,"globalID.txt"),stringsAsFactors = F,header = F)
  rownames(bcObj@global) <- bcObj@global[,1]
  bcObj@testsToRun = c(1:length(bcObj@fileNames))
  return(bcObj)
}

initISAObject <- function(isaDir,baseDir,objName)
{
  #debug
  #baseDir = "/Users/mark_enstrom/Barcode/"
  #bcDir   = "/Users/mark_enstrom/Barcode/RIS/Z09132/final_files/"
  #objName = "Z09132"
  #isaDir = "/Users/mark_enstrom/Barcode/RIS/Z08103/final_files/"
  #baseDir = "/Users/mark_enstrom/Barcode/"
  #objName = "Z08103"
  
  
  isaObj = new("BCObj")  
  isaObj@objName   = objName
  #
  # load files and labels from bc directory
  #
  df <- read.csv(paste0(isaDir,"RInputList.txt"),stringsAsFactors = FALSE)
  isaObj@fileNames = df$filename
  isaObj@labels = df$label
  isaObj@short = df$short
  #
  # read data files
  #
  files <- vector("list",length(isaObj@fileNames))
  #
  #
  #
  for (i in c(1:length(isaObj@fileNames))) {
    name = paste0(isaDir,isaObj@fileNames[i])
    td <- read.table(name,stringsAsFactors = FALSE,comment.char = '#',header = T)
    td1 <- td[which(td$oCount > 1),]
    td2 <- td[which(td1$align != "U"),]
    #message(paste0(i," td = ",nrow(td)," td1 = ",nrow(td1)," td2 = ",nrow(td2)," ",isaObj@labels[i]))
    
    files[[i]] <- td2
    rownames(files[[i]]) <- files[[i]][,1]
    isaObj@sizes[i] <- dim(files[[i]])[1]
    isaObj@labels[i] <- paste0(isaObj@labels[i]," ",dim(files[[i]])[1])
  }
  isaObj@data = files
  #
  # load global insertion list
  #
  isaObj@global = read.table(paste0(isaDir,"globalID.txt"),stringsAsFactors = F,header = F)
  rownames(isaObj@global) <- isaObj@global[,1]
  isaObj@testsToRun = c(1:length(isaObj@fileNames))
  return(isaObj)
}
