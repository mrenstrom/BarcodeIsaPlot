#******************************************************************************************
#
#   load source routines and data
#
#   
#   M Enstrom 1-1-2018
#
#******************************************************************************************
library("pheatmap")
library("RColorBrewer")
library("devtools")
library("ggplot2")
library("cluster")
library("gplots")


baseDir = "/Users/mark_enstrom/Barcode/"
plotDir = paste0(baseDir,"FinalBC_RIS/plots/")
isa8103Dir    = paste0(baseDir,"RIS/Z08103/final_files/")
isa9132Dir    = paste0(baseDir,"RIS/Z09132/final_files/")
bc8103Dir     = paste0(baseDir,"Z08103/final_files/")
bc9132Dir     = paste0(baseDir,"Z09132/final_files/")
#
# include files
#
setwd(paste0(baseDir,"FinalBC_RIS"))

source("BCObj.R")
source("LoadSourceObject.R")
source("buildColorMaps.R")
source("simFunction.R")
source("followBarcode.R")
source("finalBarplots.R")

#----------------------------------------------------------------------------------------------------
#
# init source data objects
#
#----------------------------------------------------------------------------------------------------
#
# Z08103
#
bc8103      = initBarcodeObject(bc8103Dir,baseDir,"DBS Z08103")
isa8103     = initISAObject(isa8103Dir,baseDir,"ISA Z08103")
#
# Z09132
#
bc9132      = initBarcodeObject(bc9132Dir,baseDir,"DBS Z09132")
isa9132     = initISAObject(isa9132Dir,baseDir,"ISA Z09132")
#
# load color maps
#
colArray  = buildColorMap()
heatArray = buildHeatColorMap()
#
# build permanent random colormaps for each bc and ris
#
#
#bc8103@colorMap  = buildRandomColorMap(bc8103)
#ris8103@colorMap = buildRandomColorMap(ris8103)
#isa8103@colorMap = buildRandomColorMap(isa8103)
#
#bc9132@colorMap  = buildRandomColorMap(bc9132)
#ris9132@colorMap = buildRandomColorMap(ris9132)
#isa9132@colorMap = buildRandomColorMap(isa9132)