#******************************************************************************************
#
#   defince the data storage class and provide common functions
#
#   
#   M Enstrom 1-1-2018
#
#******************************************************************************************
setClass(
  "BCObj",
  representation(
    objName="character",
    fileNames="character", 
    labels="character",
    short="character",
    sizes="numeric",
    data="list",
    filterData="list",
    global = "data.frame",
    colorMap = "matrix",
    testsToRun = "numeric"
  ),
  prototype = list(
    objName=character(),
    fileNames=character(), 
    labels=character(),
    short=character(),
    sizes=numeric(),
    data= list(),
    filterData=list(),
    global=data.frame(),
    colorMap=matrix(0,0,0),
    testsToRun = numeric())
)
