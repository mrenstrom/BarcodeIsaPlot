#******************************************************************************************
#
#  track clones - generate clone tracking plots
#
#
#   
#   M Enstrom 4-11-2018
#
#******************************************************************************************
#
#
source("InitRisBC.R")
source("followBarcode.R")

#
# make sure not to inclue "ALL" combined test cases
#

isa9132@labels

bc9132@labels
bc9132@testsToRun=c(1:26)

isa8103@labels

bc8103@labels
bc8103@testsToRun = c(1:29)





#
# these are used for the paper
#
dev.off()
clone_track(bc9132,0,doPlot=TRUE,n=25,showKey=FALSE)
clone_track(isa9132,0,doPlot=TRUE,n=25,showKey=FALSE)
clone_track(bc8103,0,doPlot=TRUE,n=25,showKey=FALSE)
clone_track(isa8103,0,doPlot=TRUE,n=25,showKey=FALSE)
#
# these are for color key
#
clone_track(bc9132,index=0,doPlot=FALSE,n=25,showKey=TRUE)
clone_track(bc9132,index=26,doPlot=TRUE,n=250,showKey=FALSE)
#
# the rest are tests
#

clone_track(bc9132,0,doPlot=TRUE,n=25,showKey=TRUE)

clone_track(bc9132,0,doPlot=TRUE,n=25,showKey=FALSE)
clone_track(isa9132,0,doPlot=TRUE,n=25,showKey=FALSE)



clone_track(bc9132,1,doPlot=TRUE,n=50)
clone_track(bc9132,2,doPlot=FALSE,n=50)
clone_track(bc9132,3,doPlot=FALSE,n=50)
clone_track(bc9132,4,doPlot=FALSE,n=50)
clone_track(bc9132,16,doPlot=FALSE,n=50)
clone_track(bc9132,22,doPlot=FALSE,n=50)
clone_track(bc9132,26,doPlot=TRUE,n=50)




clone_track(isa9132,0,doPlot=TRUE,n=50)
isaclone_track(isa9132,1,doPlot=FALSE,n=50)
clone_track(isa9132,2,doPlot=FALSE,n=50)
clone_track(isa9132,3,doPlot=FALSE,n=50)
clone_track(isa9132,23,doPlot=FALSE,n=50)



clone_track(bc8103,0,doPlot=TRUE,n=50)
clone_track(bc8103,1,doPlot=FALSE,n=50)
clone_track(bc8103,2,doPlot=FALSE,n=50)
clone_track(bc8103,3,doPlot=FALSE,n=50)
clone_track(bc8103,4,doPlot=FALSE,n=50)


clone_track(bc8103,19,doPlot=FALSE,n=50)
clone_track(bc8103,26,doPlot=FALSE,n=25)




clone_track(isa8103,0,doPlot=TRUE,n=50)
clone_track(isa8103,1,doPlot=TRUE,n=50)
clone_track(bc8103,1,doPlot=TRUE,n=50)


clone_track(isa8103,2,doPlot=FALSE,n=50)
clone_track(isa8103,3,doPlot=FALSE,n=50)

clone_track(isa8103,15,doPlot=FALSE,n=50)
clone_track(isa8103,17,doPlot=FALSE,n=50)
clone_track(isa8103,23,doPlot=FALSE,n=50)
clone_track(ris8103,0,doPlot=FALSE,n=50)
clone_track(ris8103,18,doPlot=TRUE,n=50)
#----------------------------------------------------------------------
#  Z14004
#
# global max
#
clone_track(isa4004,0,doPlot=TRUE,n=50)
#
# NK cells... and Nk
#
clone_track(isa4004,16,doPlot=TRUE,n=50)
clone_track(isa4004,23,doPlot=TRUE,n=50)
clone_track(isa4004,32,doPlot=TRUE,n=50)
clone_track(isa4004,46,doPlot=TRUE,n=50)
clone_track(isa4004,55,doPlot=TRUE,n=50)
#
# Tcells
#
clone_track(isa4004,24,doPlot=TRUE,n=50)
clone_track(isa4004,33,doPlot=TRUE,n=50)
clone_track(isa4004,47,doPlot=TRUE,n=50)
clone_track(isa4004,56,doPlot=TRUE,n=50)
#
# B cells
#
clone_track(isa4004,20,doPlot=TRUE,n=50)
clone_track(isa4004,29,doPlot=TRUE,n=50)
clone_track(isa4004,43,doPlot=TRUE,n=50)
clone_track(isa4004,52,doPlot=TRUE,n=50)
#
# mono
#
clone_track(isa4004,7,doPlot=TRUE,n=50)
clone_track(isa4004,11,doPlot=TRUE,n=50)
clone_track(isa4004,22,doPlot=TRUE,n=50)
clone_track(isa4004,31,doPlot=TRUE,n=50)
clone_track(isa4004,45,doPlot=TRUE,n=50)
clone_track(isa4004,54,doPlot=TRUE,n=50)
#
# gran
#
clone_track(isa4004,3,doPlot=TRUE,n=50)
clone_track(isa4004,5,doPlot=TRUE,n=50)
clone_track(isa4004,9,doPlot=TRUE,n=50)
clone_track(isa4004,21,doPlot=TRUE,n=50)
clone_track(isa4004,30,doPlot=TRUE,n=50)
clone_track(isa4004,44,doPlot=TRUE,n=50)
clone_track(isa4004,53,doPlot=TRUE,n=50)
#-----------------------------------------------------------------
#  Z13264
#
#
# global max
#
clone_track(isa13264,0,doPlot=TRUE,n=50)
#
# NK
#
clone_track(isa13264,10,doPlot=TRUE,n=50)
clone_track(isa13264,18,doPlot=TRUE,n=50)
clone_track(isa13264,26,doPlot=TRUE,n=50)
clone_track(isa13264,35,doPlot=TRUE,n=50)
clone_track(isa13264,45,doPlot=TRUE,n=50)
#
# mono
#
clone_track(isa13264,9,doPlot=TRUE,n=50)
clone_track(isa13264,17,doPlot=TRUE,n=50)
clone_track(isa13264,25,doPlot=TRUE,n=50)
clone_track(isa13264,34,doPlot=TRUE,n=50)
clone_track(isa13264,44,doPlot=TRUE,n=50)
#
# T cells
#
clone_track(isa13264,11,doPlot=TRUE,n=50)
clone_track(isa13264,19,doPlot=TRUE,n=50)
clone_track(isa13264,27,doPlot=TRUE,n=50)
clone_track(isa13264,36,doPlot=TRUE,n=50)
clone_track(isa13264,46,doPlot=TRUE,n=50)
#
# B cells
#
clone_track(isa13264,7,doPlot=TRUE,n=50)
clone_track(isa13264,15,doPlot=TRUE,n=50)
clone_track(isa13264,23,doPlot=TRUE,n=50)
clone_track(isa13264,32,doPlot=TRUE,n=50)
clone_track(isa13264,42,doPlot=TRUE,n=50)
#
# grans
#
clone_track(isa13264,8,doPlot=TRUE,n=50)
clone_track(isa13264,16,doPlot=TRUE,n=50)
clone_track(isa13264,24,doPlot=TRUE,n=50)
clone_track(isa13264,33,doPlot=TRUE,n=50)
clone_track(isa13264,43,doPlot=TRUE,n=50)

clone_track(isaJ02370,0,doPlot=FALSE,n=50)
clone_track(isaJ02370,9,doPlot=FALSE,n=50)
clone_track(isaJ02370,15,doPlot=FALSE,n=50)


