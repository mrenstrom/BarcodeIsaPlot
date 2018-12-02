#******************************************************************************************
#
#  build some color maps
#
#   
#   M Enstrom 4-11-2018
#
#******************************************************************************************
buildColorMap <- function() {
  rStart = 0.0
  rFinal = 0.0
  gStart = 0.0
  gFinal = 0.88
  bStart = 0.5
  bFinal = 0.94
  r = c()
  g = c()
  b = c()
  l = 399
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  rStart = 0.0
  rFinal = 1.00
  gStart = 0.88
  gFinal = 0.86
  bStart = 0.94
  bFinal = 0.0
  l = 250
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  rStart = 1.0
  rFinal = 0.94
  gStart = 0.86
  gFinal = 0.0
  bStart = 0.0
  bFinal = 0.0
  l = 250
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  rStart = 0.94
  rFinal = 0.56
  gStart = 0.0
  gFinal = 0.0
  bStart = 0.0
  bFinal = 0.0
  l = 100
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  
  
  r = round(r * 255)
  g = round(g * 255)
  b = round(b * 255)
  r
  g
  b
  
  
  
  colArray =  c()
  for (i in c(1:length(r))) {
    ri = r[i]
    gi = g[i]
    bi = b[i]
    strr = sprintf("%02x",ri)
    strg = sprintf("%02x",gi)
    strb = sprintf("%02x",bi)
    strJ1 = paste0("#",strr,strg,strb)
    colArray = c(colArray,strJ1)
  }
  #colArray = colArray[1:1000]
  length(colArray)
  colArray[1] = "#000000"
  #colArray[1000]
  return(colArray)
}
#******************************************************************************************
#
#  build index color map
#
#   
#   M Enstrom 4-11-2018
#
#******************************************************************************************
buildIndexColorMap <- function() {
  r = c()
  g = c()
  b = c()
  #
  #
  #
  rStart = 1.00
  rFinal = 0.88
  gStart = 1.00
  gFinal = 0.88
  bStart = 1.00
  bFinal = 0.88
  l = 100
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  rStart = 0.88
  rFinal = 0.88
  gStart = 0.88
  gFinal = 0.88
  bStart = 0.88
  bFinal = 0.00
  l = 100
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  rStart = 0.88
  rFinal = 0.88
  gStart = 0.88
  gFinal = 0.00
  bStart = 0.00
  bFinal = 0.00
  l = 300
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  rStart = 0.88
  rFinal = 0.33
  gStart = 0.00
  gFinal = 0.00
  bStart = 0.00
  bFinal = 0.00
  l = 300
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  r = round(r * 255)
  g = round(g * 255)
  b = round(b * 255)
  colArray =  c()
  for (i in c(1:length(r))) {
    ri = r[i]
    gi = g[i]
    bi = b[i]
    strr = sprintf("%02x",ri)
    strg = sprintf("%02x",gi)
    strb = sprintf("%02x",bi)
    strJ1 = paste0("#",strr,strg,strb)
    colArray = c(colArray,strJ1)
  }
  
  colArray = c(colArray,"#008888")
  return(colArray)
}

#******************************************************************************************
#
#  build index color map
#
#   
#   M Enstrom 4-11-2018
#
#******************************************************************************************
buildIndexColorMap2 <- function() {
  r = c()
  g = c()
  b = c()
  #
  #
  #
  rStart = 1.00
  rFinal = 0.88
  gStart = 1.00
  gFinal = 0.88
  bStart = 1.00
  bFinal = 0.88
  l = 100
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  rStart = 0.88
  rFinal = 0.88
  gStart = 0.88
  gFinal = 0.88
  bStart = 0.88
  bFinal = 0.00
  l = 100
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  rStart = 0.88
  rFinal = 0.88
  gStart = 0.88
  gFinal = 0.00
  bStart = 0.00
  bFinal = 0.00
  l = 300
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  rStart = 0.88
  rFinal = 0.50
  gStart = 0.00
  gFinal = 0.00
  bStart = 0.00
  bFinal = 0.00
  l = 300
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  r = round(r * 255)
  g = round(g * 255)
  b = round(b * 255)
  colArray =  c()
  for (i in c(1:length(r))) {
    ri = r[i]
    gi = g[i]
    bi = b[i]
    strr = sprintf("%02x",ri)
    strg = sprintf("%02x",gi)
    strb = sprintf("%02x",bi)
    strJ1 = paste0("#",strr,strg,strb)
    colArray = c(colArray,strJ1)
  }
  
  return(colArray)
}
#
# grey-yellow-red
#
grey_yellow_red_map <- function() {
  cmap = read.csv("/Volumes/SamLab/LabData/ColorMap.csv",stringsAsFactors = FALSE,quote = "")
  cmap = cmap[,2]
  cmap[1]   = "#000000"
  cmap[90]  = "#A04040"
  cmap[91]  = "#A05050"
  cmap[92]  = "#A06060"
  cmap[93]  = "#A07070"
  cmap[94]  = "#B08080"
  cmap[95]  = "#C09090"
  cmap[96]  = "#D09090"
  cmap[97]  = "#E09090"
  cmap[98]  = "#F09090"
  cmap[99]  = "#FFA0A0"
  cmap[100] = "#FFB0B0"
  
  cmap2 = cmap[1:80]
  return(cmap2)
}

buildHeatColorMap <- function() {
  rStart = 0.1
  rFinal = 1.0
  gStart = 0.0
  gFinal = 0.0
  bStart = 0.8
  bFinal = 0.0
  r = c()
  g = c()
  b = c()
  l = 100
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  
  rStart = 1.00
  rFinal = 1.00
  gStart = 0.1
  gFinal = 1.0
  bStart = 0.1
  bFinal = 0.5
  l = 400
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  
  rStart = 1.00
  rFinal = 1.00
  gStart = 1.0
  gFinal = 1.0
  bStart = 0.5
  bFinal = 1.0
  l = 100
  for (i in c(0:l)) {
    r = c(r,rStart+((i/l) * (rFinal-rStart)))
    g = c(g,gStart+((i/l) * (gFinal-gStart)))
    b = c(b,bStart+((i/l) * (bFinal-bStart)))
  }
  
  
  
  
  r = round(r * 255)
  g = round(g * 255)
  b = round(b * 255)
  
  colArray =  c()
  for (i in c(1:length(r))) {
    ri = r[i]
    gi = g[i]
    bi = b[i]
    strr = sprintf("%02x",ri)
    strg = sprintf("%02x",gi)
    strb = sprintf("%02x",bi)
    strJ1 = paste0("#",strr,strg,strb)
    colArray = c(colArray,strJ1)
  }
  colArray[1] = "#000000"
  return(colArray)
}
#------------------------------------------------------------------------------------------
#
# intensity alternating random color map
#
# n must be even
#
#------------------------------------------------------------------------------------------
randColorMapi <- function(n) {
  r = round(runif(n) * 255)
  g = round(runif(n) * 255)
  b = round(runif(n) * 255)
  i = r + g + b
  ii = order(i)
  
  indexArray = c()
  for (i in c(1:n/2)) {
    index = ii[i]
    indexArray = c(indexArray,index)
    index = ii[n-i+1]
    indexArray = c(indexArray,index)
  }
  
  colArray =  c()
  for (index in c(1:n)) {
    i = indexArray[index]
    ri = r[i]
    gi = g[i]
    bi = b[i]
    strr = sprintf("%02x",ri)
    strg = sprintf("%02x",gi)
    strb = sprintf("%02x",bi)
    strJ1 = paste0("#",strr,strg,strb)
    colArray = c(colArray,strJ1)
  }
  return(colArray)
}
#------------------------------------------------------------------------------------------
#
# random color map
#
#
#------------------------------------------------------------------------------------------
buildRandomColorMap <- function(obj) {
  n = nrow(obj@global)
  #
  # assign a random color to every entry in global barcode/refseq
  #
  r = round(runif(n) * 128)
  g = round(runif(n) * 128)
  b = round(runif(n) * 128)
  
  colArray =  c()
  for (i in c(1:n)) {
    ri = r[i]
    gi = g[i]
    bi = b[i]
    strr = sprintf("%02x",ri)
    strg = sprintf("%02x",gi)
    strb = sprintf("%02x",bi)
    strJ1 = paste0("#",strr,strg,strb)
    colArray = c(colArray,strJ1)
  }
  #
  # set top 10 global entries to bright colors
  #
  if (n >= 10) {
    colArray[1] = "#ff0000"
    colArray[2] = "#00ff00"
    colArray[3] = "#0000ff"
    colArray[4] = "#ff00ff"
    colArray[5] = "#ffff00"
    colArray[6] = "#00ffff"
    colArray[7] = "#ff80ff"
    colArray[8] = "#80ff80"
    colArray[9] = "#ff8D94"
    colArray[10]= "#f28218"
  }
  #
  # name colorMap by barcode
  #
  barcodes   = obj@global[,1]
  a <- matrix(colArray)
  names(a) <- barcodes
  #
  # change colors for barcodes found only at one datapoint to white
  #
  #or (i in c(1:nrow(obj@global))) {
  #  if (obj@global[i,6] == 1) {
  #    a[i] = "#ffffff"    
  #  }
  #}
  return(a)
}
