#******************************************************************************************
#
#  given a plasmid library of a certain size and a certain number of target cells:
#    how many cells will take up a duplicate barcode
#
#   
#   M Enstrom 11-28-18
#
#******************************************************************************************



#
# For animal Z08103, 5.76 × 106 cells were infused containing 53% GFP+ cells, 
# equivalent to 3.05 × 106 total possible clones. Assumue 1/20 repopulating cells
#
# For animal Z09132, a total of 8 × 106 cells were infused containing 29% GFP+ cells, 
# equivalent to 2.32 × 106 total possible clones. Assumue 1/20 repopulating cells
#
# We simulated the frequency of transduced repopulating cells with identical barcodes 
# by making the following assumptions: (1) the absolute library size was 1,213,714 barcodes
#
z8pop = 5.76e6 * .53 * (1/20)
z8pop
z9pop = 8.0e6 * .29 * (1/20)
z9pop
#libSize = 1213714
#libSize = 80000
libSize = length(probK562)
head(probK562)
#
# simulate each cell 
#
simDupBarcode <- function(lib,cells) {
  s = sapply(c(1:100),function(x) {
    a = sample(c(1:lib),cells,replace = T,prob = probK562)
    a = sort(a)
    a
    r = rle(a)
    r
    l2 = length(which(r$lengths >= 2))
    l3 = length(which(r$lengths >= 3))
    l4 = length(which(r$lengths >= 4))
    l5 = length(which(r$lengths >= 5))
    l6 = length(which(r$lengths >= 6))
    
    return(c(l2,l3,l4,l5,l6))  
  })
  s = t(s)
  head(s)
  #
  # s is a matrix , col = l2 l3 l4 l5
  #                 r
  ml2 = mean(s[,1])
  print(sprintf("Barcodes selected 2 times = %7.2f, f = %2.4f",ml2,ml2/cells))
  ml3 = mean(s[,2])
  print(sprintf("Barcodes selected 3 times = %7.2f, f = %2.4f",ml3,ml3/cells))
  ml4 = mean(s[,3])
  print(sprintf("Barcodes selected 4 times = %7.2f, f = %2.4f",ml4,ml4/cells))
  ml5 = mean(s[,4])
  print(sprintf("Barcodes selected 5 times = %7.2f, f = %2.4f",ml5,ml5/cells))
  ml6 = mean(s[,5])
  print(sprintf("Barcodes selected 6 times = %7.2f, f = %2.4f",ml6,ml6/cells))
}
#
# population
#
simBarcodePopulation <- function(lib,cells) {
  a = sample(c(1:lib),cells,replace = T)
  return(a)
}


simDupBarcode(libSize,z8pop)
simDupBarcode(libSize,z9pop)

#
# now for the same animals assume both are transfected with the same virus prep. How many
# of the same barcodes would transfect a cell in both animals
#
z8_bc_sel = simBarcodePopulation(libSize,z8pop)
head(z8_bc_sel)
length(z8_bc_sel)

z8_bc_sel = unique(z8_bc_sel)
length(z8_bc_sel)

z9_bc_sel = simBarcodePopulation(libSize,z9pop)
head(z9_bc_sel)
length(z9_bc_sel)
z9_bc_sel = unique(z9_bc_sel)
length(z9_bc_sel)
#
# how many barcode in both
#
bc_in_both = intersect(z8_bc_sel,z9_bc_sel)
length(bc_in_both)
#
# same simulation but randomly select transfected cells to account for how many
# clones from each animal we actually detected
#
z8_bc_sel = simBarcodePopulation(libSize,z8pop)
head(z8_bc_sel)
length(z8_bc_sel)
z8_bc_sel = sample(z8_bc_sel,41935+16000,replace = F)
length(z8_bc_sel)
z8_bc_sel = unique(z8_bc_sel)
length(z8_bc_sel)
#
# Z09132
#
z9_bc_sel = simBarcodePopulation(libSize,z9pop)
head(z9_bc_sel)
length(z9_bc_sel)
z9_bc_sel = sample(z9_bc_sel,19849+3000,replace = F)
length(z9_bc_sel)
z9_bc_sel = unique(z9_bc_sel)
length(z9_bc_sel)
#
# home many in both
#
bc_in_both = intersect(z8_bc_sel,z9_bc_sel)
length(bc_in_both)
