#******************************************************************************************
#
#  Clones that are present in plasmid & virus & Z01932 & Z08103
#
#
#   
#   M Enstrom 11-28-18
#
#******************************************************************************************

#
# load plasmid and virus data
#
nameP = "/Users/mark_enstrom/Barcode/plasmid/plasmid.txt"
nameV = "/Users/mark_enstrom/Barcode/plasmid/virus.txt"
nameV_k562 = "/Users/mark_enstrom/Barcode/virus/virus_k562.txt"
nameV_prep = "/Users/mark_enstrom/Barcode/virus/virus_prep.txt"




plasmid <- read.csv(nameP,stringsAsFactors = FALSE,comment.char = '#',header = F)
rownames(plasmid) = plasmid[,1]
plasmid[,1] = c(1:nrow(plasmid))
plasmid[,3] = log(plasmid[,2])
plasmid[,4] = "plasmid"
colnames(plasmid) <- c('rank','count','l2count','source')
head(plasmid)

virusK562 <- read.csv(nameV_k562,stringsAsFactors = FALSE,comment.char = '#',header = F)
rownames(virusK562) = virusK562[,1]
virusK562[,1] = c(1:nrow(virusK562))
virusK562[,3] = log(virusK562[,2])
virusK562[,4] = "virus_on_K562"
colnames(virusK562) <- c('rank','count','l2count','source')
head(virusK562)

virusPrep <- read.csv(nameV_prep,stringsAsFactors = FALSE,comment.char = '#',header = F)
rownames(virusPrep) = virusPrep[,1]
virusPrep[,1] = c(1:nrow(virusPrep))
virusPrep[,3] = log(virusPrep[,2])
virusPrep[,4] = "virus_prep"
colnames(virusPrep) <- c('rank','count','l2count','source')
head(virusPrep)
#
# what a distribution might look like
#
head(virusK562)
probK562 = virusK562[,2]/sum(virusK562[,2])
head(probK562)
tail(probK562)
#
#
#
plas_bc = rownames(plasmid)
virk562_bc = rownames(virusK562)
virusPrep_bc = rownames(virusPrep)
names(plas_bc) = plas_bc
names(virk562_bc) = virk562_bc
names(virusPrep_bc) = virusPrep_bc


a = intersect(plas_bc,virk562_bc)
length(a)

a = intersect(a,virusPrep_bc)
length(a)
#
# all Z09132 BC   !!! very easy cummulative clone count
#
z9_clones = c()
l = c()
for (i in seq_along(bc9132@labels)) {
  n = bc9132@data[[i]][,1]
  z9_clones = unique(c(z9_clones,n))
  print(length(z9_clones))
  l = c(l,length(z9_clones))
}
plot(l,type='l')

z8_clones = c()
l = c()
for (i in seq_along(bc8103@labels)) {
  n = bc8103@data[[i]][,1]
  z8_clones = unique(c(z8_clones,n))
  print(length(z8_clones))
  l = c(l,length(z8_clones))
}
plot(l,type='l')

b = intersect(z8_clones,z9_clones)
length(b)
x = intersect(bc8103@global[,1],bc9132@global[,1])
length(x)

y = unique(b,x)
length(y)


target = bc8103@global[,1]
head(target)

dispCommonCodes <- function(target) 
{
  print("                  Barcode   plasmid virusPrep virusK562  Z09132   Z08103")
  i = 0
  for (id in target) {
    e9 = bc9132@global[id,2]
    e8 = bc8103@global[id,2]
    ep = which(plas_bc == id)
    if (length(ep) == 0) { ep = NA}
    ev = which(virusPrep_bc== id)
    if (length(ev) == 0) { ev = NA}
    evK = which(virk562_bc== id)
    if (length(evK) == 0) { evK = NA}
    
    #if ((e9 < 100) | (e8 < 100) | (ep < 100) | (ev < 100)){
    print(sprintf("%4d-%s   %7d   %7d    %7d %6d   %6d",i,id,ep,ev,evK,e9,e8))
    #}
    i = i + 1
    if (i > 100) {
      break()
    }
  }
}

target = bc8103@global[,1]
head(target)
dispCommonCodes(target)

target = bc9132@global[,1]
head(target)
dispCommonCodes(target)


target = plas_bc
head(target)
dispCommonCodes(target)

target = virusPrep_bc
head(target)
dispCommonCodes(target)


target = virk562_bc
head(target)
dispCommonCodes(target)



omni_clones = intersect(vir_bc,plas_bc)
print(length(omni_clones))
omni_clones = intersect(omni_clones,z9_clones)
print(length(omni_clones))
omni_clones = intersect(omni_clones,z8_clones)
print(length(omni_clones))

print("Barcode   plasmid virus Z09132   Z08103")
for (i in seq_along(omni_clones)) {
  id = omni_clones[i]
  e9 = bc9132@global[id,2]
  e8 = bc8103@global[id,2]
  ep = which(plas_bc == id)
  ev = which(vir_bc == id)
  #if ((e9 < 1000) | (e8 < 1000) | (ep < 100) | (ev < 100)){
    print(sprintf("%4d-%s   %7d   %7d    %6d   %6d",i,id,ep,ev,e9,e8))
  #}
  #break()
}


timestamp()
sd = vapply(seq_along(omni_clones),FUN =  function(i){
  id = omni_clones[i]
  e9 = bc9132@global[id,2]
  e8 = bc8103@global[id,2]
  ep = plas_bc[id]
  ev = vir_bc[id]
  if ((i %% 100) == 0) {print(i)}
  return(c(id,ep,ev,e9,e8))  
}, FUN.VALUE = c(character(1),numeric(1),numeric(1),numeric(1),numeric(1)))
timestamp()

tsd = t(sd)
c1 = (tsd[,1])
c2 = (as.numeric(tsd[,2]))
c3 = (as.numeric(tsd[,3]))
c4 = (as.numeric(tsd[,4]))
c5 = (as.numeric(tsd[,5]))

df = data.frame(c2,c3,c4,c5)
colnames(df) = c("plasmid","virus","Z09132","Z08103")
rownames(df) = c1
head(df)

write.table(df,file="/Users/mark_enstrom/Barcode/FinalBC_RIS/omniClone.txt",sep = '\t')

head(df)
tail(df)

which(plas_bc == "GGGCGGTGTCTTGGCCCAGA")
which(vir_bc == "GGGCGGTGTCTTGGCCCAGA")




for (i in seq_along(omni_clones)) {
  id = omni_clones[i]
  e9 = bc9132@global[id,2]
  e8 = bc8103@global[id,2]
  ep = plas_bc[id]
  ev = vir_bc[id]
  if ((i %% 100) == 0) {print(i)}
  print(sprintf("%s   %6d %6d %6d %6d",id,ep,ev,e9,e8)) 
  break()
}

fred = as.data.frame(t(sd))
dim(fred)
head(fred)
class(fred[1,2])



barf = vapply(seq_along(omni_clones),FUN =  function(i){
  id = omni_clones[i]
  e9 = bc9132@global[id,2]
  e8 = bc8103@global[id,2]
  ep = plas_bc[id]
  ev = vir_bc[id]
  if ((i %% 100) == 0) {print(i)}
  return(c(id,ep,ev,e9,e8))  
}, FUN.VALUE = c(character(1),numeric(1),numeric(1),numeric(1),numeric(1)))
