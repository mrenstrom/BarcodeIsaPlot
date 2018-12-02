#-----------------------------------------------------------------
#
# HSC clone signature
#
#
#
# Mark Enstorm 
# 11-14-18
#
#-----------------------------------------------------------------


#-----------------------------------------------------------------
# 
# findHSC
#
# clones present in any condition in group 1 AND
# also present in groupt 2
#
#-----------------------------------------------------------------
findHSC <- function(bcobj,group1,group2) 
{
  #
  #
  clones1 = c()
  for (i in c(1:length(group1))) {
    ds = bcobj@data[[group1[i]]]
    for (j in c(1:nrow(ds))) {
      newClone = ds[j,1]
      if (newClone %in% clones1) {
        #print("already there")
      } else {
        #print("new")
        clones1 = c(clones1,newClone)
      }
    }
  }
  
  clones2 = c()
  for (i in c(1:length(group2))) {
    ds = bcobj@data[[group2[i]]]
    for (j in c(1:nrow(ds))) {
      newClone = ds[j,1]
      if (newClone %in% clones2) {
        #print("already there")
      } else {
        #print("new")
        clones2 = c(clones2,newClone)
      }
    }
  }
  
  both = intersect(clones1,clones2)
  message(length(both))
  return (both)
}

#-----------------------------------------------------------------
#
#  getTotalExpression
#
#
#
#-----------------------------------------------------------------
getTotalExpression <- function(bcobj, conditions, clones)
{
  #bcobj = bc9132
  #conditions = c(1,2,4,16,17)
  #clones = l
  exp = c()
  ld = lapply(clones, function(clone) {
    #print(clone)
    lt = lapply(conditions, function(x) {
      if (clone %in% bcobj@data[[x]][,1]) {
        d = bcobj@data[[x]][clone,4]
      } else {
        d = 0
      }
    })
    
    #print(unlist(lt))
    return(sum(unlist(lt)))
  })
  
  return(unlist(ld))
}

#
# HSC signature for this function = clone seen in short term pop(gran)
# and and other sorted cell type
#
# pick out all grans for group1
#
# pick out all non-gran sorted cell types for group2
#
bc9132@labels


l1 = findHSC(bcobj = bc9132, group1 = c(7), group2 = c(18,19,26))
l2 = findHSC(bcobj = bc9132, group1 = c(21), group2 = c(5,6,26))
hsc9 = unique(c(l1,l2))
Z09132_dbs_hsc_clones = length(hsc9)




#
# now rate these clones on extression on DPT 32,104,256,400,469
#
e = getTotalExpression(bcobj = bc9132,conditions = c(1,2,4,16,17),clones=hsc9)
head(e)
head(hsc9)
names(e) = hsc9
se = sort(e,decreasing = T)
head(se)
#
# take etop 50%
#
#clones = se[1:(length(se)/2)]
clones = se
clones = names(clones)
head(clones)
#bc9132@data[[1]][clones[1:4],4]
#
# now plot
#
g <- buildCloneBarplot(bc9132,c(1,2,4,16,17),minF = 0.01, fixedClones = clones,plotType = "Area",bkColor = "#ffffff")
pdf(paste0(plotDir,"Z09132_DBS_HSC.pdf"))
g + 
  theme_bw() +
  labs(x= "Days post transplant",y = "Frequency of HSC clones") +
  scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
  ggtitle(paste0("Z09132 DBS (",Z09132_dbs_hsc_clones,  " HSC clones)")) + 
  theme(plot.title=element_text(size='16',hjust = 0.5),
        axis.text=element_text(size=14),
        axis.title = element_text(size=14),
        legend.title = element_text(size=16))
dev.off()
#
# verify > 90%
#
s = bc9132@data[[1]][clones,5]
length(clones)
length(s)
sum(is.na(s))
s[is.na(s)] = 0
sum(s)
#
#
#
#
# ISA 9132
#
#
#
#
isa9132@labels
#
# grna/ sort
#
l1 = findHSC(bcobj = isa9132, group1 = c(9), group2 = c(13,14,15,17,18,19,21,22,23))
l2 = findHSC(bcobj = isa9132, group1 = c(16), group2 = c(17,18,19,21,22,23))
l3 = findHSC(bcobj = isa9132, group1 = c(20), group2 = c(13,14,15,21,22,23))
isa9Clones = unique(c(l1,l2,l3))
z09132_isa_hsc_clones = length(isa9Clones)
length(isa9Clones)
head(isa9Clones)
# for 32,104,256,400,469
isa9Exp = getTotalExpression(bcobj = isa9132,conditions = c(1,4,8,11,12),clones=isa9Clones)
head(isa9Exp)
length(isa9Exp)
names(isa9Exp) = isa9Clones
se = sort(isa9Exp,decreasing = T)

clones = se
clones = names(clones)
head(clones)
isa9132@data[[4]][clones[1:100],4]
g <- buildCloneBarplot(isa9132,c(1,4,8,11,12),minF = 0.01, fixedClones = clones,plotType = "Area",bkColor = "#ffffff")
pdf(paste0(plotDir,"Z09132_ISA_HSC.pdf"))
g + 
  theme_bw() +
  labs(x= "Days post transplant",y = "Frequency of HSC clones") +
  scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
  ggtitle(paste0("Z09132 ISA ( ",z09132_isa_hsc_clones," HSC clones)")) + 
  theme(plot.title=element_text(size='16',hjust = 0.5),
        axis.text=element_text(size=14),
        axis.title = element_text(size=14),
        legend.title = element_text(size=16))
dev.off()
#
#   BC8103
#
#
#
bc8103@labels
l1 = findHSC(bcobj = bc8103, group1 = c(9), group2 = c(21,22,23,29))
l2 = findHSC(bcobj = bc8103, group1 = c(24), group2 = c(7,8,29))
l = unique(c(l1,l2))
z08103_dbs_hsc_clones = length(l)
#
# now rate these clones on extression on DPT 32,104,256,400,469
#
e = getTotalExpression(bcobj = bc8103,conditions = c(1,4,6,18,20),clones=l)
names(e) = l
se = sort(e,decreasing = T)
head(se)
length(se)
#
# take etop 50%
#
clones = se[1:(length(se)/2)]
clones = se
clones = names(clones)
head(clones)
bc8103@data[[1]][clones[1:4],4]
#
# now plot
#
g <- buildCloneBarplot(bc8103,c(1,4,6,18,20),minF = 0.01, fixedClones = clones,plotType = "Area",bkColor = "#ffffff")
pdf(paste0(plotDir,"Z08103_DBS_HSC.pdf"))
g + 
  theme_bw() +
  labs(x= "Days post transplant",y = "Frequency of HSC clones") +
  scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
  ggtitle(paste0("Z08103 DBS (",z08103_dbs_hsc_clones," HSC clones)")) + 
  theme(plot.title=element_text(size='16',hjust = 0.5),
        axis.text=element_text(size=14),
        axis.title = element_text(size=14),
        legend.title = element_text(size=16))
dev.off()
#
#
# ISA 8103
#
#
#
#
#
isa8103@labels
#
# grna/ sort
#
l1 = findHSC(bcobj = isa8103, group1 = c(8), group2 = c(12,13,14,16,17,18,20,21,22))
l2 = findHSC(bcobj = isa8103, group1 = c(15), group2 = c(6,7,16,17,18,20,21,22))
l3 = findHSC(bcobj = isa8103, group1 = c(19), group2 = c(6,7,12,13,14,20,21,22))
l4 = findHSC(bcobj = isa8103, group1 = c(23), group2 = c(6,7,12,13,14,16,17,18))
isa8Clones = unique(c(l1,l2,l3,l4))
z08103_isa_hsc_clones = length(isa8Clones)
head(isa8Clones)
# for 32,104,256,400,469
isa8Exp = getTotalExpression(bcobj = isa8103,conditions = c(1,4,5,9,11),clones=isa8Clones)
head(isa8Exp)
length(isa8Exp)
names(isa8Exp) = isa8Clones
se = sort(isa8Exp,decreasing = T)

clones = se
clones = names(clones)
head(clones)
isa8103@data[[4]][clones[1:100],4]
g <- buildCloneBarplot(isa8103,c(1,4,5,9,11),minF = 0.01, fixedClones = clones,plotType = "Area",bkColor = "#ffffff")
pdf(paste0(plotDir,"Z08103_ISA_HSC.pdf"))
g + 
  theme_bw() +
  labs(x= "Days post transplant",y = "Frequency of HSC clones") +
  scale_y_continuous(breaks=c(0,.25,.5,.75,1.0),labels=c("0%","25%","50%","75%","100%")) +
  ggtitle(paste0("Z08103 ISA (",z08103_isa_hsc_clones," HSC clones)")) + 
  theme(plot.title=element_text(size='16',hjust = 0.5),
        axis.text=element_text(size=14),
        axis.title = element_text(size=14),
        legend.title = element_text(size=16))
dev.off()

