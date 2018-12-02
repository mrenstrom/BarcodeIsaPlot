#******************************************************************************************
#
#   Load and plot Hamming distance data
#
#   
#   M Enstrom 11-21-2018
#
#******************************************************************************************
#
# init clone environment
#
source("InitRisBC.R")
#******************************************************************************************
#
# plasmid hamming
#
#******************************************************************************************
ham_plas1 = read.table("/Users/mark_enstrom/Barcode/plasmid/Plasmid_Library_NA_SNP490_ACTG_NA_mp_global_ham.txt",sep=',')
ham_plas1[,3] = log(ham_plas1[,2])
colnames(ham_plas1) = c("distance","count","l2count")
head(ham_plas1,16)

pdf(paste0(plotDir,"Hamming_plasmid.pdf"))

g <- ggplot(ham_plas1,aes(x=distance,y=l2count))
g + geom_line(size=1.2) + 
  theme_bw() +
  theme(legend.position = c(0.8,0.75),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        axis.title = element_text(size=14),
        plot.title = element_text(size=16),
        legend.title = element_blank()) +
  labs(x="Hamming Distance",y="log(count)",title="Plasmid Top 100 Barcode hamming distances")# + coord_fixed(ratio=300)
dev.off()
#******************************************************************************************
#
# virus hamming
#
#******************************************************************************************
ham_virus1 = read.table("/Users/mark_enstrom/Barcode/virus/Virus_on_K562s_NA_SNP570__Virus_mp_global_ham.txt",sep=',')
ham_virus1[,3] = log(ham_virus1[,2])
colnames(ham_virus1) = c("distance","count","l2count")
head(ham_virus1)
pdf(paste0(plotDir,"Hamming_virus.pdf"))

g <- ggplot(ham_virus1,aes(x=distance,y=l2count))
g + geom_line(size=1.2) + 
  theme_bw() +
  theme(legend.position = c(0.8,0.75),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        axis.title = element_text(size=14),
        plot.title = element_text(size=16),
        legend.title = element_blank()) +
  labs(x="Hamming Distance",y="log(count)",title="Virus Top 100 Barcode hamming distances")# + coord_fixed(ratio=300)
dev.off()
#******************************************************************************************
#
# Z08103 644 grancombined hamming
#
#******************************************************************************************
ham_644 = read.table("/Users/mark_enstrom/Barcode/Z08103/hamming/644_KT21_CA_PB_Granulocytes_mp_global_ham.txt",sep=',')
ham_644[,3] = log(ham_644[,2])
colnames(ham_644) = c("distance","count","l2count")
head(ham_644)

pdf(paste0(plotDir,"Hamming_Z08103_644_Gran.pdf"))

g <- ggplot(ham_644,aes(x=distance,y=l2count))
g + geom_line(size=1.2) + 
  theme_bw() +
  theme(legend.position = c(0.8,0.75),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        axis.title = element_text(size=14),
        plot.title = element_text(size=16),
        legend.title = element_blank()) +
  labs(x="Hamming Distance",y="log(count)",title="Z08103 644 DPT Gran Top 100 Barcode hamming distances")# + coord_fixed(ratio=300)

dev.off()
#******************************************************************************************
#
#  Z09132 combined 118 gfp
#
#******************************************************************************************
ham_118 = read.table("/Users/mark_enstrom/Barcode/Z09132/hamming/118_SNP555__GFP+_PB_mp_global_ham.txt",sep=',')
ham_118[,3] = log(ham_118[,2])
colnames(ham_118) = c("distance","count","l2count")
head(ham_118)

pdf(paste0(plotDir,"Hamming_Z09132_118_GFP.pdf"))


g <- ggplot(ham_118,aes(x=distance,y=l2count))
g + geom_line(size=1.2) + 
  theme_bw() +
  theme(legend.position = c(0.8,0.75),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        axis.title = element_text(size=14),
        plot.title = element_text(size=16),
        legend.title = element_blank()) +
  labs(x="Hamming Distance",y="log(count)",title="Z09132 118 DPT GFP Top 100 Barcode hamming distances")# + coord_fixed(ratio=300)
dev.off()

#------------------------------------------------------------------------------------------------------------
#
# Z09132 118 GFP top 100
#
#
#
#
#------------------------------------------------------------------------------------------------------------
h118 = read.table("/Users/mark_enstrom/Barcode/Z09132/hamming/118_SNP555__GFP+_PB_mp_ham.txt",sep=',')
head(h118)
h118[,4] = log(h118[,3])
head(h118)
colnames(h118) = c("clone","distance","count","l2count")
head(h118)
tail(h118)
t118 = h118[which(h118[,1]<100),]

#
# good example of aes group funcion
#
g <- ggplot(t118,aes(x=distance,y=l2count))

g + geom_line(aes(color=clone,group=clone)) +
  theme_bw() +
  scale_color_gradient(low="#202020",high="#d0d0d0") +
  guides(fill=FALSE) +
  labs(x="Hamming Distance",y="log(count)",title="Z09132 118 DPT GFP Top 100 Barcode hamming distances")

pdf(paste0(plotDir,"Hamming_Z09132_118_GFP_100.pdf"))

g + geom_line(aes(color=clone,group=clone)) +
  theme_bw() +
  scale_color_gradient(low="#202020",high="#d0d0d0") +
  guides(fill=FALSE) +
  labs(x="Hamming Distance",y="log(count)",title="Z09132 118 DPT GFP Top 100 Barcode hamming distances")

dev.off()

#------------------------------------------------------------------------------------------------------------
#
# Z08103 644 gran
#
#
#
#
#------------------------------------------------------------------------------------------------------------

h644 = read.table("/Users/mark_enstrom/Barcode/Z08103/hamming/644_KT21_CA_PB_Granulocytes_mp_ham.txt",sep=',')
head(h644)
h644[,4] = log(h644[,3])
head(h644)
colnames(h644) = c("clone","distance","count","l2count")
head(h644)
tail(h644)



t644 = h644[which(h644[,1]<100),]

#
# good example of aes group funcion
#
g <- ggplot(h644,aes(x=distance,y=l2count))

g + geom_line(aes(color=clone,group=clone)) +
  theme_bw() +
  scale_color_gradient(low="#202020",high="#d0d0d0") +
  guides(fill=FALSE) +
  labs(x="Hamming Distance",y="log(count)",title="Z08103 644 DPT Gran Top 100 Barcode hamming distances")

pdf(paste0(plotDir,"Hamming_Z08103_644_gran_100.pdf"))

g + geom_line(aes(color=clone,group=clone)) +
  theme_bw() +
  scale_color_gradient(low="#202020",high="#d0d0d0") +
  guides(fill=FALSE) +
  labs(x="Hamming Distance",y="log(count)",title="Z09132 118 DPT GFP Top 100 Barcode hamming distances")

dev.off()

