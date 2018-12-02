#******************************************************************************************
#
#   load plasmid and virus barcode counts
#
#   
#   M Enstrom 11-19-2018
#
#******************************************************************************************

nameP = "/Users/mark_enstrom/Barcode/plasmid/plasmid.txt"
nameV = "/Users/mark_enstrom/Barcode/plasmid/virus.txt"

plasmid <- read.csv(nameP,stringsAsFactors = FALSE,comment.char = '#',header = F)
rownames(plasmid) = plasmid[,1]
plasmid[,1] = c(1:nrow(plasmid))
plasmid[,3] = log(plasmid[,2])
plasmid[,4] = "Plasmid"
colnames(plasmid) <- c('barcode','count','l2count','source')
head(plasmid)
dim(plasmid)

virus <- read.csv(nameV,stringsAsFactors = FALSE,comment.char = '#',header = F)
rownames(virus) = virus[,1]
virus[,1] = c(1:nrow(virus))
virus[,3] = log(virus[,2])
virus[,4] = "Virus & K562"
colnames(virus) <- c('barcode','count','l2count','source')
head(virus)
tail(virus)
dim(virus)


df <- rbind(plasmid,virus)
dim(df)
head(df)
tail(df)

pdf(paste0(plotDir,"plasm_virus_Count.pdf"))

g <- ggplot(df,aes(x=barcode,y=l2count,color=source,linetype=source))
  g + geom_line(size=1.2) + 
  theme_bw() +
  theme(legend.position = c(0.8,0.75),
          axis.text=element_text(size=14),
          legend.text = element_text(size=14),
          axis.title = element_text(size=14),
          plot.title = element_text(size=16),
          legend.title = element_blank()) +
    
  scale_color_manual(values=c("#000000", "#c0c0c0")) + 
  scale_linetype_manual(values=c("solid","dotdash")) +
  labs(x="Barcode rank",y="log(Barcode Abundance)",title="Barcode Library")# + coord_fixed(ratio=300)
dev.off()  

#
# zoom in 
#
df <- rbind(plasmid[1:100,],virus[1:100,])
pdf(paste0(plotDir,"top_100_plasmid_v.pdf"))
g <- ggplot(df,aes(x=barcode,y=count,color=source,linetype=source))
g + geom_line() + 
  theme_bw() +
  theme(legend.position = c(0.8,0.75)) +
  scale_color_manual(values=c("#000000", "#c0c0c0")) + 
  scale_linetype_manual(values=c("solid","dotdash")) +
  labs(x="Barcode rank",y="Barcode abundance",title="Top 100 Plasmid and Virus Barcodes")# + coord_fixed(ratio=300)
dev.off()
#
# haw many plasmid top-100 barcodes are in virus
#
plasNames = rownames(plasmid)
plasNames100 = plasNames[1:100]
virusNames = rownames(virus)
virusNames100 = virusNames[1:100]


plasNames100 %in% plasNames
virusNames100 %in% virusNames
plasNames100 %in% virusNames
virusNames100 %in% plasNames

plasmidMatch = plasNames100[plasNames100 %in% virusNames]
plasmidMatch
virusMatch = virusNames100[virusNames100 %in% plasNames]
virusMatch
# 
# 2 top plasmid barcodes in all of virus lib
#
s = sapply(plasmidMatch, function(bc) {
  which(virusNames == bc)
})
#
#
#
unname(unlist(s))
#[1] 635493 782042
#
# 82 top virus BC in plasmid lib, what are the ranks?
#
s = sapply(virusMatch, function(bc) {
  which(plasNames == bc)
})

sort(unname(unlist(s)))
#
# ranks of the top 82 virus barcodes in the plasmid library
#
#[1]    7356   12354   15767   17529   22829   22977   26309   33103   48951   69538   72083   78888  104787  111375  112578  120987
#[17]  122671  128402  137057  142530  143382  148131  151299  166856  169100  173797  193851  228055  238215  239102  267066  279624
#[33]  309256  330109  333908  345814  347840  350420  359913  376541  388564  422741  431653  434557  437550  442120  454684  458954
#[49]  466618  471423  478477  488266  489451  493375  500807  506934  517039  539651  573563  588756  603688  663183  680063  689893
#[65]  693551  709090  728327  737650  745912  753256  761135  771506  777108  817228  834597  895531  914367  917160 1000457 1002030
#[81] 1007764 1008641
#
# how many barcodes for Z09132
#
length(bc9132@global[,1])
#
# 21399
#
#  what are the plasmid & virus ranks of the top 100 Z09132 barcodes?
#
length(which(bc9132@global[,1] %in% plasNames))
length(which(bc9132@global[,1] %in% virusNames))
#
#  16332/21399 and 12424/21399
#
allNames = unique(c(plasNames,virusNames))
length(allNames)
#
# total = 1218230
#
length(which(bc9132@global[,1] %in% allNames))
#
# 18590/21399 in plasmid or virus
#
#
# what are the ranks of the top 100 Z09132 barcodes in plasmid lib
#
z09top100 = bc9132@global[1:100,1]
length(which(z09top100 %in% plasNames))
length(which(z09top100 %in% virusNames))
# 91 & 98
s = sapply(z09top100, function(bc) {
  which(plasNames == bc)
})
sort(unname(unlist(s)))
#[1]   9957  22300  25190  36411  38350  47530  48368  70228  83181  83321  87442  88613  94185 105367 130308 146057 153308 157717
#[19] 157929 158931 158997 159056 181986 188742 191504 195548 198888 205327 211552 212017 216141 238949 253760 296798 300080 302695
#[37] 308026 322316 322602 340576 343674 353456 365340 384832 385972 399496 401383 417686 421578 421824 422182 430641 435525 439021
#[55] 440525 457327 465667 490903 492436 498300 531752 532909 537410 551311 567345 584596 603141 604003 632649 651044 679800 684907
#[73] 692999 702577 717915 747908 749972 757091 801740 826691 831968 880249 903591 911495 927048 933359 934032 936411 947895 953772
#[91] 974665
s = sapply(z09top100, function(bc) {
  which(virusNames == bc)
})
sort(unname(unlist(s)))
#[1]  12576  16160  26689  29602  31117  31759  40366  64177  64245  64942  68092  72451 100721 106363 115985 116635 119105 123932
#[19] 126402 129662 129677 130237 130355 132125 135387 137114 137135 138922 139216 144892 145052 147979 148009 148780 151384 151496
#[37] 152559 158694 166668 169536 169953 173313 174451 174882 175427 175657 177003 181958 185088 187817 195462 195839 202549 208248
#[55] 211293 221358 226021 228712 234315 237061 242956 245150 247706 249255 252811 254435 255765 259669 266750 271389 277630 284902
#[73] 285919 288629 312545 313822 316371 340041 361188 362821 363313 380929 399927 409942 423331 424014 479972 487445 493175 528544
#[91] 570264 605113 612548 637070 703023 758902 798535 813564
z08top100 = bc8103@global[1:100,1]
length(which(z08top100 %in% plasNames))
length(which(z08top100 %in% virusNames))
# 84 & 96
s = sapply(z08top100, function(bc) {
  which(plasNames == bc)
})
sort(unname(unlist(s)))
#[1]    9957   13636   23336   30745   35469   38350   48368   52175   57617   72652   76226   95901  106888  112578  114906  120987
#[17]  128854  153308  154770  176686  181911  181978  183313  198888  236181  238864  238949  258959  260711  262597  264205  279758
#[33]  280081  290005  298692  353447  365003  366652  371884  374739  379057  384832  404510  407794  409955  430855  439021  465667
#[49]  466990  478477  480976  493375  507090  512419  530393  547896  552458  581246  583975  585424  606899  617891  625579  669961
#[65]  684688  684919  712071  717915  722702  740346  765993  766204  775143  776033  777844  801740  807147  807991  812256  829547
#[81]  856380  861706  876163 1002030

s = sapply(z08top100, function(bc) {
  which(virusNames == bc)
})
sort(unname(unlist(s)))
#
#[1]     10     26     35     39     40     70    107    123    204    369    415    575    611    694    791    831    869   1205
#[19]   1211   1231   1962   2201   2779   3842   3925   4278   4899   5495   5847   6764   6887   7372  10287  12043  12576  12665
#[37]  12765  14116  16839  19968  26637  26689  26779  30367  31117  31758  39108  40366  45789  49719  52089  62846  64177  65253
#[55]  72451  82596 111085 117978 118969 126173 131040 133773 137596 138840 141356 148009 150409 161763 190894 203644 206448 216973
#[73] 221358 230823 234192 243200 252705 255765 271389 287957 288629 360021 363030 380929 387221 406919 470152 523222 542322 556133
#[91] 600830 613624 652200 733474 835545 836349
#
s
#
# combined
#
length(plasNames)
length(virusNames)
c_lib = unique(c(plasNames,virusNames))
length(c_lib)
#
# singles
# 
singles_9 = bc9132@global[which(bc9132@global[,3]== 1),1]
length(singles_9)


match_9 = singles_9[which(singles_9 %in% c_lib)]
length(match_9)

singles_8 = bc8103@global[which(bc8103@global[,3]== 1),1]
length(singles_8)

match_8 = singles_8[which(singles_8 %in% c_lib)]
length(match_8)



singles_8 %in% c_lib


singles_8


