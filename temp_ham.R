#
# read in hamming file
#
#
ham_118 = read.table("/Users/mark_enstrom/Barcode/Z09132/hamming/118_SNP555__GFP+_PB_mp_ham.txt",sep = ',')
v0 =ham_118[which(ham_118[,1] == 0),c(2,3)]
v0[,1] = as.numeric(v0[,1])
v0[,2] = as.numeric(v0[,2])



dist = v0[,2]

t = table(dist)
t
dist = sort(dist)
r = rle(dist)
r
plot(log(r$lengths),type='l')


ham_644 = read.table("/Users/mark_enstrom/Barcode/Z08103/hamming/644_KT21_CA_PB_Granulocytes_mp_ham.txt",sep=',')
v0 =ham_644[which(ham_644[,1] == 0),c(2,3)]
v0[,1] = as.numeric(v0[,1])
v0[,2] = as.numeric(v0[,2])
head(v0)
r = rle(sort(v0[,2]))
r
counts = rep(1,20)
for (i in seq_along(r$values)) {
  print(i)  
  print(r$values[i])
  print(r$lengths[i])
  print("-----")
  counts[r$values[i]] = counts[r$values[i]]  + r$lengths[i]
}
counts



plot(log(counts),type='l')


ham_plas1 = read.table("/Users/mark_enstrom/Barcode/plasmid/Plasmid_Library_NA_SNP490_CAGT_NA_mp_ham.txt",sep=',')
