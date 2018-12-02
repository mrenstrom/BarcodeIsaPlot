isa8103@labels


#
# ISA 8103
#
for (x in c(1:24)) {
  d = isa8103@data[[x]]
  l = isa8103@labels[x]
  print(l)
  #print(nrow(d))
  #print("----------")
}
#
# BC8103
#
for (x in c(1:29)) {
  d = bc8103@data[[x]]
  l = bc8103@labels[x]
  print(l)
  #print(nrow(d))
  #print("----------")
}
#
#  combine clones for all 8 BM 314DPT samples
#
clones=c()
for (x in c(10:17)) {
  d = bc8103@data[[x]]
  clones = c(clones,d[,1])
  print(length(clones))
}
total=unique(clones)
length(total)
#
# 733 combined
#
clones=c()
for (x in c(25:28)) {
  d = bc8103@data[[x]]
  clones = c(clones,d[,1])
  print(length(clones))
}
total=unique(clones)
length(total)
#
# 9132 isa
#
isa9132@labels
for (x in c(1:24)) {
  d = isa9132@data[[x]]
  l = isa9132@labels[x]
  print(l)
  #print(nrow(d))
  #print("----------")
}
#
# BC9132 
#
bc9132@labels
for (x in c(1:26)) {
  d = bc9132@data[[x]]
  l = bc9132@labels[x]
  print(l)
  #print(nrow(d))
  #print("----------")
}
#
# 299 DPT 
#
clones=c()
for (x in c(8:15)) {
  d = bc9132@data[[x]]
  clones = c(clones,d[,1])
  print(length(clones))
}
total=unique(clones)
length(total)
#
# 701 combined
#
clones=c()
for (x in c(22:25)) {
  d = bc9132@data[[x]]
  clones = c(clones,d[,1])
  print(length(clones))
}
total=unique(clones)
length(total)
xs