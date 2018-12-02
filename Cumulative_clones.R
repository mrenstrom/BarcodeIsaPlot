#
# cumulative by day
#
#
bc8103@labels

days = list (
  list("33",c(1)) ,
  list("71",c(2)) ,
  list("145",c(3)),
  list("160",c(4)),
  list("193",c(5)),
  list("223",c(6)),
  list("281",c(7,8,9)),
  list("314",c(10:17)),
  list("368",c(18)),
  list("433",c(19)),
  list("502",c(20)),
  list("644",c(21,22,23,24)),
  list("733",c(25,26,27,28)),
  list("834",c(29))
)

isa8103@labels

isa8103_days = list (
  list("33",c(1)) ,
  list("34",c(2)) ,
  list("56",c(11)),
  list("160",c(3)),
  list("223",c(4)),
  list("308",c(5,6,7)),
  list("368",c(8)),
  list("433",c(9)),
  list("502",c(10)),
  list("602",c(12,13,14,15)),
  list("644",c(16,17,18,19)),
  list("747",c(20:24))
)
#
# !!! Much simpler way to do this uning unique
#
#
# for each day add new clones
#
cumulativeClones <- function(clones,bcObj,tests)
{
  for (i in tests) {
    message("---------")
    message(paste0("i = ",i))
    d = bcObj@data[[i]]
    added = 0
    for (i in c(1:nrow(d))) {
      bc = d[i,1]
      if (bc %in% clones) {
        #message(paste0("already in ",bc))
      } else {
        clones = c(clones,bc)
        #message(paste0("add new ",bc))
        added  = added + 1
      }
    }
    message(paste0("new clones = ",added))
    message(paste0("total = ",length(clones)))
  }
  return(clones)
}

#
# BC8103
#
clones = c()
dptCount = list(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
bc8103_y = c()
for (i in c(1:length(days))) {
  lx = unlist(days[[i]][1])
  bc8103_y = c(bc8103_y,lx)
  lv = days[[i]][2]
  lv = unlist(lv)
  #message(lx)
  #message(lv)
  #message(class(lv))
  clones = cumulativeClones(clones,bc8103,lv)
  dptCount[[i]] = length(clones)
}
bc8103_y
bc8103_dptCount = unlist(dptCount)
bc8103_dptCount
plot(bc8103_dptCount,type='l')
#
# ISA Z08103
#
clones = c()
dptCount = list(0,0,0,0,0,0,0,0,0,0,0,0)
isa8103_y = c()

for (i in c(1:length(isa8103_days))) {
  lx = unlist(isa8103_days[[i]][1])
  isa8103_y = c(isa8103_y,lx)
  lv = isa8103_days[[i]][2]
  lv = unlist(lv)
  #message(lx)
  #message(lv)
  #message(class(lv))
  clones = cumulativeClones(clones,isa8103,lv)
  dptCount[[i]] = length(clones)
}
isa8103_y
isa8103_dptCount = unlist(dptCount)
plot(isa8103_dptCount,type='l')
#
# plot results
#
i_df   = data.frame(dpt = as.numeric(isa8103_y), count = as.numeric(isa8103_dptCount), subject = "ISA Z08103")
i_df2  = data.frame(dpt = as.numeric(bc8103_y), count = as.numeric(bc8103_dptCount), subject = "DBS Z08103")

i_df3 <- rbind(i_df,i_df2)
i_df3
#
#
#
g <- ggplot(i_df3,aes(dpt,count,group=subject))

base <- g + geom_line(aes(color=subject),size=1.2) + 
  labs(x = "Days Post Transplant",y = "Unique Clones Identified", title = "Z08103 Cumulative Clones Dectected") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title = element_text(size=14),
        plot.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.title=element_blank()
  ) +
  scale_color_manual(values=c("#c0c0c0","#404040")) +
  scale_linetype_manual(values = c('dotted','solid')) +
  geom_point(aes(color=subject)) +  
  annotate("segment", x = 365, xend = 365, y = 40000, yend = 35000, color='black', size=1, alpha = 1.0,arrow=arrow()) + 
  theme(legend.position = c(0.82,0.70))


base

pdf(paste0(plotDir,"Cumulative_clone_Z08103.pdf"))

base

dev.off()
#
#BC 9132
#
bc9132@labels
bc9132_days = list (
  list("32",c(1)) ,
  list("104",c(2)) ,
  list("118",c(3)),
  list("256",c(4)),
  list("281",c(5,6,7)),
  list("299",c(8:15)),
  list("400",c(16)),
  list("469",c(17)),
  list("616",c(18:21)),
  list("701",c(22:25)),
  list("802",c(26))
)
clones = c()
dptCount = list(0,0,0,0,0,0,0,0,0,0,0)
bc9132_y = c()
for (i in c(1:length(bc9132_days))) {
  lx = unlist(bc9132_days[[i]][1])
  bc9132_y = c(bc9132_y,lx)
  lv = bc9132_days[[i]][2]
  lv = unlist(lv)
  #message(lx)
  #message(lv)
  #message(class(lv))
  clones = cumulativeClones(clones,bc9132,lv)
  dptCount[[i]] = length(clones)
}
bc9132_y
bc9132_dptCount = unlist(dptCount)
plot(bc9132_dptCount,type='l')
#
# ISA9132
#
isa9132@labels
isa9132_days = list (
  list(32,c(1)) ,
  list(59,c(2)) ,
  list(82,c(3)),
  list(104,c(4)),
  list(118,c(5)),
  list(179,c(6)),
  list(256,c(7)),
  list(281,c(8)),
  list(299,c(9)),
  list(400,c(10)),
  list(469,c(11)),
  list(579,c(12:15)),
  list(616,c(16:19)),
  list(720,c(20:23))
)
clones = c()
dptCount = list(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
isa9132_y = c()
for (i in c(1:length(isa9132_days))) {
  lday = unlist(isa9132_days[[i]][1])
  isa9132_y = c(isa9132_y,lday)
  lv = isa9132_days[[i]][2]
  lv = unlist(lv)
  #message(lx)
  #message(lv)
  #message(class(lv))
  clones = cumulativeClones(clones,isa9132,lv)
  dptCount[[i]] = length(clones)
}

isa9132_dptCount = unlist(dptCount)

df   = data.frame(dpt = as.numeric(isa9132_y), count = as.numeric(isa9132_dptCount), subject = "ISA Z09132")
df2  = data.frame(dpt = as.numeric(bc9132_y), count = as.numeric(bc9132_dptCount), subject = "DBS Z09132")

df3 <- rbind(df,df2)
df3

dfp = data.frame(x=365,y=22000,label='V',subject="ISA Z09132")

g <- ggplot(df3,aes(dpt,count,group=subject))

base <- g + geom_line(aes(color=subject),size=1.2) + 
  labs(x = "Days Post Transplant",y = "Unique Clones Identified", title = "Z09132 Cumulative Clones Dectected") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title = element_text(size=14),
        plot.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.title=element_blank()
  ) +
  scale_color_manual(values=c("#c0c0c0","#404040")) +
  scale_linetype_manual(values = c('dotted','solid')) +
  geom_point(aes(color=subject)) +  
  annotate("segment", x = 365, xend = 365, y = 20000, yend = 17500, color='black', size=1, alpha = 1.0,arrow=arrow()) + 
  theme(legend.position = c(0.8,0.75)) 


base

pdf(paste0(plotDir,"Cumulative_clone_Z09132.pdf"))
base 
dev.off()
#
# how many total
#
clones = c()
ld = lapply(c(1:26),function(f){
  message(paste0(f))
  clones <<- c(clones,bc9132@data[[f]][,1])
  message(length(clones))
  return(length(clones))
})
unlist(ld)
u = unique(clones)
length(u)

#
# which global not in this list?
#
g = bc9132@global[,1]
length(g)
length(u)







