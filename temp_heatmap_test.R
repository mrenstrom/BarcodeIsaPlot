data(mtcars)
x  <- as.matrix(mtcars)


heatmap.2(x,key = TRUE)
          

heatmap.2(x,
          key.title=NA, # no title
          key.xlab=NA  # no xlab
          #key.par=list(mgp=c(1.5, 0.5, 0),
          #             mar=c(2.5, 2.5, 1, 0))
          )

heatmap.2(x,
          key.title=NA, # no title
          key.xlab=NA  # no xlab
          #key.par=list(mgp=c(1.5, 0.5, 0),
          #             mar=c(2.5, 2.5, 1, 0))
)






heatmap.2(x,
          key.title=NA, # no title
          key.xlab=NA,  # no xlab
          key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(2.5, 2.5, 1, 0)),
          key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            return(list(
              at=parent.frame()$scale01(c(breaks[1],
                                          breaks[length(breaks)])),
              labels=c(as.character(breaks[1]),
                       as.character(breaks[length(breaks)]))
            ))
          })


heatmap.2(x,
          breaks=256,
          key.title=NA,
          key.xlab=NA,
          key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(1, 2.5, 1, 0)),
          key.xtickfun=function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("low", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("high", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          })





