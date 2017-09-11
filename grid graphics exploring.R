### https://journal.r-project.org/archive/2015-1/murrell.pdf

library(gridGraphics)
### grid.echo()    from standard graphics to grid
ggs_density(ggsObj, family="phi")  ## ggplot2
grid.ls()  ### from gridGraphics.
grid.ls()$name
upViewport(0)
downViewport("GRID.gTableParent.3038")
seekViewport("geom_density.polygon.2960")  ### no good.
vpPath(current.viewport()$name, "geom_density.polygon.2960")
seekViewport('ROOT::geom_density.polygon.2960') 
vpPath("ROOT", "GRID.gTableParent.3038")
showViewport()
############### 
library(gridBase)
library(grid)
opar <- par(no.readonly=TRUE)
# gridFIG
grid.newpage()
pushViewport(viewport(width=0.5, height=0.5))
grid.rect(gp=gpar(col="grey", lty="dashed"))
par(fig=gridFIG())
par(new=TRUE)
plot(1:10)
# multiple plots
# NOTE the use of par(mfg)
# gridOMI
par(opar)
grid.newpage()
pushViewport(viewport(width=0.5, height=0.5))
grid.rect(gp=gpar(col="grey", lty="dashed"))
par(omi=gridOMI())
par(mfrow=c(2, 2), mfg=c(1, 1), mar=c(3, 3, 1, 0))
for (i in 1:4) {
  plot(i) }
# gridPLT
par(opar)
grid.newpage()
pushViewport(viewport(width=0.5, height=0.5))
grid.rect(gp=gpar(col="grey", lwd=5))
par(plt=gridPLT())
par(new=TRUE)
plot(1:10)
# gridFIG with par(omi) set
par(opar)
grid.newpage()
par(omi=rep(1, 4))
pushViewport(viewport(width=0.5, height=0.5))
grid.rect(gp=gpar(col="grey", lwd=5))
par(fig=gridFIG())
par(new=TRUE)
plot(1:10)
# gridPLT with par(omi) set
par(opar)
grid.newpage()
par(omi=rep(1, 4))
pushViewport(viewport(width=0.5, height=0.5))
grid.rect(gp=gpar(col="grey", lwd=5))
par(plt=gridPLT())
par(new=TRUE)
plot(1:10)
# gridPAR
par(opar)
grid.newpage()
pushViewport(viewport(width=0.5, height=0.5,
                      gp=gpar(col="red", lwd=3, lty="dotted")))
grid.rect(gp=gpar(col="grey", lwd=5))
par(fig=gridFIG())
par(gridPAR())
par(new=TRUE)
plot(1:10, type="b")