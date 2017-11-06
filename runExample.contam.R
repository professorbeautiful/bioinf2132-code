runExample.contam = function(
  contam.N = 100,
  p1=0.8, sig1=1, sig2=4, mu=0,
  phi.true = c(p1=p1, sig1=sig1, sig2=sig2, mu=mu),
  savedseed=NULL,
  model.file="contaminated_normal_model.txt",
  N.CHAINS = 8,
  initList = list(
    Jminus1=rep(1,contam.N),
    mu=1,
    tau=c(1,1),
    p1=0.5),
  N.ITER = 1000,
  N.BURNIN = 10,
  N.THIN=5,
  xvalues=seq(-5,5,0.1),
  plot.make=TRUE,
  print.make=TRUE,
  plot=TRUE,
  print=TRUE
){
  require("mvbutils")
  require("R2OpenBUGS")
  dataset.contam.complete <- 
    make.dataset.contam(phi.true, 100, 
                        savedseed = savedseed,
                        plot=plot.make, print=print.make, xvalues=xvalues)
  initListAll = eval(parse(
    text="list(" %&% 
      paste(  rep("initList", N.CHAINS), collapse=",")  %&% ")"
  ))
  length(initListAll )
  
  OpenBUGS.dir=
    "/Users/Roger/.wine/drive_c/Program Files (x86)/OpenBUGS/OpenBUGS323"
  OpenBUGS.pgm=paste0(OpenBUGS.dir, "/OpenBUGS.exe")
  
  bugs.out = bugs( 
    data=list(y=dataset.contam.complete$y),
    inits=initListAll,
    parameters.to.save=cq(p1,tau,mu,J1,logTauDiff,logTauAbsDiff), 
    n.iter=N.ITER, 
    model.file=model.file,
    n.chains=N.CHAINS, 
    n.burnin=N.BURNIN, 
    n.thin=N.THIN, 
    debug=FALSE, 
    DIC=TRUE, digits=5, codaPkg=FALSE,
    OpenBUGS.pgm=OpenBUGS.pgm,
    useWINE = TRUE,
    newWINE = TRUE,
    clearWD=FALSE,  
    bugs.seed=1, summary.only=FALSE,
    over.relax = FALSE 
  )
  attr(bugs.out, "savedseed.for.data") = 
    attributes( dataset.contam.complete )$savedseed
  return(bugs.out)
}
myBugsPlot = function(.bugs.out=bugs.out,
                      component='p1') {
  attach.bugs(.bugs.out)
  PointsPerChain = .bugs.out$n.iter-.bugs.out$n.burnin
  plot(c(1,PointsPerChain), c(0,1), pch=" ",
       main=paste("chains for,", component) )
  every1000 = seq(1, PointsPerChain, length=1000)
  for(chain in 1:.bugs.out$n.chains) 
    lines(every1000, 
          .bugs.out$sims.array[every1000 , chain, component],
          col=chain+1)
  plot(density(.bugs.out$sims.list[[component]], from=0, to=1), 
       ylim=c(0,10), lwd=3,
       main=paste("density for", component, " \n",
                  .bugs.out$n.chains, " chains, BLACK=average")
  )
  for(chain in 1:.bugs.out$n.chains) 
    lines(density( 
      .bugs.out$sims.array[ , chain, component]),
      col=chain+1)
  detach.bugs()
}
bugs.out = runExample.contam()
### err:ole:CoReleaseMarshalData IMarshal::ReleaseMarshalData failed with error 0x8001011d
### does not seem to be a problem.
myBugsPlot()
attach.bugs(bugs.out)
pairs(data.frame(p1, mu, tau, logTauDiff),
      col=(logTauDiff < 0)+1)
plot(density(logTauDiff))
lines(density(logTauAbsDiff), add=T, col='red')
plot(p1, logTauDiff)
detach.bugs()