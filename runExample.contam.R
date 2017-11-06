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
	N.ITER = 10000,
	N.BURNIN = 10,
	xvalues=seq(-5,5,0.1),
	plot.make=TRUE,
	print.make=TRUE,
	plot=TRUE,
	print=TRUE
){
	require("mvbutils")
	require("R2OpenBUGS")
	dataset.contam.complete <- make.dataset.contam(phi.true, 100, 
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
		n.thin=1, 
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
	if(plot) {
		PointsPerChain = N.ITER-N.BURNIN
		plot(rep(1:PointsPerChain, N.CHAINS), bugs.out$sims.list$p1, pch=" ")
		for(chain in 1:N.CHAINS) 
			lines(1:PointsPerChain, 
			bugs.out$sims.array[ , chain, "p1"],
			col=chain+1)
		plot(density(bugs.out$sims.list$p1, from=0, to=1), ylim=c(0,10), lwd=3 )
		for(chain in 1:N.CHAINS) 
			lines(density( 
					bugs.out$sims.array[ , chain, "p1"]),
				col=chain+1)
	}
	return(bugs.out)
}

bugs.out = runExample.contam()

