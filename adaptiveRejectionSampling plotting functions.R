plot.adaptiveRejectionSampling <- function() {
	par(mfrow = c(1,2))
	plot(p.seq, glog(p.seq), type="l", ylim = c(-3, 2))
	points(Tset, glog(Tset))
	del = .6
	par(mfrow = c(1,2))
	plot(p.seq, glog(p.seq), type="l", ylim = c(-3, 2))
	points(Tset, glog(Tset))
	del = .6
## tangents
	#options(warn=-1, warning.expression=-1)
	for (i in 1:length(Tset))
		lines( Tset[i] + c(-1,1)*del, 
			glog(Tset[i]) + glogdot(Tset[i])*(c(-1,1)*del), lty=2)
	## secants
	bottom = par("usr")[3]
	lines( c(Tset[1],Tset,Tset[length(Tset)]),	
				c(bottom, glog(Tset), bottom), lty=2)
	### Now, back to original scale.
	plot(p.seq, g(p.seq), type="l", ylim = c(0, 4), xlim=c(0,1))
	points(Tset, g(Tset))
	del = .6
	## tangents
	for (i in 1:length(Tset))
		lines( Tset[i] + seq(-1,1,0.01)*del, 
			exp(glog(Tset[i]) + glogdot(Tset[i])*(seq(-1,1,0.01)*del)), lty=2)
	## secants
	bottom = par("usr")[3]
	for (i in 2:length(Tset))
		lines( Tset[i-1]*seq(0,1,0.01) + Tset[i]*seq(1,0,-0.01), 
			exp(glog(Tset[i-1])*seq(0,1,0.01) + glog(Tset[i])*seq(1,0,-0.01)), lty=2)		
}
