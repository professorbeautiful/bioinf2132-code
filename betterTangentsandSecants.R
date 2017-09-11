
####  Not completed.


drawTangentsAndSecants = function(screenNum=1, transform=function(x)x, tandel=0.6, steps=seq(-1,1,0.01), color="red"){
	screen(screenNum)
	## tangents
	for (i in 1:length(Tset))
		lines( Tset[i] + seq(-1,1,0.01)*del, 
			transform(glog(Tset[i]) + glogdot(Tset[i])*(seq(-1,1,0.01)*del)), lty=2, col=color)
	## secants
	bottom = par("usr")[3]; catn("bottom = ", bottom)
	for (i in 2:length(Tset))
		lines( Tset[i-1]*seq(0,1,0.01) + Tset[i]*seq(1,0,-0.01), 
			transform(glog(Tset[i-1])*seq(0,1,0.01) + glog(Tset[i])*seq(1,0,-0.01)), lty=2, col=color)
}		
#	segments(x0=Tset-steps*del, y0=transform(glog(Tset)-glogdot(Tset)*steps*del), 
#			  x1=Tset+steps*del, y1=transform(glog(Tset)+glogdot(Tset)*steps*del), col=color)
	## secants
#	segments(x0=Tset-steps*del, y0=transform(glog(Tset)-glogdot(Tset)*steps*del), 
#			  x1=Tset+steps*del, y1=transform(glog(Tset)+glogdot(Tset)*steps*del))
}
drawTangentsAndSecantsRight = function(del=0.6,transform){
	lines( c(Tset[1],Tset,Tset[length(Tset)]), ## secants
			c(bottom, glog(Tset), bottom), lty=2)
	screen(2)
	segments(Tset[i] + seq(-1,1,0.01)*del, 
			transform(glog(Tset[i]) + glogdot(Tset[i])*(seq(-1,1,0.01)*del)), lty=2)
	for (i in 2:length(Tset))
	segments()
		lines( Tset[i-1]*steps + Tset[i]*(1-steps)), 
			transform(glog(Tset[i-1])*steps + glog(Tset[i])*seq(1,0,-0.01)), lty=2)		
}
drawTangentsAndSecants()
