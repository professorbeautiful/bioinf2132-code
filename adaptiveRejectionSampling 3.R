options(warn=-1)
a = 2
b = 2
p.seq <- seq (0.01, 0.99, 0.01)
g = function(p,a.=a,b.=b) dbeta(p,a.,b.)   #### distribution we want to sample from
glog = function(...) log(g(...))   #### log of the  distribution we want to sample from
glogdot = function(p,a.=a,b.=b) (a.-1)/p - (b.-1)/(1-p)  ####Derivative of log of g


#########  Initialize the plot  #######
#newDevice(width=12, height=6)
split.screen(c(1,2))
screen(1)
plot(p.seq, glog(p.seq), type="l", ylim = c(-3, 2))
points(Tset, glog(Tset))
title ("log of g, the likelihood (* prior) of interest")
screen(2)
plot(p.seq, g(p.seq), type="l", ylim = c(0, 4), xlim=c(0,1))
points(Tset, g(Tset))



############## Begin sampling #####################
sampleSize = 10
redrawAt = 1:10
Tset = c(1/3, 3/4)   ### dividing points between the segments.
acceptSet = numeric(0)
sampleSet = numeric(0)
proposalSet  = numeric(0)
g.u.Set = numeric(0)
for (i in 1:sampleSize) {
	### find intersections for upper envelope
	heights = glog(Tset)		#log scale
	slopes = glogdot(Tset)		#log scale
	x.y.ydot = cbind(Tset, heights, slopes)
	K = length(Tset)
	adjacent = cbind(x.y.ydot[-1, ,drop=F], x.y.ydot[-K, ,drop=F])
	intersections = apply(adjacent, 1, function(r) {
		x1 = r[1]; y1 = r[2]; ydot1 = r[3]; 
		x2 = r[4]; y2 = r[5]; ydot2 = r[6]; 
		return( (y2 - y1 + x1*ydot1 - x2*ydot2)/(ydot1 - ydot2))
	})
	cumprobs.left = exp(heights[-1]+slopes[-1]*(intersections-Tset[-1]))
	cumprobs.right = exp(heights[-K]+slopes[-K]*(intersections-Tset[-K]))
	cumprobs.left = c(0, cumprobs.left)
	cumprobs.right = c(cumprobs.right, 0)
	p.segments = abs(cumprobs.left - cumprobs.right)
	p.segments = 	p.segments /sum(p.segments)
	catn("p.segments = " %&% p.segments) 
	segment = rmultinom(1,1,p.segments)
	u = runif(n) 
	sl.seg = slopes[segment]
	int.seg.left = c(0, intersections)[segment]
	ht.seg.left = glog(Tset[segment]) + 
			glogdot(Tset[segment])*(int.seg.left-Tset[segment])
	int.seg.right = c(intersections, 1)[segment]
	ht.seg.right = glog(Tset[segment]) + 
			glogdot(Tset[segment])*(int.seg.right-Tset[segment])
	cbind(segment, sl.seg, int.seg.right, ht.seg.right, int.seg.left, ht.seg.left)
	proposals = sl.seg^(-1) * log(
			u * exp(sl.seg*int.seg.right)
			+ (1-u) * exp(sl.seg*int.seg.left))
	if(any(is.na(proposals))) browser()
	g.u.proposals.1 = exp(ht.seg.left + (proposals-int.seg.left)*sl.seg)
	g.u.proposals.2 = exp(ht.seg.right + (proposals-int.seg.right)*sl.seg)
	cbind(proposals, g.u.proposals.1, g.u.proposals.2) 
	g.u.proposals = g.u.proposals.1
	g.u.proposals[is.na(g.u.proposals.1)] = g.u.proposals.2[is.na(g.u.proposals.1)]
	r.U.out = data.frame(proposals, g.u.proposals)
	
	proposals = r.U.out$proposals
	g.u.proposals = r.U.out$g.u.proposals
	u.n = runif(1)
	proposalSet = c(proposalSet, proposals)
	g.u.Set = c(g.u.Set, u.n*g.u.proposals)
	cat("==== ", i, "===\n")
	if(i < 20) print(g.u.proposals)
	if(i < 20) print(rbind(proposalSet, acceptSet))
	accept <- (u.n < dbeta(proposals, a, b)/g.u.proposals)
	if(!accept)
		Tset = sort(c(Tset, proposals))
	else
		sampleSet = c(sampleSet, accept)
	acceptSet = c(acceptSet, accept)
	catn(Tset)
	if(i %in% redrawAt) {
		plot.ars()
		points (proposalSet[acceptSet==1], 
			g.u.Set [acceptSet==1], col=2)
		points (proposalSet[acceptSet==0], 
			g.u.Set [acceptSet==0], pch=pch("ssquare"), col=4)
		browse()
	}
}

#plot( sort(theta[accept]), qbeta(seq(0,1,length=n.accept), a, b))

