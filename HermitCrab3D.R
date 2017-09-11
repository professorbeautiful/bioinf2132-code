###  http://www.linkedin.com/redir/redirect?url=https%3A%2F%2Fwww%2Elinkedin%2Ecom%2Fpulse%2Fhermit-crab-gets-3d-printed-designer-shell-modeled-using-jose-i-rey&urlhash=J-Ev&_t=tracking_anet
# Jose Rey:  Hermit Crab 3D 

require(rgl)
require(plot3D)
M=mesh(seq(0,2*pi, length.out = 150),
       seq(-pi/2, 5.33*pi/2, length.out=150))
u=M$x
v=M$y
x=10*1.2^v*(sin(u)^8*sin(v))
y=8*1.2^v*(sin(u)^2*cos(v))
z=10*1.2^v*(sin(u-0.4)^8*cos(u+0.4))

surface3d(x,y,z, col="chocolate", back="lines")
writeOBJ('shellR2.obj', pointRadius=0.001,
         pointShape=icosahedron3d())
movie3d(spin3d(axis = c(1, 1, 1)), duration = 20, 
        dir = getwd(), movie = "shell")
