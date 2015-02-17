function circle(x,y,r,linespec)
tc=0:0.1:2*pi+0.1;
plot(r*sin(tc)+x,r*cos(tc)+y,linespec);