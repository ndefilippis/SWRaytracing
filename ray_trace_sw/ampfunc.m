function A = ampfunc(x_,y_,xp,yp,a,sigma)


A = a*exp( -( (x_-xp).^2+(y_-yp).^2 )/(2*sigma^2));
