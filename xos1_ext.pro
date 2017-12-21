length = 9.
d0 = 0.413
d1 = 0.117
s0 = 7.e-4
f1 = 0.5
r0=d0/2.
r1=d1/2.
thick


tga = d1/2./f1

x0 = [0.,.1,.2,length*9./10.,length-.2,length-.1,length]
y0 = [d0/2.,d0/2,d0/2,(d0+d1)/4.,d1/2.+.2*tga,d1/2.+.1*tga,d1/2.]

xf = findgen(1000)/999.*length
yf = spline(x0,y0,xf)


!xtitle = 'Length, [cm]'
!ytitle = 'External radius, [cm]'


plot,xf,yf
;oplot,xf,yf


a = fltarr(3,1000)
a(0,*) = xf
b = fltarr(2,1000)
c = fltarr(2,1000)
b(0,*) = xf
b(1,*) = 0.5 * (s0 - findgen(1000)/999.*(s0 - s0*d1/d0))
c(0,*) = xf
c(1,*) = yf

waxis,'xos1.axs',a
wprf,'xos1.prf',b
wprf,'xos1.ext',c


plot,c(0,*),c(1,*)


end



