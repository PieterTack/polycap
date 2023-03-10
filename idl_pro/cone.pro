length = 9.
d0 = 0.413
d1 = 0.117
s0 = 7e-4
r0=d0/2.
r1=d1/2.


xf = findgen(1000)/999.*length
;cone
yf = r0+findgen(1000)/999.*(r1-r0) 




a = fltarr(3,1000)
a(0,*) = xf
b = fltarr(2,1000)
c = fltarr(2,1000)
b(0,*) = xf
b(1,*) = reverse(0.5 * (s0 - findgen(1000)/999.*(s0 - s0*d1/d0)))
c(0,*) = xf
c(1,*) = reverse(yf)

p = plot(b,xtitle='Length, [cm]',ytitle='Capillary radius, [cm]')
p = plot(c,xtitle='Length, [cm]',ytitle='External radius, [cm]')

waxis,'cone_conf.axs',a
wprf,'cone_conf.prf',b
wprf,'cone_conf.ext',c




end



