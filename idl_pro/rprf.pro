function rprf,name

on_ioerror,lab

if(strpos(name,'.dat') ne -1)then begin
   rdat,name
   return,0
endif


openr,unit,name,/get_lun
readf,unit,n
x=fltarr(2,n+1)
readf,unit,x
close,unit

free_lun,unit

name0 = str_sep(name,'.')
name_axis=name0(0)+'.axs'
ii = findfile(name_axis)
if(ii(0) ne '')then begin
;   spawn,'cp '+ii(0)+' axis.dat'
   print,ii(0)
   ax = raxis(ii(0))
endif else begin
   ax = fltarr(3,n+1)
   ax(0,*) = x(0,*)
   openw,1,'axis.dat'
   printf,1,long(n)
   printf,1,ax
   close,1
endelse

return,x

lab:
return,0

end


