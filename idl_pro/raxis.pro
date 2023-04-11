function raxis,name

if(name eq '')then goto,lab
on_ioerror,lab

;print,name
openr,u,name,/get_lun
readf,u,n
x=fltarr(3,n+1)
readf,u,x
close,u
free_lun,u

return,x

lab:
return,0

end


