pro wprf,name,x

on_ioerror,lab

n = size(x)
n = n(2)

openw,unit,name,/get_lun
printf,unit,n-1
printf,unit,x
close,unit

free_lun,unit

lab:
return

end


