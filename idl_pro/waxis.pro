pro waxis,name,x

if(name eq '')then name = 'axis.dat'
on_ioerror,lab

ndim = size(x)

openw,u,name,/get_lun
printf,u,ndim(2)-1
printf,u,x
close,u
free_lun,u

return
lab:
print,'Unknown I/O error...'
return

end


