function read_syn, name
 line = ''
 get_lun, luni
 openr,luni,name
 repeat begin
   readf,luni,line 
 endrep until (((strpos(line,'$SPECTRUM:') +  strpos(line,'$DATA:')) GE -1) OR eof(luni))
 if eof(luni) then return, 0
 nchan = 0
 readf,luni, nchan,nrows
 data = fltarr(nrows,nchan)
 readf, luni, data
 close, luni
 free_lun, luni
 return, data
end
