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
;=========================================

PRO plot_eff

name = dialog_pickfile(filter='*.out',/must_exist,/multiple)
colors = GET_COLOR(n_elements(name)-1,/RANDOM,DARK=400)
for i =0, n_elements(name)-1 do begin
  x=read_syn(name[i])
  if i eq 0 then p = plot(x[0,*],x[1,*]*100,xtitle='Energy [keV]',ytitle='Efficiency, %') else $
    p = plot(x[0,*],x[1,*]*100,xtitle='Energy [keV]',ytitle='Efficiency, %',/current,/overplot,color=colors[i-1])
  legname0 = strsplit(name[i],'/',/extract)
  legname = legname0[n_elements(legname0)-1]
  leg = legend(target=p,horizontal_alignment=0,shadow=0,label=legname,transparency=100,position=[0.65,0.85-0.05*i])
  endfor


end

