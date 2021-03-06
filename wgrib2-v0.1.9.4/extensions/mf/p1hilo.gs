function p1hilo(args)

i=1
chrhl=subwrd(args,i) ; i=i+1
lonhl=subwrd(args,i) ; i=i+1
lathl=subwrd(args,i) ; i=i+1
valhl=subwrd(args,i) ; i=i+1
grdhl=subwrd(args,i) ; i=i+1
lplhl=subwrd(args,i) ; i=i+1
hlfmt=subwrd(args,i) ; i=i+1

*
*	defaults to philo
*
if(hlfmt='' | hlfmt='hlfmt') ; hlfmt='%5.0f' ; endif

cvalhl=math_format(hlfmt,valhl)

*
*	plotting params
*

hlmk=1
hlmksiz=0.05
hlmkthk=0.05
hlmkcol=1

hlchoffx=0.0
hlchoffy=0.10

hlchsizl=0.10
hlchcoll=2
hlchthkl=6
hlchfntl=5

hlchsizh=0.10
hlchcolh=1
hlchthkh=6
hlchfnth=5

hlvlsiz=0.08
hlvlcol=1
hlvlthk=6

pcnth=1
pcntl=1


'q w2xy 'lonhl' 'lathl
xhl=subwrd(result,3)
yhl=subwrd(result,6)
xhlc=xhl+hlchoffx
yhlc=yhl+hlchoffy

xhlv=xhl-hlchoffx
yhlv=yhl-hlchoffy

'set string 'hlmkcol' c 'hlmkthk
'draw mark 'hlmk' 'xhl' 'yhl' 'hlmksiz


*
*	read and plot H's
*
if(chrhl = 'H')

  'set string 'hlchcolh' c 'hlchthkh
  'set strsiz 'hlchsizh
  'draw string 'xhlc' 'yhlc' `'hlchfnth%chrhl

  'set string 'hlvlcol' c 'hlvlthk
  'set strsiz 'hlvlsiz
  'draw string  'xhlv' 'yhlv' 'cvalhl
 
endif

if(chrhl = 'L')

  'set string 'hlchcoll' c 'hlchthkl
  'set strsiz 'hlchsizl
  'draw string 'xhlc' 'yhlc' `'hlchfntl%chrhl

  xhlv=xhl-hlchoffx
  yhlv=yhl-hlchoffy
  'set string 'hlvlcol' c 'hlvlthk
  'set strsiz 'hlvlsiz
  'draw string  'xhlv' 'yhlv' 'cvalhl

endif

lat1=lathl+1.0
'q w2xy 'lonhl' 'lat1
x1=subwrd(result,3)
y1=subwrd(result,6)

rc=xhl' 'yhl' 'y1

return(rc)
