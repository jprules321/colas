function test()
'reinit'
*'@ open $GADSET/model'
'open ../../pytests/data/model'
'set lon -180 180'
'set gxout grfill'
'load udxt ../re/re.udxt'
'load udxt orb.udxt'
'xts = re(ts,0.5)'
* 'set_orb dt 1'
'd orb_mask(xts,aqua,600)'
*'d orb_mask(xts,aqua)'
*'d ts'
*'orb_track aqua 30'
