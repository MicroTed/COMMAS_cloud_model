#!/bin/csh


set frame = $1
set numarg = $#argv
echo number of arguments is $numarg
echo number of arguments is $#argv
#echo RUN = ${RUN}
#exit

set i=0
foreach runname ($argv)
  #set i = $i + 1
  @ i++
  if ( $i > 1 ) then
  echo i = $i
  echo  $runname $frame $frame  $runname.${frame}.ps
  pssplit $runname $frame $frame > ${runname}.${frame}.ps
#  set OUT = out.complot.${runname}.dbzth
#  echo out = ${OUT}
  endif
end