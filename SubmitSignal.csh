#!/bin/csh -f
set INDEX = 0
set EXP = $1
set MYPATH = $2

echo $2'evtgen_exp_'$1'_BtoDDs-'${INDEX}'.mdst'
while(-e $2'evtgen_exp_'$1'_BtoDDs-'${INDEX}'.mdst')
  ./SigMC.csh $EXP $MYPATH $INDEX > scripts/SigDDs$EXP-$INDEX.csh
  chmod 755 scripts/SigDDs$EXP-${INDEX}.csh
  bsub -q s scripts/SigDDs$EXP-${INDEX}.csh
  @ INDEX = $INDEX + 1
end

