#!/bin/csh -f

set elist = "07 09 11 13 15 17 19a 19b 21 23 25 27 31 33 35 37a 37b 39 41a 41b 43 45a 45b 47 49 51 55 61 63 65"
#set elist = "07"

set MODE = "RHO"
set STR = 1
set MYPATH = "/gpfs/home/belle/vitaly/work/mixing/MCProd/gsim/mdst/BtoDDs/"

foreach l($elist)
  ./SubmitSignal.csh $l $MYPATH 
end
