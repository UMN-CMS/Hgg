#!/bin/tcsh -f
###################
#Needs one argument of the python configuration file
python $1

set Jason2="`cat nothing2.txt`"
set Jason3="`cat nothing3.txt`"

cat EffStats.C | sed "s@MYSUPERCRUST@${Jason2}@g; s@DARUNMAPCRUST@${Jason3}@g "  > nothing.C

g++ nothing.C -o nothing `root-config --cflags --libs` -L $ROOTSYS/lib -lHtml -lMinuit  -lFoam 
./nothing

rm nothing2.txt nothing3.txt nothing.C
#end of the file
