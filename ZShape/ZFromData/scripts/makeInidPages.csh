#!/bin/tcsh -f
#################################################
# This super simple script makes a basic webpage
#
#################################################
##
## Argument 1 is the directory that contains the plots 
##
set mydir=`pwd`
set PlotDir=$1
set plottype=$2

cat > ${PlotDir}/index.html <<EOF
<HTML>

<HEAD><TITLE>ZShape Result WebPages for `date`</TITLE></HEAD>
 
<BODY link="Red">
<FONT color="Black">
<h2><A name="EB"><FONT color="Black">ZShape Fit Result Web Pages for `date`</FONT></A><BR></h2>

<h3><A name="EB"><FONT color="Blue">Fit Result Histograms WP80-WP95</FONT></A><BR></h3>


EOF

foreach FromData (`/bin/ls ${PlotDir} |grep "BKGrap"| grep "png"  `)
   echo '<A HREF='${FromData}'> <img height="300" src="'${FromData}'"> </A>' >> ${PlotDir}/index.html
end

echo '<h3><A name="EB"><FONT color="Blue">Fit Result Histograms WP80-HF</FONT></A><BR></h3>' >> ${PlotDir}/index.html

foreach FromData (`/bin/ls ${PlotDir} |grep "BKGHFrap"| grep "png"  `)
   echo '<A HREF='${FromData}'> <img height="300" src="'${FromData}'"> </A>' >> ${PlotDir}/index.html
end



if (`echo $plottype | grep -c "eps"` > 0 ) then
foreach FromData (`/bin/ls ${PlotDir} |grep ".eps"`)
   ~/eps2pdf.perl ${PlotDir}/$FromData
end

mv *pdf ${PlotDir}/.

foreach FromData (`/bin/ls ${PlotDir}/.. |grep ".eps"`)
   ~/eps2pdf.perl ${PlotDir}/../$FromData
end

mv *pdf ${PlotDir}/../.

endif






#end of ZShape web page making script
