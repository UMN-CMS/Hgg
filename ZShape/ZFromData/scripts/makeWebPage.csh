#!/bin/tcsh -f
#################################################
# This super simple script makes a basic webpage
#
#################################################
##
## Argument 1 is the directory that contains the plots 
##

set PlotDir=$1

cat > ${PlotDir}/index.html <<EOF
<HTML>

<HEAD><TITLE>Electron Efficiency Web Pages for `date`</TITLE></HEAD>
 
<BODY link="Red">
<FONT color="Black">
<h2><A name="EB"><FONT color="Black">Electron Efficiency Web Pages for `date`</FONT></A><BR></h2>
<h3><A name="EB"><FONT color="Blue">Summed Efficiency Histograms</FONT></A><BR></h3>

EOF

foreach physvar ( Pt EtaDet )
   echo '<A HREF=sbs_eff_'${physvar}'.png> <img height="300" src="sbs_eff_'${physvar}'.png"> </A>' >> ${PlotDir}/index.html
end

echo '<h3><A name="EB"><FONT color="Blue">Summed Denominator Histograms</FONT></A><BR></h3>' >> ${PlotDir}/index.html

foreach physvar ( Pt EtaDet )
   echo '<A HREF=sbs_den_'${physvar}'.png> <img height="300" src="sbs_den_'${physvar}'.png"> </A>' >> ${PlotDir}/index.html
end

echo '<h3><A name="EB"><FONT color="Blue">Individual Efficiency Histograms</FONT></A><BR></h3>' >> ${PlotDir}/index.html

foreach physvar ( EtaDet Pt )
   foreach plot ( `/bin/ls ${PlotDir} |grep png | grep "eff_${physvar}_" |grep -v "_Pt_EtaDet_"` )
      echo '<A HREF='${plot}'> <img height="300" src="'${plot}'"> </A>' >> ${PlotDir}/index.html
   end
end

echo '<h3><A name="EB"><FONT color="Blue">Individual Denominator Histograms</FONT></A><BR></h3>' >> ${PlotDir}/index.html

foreach physvar ( EtaDet Pt )
   foreach plot ( `/bin/ls ${PlotDir} |grep png | grep "den_${physvar}_" |grep -v "_Pt_EtaDet_"` )
      echo '<A HREF='${plot}'> <img height="300" src="'${plot}'"> </A>' >> ${PlotDir}/index.html
   end
end

echo '<h3><A name="EB"><FONT color="Blue">2-D Efficiency Histograms</FONT></A><BR></h3>' >> ${PlotDir}/index.html

foreach plot ( `/bin/ls ${PlotDir} |grep png |grep  "eff_Pt_EtaDet_"` )
      echo '<A HREF='${plot}'> <img height="300" src="'${plot}'"> </A>' >> ${PlotDir}/index.html
   end
   
   echo '<h3><A name="EB"><FONT color="Blue">2-D Denominator Histograms</FONT></A><BR></h3>' >> ${PlotDir}/index.html

foreach plot ( `/bin/ls ${PlotDir} |grep png |grep  "den_Pt_EtaDet_"` )
      echo '<A HREF='${plot}'> <img height="300" src="'${plot}'"> </A>' >> ${PlotDir}/index.html
   end


#end of efficiency web page making script
