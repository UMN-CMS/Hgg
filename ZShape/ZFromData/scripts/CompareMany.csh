#!/bin/tcsh -f

set mydate=`date +"%b%dth%Y"`
set maindr=/data/cms004/recovery/Webpages/Electrons/2010/EffxAcc/${mydate}${1}
mkdir ${maindr}
./CompareEffAccPlots.csh ${maindr}/WP95WP95 png 1
./CompareEffAccPlots.csh ${maindr}/WP80WP80 png 2
./CompareEffAccPlots.csh ${maindr}/WP95WP95EACH png 3
./CompareEffAccPlots.csh ${maindr}/WP80WP80EACH png 4
./CompareEffAccPlots.csh ${maindr}/WP80WP95EACH png 114
./CompareEffAccPlots.csh ${maindr}/WP80WP95+WP80HFLoose png 6
./CompareEffAccPlots.csh ${maindr}/WP80WP95+WP80HFTight png 5
./CompareEffAccPlots.csh ${maindr}/WP80WP95EACH+WP80HFLoose png 14

./CompareEffAccPlots.csh ${maindr}/WP80NTTight png 160
./CompareEffAccPlots.csh ${maindr}/WP80WP95EACH+WP80HFLoose+WP80NTTight png 161
./CompareEffAccPlots.csh ${maindr}/WP80WP95+WP80HFLoose+WP80NTTight png 16




#./CompareEffAccPlots.csh ${maindr}/WP95WP95+WP80HFLoose eps 6
#./CompareEffAccPlots.csh ${maindr}/WP95WP95+WP80HFLoose C 6
#./CompareEffAccPlots.csh ${maindr}/WP95WP95+WP80HFLoose root 6
#./CompareEffAccPlots.csh ${maindr}/WP95WP95EACH+WP80HFLoose eps 14
#./CompareEffAccPlots.csh ${maindr}/WP95WP95EACH+WP80HFLoose C 14
#./CompareEffAccPlots.csh ${maindr}/WP95WP95EACH+WP80HFLoose root 14

./CompareEffAccPlots.csh ${maindr}/WP80WP95+WP80HFLoose eps 6
./CompareEffAccPlots.csh ${maindr}/WP80WP95+WP80HFLoose C 6
./CompareEffAccPlots.csh ${maindr}/WP80WP95+WP80HFLoose root 6
./CompareEffAccPlots.csh ${maindr}/WP80WP95+WP80HFLoose+WP80NTTight eps 16


echo "Done with all the default directories and plots"

#end of the file
