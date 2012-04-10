# -*- coding: utf-8 -*-
import sys


B0=['TH1F','MYH','EffStatsError']

B1={'RedCirc':['kRed','kFullCircle','0.9','3'], 
    'BlueTri':['kBlue','kFullTriangleDown','1.1','2'], 
    'GrnSqur':['kGreen','kFullSquare','0.6','2'], 
    'BlkStar':['kBlack','kFullStar','1.6','3'],
    'CynCirc':['kCyan','kFullCircle','0.9','3'], 
    'MagnTri':['kMagenta','kOpenTriangleUp','1.3','2'],
    'OrngTri':['kOrange','kOpenTriangleUp','1.3','3'],
    'BlkCirc':['kBlack','kFullCircle','1.0','2'],}
B2={'WP95.root':['WP95','CynCirc'],
    'WP80.root':['WP80','RedCirc'],
    'HFEID.root':['HF Electron ID','GrnSqur'],
    'HLT.root':['HLT','MagnTri'],
    'GSF.root':['GSF Track Match','BlueTri'],
    'SC.root':['Supercluster','OrngTri'],
    'ZSummedQuadHist':['Sum in Quad','BlkCirc']
    }
#G ratios, bayes ratios,  doscale,doIndividualFITS,
G=['false',      'false',   'false','false','']

LEGEND=['0.45','0.64','0.65','0.95','l','l','Error From']
#LEGEND=['0.33','0.24','0.68','0.55','p','EB-EE Normalized Rapidity']

XY=['-3.6','3.59','0.00015','0.1']
#XY=['-2.5','2.5','0','3.0']
PLOT=[' HIST',' SAME HIST']




#HERE IS WHERE I ACTUALLY MADE THE FILES TO DO THE PLOTTING
###############THEN THE OUTPUT IS USED BY A SCRIPT to compile with the macro template.

file = open('nothing2.txt','w')
sys.stdout=file

print 'myTClass <'+B0[0]+'* > mc("'+B0[1]+'");'
print 'std::string mystr  = "'+B0[2]+'";'
print 'std::string myops  = "'+PLOT[0]+'";'
print 'std::string my2ops = "'+PLOT[1]+'";'


for myt in B1: 
    print 'myCust_t '+myt+'; '+myt+'.color = '+B1[myt][0]+';   '+myt+'.marker = '+B1[myt][1]+';  '+myt+'.thickness = '+B1[myt][2]+'; '+myt+'.linethick = '+B1[myt][3]+';'

print 'mc.setXYLimits('+XY[0]+','+XY[1]+','+XY[2]+','+XY[3]+');'
print 'mc.setPlotType('+G[0]+','+G[1]+','+G[2]+');'
print 'mc.setSTlab(mystr);'
print 'mc.setIndiPlots('+G[3]+');'

file.close()

file = open('nothing3.txt','w')
sys.stdout=file
keys = B2.keys()

for myfile in sorted(keys, reverse = True ):
    print 'myRunMap["'+myfile+'"]["'+B2[myfile][0]+'"]='+B2[myfile][1]+';'

print 'TLegend *myLeg = new TLegend('+LEGEND[0]+','+LEGEND[1]+','+LEGEND[2]+','+LEGEND[3]+',"'+LEGEND[6]+'");'
print 'std::string legops = "'+LEGEND[5]+'";'
print 'std::string legfops = "'+LEGEND[4]+'";'


sys.stdout = sys.__stdout__   

#sorted(B2, key=lambda B: B[0], reverse = True)

sorted(B2, key=lambda B: B[0])

keys = B2.keys()

for myfile in sorted(keys, reverse = True ):
    print 'myRunMap["'+myfile+'"]["'+B2[myfile][0]+'"]='+B2[myfile][1]+';'


