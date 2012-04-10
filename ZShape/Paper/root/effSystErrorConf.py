# -*- coding: utf-8 -*-
#file = open('nothing2.txt')
import sys


B0=['TH1F','MYNA','EffSyst']
#B0=['TH1F','ZFromData/ECAL95EB-ECAL95EE/C08-m(70,110)/e1_P_t','EB_ET_Shape_Changes_EB-EE']
#B0=['TH1F','ZFromData/ECAL95EB-ECAL95EE/C07-HLT-GSF/e2_P_t','EE_ET_Shape_Changes_EB-EE']
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
    'SC.root':['Supercluster','OrngTri']#,
    #'ZSummedQuadHist':['Sum in Quad','BlkCirc']
    }
#G ratios, bayes ratios,  doscale,doIndividualFITS,
G=['false',      'false',   'false','false','']

LEGEND=['0.45','0.64','0.65','0.95','l','l','Error From']
#LEGEND=['0.33','0.24','0.68','0.55','p','EB-EE Normalized Rapidity']

XY=['-3.6','3.59','0.000005','0.08']
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




B0=['TH1F','mcEff/ECAL95EB-ECAL95EB/C07-HLT-GSF/Z0_mass','EB-EB_Mass_a_term']
#B0=['TH1F','ZFromData/ECAL95EB-ECAL95EE/C08-m(70,110)/e1_P_t','EB_ET_Shape_Changes_EB-EE']
#B0=['TH1F','ZFromData/ECAL95EB-ECAL95EE/C07-HLT-GSF/e2_P_t','EE_ET_Shape_Changes_EB-EE']
B1={'RedCirc':['kRed','kFullCircle','0.9','2'], 
    'BlueTri':['kBlue','kFullTriangleDown','1.1','2'], 
    'GrnSqur':['kGreen','kFullSquare','0.6','2'], 
    'BlkStar':['kBlack','kFullStar','1.6','3'],
    'MagnTri':['kMagenta','kOpenTriangleUp','1.3','2'],
    'BlkCirc':['kBlack','kFullCircle','1.0','2'],}
B2={'AZFD.4.5.NEW.root':['4.5 pb^{-1} Data','BlkCirc'],
    'EcalScanEBp0_0_010.root':['a = 0.01','BlkStar'],
    'EcalScanEBp0_0_020.root':['a = 0.02','BlueTri'],
    'EcalScanEBp0_0_030.root':['a = 0.03 #leftarrow','RedCirc'],
    'EcalScanEBp0_0_040.root':['a = 0.04','GrnSqur'],
    'EcalScanEBp0_0_050.root':['a = 0.05','MagnTri']
    }
#G ratios, bayes ratios,  doscale,
G=['false',      'false',   'true','','']

LEGEND=['0.78','0.64','0.98','0.95','p','l','MC EB-EB Zmass']
#LEGEND=['0.33','0.24','0.68','0.55','p','EB-EE Normalized Rapidity']

XY=['50','140','0','220']
#XY=['-2.5','2.5','0','3.0']
PLOT=[' p',' SAME HIST']






#ZFromData/ECAL80-HFTight-15/C08-m(70,110)/Z0_Y
#ZFromData/ECAL95-ECAL95-15/C08-m(70,110)/Z0_Y

#Z0_mass; Z0_Y, e1_P_t; e2_P_t;
#B0=['TH1F','ZFromData/ECAL95EE-ECAL95EE/C08-m(70,110)/Z0_Y','Rapidity_Changes_EE-EE']
B0=['TH1F','ZFromData/ECAL95EB-ECAL95EE/C08-m(70,110)/Z0_Y','Ratios_Rapidity_Changes_EB-EE_Zoomed2']
#B0=['TH1F','ZFromData/ECAL95EB-ECAL95EE/C08-m(70,110)/e1_P_t','EB_ET_Shape_Changes_EB-EE']
#B0=['TH1F','ZFromData/ECAL95EB-ECAL95EE/C07-HLT-GSF)/e2_P_t','EE_ET_Shape_Changes_EB-EE']
B1={'RedCirc':['kRed','kFullCircle','0.9','2'], 
    'BlueTri':['kBlue','kFullTriangleDown','1.1','2'], 
    'GrnSqur':['kGreen','kFullSquare','0.6','2'], 
    'BlkStar':['kBlack','kFullStar','1.6','3'],
    'MagnTri':['kMagenta','kOpenTriangleUp','1.3','2']}
B2={'DZFDRun2010B147251-148000.root':['6-11.5 pb^{-1}/0-3 pb^{-1} 147251-148000','BlkStar'],
    'BZFDRun2010B146241-147250.root':['3-6 pb^{-1}/0-3 pb^{-1}  146241-147250','BlueTri'],
    'AZFDRun2010ASept17th.root':['0-3 pb^{-1}/0-3 pb^{-1} Sept17thReco','RedCirc']
    }
#G ratios, bayes ratios,  doscale,
G=['true',      'false',   'true','','']

LEGEND=['0.63','0.74','0.98','0.95','p','EB-EE Normalized Rapidity']
#LEGEND=['0.33','0.24','0.68','0.55','p','EB-EE Normalized Rapidity']

#XY=['50','120','0','160']
XY=['-2.5','2.5','0','3.0']
PLOT=[' p',' SAME p']










B0=['TH1F','ZFromData/ECAL95EE-ECAL95EE/C07-HLT-GSF/Z0_mass','Normed_Mass_Changes_EE-EE']
#B0=['TH1F','ZFromData/ECAL95EB-ECAL95EE/C08-m(70,110)/e1_P_t','EB_ET_Shape_Changes_EB-EE']
#B0=['TH1F','ZFromData/ECAL95EB-ECAL95EE/C07-HLT-GSF)/e2_P_t','EE_ET_Shape_Changes_EB-EE']
B1={'RedCirc':['kRed','kFullCircle','0.9','2'], 
    'BlueTri':['kBlue','kFullTriangleDown','1.1','2'], 
    'GrnSqur':['kGreen','kFullSquare','0.6','2'], 
    'BlkStar':['kBlack','kFullStar','1.6','3'],
    'MagnTri':['kMagenta','kOpenTriangleUp','1.3','2']}
B2={'AZFDRun2010B147251-148000.root':['6-11.5 pb^{-1} 147251-148000','BlkStar'],
    'BZFDRun2010B146241-147250.root':['3-6 pb^{-1}  146241-147250','BlueTri'],
    'DZFDRun2010ASept17th.root':['0-3 pb^{-1} Sept17thReco','RedCirc']
    }
LEGEND=['0.6','0.7','0.98','0.95','pl','EE-EE Normalized Mass']
XY=['50','120','0','57']






















B0=['TH1F','mcEff/ECAL95EB-ECAL95EE/C08-m(70,110)/e1_P_t','MCEB_ET_Changes_EB-EE']
B1={'RedCirc':['kRed','kFullCircle','1','4'], 
    'BlueTri':['kBlue','kFullTriangleDown','0.9','3'], 
    'GrnSqur':['kGreen','kFullSquare','0.6','2'], 
    'BlkStar':['kBlack','kFullStar','1.1','2'],
    'MagnTri':['kMagenta','kOpenTriangleUp','1.3','2']}
B2={'MCHist.root':['NEW EffxAcc','BlkStar'],
    'MCHist4thZSmear.root':['5.3 invpb EffxAcc','RedCirc']
    }



#B0=['TH1F','ZFromData/ECAL95EE-ECAL95EE/C08-m(70,110)/Z0_Y','Rapidity_Changes_EE-EE']
#B0=['TH1F','ZFromData/ECAL95EB-ECAL95EB/C08-m(70,110)/Z0_mass','Mass_Changes_EB-EB']
B0=['TH1F','ZFromData/ECAL95EB-ECAL95EE/C08-m(70,110)/e1_P_t','EE_ET_Shape_Changes_EB-EE']
B1={'RedCirc':['kRed','kFullCircle','1','4'], 
    'BlueTri':['kBlue','kFullTriangleDown','0.9','3'], 
    'GrnSqur':['kGreen','kFullSquare','0.6','2'], 
    'BlkStar':['kBlack','kFullStar','1.1','2'],
    'MagnTri':['kMagenta','kOpenTriangleUp','1.3','2']}
B2={'AZFDRun2010B147251-148000.root':['147251-14800','BlkStar'],
    'BZFDRun2010B146751-147500.root':['146751-147250','GrnSqur'],
    'CZFDRun2010B146241-146750.root':['146241-146750','BlueTri'],
    'DZFDRun2010ASept17th.root':['Sept17thReco','RedCirc']
    }


LEGEND=['0.01','0.4','0.12','0.99']



B0=['TH1F','NewFullMCZFD.root','NEWZeeECALHFRelTGraph']
B1={'RedCirc':['kRed','kFullCircle','1','4'], 
    'BlueTri':['kBlue','kFullTriangleDown','0.9','3'], 
    'GrnSqur':['kGreen','kFullSquare','0.6','2'], 
    'BlkStar':['kBlack','kFullStar','1.1','2'],
    'MagnTri':['kMagenta','kOpenTriangleUp','1.3','2']}
B2={
    'ZFromData/ECAL80-HFTight-15/C08-m(70,110)/Z0_Y':['E_T 15','RedCirc'],
    'ZFromData/ECAL80-HFTight-17/C08-m(70,110)/Z0_Y':['E_T 17','BlueTri'],
    'ZFromData/ECAL80-HFTight-20/C08-m(70,110)/Z0_Y':['E_T 20','GrnSqur'],
    'ZFromData/ECAL80-HFTight-25/C08-m(70,110)/Z0_Y':['E_T 25','BlkStar'],
    'ZFromData/ECAL80-HFTight-30/C08-m(70,110)/Z0_Y':['E_T 30','MagnTri']
    }

B0=['TProfile','SM_timing','Laser_Timing_Profile_Before']
#B0=['TH1F','SM_614/SM_614_timingAll_1','Laser_Timing_Hist4Runs614']
B1={'RedCirc':['kRed','kFullCircle','1.0','2'], 
    'BlueTri':['kBlue','kFullTriangleUp','1.0','2'], 
    'GrnSqur':['kGreen','kFullSquare','1.0','2'], 
    'BlkStar':['kBlack','kFullStar','1.3','2']}
B2={
    #'Laser_133463.root':['Run 133463','RedCirc'],
    #'Laser_134648.root':['Run 134648','GrnSqur'],
    'Laser_134432.root':['Run 134432','GrnSqur'],
    'Laser_134438.root':['Run 134438','BlkStar'],
    #'Laser_133877.root':['Run 133877','BlueTri'],
    'Laser_132658.root':['Run 132658','RedCirc'],
    'Laser_132659.root':['Run 132659','BlueTri']

    }


B0=['TProfile','SM_timing','Laser_Timing_PTEST']
#B0=['TH1F','SM_614/SM_614_timingAll_1','Laser_Timing_Hist4Runs614']
B1={'RedCirc':['kRed','kFullCircle','1.0','2'], 
    'BlueTri':['kBlue','kFullTriangleDown','1.0','2'], 
    'GrnSqur':['kGreen','kFullSquare','1.0','2'], 
    'BlkStar':['kBlack','kFullStar','1.3','2'],
    'MagnTri':['kMagenta','kOpenTriangleUp','1.0','2']}
B2={
    'Laser_133463.root':['Run 133463','RedCirc'],
    #'Laser_134648.root':['Run 134648','GrnSqur'],
    'Laser_133877.root':['Run 133877','BlueTri'],
    'Laser_133511.root':['Run 133511','GrnSqur'],
    'Laser_133744.root':['Run 133744','BlkStar'],
    'Laser_134283.root':['Run 134283','BlkStar'],
    'Laser_134285.root':['Run 134285','BlkStar'],
    'Laser_134286.root':['Run 134286','BlkStar'],
    'Laser_134289.root':['Run 134289','BlkStar'],
    'Laser_134290.root':['Run 134290','BlkStar'],
    'Laser_134291.root':['Run 134291','BlkStar'],
    'Laser_134293.root':['Run 134293','BlkStar'],
    'Laser_134298.root':['Run 134298','BlkStar'],
    'Laser_134313.root':['Run 134313','BlkStar'],
    'Laser_134314.root':['Run 134314','BlkStar'],
    'Laser_134316.root':['Run 134316','BlkStar'],
    'Laser_134318.root':['Run 134318','BlkStar'],
    'Laser_134327.root':['Run 134327','BlkStar'],
    'Laser_134349.root':['Run 134349','BlkStar'],
    'Laser_134357.root':['Run 134357','BlkStar'],
    'Laser_134365.root':['Run 134365','BlkStar'],
    'Laser_134371.root':['Run 134371','BlkStar'],
    'Laser_134444.root':['Run 134444','MagnTri'],
    'Laser_134447.root':['Run 134447','MagnTri'],
    'Laser_134449.root':['Run 134449','MagnTri'],
    'Laser_134450.root':['Run 134450','MagnTri'],
    'Laser_134452.root':['Run 134452','MagnTri'],
    'Laser_134453.root':['Run 134453','MagnTri'],
    'Laser_134456.root':['Run 134456','MagnTri'],
    'Laser_134457.root':['Run 134457','MagnTri'],
    'Laser_134458.root':['Run 134458','MagnTri'],
    'Laser_134459.root':['Run 134459','MagnTri'],
    'Laser_134463.root':['Run 134463','MagnTri'],
    'Laser_134466.root':['Run 134466','MagnTri'],
    'Laser_134469.root':['Run 134469','MagnTri'],
    'Laser_134470.root':['Run 134470','MagnTri'],
    'Laser_134472.root':['Run 134472','MagnTri'],
    'Laser_134473.root':['Run 134473','MagnTri']
    }

#import os
#os.system('cat LaserMay2010.C | sed "" | sed "" > nothing.C')
#I should do this in python, but I will do it in shell as I don't have pyroot installed.
