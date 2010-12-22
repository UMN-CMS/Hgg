
#include <math.h>
 
#include "ZShape/EffAcc/interface/EfficiencyPlots.h"
#include "ZShape/Base/interface/EfficiencyStatistics.h"
#include <TGraphAsymmErrors.h>
#include <time.h>
#include <TText.h>
#include <TROOT.h>
#include <TPaveLabel.h>
#include <TStyle.h>


EfficiencyPlots::EfficiencyPlots(char *txtFile, char *cutType, int binStart, char *binsFile, int effDom) :
  txtFile_(txtFile),
  cutType_(cutType),
  binStart_(binStart),
  binsFile_(binsFile),
  effDom_(effDom)
{
   trueeff = "Supercluster Creation";
   if (strstr(cutType_,"Gsf") ) trueeff = "GSF Track Matching";
   else if (strstr(cutType_,"Iso") ) {trueeff = "Isolation WP95"; if ( strstr(cutType_,"80")) trueeff = "Isolation WP80"; }
   else if (strstr(cutType_,"ElectronId") ) {trueeff = "Electron Identification WP95";if (strstr(cutType_,"80")) trueeff = "Electron Identification WP80";}
   else if (strstr(cutType_,"HLT") ) trueeff = "Trigger Efficiency";
   else if (strstr(cutType_,"HFEIdT") ) trueeff = "HF Tight Electron Identification";
   else if (strstr(cutType_,"HF") ) trueeff = "HF Loose Electron Identification";
   else if (strstr(cutType_,"WP95") ) trueeff = "WP95";
   else if (strstr(cutType_,"WP90") ) trueeff = "WP90";
   else if (strstr(cutType_,"WP85") ) trueeff = "WP85";
   else if (strstr(cutType_,"WP80") ) trueeff = "WP80";
   else if (strstr(cutType_,"WP70") ) trueeff = "WP70";
   else if (strstr(cutType_,"WP60") ) trueeff = "WP60";
   else if (strstr(cutType_,"NTEIdT") ) trueeff = "NT Tight Electron Identification";
   else if (strstr(cutType_,"NT") ) trueeff = "NT Loose Electron Identification";

   //File we initialize the strings for the efficiencies
   Effs_[0] = "GsfTrack"; 
   Effs_[1] = "Isolation";
   Effs_[2] = "ElectronId";
   Effs_[3] = "HLT";
   Effs_[4] = "HFElectronId";

  //Next we initialize the colors
   Cols_[0] = kBlack;
   Cols_[1] = kRed;
   Cols_[2] = kMagenta;
   Cols_[3] = kBlue;
   Cols_[4] = kCyan+1;
   Cols_[5] = kGreen;
   Cols_[6] = kYellow+2;

   Marker_[0] = kOpenCircle;
   Marker_[1] = kOpenTriangleUp;
   Marker_[2] = kOpenSquare;
   Marker_[3] = kOpenStar;
   Marker_[4] = kOpenDiamond;
   Marker_[5] = kOpenCross;
   Marker_[6] = kPlus;

   // TStyle STuff
   mystyle_ = new TStyle("mystyle","mystyle");
   // For the axis:
// For the canvas:
  mystyle_->SetCanvasBorderMode(0);
  mystyle_->SetCanvasColor(kWhite);
  mystyle_->SetCanvasDefH(600); //Height of canvas
  mystyle_->SetCanvasDefW(600); //Width of canvas
  mystyle_->SetCanvasDefX(0);   //POsition on screen
  mystyle_->SetCanvasDefY(0);

// For the Pad:
  mystyle_->SetPadBorderMode(0);
  // mystyle_->SetPadBorderSize(Width_t size = 1);
  mystyle_->SetPadColor(kWhite);
  mystyle_->SetPadGridX(false);
  mystyle_->SetPadGridY(false);
  mystyle_->SetGridColor(0);
  mystyle_->SetGridStyle(3);
  mystyle_->SetGridWidth(1);

// For the frame:
  mystyle_->SetFrameBorderMode(0);
  mystyle_->SetFrameBorderSize(1);
  mystyle_->SetFrameFillColor(0);
  mystyle_->SetFrameFillStyle(0);
  mystyle_->SetFrameLineColor(1);
  mystyle_->SetFrameLineStyle(1);
  mystyle_->SetFrameLineWidth(1);

  mystyle_->SetAxisColor(1, "XYZ");
  mystyle_->SetStripDecimals(kTRUE);
  mystyle_->SetTickLength(0.03, "XYZ");
  mystyle_->SetNdivisions(510, "XYZ");
  mystyle_->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  mystyle_->SetPadTickY(1);
// For the axis titles:

  mystyle_->SetTitleColor(1, "XYZ");
  mystyle_->SetTitleFont(42, "XYZ");
  mystyle_->SetTitleSize(0.06, "XYZ");
  // mystyle_->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // mystyle_->SetTitleYSize(Float_t size = 0.02);
  mystyle_->SetTitleXOffset(0.9);
  mystyle_->SetTitleYOffset(1.25);
  // mystyle_->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
// Margins:
  mystyle_->SetPadTopMargin(0.05);
  mystyle_->SetPadBottomMargin(0.13);
  mystyle_->SetPadLeftMargin(0.16);
  mystyle_->SetPadRightMargin(0.02);
// For the statistics box:
  mystyle_->SetOptFile(0);
  mystyle_->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  mystyle_->SetStatColor(kWhite);
  mystyle_->SetStatFont(42);
  mystyle_->SetStatFontSize(0.025);
  mystyle_->SetStatTextColor(1);
  mystyle_->SetStatFormat("6.4g");
  mystyle_->SetStatBorderSize(1);
  mystyle_->SetStatH(0.1);
  mystyle_->SetStatW(0.15);
  mystyle_->SetPaperSize(20.,20.);

// For the Global title:

  mystyle_->SetOptTitle(0);
  mystyle_->SetTitleFont(42);
  mystyle_->SetTitleColor(1);
  mystyle_->SetTitleTextColor(1);
  mystyle_->SetTitleFillColor(10);
  mystyle_->SetTitleFontSize(0.05);

  mystyle_->SetLabelColor(1, "XYZ");
  mystyle_->SetLabelFont(42, "XYZ");
  mystyle_->SetLabelOffset(0.007, "XYZ");
  mystyle_->SetLabelSize(0.05, "XYZ");

  mystyle_->cd();


   //Cout some of the information for the record

   for (iteffstr it = Effs_.begin(); it != Effs_.end(); it++)
   {
      int type = (*it).first;
      // rtype = (*it).second;
      std::cout <<" Working on " << type << " " << (*it).second <<  std::endl;
   }
   
   efftable_ = new EffTableLoader(txtFile_); //This will automatically load in the txt file.
   binsTable_ = new EffTableLoader(binsFile_);
   numbins_ = efftable_->size();

    //Now it is time to get the TIME
   time_t rawtime;
   struct tm * timeinfo;
   time (&rawtime);
   timeinfo = localtime (&rawtime);
   time_ = asctime (timeinfo);
   std::cout << " The time is " << time_ << std::endl;

   //gROOT->SetStyle("Plain");

   physTypeVar_ = ( effDom ) ? ("P_{T}") : ("#eta_{d}");
   physTypeDom_ = ( effDom ) ? ("#eta_{d}") : ("P_{T}");

   physUnitsVar_ = ( effDom ) ? ("[GeV/c]") : ("");
   physUnitsDom_ = ( effDom ) ? ("") : ("[GeV/c]");

   vareff_ = (effDom) ? 0 : 1;
}

void
EfficiencyPlots::makeHistos(void)
{
        mystyle_->cd();
	float min=-50.;
	int binnumber = 0;
        std::cout << " size of table " << binsTable_->size() << std::endl;
	for (int b=0; b<binsTable_->size(); ++b) 
	  {
           std::vector<std::pair<float, float> > BinInfo = binsTable_->GetCellInfo(b);
	   std::vector<float> BinVals = binsTable_->correctionEff(b);
	   std::cout << " inside b " << b << std::endl;
	   TGraphAsymmErrors *mytgraph = new TGraphAsymmErrors(int (BinVals.size()) -1);
	   mytgraph->SetName(Form("sbs_eff_%s_bin_%d",cutType_, b));
	   mytgraph->SetTitle(Form("sbs_eff_%s_bin_%d;%s %s;%s Efficiency",cutType_, b, physTypeDom_.c_str() , physUnitsDom_.c_str(), cutType_));
           mytgraph->GetXaxis()->CenterTitle(kTRUE);
           mytgraph->GetYaxis()->CenterTitle(kTRUE);
           std::cout << " size of bins " << BinVals.size() << std::endl;
           for (int r=1; r< int (BinVals.size()); ++r) 
	   {
               std::cout << " inside r " << r << std::endl;
	      float mybin = BinVals[r];
	      int mybinint = int (mybin);
              std::vector<std::pair<float, float> > EffInfo = efftable_->GetCellInfo(mybinint);
	      std::vector<float> EffVals = efftable_->correctionEff(mybinint);
	      //effhists[b]->Fill(r, EffVals[0] );
	      float xcenter = 0.5*(EffInfo[effDom_].first + EffInfo[effDom_].second );
	      float xerror = fabs(EffInfo[effDom_].second - xcenter);
	      double eff = double((EffVals[1]));
	      int den = int ((EffVals[4]));
	      EfficiencyStatistics effstats(eff,den);
              double sigmaminus = fabs(effstats.sigma(-1) - eff);
              double sigmaplus  = fabs(effstats.sigma(1) - eff);
	      mytgraph->SetPoint(r-1,xcenter,eff);
	      mytgraph->SetPointError(r-1,xerror,xerror,sigmaminus, sigmaplus);

              std::cout << " Eff is " << eff << " den is " << den << " stat- " << effstats.sigma(-1) << " stat+ " << effstats.sigma(1) << std::endl;
              std::cout << " xcetner " << xcenter << " xerror " << xerror << std::endl;
	   }
           std::cout << " done with r loop" << std::endl;
	   effgraphs_.push_back(mytgraph);
       }
}



void
EfficiencyPlots::printIndividualHistos(const char *ftype, bool withcolor)
{
        mystyle_->cd();
	for (int b=0; b<binsTable_->size(); ++b) 
	{
           std::vector<std::pair<float, float> > BinInfo = binsTable_->GetCellInfo(b);
	   TCanvas *tempCan = new TCanvas("tempcan","tempcan");
           TLegend *myLeg = new TLegend(.67,.12,.88,.38);
           if (withcolor) effgraphs_[b]->SetLineColor(Cols_[b]);
           if (withcolor) effgraphs_[b]->SetMarkerColor(Cols_[b]);
           effgraphs_[b]->SetMarkerStyle(Marker_[b]);
           effgraphs_[b]->GetXaxis()->CenterTitle(kTRUE);
           effgraphs_[b]->GetYaxis()->CenterTitle(kTRUE);
           effgraphs_[b]->Draw("AP");
           myLeg->AddEntry(effgraphs_[b],Form("%s %.2f to %.2f %s",
                            physTypeVar_.c_str(), BinInfo[vareff_].first,BinInfo[vareff_].second,physUnitsVar_.c_str()),"lp");
           myLeg->SetFillColor(kWhite);
           if(binsTable_->size() > 1) myLeg->Draw();
           TText *tlabel = new TText();
           tlabel-> SetNDC();
           tlabel -> SetTextFont(1);
           tlabel -> SetTextColor(1);
           tlabel -> SetTextSize(0.03);
           tlabel -> SetTextAlign(22);
           tlabel -> SetTextAngle(0);
           tlabel -> DrawText(0.5, 0.17, Form("PRELIMINARY: %s",time_));
           

           tempCan->Print(Form("sbs_eff_%s_%d.%s",cutType_,b,ftype));
	} 

}


void
EfficiencyPlots::printSumHistos(const char *ftype)
{
        mystyle_->cd();
	TCanvas *tempCan = new TCanvas("tempcan","tempcan");
        TLegend *myLeg = new TLegend(.42,.13,.73,.39);
	for (int b=0; b<binsTable_->size(); ++b) 
	{
           std::vector<std::pair<float, float> > BinInfo = binsTable_->GetCellInfo(b);
	   effgraphs_[b]->GetXaxis()->CenterTitle(kTRUE);
           effgraphs_[b]->GetYaxis()->CenterTitle(kTRUE);
           effgraphs_[b]->SetLineColor(Cols_[b]);
           effgraphs_[b]->SetMarkerColor(Cols_[b]);
           effgraphs_[b]->SetMarkerStyle(Marker_[b]);
           if ( b == 0) 
           {
               effgraphs_[b]->Draw("AP");
               effgraphs_[b]->SetTitle(Form(";%s %s;Single Electron Efficiency",physTypeDom_.c_str() , physUnitsDom_.c_str()));
               double min = effgraphs_[b]->GetYaxis()->GetXmin();
               double max = effgraphs_[b]->GetYaxis()->GetXmax();
	       double dmin = min-0.15;
	       double dmax = 1.02;
	       if ( max - min < 0.05 ) {dmin = min - 0.02; dmax = 1.01;}
               if (max > 0.98 || binsTable_->size() > 1 ) effgraphs_[b]->SetMaximum(dmax);
               if (binsTable_->size() > 1) effgraphs_[b]->SetMinimum(dmin);
               else if (min < 0.2) effgraphs_[b]->SetMinimum(0.70);
           }
           else effgraphs_[b]->Draw("P");
           myLeg->AddEntry(effgraphs_[b],Form("%s %.2f to %.2f %s",
                            physTypeVar_.c_str(), BinInfo[vareff_].first,BinInfo[vareff_].second,physUnitsVar_.c_str()),"lp");
	}
        
	myLeg->SetFillColor(kWhite);
        if(binsTable_->size() > 1) myLeg->Draw();
		TText *plabel = new TText();
        plabel-> SetNDC();
        plabel -> SetTextFont(1);
        plabel -> SetTextColor(1);
        plabel -> SetTextSize(0.04);
        plabel -> SetTextAlign(22);
        plabel -> SetTextAngle(0);
        plabel -> DrawText(0.29, 0.20, "PRELIMINARY");
        TText *tlabel = new TText();
        tlabel-> SetNDC();
        tlabel -> SetTextFont(1);
        tlabel -> SetTextColor(1);
        tlabel -> SetTextSize(0.02);
        tlabel -> SetTextAlign(22);
        tlabel -> SetTextAngle(0);
        tlabel -> DrawText(0.29, 0.17, Form("%s",time_));

        TPaveLabel *tpave = new TPaveLabel(0.50,0.89,0.98,0.95,Form("%s",trueeff.c_str()),"NDC");
        tpave->SetTextSize(0.6);
        tpave->SetFillColor(kWhite);
        tpave->Draw();
        tempCan->Print(Form("sbs_eff_%s.%s",cutType_,ftype));
   
 
}

EfficiencyPlots::~EfficiencyPlots()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

