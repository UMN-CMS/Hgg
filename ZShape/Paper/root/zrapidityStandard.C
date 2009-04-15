#include "TText.h"
#include "TDatime.h"

void zrap_Prelim(double px, double py, double tx=-1, double ty=-1) {
  const char* time_;
  TDatime mytime;
  time_ = mytime.AsString();

  TText *plabel = new TText();
  plabel-> SetNDC();
  plabel -> SetTextFont(1);
  plabel -> SetTextColor(1);
  plabel -> SetTextSize(0.04);
  plabel -> SetTextAlign(22);
  plabel -> SetTextAngle(0);

  
  //Then for each plot, pick a nice spot and draw
  plabel -> DrawText(px, py, "PRELIMINARY");
  if (tx>=0 && ty>=0) {
    TText *tlabel = new TText();
    tlabel-> SetNDC();
    tlabel -> SetTextFont(1);
    tlabel -> SetTextColor(1);
    tlabel -> SetTextSize(0.02);
    tlabel -> SetTextAlign(22);
    tlabel -> SetTextAngle(0);
    
    tlabel -> DrawText(tx, ty, Form("%s",time_));
  }
}
