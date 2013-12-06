#include "TLatex.h"

TLatex * CMSPrel (Float_t Lumi, TString _decaychannel) {

  TLatex * latexLabel = new TLatex ();
  latexLabel->SetTextSize (0.0275);
  latexLabel->SetTextFont (42);
  latexLabel->SetLineWidth (2);
  latexLabel->SetNDC ();

//  latexLabel->DrawLatex (0.25, 0.35, "CMS Preliminary");
  latexLabel->DrawLatex (0.25, 0.30, Form("#sqrt{s} = 8 TeV  #int Ldt = %.1f fb^{-1}", Lumi));
  latexLabel->DrawLatex (0.25, 0.25, "anti-k_{T} (R = 0.5) PF Jets > 30 GeV");
  latexLabel->DrawLatex (0.25, 0.08, _decaychannel);

  return latexLabel;
}

TLatex * CMSPrel (Float_t Lumi, TString _decaychannel, double x, double y) {

  TLatex *latexLabel = new TLatex ();
  latexLabel->SetTextSize (0.0275);
  latexLabel->SetTextFont (42);
  latexLabel->SetLineWidth (2);
  latexLabel->SetNDC ();

//  latexLabel->DrawLatex (x, y, Form("CMS Preliminary            #sqrt{s} = 8 TeV #int Ldt = %.1f fb^{-1}", Lumi));
  latexLabel->DrawLatex (x, y, Form("                           #sqrt{s} = 8 TeV #int Ldt = %.1f fb^{-1}", Lumi));
  latexLabel->DrawLatex (x, y - 0.12, _decaychannel);

  return latexLabel;
}
