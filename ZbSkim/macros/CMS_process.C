#include "TPad.h"
#include "TLatex.h"

void CMS_process (TPad* pad, TString decaychannel, double x, double y) {

  TLatex *latexLabel = new TLatex();
  latexLabel->SetTextFont(43);
  latexLabel->SetTextSize(16);
  latexLabel->SetLineWidth(2);
  latexLabel->SetNDC();

  latexLabel->DrawLatex(x,y,decaychannel);
  latexLabel->DrawLatex(x,y-0.04,"anti-k_{T} (R = 0.5) jets");
  latexLabel->DrawLatex(x,y-0.08,"p_{T}^{jet} > 30 GeV, |#eta^{jet}| < 2.4");

  pad->cd();
  latexLabel->Draw("same");

  return;
}

