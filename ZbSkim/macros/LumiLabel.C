
TLatex * CMSPrel2 (Float_t Lumi, TString _decaychannel, int pos, double x, double y) {

  _decaychannel = _decaychannel + "";

  TLatex *latexLabel = new TLatex();
  latexLabel->SetTextFont(43);
  latexLabel->SetTextSize(16);
  latexLabel->SetLineWidth(2);
  latexLabel->SetNDC();

  latexLabel->DrawLatex(x,y,_decaychannel);
  latexLabel->DrawLatex(x,y-0.04,"anti-k_{T} (R = 0.5) jets");
  latexLabel->DrawLatex(x,y-0.08,"p_{T}^{jet} > 30 GeV, |#eta^{jet}| < 2.4");

  TLatex *latexLabel2 = new TLatex();
  latexLabel2->SetTextFont(43);
  latexLabel2->SetLineWidth(3);
  latexLabel2->SetNDC();

  latexLabel2->SetTextSize(20);
  latexLabel2->DrawLatex(0.735,0.945,Form("%.1f fb^{-1} (8 TeV)", Lumi));

  TLatex *latexLabel3 = new TLatex();
  latexLabel3->SetTextFont(43);
  latexLabel3->SetLineWidth(3);
  latexLabel3->SetNDC();

  latexLabel3->SetTextSize(20);
  if (pos==0)  latexLabel3->DrawLatex(0.125,0.89,"CMS Preliminary");
  if (pos==1)  latexLabel3->DrawLatex(0.1,0.94,"CMS Preliminary");

  return latexLabel;
}

TLatex * CMSPrel2New (TString _decaychannel, double x, double y) {

  _decaychannel = _decaychannel + "";

  TLatex *latexLabel = new TLatex();
  latexLabel->SetTextFont(43);
  latexLabel->SetTextSize(16);
  latexLabel->SetLineWidth(2);
  latexLabel->SetNDC();

  latexLabel->DrawLatex(x,y,_decaychannel);
  latexLabel->DrawLatex(x,y-0.04,"anti-k_{T} (R = 0.5) jets");
  latexLabel->DrawLatex(x,y-0.08,"p_{T}^{jet} > 30 GeV, |#eta^{jet}| < 2.4");

  return latexLabel;
}

