TLatex * CMSPrel (Float_t Lumi, TString _decaychannel) {

  Lumi = 1.0*Lumi;
  _decaychannel = _decaychannel + "";

  TLatex * latexLabel = new TLatex ();
  latexLabel->SetTextSize (0.0275);
  latexLabel->SetTextFont (42);
  latexLabel->SetLineWidth (2);
  latexLabel->SetNDC ();

  latexLabel->DrawLatex (0.25, 0.35, "CMS Preliminary");
  //latexLabel->DrawLatex (0.25, 0.30, Form("#sqrt{s} = 8 TeV  #int Ldt = %.1f fb^{-1}", Lumi));
  latexLabel->DrawLatex (0.25, 0.30, "19.7 fb^{-1} (8 TeV)");
  latexLabel->DrawLatex (0.25, 0.25, "anti-k_{T} (R = 0.5) PF Jets > 30 GeV");
  latexLabel->DrawLatex (0.25, 0.08, _decaychannel);

  return latexLabel;
}

TLatex * CMSPrel (Float_t Lumi, TString _decaychannel, double x, double y) {

  Lumi = 1.0*Lumi;
  _decaychannel = _decaychannel + "";
  x = 1*x;
  y = 1*y;

  TLatex *latexLabel = new TLatex ();
  latexLabel->SetTextSize (0.0275);
  latexLabel->SetTextFont (42);
  latexLabel->SetLineWidth (2);
  latexLabel->SetNDC ();

  latexLabel->DrawLatex (x, y, Form("CMS Preliminary            #sqrt{s} = 8 TeV #int Ldt = %.1f fb^{-1}", Lumi));
  latexLabel->DrawLatex (x, y - 0.12, _decaychannel);

  return latexLabel;
}

TLatex * CMSPrel2 (Float_t Lumi, TString _decaychannel, int pos, double x, double y) {

  Lumi = 1.0*Lumi;
  _decaychannel = _decaychannel + "";
  x = 1*x;
  y = 1*y;

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
  latexLabel2->DrawLatex(0.735,0.945,"19.7 fb^{-1} (8 TeV)");

  TLatex *latexLabel3 = new TLatex();
  latexLabel3->SetTextFont(43);
  latexLabel3->SetLineWidth(3);
  latexLabel3->SetNDC();

  latexLabel3->SetTextSize(20);
  if (pos==0)  latexLabel3->DrawLatex(0.125,0.89,"CMS Preliminary");
  if (pos==1)  latexLabel3->DrawLatex(0.1,0.94,"CMS Preliminary");

  return latexLabel;
}

TLatex * CMSPrel3 (Float_t Lumi, TString _decaychannel, int pos, double x, double y) {

  Lumi = 1.0*Lumi;
  _decaychannel = _decaychannel + "";
  x = 1*x;
  y = 1*y;

  TLatex *latexLabel = new TLatex();

  latexLabel->SetTextFont(43);
  latexLabel->SetLineWidth(3);
  latexLabel->SetNDC();

  latexLabel->SetTextSize(20);
  latexLabel->DrawLatex(0.735,0.945,"19.7 fb^{-1} (8 TeV)");

  TLatex *latexLabel2 = new TLatex();
  latexLabel2->SetTextFont(43);
  latexLabel2->SetLineWidth(3);
  latexLabel2->SetNDC();

  latexLabel2->SetTextSize(20);
  if (pos==0)  latexLabel2->DrawLatex(0.125,0.89,"CMS Preliminary");
  if (pos==1)  latexLabel2->DrawLatex(0.1,0.94,"CMS Preliminary");

  return latexLabel;
}


TLatex * CMSFinal (Float_t Lumi, TString _decaychannel, int pos, double x, double y) {

  Lumi = 1.0*Lumi;
  _decaychannel = _decaychannel + "";
  x = 1*x;
  y = 1*y;

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
  latexLabel2->DrawLatex(0.735,0.945,"19.7 fb^{-1} (8 TeV)");
 
  TLatex *latexLabel3 = new TLatex();
  latexLabel3->SetTextFont(43);
  latexLabel3->SetLineWidth(3);
  latexLabel3->SetNDC();

  latexLabel3->SetTextSize(20);
  if (pos==0)  latexLabel3->DrawLatex(0.125,0.89,"CMS");
  if (pos==1)  latexLabel3->DrawLatex(0.1,0.94,"CMS"); 

  return latexLabel;
}

TLatex * CMSFinal2 (Float_t Lumi, TString _decaychannel, int pos, double x, double y) {

  Lumi = 1.0*Lumi;
  _decaychannel = _decaychannel + "";
  x = 1*x;
  y = 1*y;

  TLatex *latexLabel = new TLatex();

  latexLabel->SetTextFont(43);
  latexLabel->SetLineWidth(3);
  latexLabel->SetNDC();

  latexLabel->SetTextSize(20);
  latexLabel->DrawLatex(0.735,0.945,"19.7 fb^{-1} (8 TeV)");

  TLatex *latexLabel2 = new TLatex();
  latexLabel2->SetTextFont(43);
  latexLabel2->SetLineWidth(3);
  latexLabel2->SetNDC();

  latexLabel2->SetTextSize(20);
  if (pos==0)  latexLabel2->DrawLatex(0.125,0.89,"CMS");
  if (pos==1)  latexLabel2->DrawLatex(0.1,0.94,"CMS");  

  return latexLabel;
}

TLatex * CMSPrel2New (Float_t Lumi, TString _decaychannel, double x, double y) {

  Lumi = 1.0*Lumi;
  _decaychannel = _decaychannel + "";
  x = 1*x;
  y = 1*y;

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
