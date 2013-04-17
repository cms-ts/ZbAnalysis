TLatex *
CMSPrel (Float_t Lumi, TString _decaychannel)
{

  TLatex *latexLabel = new TLatex ();
  latexLabel->SetTextSize (0.0275);
  latexLabel->SetTextFont (42);
  latexLabel->SetLineWidth (2);
  latexLabel->SetNDC ();

  latexLabel->DrawLatex (0.25, 0.35, "CMS Preliminary");
  latexLabel->DrawLatex (0.25, 0.30,
			 "#sqrt{s} = 8 TeV  #int Ldt = 19.3 fb^{-1}");
  latexLabel->DrawLatex (0.25, 0.25, "anti-k_{T} (R = 0.5) PF Jets > 30 GeV");
  //latexLabel->DrawLatex(0.25,0.18,(TString)Form("#int Ldt = %.3f fb^{-1}",Lumi));  
  //latexLabel->DrawLatex(0.25,0.18,"#int Ldt = 4.890 fb^{-1}");  
  //latexLabel->DrawLatex(0.25,0.13,"Z#rightarrow ee channel");
  latexLabel->DrawLatex (0.25, 0.08, _decaychannel);


  //return latexLabel;
}

TLatex *
CMSPrel (Float_t Lumi, TString _decaychannel, double x, double y)
{

  TLatex *latexLabel = new TLatex ();
  latexLabel->SetTextSize (0.0275);
  latexLabel->SetTextFont (42);
  latexLabel->SetLineWidth (2);
  latexLabel->SetNDC ();

  latexLabel->DrawLatex (x, y, "CMS Preliminary            #sqrt{s} = 8 TeV #int Ldt = 19.3 fb^{-1}");
  //latexLabel->DrawLatex (x, y - 0.04,
  //			 "#sqrt{s} = 8 TeV #int Ldt = 19.3 fb^{-1}");
  //latexLabel->DrawLatex (x, y - 0.08,
  //			 "anti-k_{T} (R = 0.5) PF Jets > 30 GeV ");
  //latexLabel->DrawLatex(x,0.18,(TString)Form("#int Ldt = %.3f fb^{-1}",Lumi));
  //latexLabel->DrawLatex(x,y-0.18,"Z#rightarrow ee channel");
  latexLabel->DrawLatex (x, y - 0.12, _decaychannel);


  //return latexLabel;
}
