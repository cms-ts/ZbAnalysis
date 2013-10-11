#include "LumiLabel.C"
#include "LumiInfo_v11.h"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";

TH1F* h_data = 0;
TH1F* h_data_fit = 0;

void fcn(int& npar, double* gin, double& fun, double* par, int iflag) {
  double chisq = 0.0;
  for (int i=1; i<=h_data->GetNbinsX(); i++) {
    double xn = h_data->GetBinContent(i);
    double xd = h_data->GetBinError(i)**2;
    if (npar>0) {
      xn = xn - par[0]*h_data_fit->GetBinContent(i);
      xd = xd + (par[0]*h_data_fit->GetBinError(i))**2;
    }
    if (xd!=0) chisq = chisq + (xn*xn)/xd;
  }
  fun = chisq;
}

void DataMCComp5(string& title="", int plot=0, int ilepton=1, int doFit=0) {

//int useFitResults=0; // use MC predictions for c_t
int useFitResults=1;  // use fit results for c_t

      double c_t=1.0;

      if (ilepton==1 && useFitResults) c_t = 1.014;
      if (ilepton==2 && useFitResults) c_t = 1.024;

      double Lumi2012;
      
      if (ilepton==1) Lumi2012 = Lumi2012_ele;
      if (ilepton==2) Lumi2012 = Lumi2012_muon;

      double norm1 = ((Lumi2012 * Xsec_dy) / Ngen_dy);
      double norm2 = ((Lumi2012 * Xsec_tt) / Ngen_tt);
      double norm3 = ((Lumi2012 * Xsec_zz) / Ngen_zz);
      double norm4 = ((Lumi2012 * Xsec_wz) / Ngen_wz);
      double norm5 = ((Lumi2012 * Xsec_qcd) / Ngen_qcd);
      double norm6 = ((Lumi2012 * Xsec_ww) / Ngen_ww);
      double norm7 = ((Lumi2012 * Xsec_wj) / Ngen_wj);

      double norm1_fit = ((Lumi2012_ele_muon * Xsec_dy) / Ngen_dy);
      double norm2_fit = ((Lumi2012_ele_muon * Xsec_tt) / Ngen_tt);
      double norm3_fit = ((Lumi2012_ele_muon * Xsec_zz) / Ngen_zz);
      double norm4_fit = ((Lumi2012_ele_muon * Xsec_wz) / Ngen_wz);
      double norm5_fit = ((Lumi2012_ele_muon * Xsec_qcd) / Ngen_qcd);
      double norm6_fit = ((Lumi2012_ele_muon * Xsec_ww) / Ngen_ww);
      double norm7_fit = ((Lumi2012_ele_muon * Xsec_wj) / Ngen_wj);

      if (title.empty()) title = "w_jetmultiplicity";

      if (ilepton==1) {
        if (title.find("muon")!=string::npos) return;
	if (title.find("mm")!=string::npos) return;
      }
      if (ilepton==2) {
	if (title.find("ele")!=string::npos) return;
	if (title.find("ee")!=string::npos) return;
      }

      TFile *data;
      if (ilepton==1) data = TFile::Open((path + "/" + version + "/" + "DoubleElectron_2012_merge.root").c_str());
      if (ilepton==2) data = TFile::Open((path + "/" + version + "/" + "DoubleMu_2012_merge.root").c_str());

      TFile *data_fit;
      data_fit = TFile::Open((path + "/" + version + "/" + "MuEG_2012_merge.root").c_str());

      TFile *mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL.root").c_str());
      TFile *mc2 = TFile::Open((path + "/" + version + "/" + "TTbar.root").c_str());
      TFile *mc3 = TFile::Open((path + "/" + version + "/" + "ZZ.root").c_str());
      TFile *mc4 = TFile::Open((path + "/" + version + "/" + "WZ.root").c_str());
//    TFile *mc5 = TFile::Open((path + "/" + version + "/" + "QCD.root").c_str());
      TFile *mc6 = TFile::Open((path + "/" + version + "/" + "WW.root").c_str());
      TFile *mc7 = TFile::Open((path + "/" + version + "/" + "Wj.root").c_str());

      string title_fit = title;

      if (title_fit.find("ee")!=string::npos) {
        title_fit.replace(title_fit.find("ee")+1, 1, "m");
      }
      if (title_fit.find("mm")!=string::npos) {
        title_fit.replace(title_fit.find("mm"), 1, "e");
      }

      if (ilepton==1) data->cd("demoEle");
      if (ilepton==2) data->cd("demoMuo");
      h_data = (TH1F*)gDirectory->Get(title.c_str());
      data_fit->cd("demoEleMuo");
      h_data_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc1->cd("demoEle");
      if (ilepton==2) mc1->cd("demoMuo");
      TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
      mc1->cd("demoEleMuo");
      TH1F* h_mc1_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc2->cd("demoEle");
      if (ilepton==2) mc2->cd("demoMuo");
      TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
      mc2->cd("demoEleMuo");
      TH1F* h_mc2_fit = (TH1F*)gDirectory->Get(title_fit.c_str());
      
      if (ilepton==1) mc3->cd("demoEle");
      if (ilepton==2) mc3->cd("demoMuo");
      TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());
      mc3->cd("demoEleMuo");
      TH1F* h_mc3_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc4->cd("demoEle");
      if (ilepton==2) mc4->cd("demoMuo");
      TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());
      mc4->cd("demoEleMuo");
      TH1F* h_mc4_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

//    if (ilepton==1) mc5->cd("demoEle");
//    if (ilepton==2) mc5->cd("demoMuo");
//    TH1F* h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
//    mc5->cd("demoEleMuo");
//    TH1F* h_mc5_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc6->cd("demoEle");
      if (ilepton==2) mc6->cd("demoMuo");
      TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());
      mc6->cd("demoEleMuo");
      TH1F* h_mc6_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc7->cd("demoEle");
      if (ilepton==2) mc7->cd("demoMuo");
      TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());
      mc7->cd("demoEleMuo");
      TH1F* h_mc7_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      h_data->Sumw2();
      h_data_fit->Sumw2();

      h_mc1->Sumw2();
      h_mc2->Sumw2();
      h_mc3->Sumw2();
      h_mc4->Sumw2();
//    h_mc5->Sumw2();
      h_mc6->Sumw2();
      h_mc7->Sumw2();
      
      h_mc1_fit->Sumw2();
      h_mc2_fit->Sumw2();
      h_mc3_fit->Sumw2();
      h_mc4_fit->Sumw2();
//    h_mc5_fit->Sumw2();
      h_mc6_fit->Sumw2();
      h_mc7_fit->Sumw2();

      h_mc1->Scale(norm1);
      h_mc2->Scale(norm2*c_t);
      h_mc3->Scale(norm3);
      h_mc4->Scale(norm4);
//    h_mc5->Scale(norm5);
      h_mc6->Scale(norm6);
      h_mc7->Scale(norm7);
      
      h_mc1_fit->Scale(norm1_fit);
      h_mc2_fit->Scale(norm2_fit*c_t);
      h_mc3_fit->Scale(norm3_fit);
      h_mc4_fit->Scale(norm4_fit);
//    h_mc5_fit->Scale(norm5_fit);
      h_mc6_fit->Scale(norm6_fit);
      h_mc7_fit->Scale(norm7_fit);

      h_data->Add(h_mc7, -1.);
      h_data->Add(h_mc6, -1.);
//    h_data->Add(h_mc5, -1.);
      h_data->Add(h_mc4, -1.);
      h_data->Add(h_mc3, -1.);
      h_data->Add(h_mc1, -1.);

      h_data_fit->Add(h_mc7_fit, -1.);
      h_data_fit->Add(h_mc6_fit, -1.);
//    h_data_fit->Add(h_mc5_fit, -1.);
      h_data_fit->Add(h_mc4_fit, -1.);
      h_data_fit->Add(h_mc3_fit, -1.);
      h_data_fit->Add(h_mc1_fit, -1.);

      TH1F* h_data_fit_raw = (TH1F*)h_data_fit->Clone();

      h_data_fit->Scale(Lumi2012/Lumi2012_ele_muon);

      if (doFit==1) {
        for (int i=0; i<=h_data_fit->GetNbinsX()+1; i++) {
	  bool skip = false;
          if (title.find("w_mass")!=string::npos) {
            if (h_data_fit->GetXaxis()->GetBinCenter(i)>85 && h_data_fit->GetXaxis()->GetBinCenter(i)<100) {
	      skip = true;
            }
	  }
          if (title=="w_MET") {
            if (h_data_fit->GetXaxis()->GetBinCenter(i)<125) {
	      skip = true;
            }
	  }
          if (title=="w_MET_sign") {
            if (h_data_fit->GetXaxis()->GetBinCenter(i)<50) {
	      skip = true;
            }
	  }
          if (title=="w_MET_b") {
            if (h_data_fit->GetXaxis()->GetBinCenter(i)<90) {
	      skip = true;
            }
          }
          if (title=="w_MET_sign_b") {
            if (h_data_fit->GetXaxis()->GetBinCenter(i)<30) {
	      skip = true;
            }
          }
	  if (skip) {
  	    h_data->SetBinContent(i, 0);
	    h_data->SetBinError(i, 0);
	    h_data_fit->SetBinContent(i, 0);
	    h_data_fit->SetBinError(i, 0);
	    h_mc2->SetBinContent(i, 0);
	    h_mc2->SetBinError(i, 0);
	  }
        }
	TVirtualFitter* fitter = TVirtualFitter::Fitter(0, 1);
	fitter->SetFCN(fcn);
	double arglist[1] = {-1.0};
	fitter->ExecuteCommand("SET PRINT", arglist, 1);
	fitter->SetParameter(0, "c(t)", 0.5, 0.1, 0.0, 1.0);
	fitter->ExecuteCommand("MIGRAD", arglist, 0);
	h_data_fit->Scale(fitter->GetParameter(0));
      }

      TCanvas* c1;
      if (doFit) {
        c1 = new TCanvas("c", "c", 800, 600);
        c1->cd();

	TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
	pad1->SetBottomMargin(0.001);
	pad1->Draw();
	pad1->cd();

        h_data->SetMinimum(0);
        h_data->SetMaximum(-1111);

        h_data->Draw("EP");
        h_data_fit->Draw("HISTSAME");
        h_data->SetMarkerColor(kBlack);
        h_data_fit->SetMarkerColor(kRed);
        h_data_fit->SetLineColor(kRed);
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize (1.0);
        h_data->SetTitle("");
        //h_data->SetStats(0);

        h_mc2->Draw("HISTSAME");
        h_mc2->SetLineColor(kBlue);

	pad1->Update();
	c1->Update();

	c1->cd();

	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();

	TH1F* h_ratio = h_data->Clone("h_ratio");
	h_ratio->Divide(h_data_fit);

	h_ratio->SetTitle("");
        h_ratio->SetStats(0);

	h_ratio->GetXaxis()->SetTitleOffset(0.9);
	h_ratio->GetXaxis()->SetTitleSize(0.1);
        h_ratio->GetXaxis()->SetLabelFont(42);
        h_ratio->GetXaxis()->SetLabelSize(0.08);
        h_ratio->GetXaxis()->SetTitleFont(42);
        h_ratio->GetYaxis()->SetTitle("ratio");
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.09);
        h_ratio->GetYaxis()->SetLabelSize(0.08);
        h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
        h_ratio->GetYaxis()->SetTitleOffset(0.4);

        h_ratio->SetMarkerStyle(20);
        h_ratio->SetMarkerColor(kBlack);

        h_ratio->Draw("EPX");

        TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
        OLine->SetLineColor(kGreen);
        OLine->SetLineWidth(1);
        OLine->Draw();

        c1->cd();

        TLatex *latexLabel = CMSPrel(Lumi2012/1000.,"",0.15,0.94);
        latexLabel->Draw("same");

        TLatex *fitLabel = new TLatex();
        fitLabel->SetTextSize(0.0275);
        fitLabel->SetTextFont(42);
        fitLabel->SetLineWidth(2);
        fitLabel->SetNDC();
        char buff[100];
        sprintf(buff, "c_{t} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
        fitLabel->DrawLatex(0.68, 0.68, buff);
        if (ilepton==1) sprintf(buff, "I_{ee} = %5.1f", h_data->Integral());
        if (ilepton==2) sprintf(buff, "I_{#mu#mu} = %5.1f", h_data->Integral());
        fitLabel->DrawLatex(0.68, 0.63, buff);
        sprintf(buff, "I_{e#mu} = %5.1f #pm %5.1f", h_data_fit->Integral(), h_data_fit->Integral()*fitter->GetParError(0)/fitter->GetParameter(0));
        fitLabel->DrawLatex(0.68, 0.58, buff);
        sprintf(buff, "I_{t#bar{t}}  = %5.1f", h_mc2->Integral());
        fitLabel->DrawLatex(0.68, 0.53, buff);
      }

      if (plot) {
        if (doFit) title = title + "_doFit";
        if (ilepton==1) {
          gSystem->mkdir((path + "/electrons/" + version + "/ttbar_sub/").c_str(), kTRUE);
          if (c1) c1->SaveAs((path + "/electrons/" + version + "/ttbar_sub/" + title + ".pdf").c_str());
	  if (doFit==0) {
	    TFile f((path + "/electrons/" + version + "/ttbar_sub/" + title + ".root").c_str(),"RECREATE");
            h_data_fit_raw->Write(title.c_str());
            f.Close();
	  }
        }
        if (ilepton==2) {
          gSystem->mkdir((path + "/muons/" + version + "/ttbar_sub/").c_str(), kTRUE);
          if (c1) c1->SaveAs((path + "/muons/" + version + "/ttbar_sub/" + title + ".pdf").c_str());
	  if (doFit==0) {
	    TFile f((path + "/muons/" + version + "/ttbar_sub/" + title + ".root").c_str(),"RECREATE");
            h_data_fit_raw->Write(title.c_str());
            f.Close();
	  }
        }
      }
}
