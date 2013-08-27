#include "LumiLabel.C"
#include "LumiInfo_v10.h"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";

TH1F* h_data = 0;
TH1F* h_data_fit = 0;

double func(double* x, double* p) {
  double val = 0.0;
  int i = h_data->GetXaxis()->FindBin(x[0]);
  if (h_data_fit) val = val + p[0]*h_data_fit->GetBinContent(i);
  return val;
}

void DataMCComp5(string& title="", int plot=0, int ilepton=1, int doFit=0) {

      double Lumi2012;
      
      if (ilepton==1) Lumi2012 = Lumi2012_ele;
      if (ilepton==2) Lumi2012 = Lumi2012_muon;

      double norm1 = ( (Lumi2012 * Xsec_dy) / Ngen_dy);
      double norm3 = ( (Lumi2012 * Xsec_zz) / Ngen_zz);
      double norm4 = ( (Lumi2012 * Xsec_wz) / Ngen_wz);
      double norm5 = ( (Lumi2012 * Xsec_qcd) / Ngen_qcd);
      double norm6 = ( (Lumi2012 * Xsec_ww) / Ngen_ww);
      double norm7 = ( (Lumi2012 * Xsec_wj) / Ngen_wj);

      if (title.empty()) title = "w_jetmultiplicity";

      if (ilepton==1) {
        if (title.find("muon")!=string::npos) return;
	if (title.find("mm")!=string::npos) return;
      }
      if (ilepton==2) {
	if (title.find("ele")!=string::npos) return;
	if (title.find("ee")!=string::npos) return;
      }

      TFile *data_ee = TFile::Open((path + "/" + version + "/" + "DoubleElectron_2012_merge.root").c_str());
      TFile *data_mm = TFile::Open((path + "/" + version + "/" + "DoubleMu_2012_merge.root").c_str());

      if (ilepton==1) data_ee->cd("demoEle");
      if (ilepton==2) data_mm->cd("demoMuo");
      h_data = (TH1F*)gDirectory->Get(title.c_str());

      string title_em = title;

      if (title_em.find("ee")!=string::npos) {
        title_em.replace(title_em.find("ee")+1, 1, "m");
      }
      if (title_em.find("mm")!=string::npos) {
        title_em.replace(title_em.find("mm"), 1, "e");
      }

      TFile *data_em = TFile::Open((path + "/" + version + "/" + "MuEG_2012_merge.root").c_str());

      TFile *mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL.root").c_str());
      TFile *mc3 = TFile::Open((path + "/" + version + "/" + "ZZ.root").c_str());
      TFile *mc4 = TFile::Open((path + "/" + version + "/" + "WZ.root").c_str());
//    TFile *mc5 = TFile::Open((path + "/" + version + "/" + "QCD.root").c_str());
      TFile *mc6 = TFile::Open((path + "/" + version + "/" + "WW.root").c_str());
      TFile *mc7 = TFile::Open((path + "/" + version + "/" + "Wj.root").c_str());

      data_em->cd("demoEleMuo");
      h_data_fit = (TH1F*)gDirectory->Get(title_em.c_str());

      h_data->Sumw2();
      h_data_fit->Sumw2();

      if (!doFit) {
        mc1->cd("demoEleMuo");
        TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());

        mc3->cd("demoEleMuo");
        TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());

        mc4->cd("demoEleMuo");
        TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());

//      mc5->cd("demoEleMuo");
//      TH1F* h_mc5 = (TH1F*)gDirectory->Get(title.c_str());

        mc6->cd("demoEleMuo");
        TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());

        mc7->cd("demoEleMuo");
        TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());

        h_mc1->Sumw2();
        h_mc3->Sumw2();
        h_mc4->Sumw2();
//      h_mc5->Sumw2();
        h_mc6->Sumw2();
        h_mc7->Sumw2();

        h_mc1->Scale(norm1);
        h_mc3->Scale(norm3);
        h_mc4->Scale(norm4);
//      h_mc5->Scale(norm5);
        h_mc6->Scale(norm6);
        h_mc7->Scale(norm7);

        h_data_fit->Add(h_mc7, -1.);
        h_data_fit->Add(h_mc6, -1.);
//      h_data_fit->Add(h_mc5, -1.);
        h_data_fit->Add(h_mc4, -1.);
        h_data_fit->Add(h_mc3, -1.);
      }

      TCanvas* c1;
      TF1 *f1 = new TF1("f1", func, 60, 120.00, 1);
      if (doFit) {
        for (int i=0; i<=h_data_fit->GetNbinsX()+1; i++) {
          if (h_data_fit->GetXaxis()->GetBinCenter(i)>=70 && h_data_fit->GetXaxis()->GetBinCenter(i)<=110) {
  	    h_data->SetBinContent(i, 0);
	    h_data->SetBinError(i, 0);
	    h_data_fit->SetBinContent(i, 0);
	    h_data_fit->SetBinError(i, 0);
          }
        }
        f1->SetParameter(0,0.5);
        f1->SetParNames("c(t)", "dummy", "dummy");
        h_data->Fit("f1", "Q0");
        h_data_fit->Scale(f1->GetParameter(0));

        c1 = new TCanvas("c", "c", 800, 600);
        c1->cd();

        h_data->SetMinimum(-1111);
        h_data->SetMaximum(-1111);

        h_data->Draw("EP");
        h_data_fit->Draw("HISTSAMES");
        h_data->SetMarkerColor(kBlack);
        h_data_fit->SetMarkerColor(kRed);
        h_data_fit->SetLineColor(kRed);
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize (1.0);
        h_data->SetTitle("");
        //h_data->SetStats(0);

        TLatex *latexLabel = CMSPrel(Lumi2012/1000.,"",0.15,0.94);
        latexLabel->Draw("same");
        TLatex *fitLabel = new TLatex();
        fitLabel->SetTextSize(0.0275);
        fitLabel->SetTextFont(42);
        fitLabel->SetLineWidth(2);
        fitLabel->SetNDC();
        char buff[100];
        sprintf(buff, "c_{t} = %5.3f #pm %5.3f", f1->GetParameter(0), f1->GetParError(0));
        fitLabel->DrawLatex(0.68, 0.58, buff);
        fitLabel->Draw("same");
      }

      if (plot) {
        if (doFit) title = title + "_doFit";
        if (ilepton==1) {
          gSystem->mkdir((path + "/electrons/" + version + "/ttbar_sub/").c_str(), kTRUE);
          if (c1) c1->SaveAs((path + "/electrons/" + version + "/ttbar_sub/" + title + ".pdf").c_str());
	  TFile f((path + "/electrons/" + version + "/ttbar_sub/" + title + ".root").c_str(),"RECREATE");
          h_data_fit->Write(title.c_str());
          f.Close();
        }
        if (ilepton==2) {
          gSystem->mkdir((path + "/muons/" + version + "/ttbar_sub/").c_str(), kTRUE);
          if (c1) c1->SaveAs((path + "/muons/" + version + "/ttbar_sub/" + title + ".pdf").c_str());
	  TFile f((path + "/muons/" + version + "/ttbar_sub/" + title + ".root").c_str(),"RECREATE");
          h_data_fit->Write(title.c_str());
          f.Close();
        }
      }
}
