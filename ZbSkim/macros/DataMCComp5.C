#include "LumiLabel.C"
#include "LumiInfo_v09.h"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";

TH1F* h_data = 0;
TH1F* h_data_fit = 0;

double func(double* x, double* p) {
  double val = 0.0;
  int i = h_data->GetXaxis()->FindBin(x[0]);
  if (h_data_fit) val = val + p[0]*h_data_fit->GetBinContent(i);
  return val;
}

void DataMCComp5(string& title="", int plot=0, int ilepton=1) {

      if (ilepton<1 || ilepton>2) {
        ilepton = 1 + ilepton % 2;
      }

      if (ilepton==1) {
        if (title.find("muon")!=string::npos) return;
	if (title.find("mm")!=string::npos) return;
      }
      if (ilepton==2) {
	if (title.find("ele")!=string::npos) return;
	if (title.find("ee")!=string::npos) return;
      }

      double Lumi2012;

      if (ilepton==1) Lumi2012 = Lumi2012_ele;
      if (ilepton==2) Lumi2012 = Lumi2012_muon;

      if (title.empty()) title = "w_jetmultiplicity";

      TFile *data_ee = TFile::Open((path + "/" + version + "/" + "DoubleElectron_2012_merge.root").c_str());
      TFile *data_mm = TFile::Open((path + "/" + version + "/" + "DoubleMu_2012_merge.root").c_str());

      if (ilepton==1) data_ee->cd("demoEle");
      if (ilepton==2) data_mm->cd("demoMuo");
      h_data = (TH1F*)gDirectory->Get(title.c_str());

      string title_em = title;

      if (title_em.find("_ee_")!=string::npos) {
        title_em.replace(title_em.find("_ee_")+2, 1, "m");
      }
      if (title_em.find("_mm_")!=string::npos) {
        title_em.replace(title_em.find("_mm_")+1, 1, "e");
      }

      TFile *data_em = TFile::Open((path + "/" + version + "/" + "MuEG_2012_merge.root").c_str());

      data_em->cd("demoEleMuo");
      h_data_fit = (TH1F*)gDirectory->Get(title_em.c_str());

      h_data->Sumw2();
      h_data_fit->Sumw2();

      h_data->Scale(1./Lumi2012);
      h_data_fit->Scale(1./Lumi2012);

      TF1 *f1 = new TF1("f1", func, 60, 120.00, 1);
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

      TCanvas* c1 = new TCanvas("c", "c", 800, 600);
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

      c1->cd();

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

      if (plot) {
        if (ilepton==1) {
          gSystem->mkdir((path + "electrons/" + version + "/ttbar_sub/").c_str(), kTRUE);
          c1->SaveAs((path + "electrons/" + version + "/ttbar_sub/" + title + ".pdf").c_str());
        }
        if (ilepton==2) {
          gSystem->mkdir((path + "muons/" + version + "/ttbar_sub/").c_str(), kTRUE);
          c1->SaveAs((path + "muons/" + version + "/ttbar_sub/" + title + ".pdf").c_str());
        }
      }
}
