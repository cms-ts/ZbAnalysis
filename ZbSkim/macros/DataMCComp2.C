#include "LumiLabel.C"
#include "LumiInfo_v09.h"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/" + version + "/";

void DataMCComp2(string& title="", int plot=0, int ilepton=1) {

if (ilepton<1 || ilepton>2) {
  ilepton = 1 + ilepton % 2;
}

	double c_b;
	double c_c; 
	double c_l;

	if (ilepton==1) {
	  c_b = 0.798;
	  c_c = 2.118; 
	  c_l = 1.354;
	}
	if (ilepton==2) {
	  c_b = 0.793;
	  c_c = 1.925; 
	  c_l = 1.611;
	}

	double Lumi2012;

	if (ilepton==1) Lumi2012 = Lumi2012_ele;
	if (ilepton==2) Lumi2012 = Lumi2012_muon;

	double norm1 = ( (Lumi2012 * Xsec_dy ) / Ngen_dy);
	double norm2 = ( (Lumi2012 * Xsec_tt) / Ngen_tt);
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

	if (ilepton==1)
	  TFile *data = TFile::Open((path + "DoubleElectron_2012_merge.root").c_str());
	if (ilepton==2)
	  TFile *data = TFile::Open((path + "DoubleMu_2012_merge.root").c_str());

	TFile *mc1 = TFile::Open((path + "DYJetsToLL.root").c_str());
	TFile *mc2 = TFile::Open((path + "TTbar.root").c_str());
	TFile *mc3 = TFile::Open((path + "ZZ.root").c_str());
	TFile *mc4 = TFile::Open((path + "WZ.root").c_str());
//	TFile *mc5 = TFile::Open((path + "QCD.root").c_str());
	TFile *mc6 = TFile::Open((path + "WW.root").c_str());
	TFile *mc7 = TFile::Open((path + "Wj.root").c_str());

	string title_b = title;

	if (title.find("_bjet_")!=string::npos) {
	  title.erase(title.find("_bjet_")+1, 1);
	} else {
	  title_b = title + "_b";
        }

	if (ilepton==1) data->cd("demo_ee");
	if (ilepton==2) data->cd("demo_mm");
	TH1F* h_data = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_data_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc1->cd("demo_ee");
	if (ilepton==2) mc1->cd("demo_mm");
	TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc1_b = (TH1F*)gDirectory->Get(title_b.c_str());
	TH1F* h_mc1b_b = (TH1F*)gDirectory->Get(("b"+title_b.substr(1)).c_str());
	TH1F* h_mc1c_b = (TH1F*)gDirectory->Get(("c"+title_b.substr(1)).c_str());

	if (ilepton==1) mc2->cd("demo_ee");
	if (ilepton==2) mc2->cd("demo_mm");
	TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc2_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc3->cd("demo_ee");
	if (ilepton==2) mc3->cd("demo_mm");
	TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc3_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc4->cd("demo_ee");
	if (ilepton==2) mc4->cd("demo_mm");
	TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc4_b = (TH1F*)gDirectory->Get(title_b.c_str());

//	if (ilepton==1) mc5->cd("demo_ee");
//	if (ilepton==2) mc5->cd("demo_mm");
//	TH1F* h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
//	TH1F* h_mc5_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc6->cd("demo_ee");
	if (ilepton==2) mc6->cd("demo_mm");
	TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc6_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc7->cd("demo_ee");
	if (ilepton==2) mc7->cd("demo_mm");
	TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc7_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (title_b.find("_jet_")!=string::npos) {
	  title.insert(title.find("_jet_")+1, "b");
	}

	h_data -> Sumw2();

	h_mc1 -> Sumw2();
	h_mc2 -> Sumw2();
	h_mc3 -> Sumw2();
	h_mc4 -> Sumw2();
//	h_mc5 -> Sumw2();
	h_mc6 -> Sumw2();
	h_mc7 -> Sumw2();

	h_mc1_b -> Sumw2();
	if (h_mc1b_b) h_mc1b_b -> Sumw2();
	if (h_mc1c_b) h_mc1c_b -> Sumw2();
	h_mc2_b -> Sumw2();
	h_mc3_b -> Sumw2();
	h_mc4_b -> Sumw2();
//	h_mc5_b -> Sumw2();
	h_mc6_b -> Sumw2();
	h_mc7_b -> Sumw2();

	h_mc1->Scale(norm1);
	h_mc2->Scale(norm2);
	h_mc3->Scale(norm3);
	h_mc4->Scale(norm4);
//	h_mc5->Scale(norm5);
	h_mc6->Scale(norm6);
	h_mc7->Scale(norm7);

	h_mc1_b->Scale(norm1);
	if (h_mc1b_b) h_mc1b_b->Scale(norm1);
	if (h_mc1c_b) h_mc1c_b->Scale(norm1);
	h_mc2_b->Scale(norm2);
	h_mc3_b->Scale(norm3);
	h_mc4_b->Scale(norm4);
//	h_mc5_b->Scale(norm5);
	h_mc6_b->Scale(norm6);
	h_mc7_b->Scale(norm7);

	h_data->Add(h_mc7, -1.);
	h_data->Add(h_mc6, -1.);
//	h_data->Add(h_mc5, -1.);
	h_data->Add(h_mc4, -1.);
	h_data->Add(h_mc3, -1.);
	h_data->Add(h_mc2, -1.);

	h_data_b->Add(h_mc7_b, -1.);
	h_data_b->Add(h_mc6_b, -1.);
//	h_data_b->Add(h_mc5_b, -1.);
	h_data_b->Add(h_mc4_b, -1.);
	h_data_b->Add(h_mc3_b, -1.);
	h_data_b->Add(h_mc2_b, -1.);

	TH1F *h_mc1uds_b = h_mc1_b->Clone("h_mc1uds_b");
	if (h_mc1b_b) h_mc1uds_b->Add(h_mc1b_b, -1);
	if (h_mc1c_b) h_mc1uds_b->Add(h_mc1c_b, -1);
	for (int i=0; i<=h_mc1uds_b->GetNbinsX()+1; i++) {
	  float e = h_mc1uds_b->GetBinError(i)**2;
	  if (h_mc1b_b) e = e - h_mc1b_b->GetBinError(i)**2;
	  if (h_mc1c_b) e = e - h_mc1c_b->GetBinError(i)**2;
	  h_mc1uds_b->SetBinError(i, TMath::Sqrt(e));
	}

	h_mc1uds_b->Scale(c_l);
	h_mc1b_b->Scale(c_b);
	h_mc1c_b->Scale(c_c);

	h_data_b->Add(h_mc1c_b, -1.);
	h_data_b->Add(h_mc1uds_b, -1.);

	h_data_b->Divide(h_data);
	h_data_b->Scale(100.);

	h_mc1b_b->Divide(h_mc1);
	h_mc1b_b->Scale(100.);

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
	c1->cd();

	h_data_b->SetTitle("");
	if (title=="w_first_jet_pt") {
	  h_data_b->GetXaxis ()->SetTitle("leading jet p_{T} [GeV/c]");
	  h_data_b->GetXaxis()->SetRangeUser(0, 200);
	  h_data_b->GetYaxis()->SetRangeUser(0, 4);
	} else if (title=="w_first_jet_eta") {
	  h_data_b->GetXaxis ()->SetTitle("leading jet #eta");
	  h_data_b->GetYaxis()->SetRangeUser(0, 4);
	} else if (title=="w_first_bjet_pt") {
	  h_data_b->GetXaxis ()->SetTitle("leading b-jet p_{T} [GeV/c]");
	  h_data_b->GetXaxis()->SetRangeUser(0, 200);
	  h_data_b->GetYaxis()->SetRangeUser(0, 4);
	} else if (title=="w_first_bjet_eta") {
	  h_data_b->GetXaxis ()->SetTitle("leading b-jet #eta");
	  h_data_b->GetYaxis()->SetRangeUser(0, 4);
	} else if (title=="w_pt_Z_ee"||title=="w_pt_Z_mm") {
	  h_data_b->GetXaxis ()->SetTitle("Z boson p_{T} [GeV/c]");
	  h_data_b->GetXaxis()->SetRangeUser(0, 200);
	  h_data_b->GetYaxis()->SetRangeUser(0, 5);
	} else if (title=="w_Ht") {
	  h_data_b->GetXaxis ()->SetTitle("H_{T} [GeV/c]");
	  h_data_b->GetXaxis()->SetRangeUser(0, 250);
	  h_data_b->GetYaxis()->SetRangeUser(0, 8);
	}

	h_data_b->Draw("EP");
	h_data_b->GetYaxis()->SetTitle("#sigma_{Z+b-jets}/#sigma_{Z+jets} [%]");
	h_data_b->GetYaxis()->SetTitleOffset(1.2);
	h_data_b->GetXaxis()->SetTitleOffset(1.3);
	h_data_b->SetMarkerColor(kBlack);
	h_data_b->SetLineColor(kBlack);
	h_data_b->SetMarkerStyle(20);
	h_data_b->SetMarkerSize (1.0);
	h_data_b->SetStats(0);

	h_mc1b_b->SetLineColor(kRed);
	h_mc1b_b->SetMarkerColor(kRed);
	h_mc1b_b->SetMarkerStyle(20);
	h_mc1b_b->SetMarkerSize (1.0);
	h_mc1b_b->Draw("HISTSAME");
	h_mc1b_b->SetStats(0);

	leg = new TLegend(0.58, 0.08, 0.78, 0.38);
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (ilepton==1) leg->AddEntry(h_data_b,"Z(#rightarrow ee)+jets DATA","p");
	if (ilepton==2) leg->AddEntry(h_data_b,"Z(#rightarrow #mu#mu)+jets DATA","p");
	leg->AddEntry(h_mc1b_b,"Z(#rightarrow ee)+jets MC","l");

	leg->Draw();

	c1->Update();

	c1->cd();

 	TLatex *latexLabel = CMSPrel(Lumi2012/1000.,"",0.15,0.94);
	latexLabel->Draw("same");

	if (plot) {
	  if (ilepton==1) {
	    gSystem->mkdir(("electrons/" + version).c_str());
	    c1->SaveAs(("electrons/" + version + "/" + title + "_ratio" + ".pdf").c_str());
	  }
	  if (ilepton==2) {
	    gSystem->mkdir(("muons/" + version).c_str());
	    c1->SaveAs(("muons/" + version + "/" + title + "_ratio" + ".pdf").c_str());
	  }
	}
}

