#include <iostream>
#include <sstream>
#include <string.h>

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TLine.h"
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TTree.h"

#include "LumiLabel.C"
#include "LumiInfo_v06.h"

string path = "/gpfs/cms/users/candelis/work/Zb/data/" + version + "/";

void DataMCComp2(string& title="", int plot=0, int ilepton=1) {

if (ilepton<1 || ilepton>2) {
  ilepton = 1 + ilepton % 2;
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

	if (ilepton==1) data->cd("demo_ee");
	if (ilepton==2) data->cd("demo_mm");
	TH1F* h_data = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_data_b = (TH1F*)gDirectory->Get((title + "_b").c_str());

	if (ilepton==1) mc1->cd("demo_ee");
	if (ilepton==2) mc1->cd("demo_mm");
	TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc1_b = (TH1F*)gDirectory->Get((title + "_b").c_str());

	if (ilepton==1) mc2->cd("demo_ee");
	if (ilepton==2) mc2->cd("demo_mm");
	TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc2_b = (TH1F*)gDirectory->Get((title + "_b").c_str());

	if (ilepton==1) mc3->cd("demo_ee");
	if (ilepton==2) mc3->cd("demo_mm");
	TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc3_b = (TH1F*)gDirectory->Get((title + "_b").c_str());

	if (ilepton==1) mc4->cd("demo_ee");
	if (ilepton==2) mc4->cd("demo_mm");
	TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc4_b = (TH1F*)gDirectory->Get((title + "_b").c_str());

//	if (ilepton==1) mc5->cd("demo_ee");
//	if (ilepton==2) mc5->cd("demo_mm");
//	TH1F* h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
//	TH1F* h_mc5_b = (TH1F*)gDirectory->Get((title + "_b").c_str());

	if (ilepton==1) mc6->cd("demo_ee");
	if (ilepton==2) mc6->cd("demo_mm");
	TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc6_b = (TH1F*)gDirectory->Get((title + "_b").c_str());

	if (ilepton==1) mc7->cd("demo_ee");
	if (ilepton==2) mc7->cd("demo_mm");
	TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc7_b = (TH1F*)gDirectory->Get((title + "_b").c_str());

	h_data -> Sumw2();

	h_mc1 -> Sumw2();
	h_mc1 -> SetLineColor(kBlack);
	h_mc1 -> SetFillColor(kYellow-4);
	//h_mc1 -> SetFillStyle(3004);

	h_mc2 -> Sumw2();
	h_mc2 -> SetLineColor(kBlack);
	h_mc2 -> SetFillColor(kBlue);
	//h_mc2 -> SetFillStyle(3004);

	h_mc3 -> Sumw2();
	h_mc3 -> SetLineColor(kBlack);
	h_mc3 -> SetFillColor(kGray+2);
	//h_mc3 -> SetFillStyle(3004);

	h_mc4 -> Sumw2();
	h_mc4 -> SetLineColor(kBlack);
	h_mc4 -> SetFillColor(kGray+3);
	//h_mc4 -> SetFillStyle(3004);

//	h_mc5 -> Sumw2();
//	h_mc5 -> SetLineColor(kBlack);
//	h_mc5 -> SetFillColor(kGray+3);
//	//h_mc5 -> SetFillStyle(3004);

	h_mc6 -> Sumw2();
	h_mc6 -> SetLineColor(kBlack);
	h_mc6 -> SetFillColor(kRed+2);
	//h_mc6 -> SetFillStyle(3004);

	h_mc7 -> Sumw2();
	h_mc7 -> SetLineColor(kBlack);
	h_mc7 -> SetFillColor(kGray);
	//h_mc7 -> SetFillStyle(3004);

	h_mc1->Scale(norm1);
	h_mc2->Scale(norm2);
	h_mc3->Scale(norm3);
	h_mc4->Scale(norm4);
//	h_mc5->Scale(norm5);
	h_mc6->Scale(norm6);
	h_mc7->Scale(norm7);

	TH1F *ht = h_data->Clone("ht");
	ht->Add(h_mc7, -1);
	ht->Add(h_mc6, -1);
//	ht->Add(h_mc5, -1);
	ht->Add(h_mc4, -1);
	ht->Add(h_mc3, -1);
	ht->Add(h_mc2, -1);

	TH1F *ht_b = h_data_b->Clone("ht_b");
	ht_b->Add(h_mc7_b, -1);
	ht_b->Add(h_mc6_b, -1);
//	ht_b->Add(h_mc5_b, -1);
	ht_b->Add(h_mc4_b, -1);
	ht_b->Add(h_mc3_b, -1);
	ht_b->Add(h_mc2_b, -1);

	ht_b->Divide(ht);

	h_mc1_b->Divide(h_mc1);

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
	c1->cd();

	ht_b->Draw("EPX0");
	ht_b->GetYaxis()->SetTitle("Events");
	ht_b->SetMarkerColor(kBlack);
	ht_b->SetMarkerStyle(20);
	ht_b->SetMarkerSize (1.0);
	ht_b->GetXaxis()->SetRangeUser(0, 100);
	ht_b->GetYaxis()->SetRangeUser(0, 0.02);

	h_mc1_b->SetLineColor(kRed);
	h_mc1_b->Draw("LSAME");

	leg = new TLegend(0.68, 0.58, 0.88, 0.88);
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (ilepton==1) leg->AddEntry(h_data,"Z(#rightarrow ee)+jets","p");
	if (ilepton==2) leg->AddEntry(h_data,"Z(#rightarrow #mu#mu)+jets","p");

	leg->Draw();

	c1->Update();

	c1->cd();

 	TLatex *latexLabel=CMSPrel(19.3,"",0.1,0.94); // make fancy label
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

