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

TH1F* h_data = 0;
TH1F* h_mc1 = 0;
TH1F* h_mc1b = 0;
TH1F* h_mc1c = 0;

double func(double* x, double* p) {
  int i = h_mc1->GetXaxis()->FindBin(x[0]);
  float c_mc1 = h_mc1->GetBinContent(i);
  float c_mc1b = h_mc1b->GetBinContent(i);
  float c_mc1c = h_mc1c->GetBinContent(i);
  return (1.0-(1.0-p[0])-(1.0-p[1]))*c_mc1 + p[0]*c_mc1b + p[1]*c_mc1c;
}

void DataMCComp(string& title="", int plot=0, int ilepton=1, int doBkg=0, int doFit=0) {

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
	  if (title=="first_muon_pt") return;
	  if (title=="first_muon_eta") return;
	  if (title=="w_mm_inv") return;
	  if (title=="b_mm_inv") return;
	  if (title=="w_delta_phi_mm") return;
	  if (title=="Z_pt_mm") return;
	  if (title=="Z_pt_mm_b") return;
	}
	if (ilepton==2) {
	  if (title=="first_ele_pt") return;
	  if (title=="first_ele_eta") return;
	  if (title=="w_ee_inv") return;
	  if (title=="b_ee_inv") return;
	  if (title=="w_delta_phi_ee") return;
	  if (title=="Z_pt_ee") return;
	  if (title=="Z_pt_ee_b") return;
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
	h_data = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc1->cd("demo_ee");
	if (ilepton==2) mc1->cd("demo_mm");
	h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	if (title == "w_secondvtx_N") {
	  h_mc1b = (TH1F*)gDirectory->Get("bquarks");
	}
	if (title == "SVTX_mass_jet") {
	  h_mc1b = (TH1F*)gDirectory->Get("SVTX_mass_jet_b");
	  h_mc1c = (TH1F*)gDirectory->Get("SVTX_mass_jet_c");
	}
	if (title == "SVTX_mass_trk") {
	  h_mc1b = (TH1F*)gDirectory->Get("SVTX_mass_trk_b");
	  h_mc1c = (TH1F*)gDirectory->Get("SVTX_mass_trk_c");
	}
	if (title == "SVTX_mass") {
	  h_mc1b = (TH1F*)gDirectory->Get("SVTX_mass_b");
	  h_mc1c = (TH1F*)gDirectory->Get("SVTX_mass_c");
	}
	if (title == "w_first_jet_pt") {
	  h_mc1b = (TH1F*)gDirectory->Get("b_first_jet_pt");
	}
	if (title == "w_first_jet_eta") {
	  h_mc1b = (TH1F*)gDirectory->Get("b_first_jet_eta");
	}
	if (title == "Z_pt_ee") {
	  h_mc1b = (TH1F*)gDirectory->Get("b_pt_Z_ee");
	}
	if (title == "Z_pt_mm") {
	  h_mc1b = (TH1F*)gDirectory->Get("b_pt_Z_mm");
	}
	if (title == "w_mm_inv") {
	  h_mc1b = (TH1F*)gDirectory->Get("b_invMass_mm");
	}
	if (title == "w_ee_inv") {
	  h_mc1b = (TH1F*)gDirectory->Get("b_invMass_ee");
	}

	if (doFit && (h_mc1b==0 || h_mc1c==0)) return;

	if (ilepton==1) mc2->cd("demo_ee");
	if (ilepton==2) mc2->cd("demo_mm");
	TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc3->cd("demo_ee");
	if (ilepton==2) mc3->cd("demo_mm");
	TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc4->cd("demo_ee");
	if (ilepton==2) mc4->cd("demo_mm");
	TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());

//	if (ilepton==1) mc5->cd("demo_ee");
//	if (ilepton==2) mc5->cd("demo_mm");
//	TH1F* h_mc5 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc6->cd("demo_ee");
	if (ilepton==2) mc6->cd("demo_mm");
	TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc7->cd("demo_ee");
	if (ilepton==2) mc7->cd("demo_mm");
	TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());

	h_data -> Sumw2();

	h_mc1 -> Sumw2();
	h_mc1 -> SetLineColor(kBlack);
	h_mc1 -> SetFillColor(kYellow-4);
	//h_mc1 -> SetFillStyle(3004);

	if (h_mc1b) {
	  h_mc1b -> Sumw2();
	  h_mc1b->SetLineColor(kBlack);
	  h_mc1b->SetFillColor(kYellow-4);
	  h_mc1b->SetFillStyle(3254);
	}
	if (h_mc1c) {
	  h_mc1c -> Sumw2();
	  h_mc1c->SetLineColor(kBlack);
	  h_mc1c->SetFillColor(kYellow-4);
	  h_mc1c->SetFillStyle(3245);
	}

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
	if (h_mc1b) h_mc1b->Scale(norm1);
	if (h_mc1c) h_mc1c->Scale(norm1);
	h_mc2->Scale(norm2);
	h_mc3->Scale(norm3);
	h_mc4->Scale(norm4);
//	h_mc5->Scale(norm5);
	h_mc6->Scale(norm6);
	h_mc7->Scale(norm7);

	if (doBkg) {
	  h_data->Add(h_mc7, -1.);
	  h_data->Add(h_mc6, -1.);
//	  h_data->Add(h_mc5, -1.);
	  h_data->Add(h_mc4, -1.);
	  h_data->Add(h_mc3, -1.);
	  h_data->Add(h_mc2, -1.);
	}

	for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
	  if (h_mc1b) h_mc1b->SetBinError(i, 0.);
	  if (h_mc1c) h_mc1c->SetBinError(i, 0.);
	}
	if (h_mc1b) h_mc1->Add(h_mc1b, -1.);
	if (h_mc1c) h_mc1->Add(h_mc1c, -1.);

	TH1F *h_data_fit = h_data->Clone("h_data_fit");
	TF1 *f1 = new TF1("f1", func, 0.00, 1.00, 2);
	if (doFit) {
	  if (!doBkg) {
	    h_data_fit->Add(h_mc7, -1.);
	    h_data_fit->Add(h_mc6, -1.);
//	    h_data_fit->Add(h_mc5, -1.);
	    h_data_fit->Add(h_mc4, -1.);
	    h_data_fit->Add(h_mc3, -1.);
	    h_data_fit->Add(h_mc2, -1.);
	  }
	  for (int i=0; i<=h_data_fit->GetNbinsX()+1; i++) {
	    float e = h_data_fit->GetBinError(i)**2 + h_mc1->GetBinError(i)**2;
	    h_data_fit->SetBinError(i, TMath::Sqrt(e));
	  }
	  f1->SetParameters(1.0, 1.0);
	  f1->SetParNames("f_b", "f_c");
	  h_data_fit->Fit("f1");
	  h_mc1->Scale(1.0-(1.0-f1->GetParameter(0))-(1.0-f1->GetParameter(1)));
	  h_mc1b->Scale(f1->GetParameter(0));
	  h_mc1c->Scale(f1->GetParameter(1));
	}

	TH1F *ht = h_mc1->Clone("ht");
	ht->Reset();
	if (!doBkg) {
	  ht->Add(h_mc7);
	  ht->Add(h_mc6);
//	  ht->Add(h_mc5);
	  ht->Add(h_mc4);
	  ht->Add(h_mc3);
	  ht->Add(h_mc2);
	}
	if (h_mc1b) ht->Add(h_mc1b);
	if (h_mc1c) ht->Add(h_mc1c);
	ht->Add(h_mc1);

	THStack *hs = new THStack("hs","");
	if (!doBkg) {
//	  hs->Add(h_mc5);
	  hs->Add(h_mc6);
	  hs->Add(h_mc7);
	  hs->Add(h_mc4);
	  hs->Add(h_mc3);
	  hs->Add(h_mc2);
	}
	if (h_mc1b) hs->Add(h_mc1b);
	if (h_mc1c) hs->Add(h_mc1c);
	hs->Add(h_mc1);

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
	c1->cd();

	TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
	pad1->SetBottomMargin(0.001);
	pad1->Draw();
	pad1->cd();
	pad1->SetLogy();

	TH1F *h3=h_data->DrawCopy();
	h3->SetMinimum(-100);

	hs->Draw("HIST");
	hs->GetYaxis()->SetTitle("Events");
 	hs->GetXaxis()->SetLabelSize(0.08);
	hs->GetXaxis()->SetTitleOffset (0.7);
	hs->SetMinimum(8);

	h_data->Draw("EPX0SAMES");
	h_data->SetMarkerColor(kBlack);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize (1.0);
	//h_data->SetMaximum(1.2);

	leg = new TLegend(0.68, 0.58, 0.88, 0.88);
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (ilepton==1) leg->AddEntry(h_data,"Z(#rightarrow ee)+jets","p");
	if (ilepton==2) leg->AddEntry(h_data,"Z(#rightarrow #mu#mu)+jets","p");

	leg->AddEntry(h_mc1,"Z+jets","f");
	if (h_mc1b) leg->AddEntry(h_mc1b,"Z+b-jets","f");
	if (h_mc1c) leg->AddEntry(h_mc1c,"Z+c-jets","f");
	if (!doBkg) {
	  leg->AddEntry(h_mc2,"t#bar{t}","f");
	  leg->AddEntry(h_mc3,"ZZ","f");
	  leg->AddEntry(h_mc4,"WZ","f");
	  leg->AddEntry(h_mc7,"W+jets", "f");
	  leg->AddEntry(h_mc6,"WW","f");
//	  leg->AddEntry(h_mc5,"QCD","f");
	}
	leg->Draw();

	pad1->Update();
	c1->Update();

	c1->cd();

	TH1F *h_ratio = h_data->Clone("h_ratio");

	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();
	h_ratio->SetStats(0);
	h_ratio->SetTitle("");
	h_ratio->SetStats(0);

	if (title=="w_jetmultiplicity") {
	  h_ratio->GetXaxis ()->SetTitle("Number of jets");
	} else if (title=="w_mm_inv") {
	  h_ratio->GetXaxis ()->SetTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
	} else if (title=="w_ee_inv") {
	  h_ratio->GetXaxis ()->SetTitle("e^{+}e^{-} invariant mass [GeV/c^{2}]");
	} else if (title=="first_muon_pt") {
	  h_ratio->GetXaxis ()->SetTitle("muon p_{T} [GeV/c]");
	} else if (title=="first_ele_pt") {
	  h_ratio->GetXaxis ()->SetTitle("electron p_{T} [GeV/c]");
	} else if (title=="first_jet_pt") {
	  h_ratio->GetXaxis ()->SetTitle("jet p_{T} [GeV/c]");
	} else if (title=="w_secondvtx_N") {
	  h_ratio->GetXaxis ()->SetTitle("CSV discriminator");
	} else if (title=="w_recoVTX") {
	  h_ratio->GetXaxis ()->SetTitle("Number of offline vertices");
	} else if (title=="w_muon_pt") {
	  h_ratio->GetXaxis ()->SetTitle("muon p_{T} [GeV/c]");
	} else if (title=="w_MET_sign") {
	  h_ratio->GetXaxis ()->SetTitle("MET [GeV/c]");
	} else if (title=="w_MET") {
	  h_ratio->GetXaxis ()->SetTitle("MET Significance [GeV/c]");
	} else if (title=="Z_pt_ee") {
	  h_ratio->GetXaxis ()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="Z_pt_mm") {
	  h_ratio->GetXaxis ()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_bjetmultiplicity") {
	  h_ratio->GetXaxis ()->SetTitle("b quark jets multiplicity");
	} else if (title=="w_first_jet_pt_b") {
	  h_ratio->GetXaxis ()->SetTitle("leading b quark p_{T} [GeV/c]");
	} else if (title=="w_first_jet_eta_b") {
	  h_ratio->GetXaxis ()->SetTitle("leading b quark #eta");
	} else if (title=="Z_pt_mm_b") {
	  h_ratio->GetXaxis ()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="Z_pt_ee_b") {
	  h_ratio->GetXaxis ()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_delta_phi_ee") {
	  h_ratio->GetXaxis ()->SetTitle("#Delta#phi(bZ) [rad]");
	} else if (title=="w_delta_phi_mm") {
	  h_ratio->GetXaxis ()->SetTitle("#Delta#phi(bZ) [rad]");
	} else if (title=="b_ee_inv"||title=="b_mm_inv") {
	  h_ratio->GetXaxis ()->SetTitle("Z mass + 1 b quark [GeV/c^{2}]");
	} else if (title=="SVTX_mass_jet"||title=="SVTX_mass_trk"||title=="SVTX_mass") {
	  h_ratio->GetXaxis ()->SetTitle("SV mass [GeV/c^{2}]");
	}

	h_ratio->GetXaxis()->SetTitleOffset(0.9);
 	h_ratio->GetXaxis()->SetTitleSize(0.1);
	h_ratio->GetXaxis()->SetLabelFont(42);
 	h_ratio->GetXaxis()->SetLabelSize(0.08);
	h_ratio->GetXaxis()->SetTitleFont(42);
	h_ratio->GetYaxis()->SetTitle("Data/MC");
	h_ratio->GetYaxis()->SetNdivisions(505);
	h_ratio->GetYaxis()->SetTitleSize(0.09);
	h_ratio->GetYaxis()->SetLabelSize(0.08);
	h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
	h_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_ratio->Divide(ht);
	h_ratio->SetMarkerStyle(20);
	h_ratio->Draw("EPX0");

	TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
	OLine->SetLineColor(kRed);
	OLine->SetLineStyle(2);
	OLine->Draw();

	c1->cd();

 	TLatex *latexLabel=CMSPrel(19.3,"",0.1,0.94); // make fancy label
	latexLabel->Draw("same");

	if (doFit) {
	  TLatex *fitLabel = new TLatex();
	  fitLabel->SetTextSize(0.0275);
	  fitLabel->SetTextFont(42);
	  fitLabel->SetLineWidth(2);
	  fitLabel->SetNDC();
	  char buff[100];
	  sprintf(buff, "f_{b} = %5.3f #pm %5.3f", f1->GetParameter(0), f1->GetParError(0));
	  fitLabel->DrawLatex(0.68, 0.53, buff);
	  sprintf(buff, "f_{c} = %5.3f #pm %5.3f", f1->GetParameter(1), f1->GetParError(1));
	  fitLabel->DrawLatex(0.68, 0.48, buff);
	  fitLabel->Draw("same");
	}

	if (plot) {
	  if (doBkg) title = title + "_doBkg";
	  if (doFit) title = title + "_doFit";
	  if (ilepton==1) {
	    gSystem->mkdir(("electrons/" + version).c_str());
	    c1->SaveAs(("electrons/" + version + "/" + title + ".pdf").c_str());
	  }
	  if (ilepton==2) {
	    gSystem->mkdir(("muons/" + version).c_str());
	    c1->SaveAs(("muons/" + version + "/" + title + ".pdf").c_str());
	  }
	}
}

