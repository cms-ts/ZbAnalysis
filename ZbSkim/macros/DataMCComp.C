#include "LumiLabel.C"
#include "LumiInfo_v10.h"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";

TH1F* h_data_fit = 0;
TH1F* h_mc_fit0 = 0;
TH1F* h_mc_fit1 = 0;
TH1F* h_mc_fit2 = 0;

double func(double* x, double* p) {
  double val = 0.0;
  int i = h_data_fit->GetXaxis()->FindBin(x[0]);
  if (h_mc_fit0) val = val + p[0]*h_mc_fit0->GetBinContent(i);
  if (h_mc_fit1) val = val + p[1]*h_mc_fit1->GetBinContent(i);
  if (h_mc_fit2) val = val + p[2]*h_mc_fit2->GetBinContent(i);
  return val;
}

void DataMCComp(string& title="", int plot=0, int ilepton=1, int doBkg=0, int doFit=0) {

//int useEleMuo = 0;
int useEleMuo = 1;

	double c_t=1.0;

	if (ilepton==1 && doFit==3) c_t=1.062;
	if (ilepton==2 && doFit==3) c_t=1.013;
	if (ilepton==3 && doFit==3) c_t=0.989;

	double a1_t=1.0;
	double a2_t=1.0;

	if (ilepton==1 && useEleMuo) {
	  a1_t=0.451;
	  a2_t=0.430;
	}
	if (ilepton==2 && useEleMuo) {
	  a1_t=0.573;
	  a2_t=0.553;
	}

	double Lumi2012;

	if (ilepton==1) Lumi2012 = Lumi2012_ele;
	if (ilepton==2) Lumi2012 = Lumi2012_muon;
	if (ilepton==3) Lumi2012 = Lumi2012_ele_muon;

	double norm1 = ( (Lumi2012 * Xsec_dy) / Ngen_dy);
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
	  TFile *data = TFile::Open((path + "/" + version + "/" + "DoubleElectron_2012_merge.root").c_str());
	if (ilepton==2)
	  TFile *data = TFile::Open((path + "/" + version + "/" + "DoubleMu_2012_merge.root").c_str());
	if (ilepton==3)
	  TFile *data = TFile::Open((path + "/" + version + "/" + "MuEG_2012_merge.root").c_str());

	TFile *mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL.root").c_str());
	TFile *mc2 = TFile::Open((path + "/" + version + "/" + "TTbar.root").c_str());
	TFile *mc3 = TFile::Open((path + "/" + version + "/" + "ZZ.root").c_str());
	TFile *mc4 = TFile::Open((path + "/" + version + "/" + "WZ.root").c_str());
//	TFile *mc5 = TFile::Open((path + "/" + version + "/" + "QCD.root").c_str());
	TFile *mc6 = TFile::Open((path + "/" + version + "/" + "WW.root").c_str());
	TFile *mc7 = TFile::Open((path + "/" + version + "/" + "Wj.root").c_str());

	if (ilepton==1) data->cd("demoEle");
	if (ilepton==2) data->cd("demoMuo");
	if (ilepton==3) data->cd("demoEleMuo");
	TH1F* h_data = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc1->cd("demoEle");
	if (ilepton==2) mc1->cd("demoMuo");
	if (ilepton==3) mc1->cd("demoEleMuo");
	TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc1b = (TH1F*)gDirectory->Get(("b"+title.substr(1)).c_str());
	TH1F* h_mc1c = (TH1F*)gDirectory->Get(("c"+title.substr(1)).c_str());

	if (ilepton==1) mc2->cd("demoEle");
	if (ilepton==2) mc2->cd("demoMuo");
	if (ilepton==3) mc2->cd("demoEleMuo");
	TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());

	if (useEleMuo) {
	  if (ilepton==1) {
	    mc2 = TFile::Open((path + "/electrons/" + version + "/ttbar_sub/" + title + ".root").c_str());
	    h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	  }
	  if (ilepton==2) {
	    mc2 = TFile::Open((path + "/muons/" + version + "/ttbar_sub/" + title + ".root").c_str());
	    h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	  }
	}

	if (ilepton==1) mc3->cd("demoEle");
	if (ilepton==2) mc3->cd("demoMuo");
	if (ilepton==3) mc3->cd("demoEleMuo");
	TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc4->cd("demoEle");
	if (ilepton==2) mc4->cd("demoMuo");
	if (ilepton==3) mc4->cd("demoEleMuo");
	TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());

//	if (ilepton==1) mc5->cd("demoEle");
//	if (ilepton==2) mc5->cd("demoMuo");
//	if (ilepton==3) mc5->cd("demoEleMuo");
//	TH1F* h_mc5 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc6->cd("demoEle");
	if (ilepton==2) mc6->cd("demoMuo");
	if (ilepton==3) mc6->cd("demoEleMuo");
	TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc7->cd("demoEle");
	if (ilepton==2) mc7->cd("demoMuo");
	if (ilepton==3) mc7->cd("demoEleMuo");
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
	  h_mc1c->SetFillColor(kOrange);
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
	h_mc2->Scale(norm2*c_t);
	if (useEleMuo) {
	  if (title.find("_b")==string::npos) {
	    h_mc2->Scale(a1_t*(Lumi2012/Lumi2012_ele_muon)/(norm2*c_t));
	  } else {
	    h_mc2->Scale(a2_t*(Lumi2012/Lumi2012_ele_muon)/(norm2*c_t));
	  }
	}
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

	if (h_mc1b) h_mc1->Add(h_mc1b, -1.);
	if (h_mc1c) h_mc1->Add(h_mc1c, -1.);
	for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
	  float e = h_mc1->GetBinError(i)**2;
	  if (h_mc1b) e = e - h_mc1b->GetBinError(i)**2;
	  if (h_mc1c) e = e - h_mc1c->GetBinError(i)**2;
	  h_mc1->SetBinError(i, TMath::Sqrt(e));
	}

	TF1 *f1 = new TF1("f1", func, 0.00, 100.00, 3);
	if (doFit==1) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    h_data_fit->Add(h_mc7, -1.);
	    h_data_fit->Add(h_mc6, -1.);
//	    h_data_fit->Add(h_mc5, -1.);
	    h_data_fit->Add(h_mc4, -1.);
	    h_data_fit->Add(h_mc3, -1.);
	  }
	  h_data_fit->Add(h_mc1, -1.);
	  if (h_mc1b) h_data_fit->Add(h_mc1b, -1.);
	  if (h_mc1c) h_data_fit->Add(h_mc1c, -1.);
	  h_mc_fit0 = h_mc2;
	  for (int i=0; i<=h_data_fit->GetNbinsX()+1; i++) {
	    float e = h_data_fit->GetBinError(i)**2;
	    e = e + h_mc_fit0->GetBinError(i)**2;
	    if (title=="w_MET" && h_data_fit->GetXaxis()->GetBinCenter(i) < 125.) e = 1.e10;
	    if (title=="w_MET_b" && h_data_fit->GetXaxis()->GetBinCenter(i) < 90.) e = 1.e10;
	    if (title=="w_MET_sign" && h_data_fit->GetXaxis()->GetBinCenter(i) < 50.) e = 1.e10;
	    if (title=="w_MET_sign_b" && h_data_fit->GetXaxis()->GetBinCenter(i) < 30.) e = 1.e10;
	    h_data_fit->SetBinError(i, TMath::Sqrt(e));
	  }
	  f1->SetParameters(1.0, 0.0, 0.0);
	  f1->SetParNames("c(t)", "dummy", "dummy");
	  f1->FixParameter(1, 0.0);
	  f1->FixParameter(2, 0.0);
	  h_data_fit->Fit("f1", "Q0");
	  h_mc_fit0->Scale(f1->GetParameter(0));
	}
	if (doFit==2) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    h_data_fit->Add(h_mc7, -1.);
	    h_data_fit->Add(h_mc6, -1.);
//	    h_data_fit->Add(h_mc5, -1.);
	    h_data_fit->Add(h_mc4, -1.);
	    h_data_fit->Add(h_mc3, -1.);
	  }
	  h_mc_fit0 = h_mc1;
	  if (h_mc1b) h_mc_fit0->Add(h_mc1b, 1.);
	  if (h_mc1c) h_mc_fit0->Add(h_mc1c, 1.);
	  h_mc_fit1 = h_mc2;
	  for (int i=0; i<=h_data_fit->GetNbinsX()+1; i++) {
	    float e = h_data_fit->GetBinError(i)**2;
	    e = e + h_mc_fit0->GetBinError(i)**2;
	    e = e + h_mc_fit1->GetBinError(i)**2;
	    h_data_fit->SetBinError(i, TMath::Sqrt(e));
	  }
	  f1->SetParameters(1.0, 1.0, 0.0);
	  f1->SetParNames("c(Z+jets)", "c(t)", "dummy");
	  f1->FixParameter(2, 0.0);
	  h_data_fit->Fit("f1", "Q0");
	  if (h_mc1b) h_mc_fit0->Add(h_mc1b, -1.);
	  if (h_mc1c) h_mc_fit0->Add(h_mc1c, -1.);
	  h_mc_fit0->Scale(f1->GetParameter(0));
	  if (h_mc1b) h_mc1b->Scale(f1->GetParameter(0));
	  if (h_mc1c) h_mc1c->Scale(f1->GetParameter(0));
	  h_mc_fit1->Scale(f1->GetParameter(1));
	}
	if (doFit==3) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    h_data_fit->Add(h_mc7, -1.);
	    h_data_fit->Add(h_mc6, -1.);
//	    h_data_fit->Add(h_mc5, -1.);
	    h_data_fit->Add(h_mc4, -1.);
	    h_data_fit->Add(h_mc3, -1.);
	    h_data_fit->Add(h_mc2, -1.);
	  }
	  h_mc_fit0 = h_mc1;
	  h_mc_fit1 = h_mc1b;
	  h_mc_fit2 = h_mc1c;
	  for (int i=0; i<=h_data_fit->GetNbinsX()+1; i++) {
	    float e = h_data_fit->GetBinError(i)**2;
	    e = e + h_mc_fit0->GetBinError(i)**2;
	    e = e + h_mc_fit1->GetBinError(i)**2;
	    e = e + h_mc_fit2->GetBinError(i)**2;
	    if (title=="w_SVTX_mass" && h_data_fit->GetXaxis()->GetBinCenter(i) < 0.2) e = 1.e10;
	    h_data_fit->SetBinError(i, TMath::Sqrt(e));
	  }
	  f1->SetParameters(1.0, 1.0, 1.0);
	  f1->SetParNames("c(uds)", "c(b)", "c(c)");
	  h_data_fit->Fit("f1", "Q0");
	  h_mc_fit0->Scale(f1->GetParameter(0));
	  h_mc_fit1->Scale(f1->GetParameter(1));
	  h_mc_fit2->Scale(f1->GetParameter(2));
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

	TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
	pad1->SetBottomMargin(0.001);
	pad1->Draw();
	pad1->cd();
	pad1->SetLogy();
	if (title.find("MET")!=string::npos) pad1->SetLogy(0);

	hs->Draw("HIST");
	hs->GetYaxis()->SetTitle("Events");
 	hs->GetXaxis()->SetLabelSize(0.08);
	hs->GetXaxis()->SetTitleOffset(0.7);
	hs->SetMinimum(8);

	h_data->Draw("EPX0SAMES");
	h_data->SetMarkerColor(kBlack);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize (1.0);
	//h_data->SetStats(0);

	TLegend *leg;
	if (doBkg) {
	  if (h_mc1c && h_mc1b) {
	    leg = new TLegend(0.62, 0.747, 0.88, 0.88);
	  } else if (h_mc1c || h_mc1b) {
	    leg = new TLegend(0.62, 0.780, 0.88, 0.88);
	  } else {
	    leg = new TLegend(0.62, 0.813, 0.88, 0.88);
	  }
	} else {
	  if (h_mc1c && h_mc1b) {
	    leg = new TLegend(0.62, 0.580, 0.88, 0.88);
	  } else if (h_mc1c || h_mc1b) {
	    leg = new TLegend(0.62, 0.613, 0.88, 0.88);
	  } else {
	    leg = new TLegend(0.62, 0.647, 0.88, 0.88);
	  }
	}
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (ilepton==1) leg->AddEntry(h_data,"Z(#rightarrow ee)+jets","p");
	if (ilepton==2) leg->AddEntry(h_data,"Z(#rightarrow #mu#mu)+jets","p");
	if (ilepton==2) leg->AddEntry(h_data,"Z(#rightarrow e#mu)+jets","p");

	leg->AddEntry(h_mc1,"Z+jets","f");
	if (h_mc1c) leg->AddEntry(h_mc1c,"Z+c-jets","f");
	if (h_mc1b) leg->AddEntry(h_mc1b,"Z+b-jets","f");
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
	h_ratio->SetTitle("");
	h_ratio->SetStats(0);

	if (title=="w_jetmultiplicity") {
	  h_ratio->GetXaxis ()->SetTitle("jet multiplicity");
	} else if (title=="w_mass_ee"||title=="w_mass_mm") {
	  h_ratio->GetXaxis ()->SetTitle("invariant mass [GeV/c^{2}]");
	} else if (title=="w_first_muon_pt") {
	  h_ratio->GetXaxis ()->SetTitle("muon p_{T} [GeV/c]");
	} else if (title=="w_first_ele_pt") {
	  h_ratio->GetXaxis ()->SetTitle("electron p_{T} [GeV/c]");
	} else if (title=="w_first_jet_pt") {
	  h_ratio->GetXaxis ()->SetTitle("jet p_{T} [GeV/c]");
	} else if (title=="w_secondvtx_N") {
	  h_ratio->GetXaxis ()->SetTitle("CSV discriminator");
	} else if (title=="w_recoVTX") {
	  h_ratio->GetXaxis ()->SetTitle("Number of offline vertices");
	} else if (title=="w_muon_pt") {
	  h_ratio->GetXaxis ()->SetTitle("muon p_{T} [GeV/c]");
	} else if (title=="w_MET") {
	  h_ratio->GetXaxis ()->SetTitle("MET [GeV/c]");
	} else if (title=="w_MET_sign") {
	  h_ratio->GetXaxis ()->SetTitle("MET Significance [GeV/c]");
	} else if (title=="w_Ht") {
	  h_ratio->GetXaxis ()->SetTitle("H_{T} [GeV/c]");
	} else if (title=="w_pt_Z_ee"||title=="w_pt_Z_mm") {
	  h_ratio->GetXaxis ()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_bjetmultiplicity") {
	  h_ratio->GetXaxis ()->SetTitle("b-jet multiplicity");
	} else if (title=="w_first_jet_pt_b") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet p_{T} [GeV/c]");
	} else if (title=="w_first_jet_eta_b") {
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet #eta");
	} else if (title=="w_second_jet_pt_b") {
	  h_ratio->GetXaxis ()->SetTitle("subleading b-jet p_{T} [GeV/c]");
	} else if (title=="w_second_jet_eta_b") {
	  h_ratio->GetXaxis ()->SetTitle("subleading b-jet #eta");
	} else if (title=="w_third_jet_pt_b") {
	  h_ratio->GetXaxis ()->SetTitle("subsubleading b-jet p_{T} [GeV/c]");
	} else if (title=="w_third_jet_eta_b") {
	  h_ratio->GetXaxis ()->SetTitle("subsubleading b-jet #eta");
	} else if (title=="w_mass_ee_b"||title=="w_mm_mass_b") {
	  h_ratio->GetXaxis ()->SetTitle("Z mass + (#geq 1 b-jet) [GeV/c^{2}]");
	} else if (title=="w_pt_Z_ee_b"||title=="w_pt_Z_mm_b") {
	  h_ratio->GetXaxis ()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_delta_phi_ee"||title=="w_delta_phi_mm") {
	  h_ratio->GetXaxis ()->SetTitle("#Delta#phi(jZ) [rad]");
	} else if (title=="w_delta_phi_ee_b"||title=="w_delta_phi_mm_b") {
	  h_ratio->GetXaxis ()->SetTitle("#Delta#phi(bZ) [rad]");
	} else if (title=="SVTX_mass_jet"||title=="SVTX_mass_trk"||title=="SVTX_mass") {
	  h_ratio->GetXaxis ()->SetTitle("SV mass [GeV/c^{2}]");
	} else if (title=="w_BJP"||title=="w_BPJ") {
	  h_ratio->GetXaxis ()->SetTitle("JP Discriminator");
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
	h_ratio->GetYaxis()->SetTitleOffset(0.4);
	h_ratio->Divide(ht);
	h_ratio->SetMarkerStyle(20);
	h_ratio->Draw("EPX0");

	TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
	OLine->SetLineColor(kRed);
	OLine->SetLineWidth(2);
	OLine->Draw();

	c1->cd();

 	TLatex *latexLabel = CMSPrel(Lumi2012/1000.,"",0.15,0.94);
	latexLabel->Draw("same");

	if (doFit) {
	  TLatex *fitLabel = new TLatex();
	  fitLabel->SetTextSize(0.0275);
	  fitLabel->SetTextFont(42);
	  fitLabel->SetLineWidth(2);
	  fitLabel->SetNDC();
	  char buff[100];
	  if (doFit==1) {
	    sprintf(buff, "c_{t} = %5.3f #pm %5.3f", f1->GetParameter(0), f1->GetParError(0));
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	  }
	  if (doFit==2) {
	    sprintf(buff, "c_{Z+jets} = %5.3f #pm %5.3f", f1->GetParameter(0), f1->GetParError(0));
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	    sprintf(buff, "c_{t} = %5.3f #pm %5.3f", f1->GetParameter(1), f1->GetParError(1));
	    fitLabel->DrawLatex(0.68, 0.43, buff);
	  }
	  if (doFit==3) {
	    sprintf(buff, "c_{uds} = %5.3f #pm %5.3f", f1->GetParameter(0), f1->GetParError(0));
	    fitLabel->DrawLatex(0.38, 0.48, buff);
	    sprintf(buff, "c_{b}   = %5.3f #pm %5.3f", f1->GetParameter(1), f1->GetParError(1));
	    fitLabel->DrawLatex(0.38, 0.43, buff);
	    sprintf(buff, "c_{c}   = %5.3f #pm %5.3f", f1->GetParameter(2), f1->GetParError(2));
	    fitLabel->DrawLatex(0.38, 0.38, buff);
	    float f_uds = 100*h_mc_fit0->Integral()/(h_mc_fit0->Integral()+h_mc_fit1->Integral()+h_mc_fit2->Integral());
	    float ef_uds = f_uds*(f1->GetParError(0)/f1->GetParameter(0));
	    sprintf(buff, "f_{uds} = %4.1f #pm %3.1f %%", f_uds, ef_uds);
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	    float f_b = 100*h_mc_fit1->Integral()/(h_mc_fit0->Integral()+h_mc_fit1->Integral()+h_mc_fit2->Integral());
	    float ef_b = f_b*(f1->GetParError(1)/f1->GetParameter(1));
	    sprintf(buff, "f_{b}   = %4.1f #pm %3.1f %%", f_b, ef_b);
	    fitLabel->DrawLatex(0.68, 0.43, buff);
	    float f_c = 100*h_mc_fit2->Integral()/(h_mc_fit0->Integral()+h_mc_fit1->Integral()+h_mc_fit2->Integral());
	    float ef_c = f_c*(f1->GetParError(2)/f1->GetParameter(2));
	    sprintf(buff, "f_{c}   = %4.1f #pm %3.1f %%", f_c, ef_c);
	    fitLabel->DrawLatex(0.68, 0.38, buff);
	  }
	  fitLabel->Draw("same");
	}

	if (plot) {
	  if (doBkg) title = title + "_doBkg";
	  if (doFit) title = title + "_doFit";
	  if (ilepton==1) {
	    gSystem->mkdir((path + "/electrons/" + version + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electrons/" + version + "/distributions/" + title + ".pdf").c_str());
	  }
	  if (ilepton==2) {
	    gSystem->mkdir((path + "/muons/" + version + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/muons/" + version + "/distributions/" + title + ".pdf").c_str());
	  }
	  if (ilepton==3) {
	    gSystem->mkdir((path + "/electrons+muons/" + version + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electrons+muons/" + version + "/distributions/" + title + ".pdf").c_str());
	  }
	}
}

