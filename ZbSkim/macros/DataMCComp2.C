#include "LumiLabel.C"
#include "LumiInfo_v11.h"

#include "fixrange.C"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";

void DataMCComp2(string& title="", int plot=0, int ilepton=1, int isratio=1, int unfold=0) {

//int useFitResults=0; // use MC predictions for c_b, c_c, c_uds, c_t
int useFitResults=1;  // use fit results for c_b, c_c, c_uds, c_t

//int useEleMuo = 0; // use MC or fit results for c_t
int useEleMuo = 1; // use e-mu fit results for c_t

	if (gROOT->GetVersionInt() >= 53401) {
	  gROOT->GetColor(kRed)->SetAlpha(0.5);
	  gROOT->GetColor(kGreen+2)->SetAlpha(0.5);
	  gROOT->GetColor(kMagenta-6)->SetAlpha(0.5);
	  gROOT->GetColor(kBlue-4)->SetAlpha(0.5);
	}

	/* purity */

	double c_b=1.0;
	double ec_b=0.0;
	double c_c=1.0;
	double ec_c=0.0;
	double c_uds=1.0;
	double ec_uds=0.0;

	/* top */

	double c1_t=1.0;
	double ec1_t=0.0;
	double c2_t=c1_t;
	double ec2_t=ec1_t;

	if (ilepton==1) {
	  if (useFitResults) {
	    c_b   = 0.709;
	    ec_b   = 0.033;
	    c_c   = 1.665;
	    ec_c   = 0.198;
	    c_uds = 1.632;
	    ec_uds = 0.382;
	    c1_t   = 1.080;
	    ec1_t  = 0.022;
	    c2_t   = 0.924;
	    ec2_t  = 0.020;
	    if (useEleMuo) {
	      c1_t  = 0.457;
	      ec1_t = 0.008;
	      c2_t  = 0.438;
	      ec2_t = 0.010;
	    }
	  }
	}

	if (ilepton==2) {
	  if (useFitResults) {
	    c_b   = 0.686;
	    ec_b   = 0.026;
	    c_c   = 1.617;
	    ec_c   = 0.146;
	    c_uds = 1.594;
	    ec_uds = 0.244;
	    c1_t   = 1.028;
	    ec1_t  = 0.019;
	    c2_t  = 0.891;
	    ec2_t = 0.017;
	    if (useEleMuo) {
	      c1_t  = 0.580;
	      ec1_t = 0.010;
	      c2_t  = 0.560;
	      ec2_t = 0.011;
	    }
	  }
	}

	double Lumi2012;

	if (ilepton==1) Lumi2012 = Lumi2012_ele;
	if (ilepton==2) Lumi2012 = Lumi2012_muon;

	double norm1 = ((Lumi2012 * Xsec_dy) / Ngen_dy);
	double norm1_1 = ((Lumi2012 * Xsec_dy_1) / Ngen_dy_1);
	double norm1_2;
	if (ilepton==1) norm1_2 = ((Lumi2012 * Xsec_dy_2) / Ngen_dy_2_ee);
	if (ilepton==2) norm1_2 = ((Lumi2012 * Xsec_dy_2) / Ngen_dy_2_mm);
	double norm2 = ((Lumi2012 * Xsec_tt) / Ngen_tt);
	if (useEleMuo) norm2 = (Lumi2012 / Lumi2012_ele_muon);
	double norm3 = ((Lumi2012 * Xsec_zz) / Ngen_zz);
	double norm4 = ((Lumi2012 * Xsec_wz) / Ngen_wz);
	double norm5 = ((Lumi2012 * Xsec_qcd) / Ngen_qcd);
	double norm6 = ((Lumi2012 * Xsec_ww) / Ngen_ww);
	double norm7 = ((Lumi2012 * Xsec_wj) / Ngen_wj);

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

	TFile *mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL.root").c_str());
	TFile *mcg = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	TFile *mcg1 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());
	TFile *mcg2;
	if (ilepton==1) mcg2 = TFile::Open((path + "/" + version + "/" + "DYToEE_powheg_gen.root").c_str());
	if (ilepton==2) mcg2 = TFile::Open((path + "/" + version + "/" + "DYToMuMu_powheg_gen.root").c_str());
	TFile *mc2 = TFile::Open((path + "/" + version + "/" + "TTbar.root").c_str());
	TFile *mc3 = TFile::Open((path + "/" + version + "/" + "ZZ.root").c_str());
	TFile *mc4 = TFile::Open((path + "/" + version + "/" + "WZ.root").c_str());
//	TFile *mc5 = TFile::Open((path + "/" + version + "/" + "QCD.root").c_str());
	TFile *mc6 = TFile::Open((path + "/" + version + "/" + "WW.root").c_str());
	TFile *mc7 = TFile::Open((path + "/" + version + "/" + "Wj.root").c_str());

	string title_b = title;

	if (title.find("_bjet_")!=string::npos) {
	  title.erase(title.find("_bjet_")+1, 1);
	} else {
	  title_b = title + "_b";
        }

	if (ilepton==1) data->cd("demoEle");
	if (ilepton==2) data->cd("demoMuo");
	TH1F* h_data;
	TH1F* h_data_b;
	if (unfold==0) {
	  h_data = (TH1F*)gDirectory->Get(title.c_str());
	  h_data_b = (TH1F*)gDirectory->Get(title_b.c_str());
	}
	if (unfold==1) {
          if (ilepton==1) {
	    TFile f((path + "/electrons/" + version + "/unfolding/" + title + "_unfolding.root").c_str());
	    TFile f_b((path + "/electrons/" + version + "/unfolding/" + title_b + "_unfolding.root").c_str());
	    h_data = (TH1F*)f.Get(title.c_str())->Clone();
	    h_data_b = (TH1F*)f_b.Get(title_b.c_str())->Clone();
	    h_data->SetDirectory(0);
	    h_data_b->SetDirectory(0);
	    f.Close();
	    f_b.Close();
          }
          if (ilepton==2) {
	    TFile f((path + "/muons/" + version + "/unfolding/" + title + "_unfolding.root").c_str());
	    TFile f_b((path + "/muons/" + version + "/unfolding/" + title_b + "_unfolding.root").c_str());
	    h_data = (TH1F*)f.Get(title.c_str())->Clone();
	    h_data_b = (TH1F*)f_b.Get(title_b.c_str())->Clone();
	    h_data->SetDirectory(0);
	    h_data_b->SetDirectory(0);
	    f.Close();
	    f_b.Close();
          }
	  h_data->SetStats(0);
	  h_data_b->SetStats(0);
	}

	if (ilepton==1) mc1->cd("demoEle");
	if (ilepton==2) mc1->cd("demoMuo");
	TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc1_b = (TH1F*)gDirectory->Get(title_b.c_str());
	TH1F* h_mc1b_b = (TH1F*)gDirectory->Get(("b"+title_b.substr(1)).c_str());
	TH1F* h_mc1c_b = (TH1F*)gDirectory->Get(("c"+title_b.substr(1)).c_str());
	TH1F* h_mc1t_b = (TH1F*)gDirectory->Get(("t"+title_b.substr(1)).c_str());

	if (ilepton==1) mcg->cd("demoEleGen");
	if (ilepton==2) mcg->cd("demoMuoGen");
	TH1F* h_mcg = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mcg_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mcg1->cd("demoEleGen");
	if (ilepton==2) mcg1->cd("demoMuoGen");
	TH1F* h_mcg1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mcg1_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mcg2->cd("demoEleGen");
	if (ilepton==2) mcg2->cd("demoMuoGen");
	TH1F* h_mcg2 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mcg2_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc2->cd("demoEle");
	if (ilepton==2) mc2->cd("demoMuo");
	TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc2_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (useEleMuo) {
	  if (ilepton==1) {
	    mc2 = TFile::Open((path + "/electrons/" + version + "/ttbar_sub/" + title + ".root").c_str());
	    h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	    mc2 = TFile::Open((path + "/electrons/" + version + "/ttbar_sub/" + title_b + ".root").c_str());
	    h_mc2_b = (TH1F*)gDirectory->Get(title_b.c_str());
	  }
	  if (ilepton==2) {
	    mc2 = TFile::Open((path + "/muons/" + version + "/ttbar_sub/" + title + ".root").c_str());
	    h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	    mc2 = TFile::Open((path + "/muons/" + version + "/ttbar_sub/" + title_b + ".root").c_str());
	    h_mc2_b = (TH1F*)gDirectory->Get(title_b.c_str());
	  }
	}

	if (ilepton==1) mc3->cd("demoEle");
	if (ilepton==2) mc3->cd("demoMuo");
	TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc3_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc4->cd("demoEle");
	if (ilepton==2) mc4->cd("demoMuo");
	TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc4_b = (TH1F*)gDirectory->Get(title_b.c_str());

//	if (ilepton==1) mc5->cd("demoEle");
//	if (ilepton==2) mc5->cd("demoMuo");
//	TH1F* h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
//	TH1F* h_mc5_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc6->cd("demoEle");
	if (ilepton==2) mc6->cd("demoMuo");
	TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc6_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc7->cd("demoEle");
	if (ilepton==2) mc7->cd("demoMuo");
	TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc7_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (unfold==0) {
	  h_data->Sumw2();
	  h_data_b->Sumw2();
	}

	h_mc1->Sumw2();
	h_mcg->Sumw2();
	h_mcg1->Sumw2();
	h_mcg2->Sumw2();
	h_mc2->Sumw2();
	h_mc3->Sumw2();
	h_mc4->Sumw2();
//	h_mc5->Sumw2();
	h_mc6->Sumw2();
	h_mc7->Sumw2();

	h_mc1_b->Sumw2();
	if (h_mc1b_b) h_mc1b_b->Sumw2();
	if (h_mc1c_b) h_mc1c_b->Sumw2();
	if (h_mc1t_b) h_mc1t_b->Sumw2();
	h_mcg_b->Sumw2();
	h_mcg1_b->Sumw2();
	h_mcg2_b->Sumw2();
	h_mc2_b->Sumw2();
	h_mc3_b->Sumw2();
	h_mc4_b->Sumw2();
//	h_mc5_b->Sumw2();
	h_mc6_b->Sumw2();
	h_mc7_b->Sumw2();

	h_mc1->Scale(norm1);
	h_mcg->Scale(norm1);
	h_mcg1->Scale(norm1_1);
	h_mcg2->Scale(norm1_2);
	h_mc2->Scale(norm2*c1_t);
	for (int i=0; i<=h_mc2->GetNbinsX()+1; i++) {
	  float e = h_mc2->GetBinError(i)**2;
	  e = e + (h_mc2->GetBinContent(i)*(ec1_t/c1_t))**2;
	  h_mc2->SetBinError(i, TMath::Sqrt(e));
	}
	h_mc3->Scale(norm3);
	h_mc4->Scale(norm4);
//	h_mc5->Scale(norm5);
	h_mc6->Scale(norm6);
	h_mc7->Scale(norm7);

	h_mc1_b->Scale(norm1);
	if (h_mc1b_b) h_mc1b_b->Scale(norm1);
	if (h_mc1c_b) h_mc1c_b->Scale(norm1);
	if (h_mc1t_b) h_mc1t_b->Scale(norm1);
	h_mcg_b->Scale(norm1);
	h_mcg1_b->Scale(norm1_1);
	h_mcg2_b->Scale(norm1_2);
	h_mc2_b->Scale(norm2*c2_t);
	for (int i=0; i<=h_mc2_b->GetNbinsX()+1; i++) {
	  float e = h_mc2_b->GetBinError(i)**2;
	  e = e + (h_mc2_b->GetBinContent(i)*(ec2_t/c2_t))**2;
	  h_mc2_b->SetBinError(i, TMath::Sqrt(e));
	}
	h_mc3_b->Scale(norm3);
	h_mc4_b->Scale(norm4);
//	h_mc5_b->Scale(norm5);
	h_mc6_b->Scale(norm6);
	h_mc7_b->Scale(norm7);

	if (unfold==0) {
	  h_data->Add(h_mc7, -1.);
	  h_data->Add(h_mc6, -1.);
//	  h_data->Add(h_mc5, -1.);
	  h_data->Add(h_mc4, -1.);
	  h_data->Add(h_mc3, -1.);
	  h_data->Add(h_mc2, -1.);

	  h_data_b->Add(h_mc7_b, -1.);
	  h_data_b->Add(h_mc6_b, -1.);
//	  h_data_b->Add(h_mc5_b, -1.);
	  h_data_b->Add(h_mc4_b, -1.);
	  h_data_b->Add(h_mc3_b, -1.);
	  h_data_b->Add(h_mc2_b, -1.);
	}

	TH1F *h_mc1uds_b = h_mc1_b->Clone("h_mc1uds_b");
	if (h_mc1b_b) h_mc1uds_b->Add(h_mc1b_b, -1);
	if (h_mc1c_b) h_mc1uds_b->Add(h_mc1c_b, -1);
	if (h_mc1t_b) h_mc1uds_b->Add(h_mc1t_b, -1);
	for (int i=0; i<=h_mc1uds_b->GetNbinsX()+1; i++) {
	  float e = h_mc1uds_b->GetBinError(i)**2;
	  if (h_mc1b_b) e = e - h_mc1b_b->GetBinError(i)**2;
	  if (h_mc1c_b) e = e - h_mc1c_b->GetBinError(i)**2;
	  if (h_mc1t_b) e = e - h_mc1t_b->GetBinError(i)**2;
	  h_mc1uds_b->SetBinError(i, TMath::Sqrt(e));
	}

	if (h_mc1uds_b) {
	  h_mc1uds_b->Scale(c_uds);
	  for (int i=0; i<=h_mc1uds_b->GetNbinsX()+1; i++) {
	    float e = h_mc1uds_b->GetBinError(i)**2;
	    e = e + (h_mc1uds_b->GetBinContent(i)*(ec_uds/c_uds))**2;
	    h_mc1uds_b->SetBinError(i, TMath::Sqrt(e));
	  }
	}
	if (h_mc1b_b) {
	  h_mc1b_b->Scale(c_b);
	  for (int i=0; i<=h_mc1b_b->GetNbinsX()+1; i++) {
	    float e = h_mc1b_b->GetBinError(i)**2;
	    e = e + (h_mc1b_b->GetBinContent(i)*(ec_b/c_b))**2;
	    h_mc1b_b->SetBinError(i, TMath::Sqrt(e));
	  }
	}
	if (h_mc1c_b) {
	  h_mc1c_b->Scale(c_c);
	  for (int i=0; i<=h_mc1c_b->GetNbinsX()+1; i++) {
	    float e = h_mc1c_b->GetBinError(i)**2;
	    e = e + (h_mc1c_b->GetBinContent(i)*(ec_c/c_c))**2;
	    h_mc1c_b->SetBinError(i, TMath::Sqrt(e));
	  }
	}

	if (unfold==0) {
	  h_data_b->Add(h_mc1c_b, -1.);
	  h_data_b->Add(h_mc1uds_b, -1.);
	}

	TH1F *h_data_raw;
	TH1F *h_data_b_raw;
	if (unfold==0) {
	  h_data_raw = (TH1F*)h_data->Clone();
	  h_data_b_raw = (TH1F*)h_data_b->Clone();
	}

        if (ilepton==1) {
	  TFile f((path + "/electrons/" + version + "/efficiency/" + string(h_data->GetName()) + "_efficiency.root").c_str());
	  TFile f_b((path + "/electrons/" + version + "/efficiency/" + string(h_data_b->GetName()) + "_efficiency.root").c_str());
	  TH1F* h = (TH1F*)f.Get(h_data->GetName())->Clone();
	  TH1F* h_b = (TH1F*)f_b.Get(h_data_b->GetName())->Clone();
	  h->SetDirectory(0);
	  h_b->SetDirectory(0);
	  f.Close();
	  f_b.Close();
	  if (unfold==0) {
	    h_data->Divide(h);
	    h_data_b->Divide(h_b);
	  }
	  h_mc1->Divide(h);
	  h_mc1b_b->Divide(h_b);
        }
	if (ilepton==2) {
	  TFile f((path + "/muons/" + version + "/efficiency/" + string(h_data->GetName()) + "_efficiency.root").c_str());
	  TFile f_b((path + "/muons/" + version + "/efficiency/" + string(h_data_b->GetName()) + "_efficiency.root").c_str());
	  TH1F* h = (TH1F*)f.Get(h_data->GetName())->Clone();
	  TH1F* h_b = (TH1F*)f_b.Get(h_data_b->GetName())->Clone();
	  h->SetDirectory(0);
	  h_b->SetDirectory(0);
	  f.Close();
	  f_b.Close();
	  if (unfold==0) {
	    h_data->Divide(h);
	    h_data_b->Divide(h_b);
	  }
	  h_mc1->Divide(h);
	  h_mc1b_b->Divide(h_b);
        }

	h_data->Scale(1./Lumi2012, "width");
	h_data_b->Scale(1./Lumi2012, "width");
	h_mc1->Scale(1./Lumi2012, "width");
	h_mc1b_b->Scale(1./Lumi2012, "width");
	if (isratio==1) {
	  h_data_b->Divide(h_data);
	  h_data_b->Scale(100.);
	  h_mc1b_b->Divide(h_mc1);
	  h_mc1b_b->Scale(100.);
	}

	h_mcg->Scale(1./Lumi2012, "width");
	h_mcg1->Scale(1./Lumi2012, "width");
	h_mcg2->Scale(1./Lumi2012, "width");
	h_mcg_b->Scale(1./Lumi2012, "width");
	h_mcg1_b->Scale(1./Lumi2012, "width");
	h_mcg2_b->Scale(1./Lumi2012, "width");
	if (isratio==1) {
	  h_mcg_b->Divide(h_mcg);
	  h_mcg1_b->Divide(h_mcg1);
	  h_mcg2_b->Divide(h_mcg2);
	  h_mcg_b->Scale(100.);
	  h_mcg1_b->Scale(100.);
	  h_mcg2_b->Scale(100.);
	}

	h_data = fixrange(h_data);
	h_data_b = fixrange(h_data_b);
	h_mc1 = fixrange(h_mc1);
	h_mc1b_b = fixrange(h_mc1b_b);
	h_mcg = fixrange(h_mcg);
	h_mcg_b = fixrange(h_mcg_b);
	h_mcg1 = fixrange(h_mcg1);
	h_mcg1_b = fixrange(h_mcg1_b);
	h_mcg2 = fixrange(h_mcg2);
	h_mcg2_b = fixrange(h_mcg2_b);

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
	c1->cd();

	h_mc1b_b->SetTitle("");
	if (isratio==1) {
	  h_mc1b_b->GetYaxis()->SetTitle("#sigma(Z+b) / #sigma(Z+j) [%]");
	} else {
	  h_mc1b_b->GetYaxis()->SetTitle("#sigma [pb]");
	}

	TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
	pad1->SetTopMargin(0.115);
	pad1->SetBottomMargin(0.0001);
	pad1->Draw();
	pad1->cd();

	h_mc1b_b->SetLineColor(kRed);
	h_mc1b_b->SetLineWidth(2);
	h_mc1b_b->SetMarkerColor(kRed);
	h_mc1b_b->SetFillColor(kRed);
	h_mc1b_b->SetStats(0);
	if (isratio==1) {
	  h_mc1b_b->Draw("E5");
	}

	h_mcg_b->SetLineColor(kGreen+2);
	h_mcg_b->SetLineWidth(2);
	h_mcg_b->SetFillColor(kGreen+2);
	h_mcg_b->SetMarkerColor(kGreen+2);
	if (isratio==1) {
	  h_mcg_b->Draw("E5SAME");
	}

	h_mcg1_b->SetLineColor(kMagenta-6);
	h_mcg1_b->SetLineWidth(2);
	h_mcg1_b->SetFillColor(kMagenta-6);
	h_mcg1_b->SetMarkerColor(kMagenta-6);
	if (isratio==1) {
	  h_mcg1_b->Draw("E5SAME");
	}

	h_mcg2_b->SetLineColor(kBlue-4);
	h_mcg2_b->SetLineWidth(2);
	h_mcg2_b->SetFillColor(kBlue-4);
	h_mcg2_b->SetMarkerColor(kBlue-4);
	if (isratio==1) {
	  h_mcg2_b->Draw("E5SAME");
	}

	if (isratio==1) {
	  h_data_b->GetYaxis()->SetTitle("#sigma_{Z+b-jets}/#sigma_{Z+jets} [%]");
	}
	h_data_b->GetYaxis()->SetTitleOffset(1.2);
	h_data_b->GetXaxis()->SetTitleOffset(1.3);
	h_data_b->SetMarkerColor(kBlack);
	h_data_b->SetLineColor(kBlack);
	h_data_b->SetMarkerStyle(24);
	h_data_b->SetMarkerSize(0.7);
	h_data_b->SetStats(0);
	if (isratio==1) {
	  h_data_b->Draw("EPX0SAME");
	}

	TLegend *leg = new TLegend(0.62, 0.580, 0.88, 0.88);
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (isratio==0) {
	  pad1->SetLogy();

	  h_mc1b_b->SetMaximum(4*h_data->GetMaximum());
	  h_mc1b_b->SetMinimum(TMath::Max(0.000002,0.25*h_data_b->GetBinContent(h_data_b->GetMinimumBin())));

	  h_mc1b_b->Draw("E5");
	  TH1F* tmp1 = h_mc1b_b->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp1->GetMinimum()==0) tmp1->GetXaxis()->SetRangeUser(0, tmp1->GetBinCenter(tmp1->GetMinimumBin()-1));
	  }
	  tmp1->SetFillColor(0);
	  tmp1->DrawClone("HISTLSAME");

	  h_mcg_b->Draw("E5SAME");
	  TH1F* tmp2 = h_mcg_b->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp2->GetMinimum()==0) tmp2->GetXaxis()->SetRangeUser(0, tmp2->GetBinCenter(tmp2->GetMinimumBin()-1));
	  }
	  tmp2->SetFillColor(0);
	  tmp2->DrawClone("HISTLSAME");

	  h_mcg1_b->Draw("E5SAME");
	  TH1F* tmp2_1 = h_mcg1_b->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp2_1->GetMinimum()==0) tmp2_1->GetXaxis()->SetRangeUser(0, tmp2_1->GetBinCenter(tmp2_1->GetMinimumBin()-1));
	  }
	  tmp2_1->SetFillColor(0);
	  tmp2_1->DrawClone("HISTLSAME");

	  h_mcg2_b->Draw("E5SAME");
	  TH1F* tmp2_2 = h_mcg2_b->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp2_2->GetMinimum()==0) tmp2_2->GetXaxis()->SetRangeUser(0, tmp2_2->GetBinCenter(tmp2_2->GetMinimumBin()-1));
	  }
	  tmp2_2->SetFillColor(0);
	  tmp2_2->DrawClone("HISTLSAME");

	  h_data_b->Draw("EPX0SAME");

	  h_mc1->SetLineColor(kRed);
	  h_mc1->SetLineWidth(2);
	  h_mc1->SetMarkerColor(kRed);
	  h_mc1->SetFillColor(kRed);
	  h_mc1->Draw("E5SAME");
	  TH1F* tmp3 = h_mc1->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp3->GetMinimum()==0) tmp3->GetXaxis()->SetRangeUser(0, tmp3->GetBinCenter(tmp3->GetMinimumBin()-1));
	  }
	  tmp3->SetFillColor(0);
	  tmp3->DrawClone("HISTLSAME");

	  h_mcg->SetLineColor(kGreen+2);
	  h_mcg->SetLineWidth(2);
	  h_mcg->SetMarkerColor(kGreen+2);
	  h_mcg->SetFillColor(kGreen+2);
	  h_mcg->Draw("E5SAME");
	  TH1F* tmp4 = h_mcg->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp4->GetMinimum()==0) tmp4->GetXaxis()->SetRangeUser(0, tmp4->GetBinCenter(tmp4->GetMinimumBin()-1));
	  }
	  tmp4->SetFillColor(0);
	  tmp4->DrawClone("HISTLSAME");

	  h_mcg1->SetLineColor(kMagenta-6);
	  h_mcg1->SetLineWidth(2);
	  h_mcg1->SetMarkerColor(kMagenta-6);
	  h_mcg1->SetFillColor(kMagenta-6);
	  h_mcg1->Draw("E5SAME");
	  TH1F* tmp4_1 = h_mcg1->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp4_1->GetMinimum()==0) tmp4_1->GetXaxis()->SetRangeUser(0, tmp4_1->GetBinCenter(tmp4_1->GetMinimumBin()-1));
	  }
	  tmp4_1->SetFillColor(0);
	  tmp4_1->DrawClone("HISTLSAME");

	  h_mcg2->SetLineColor(kBlue-4);
	  h_mcg2->SetLineWidth(2);
	  h_mcg2->SetMarkerColor(kBlue-4);
	  h_mcg2->SetFillColor(kBlue-4);
	  h_mcg2->Draw("E5SAME");
	  TH1F* tmp4_2 = h_mcg2->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp4_2->GetMinimum()==0) tmp4_2->GetXaxis()->SetRangeUser(0, tmp4_2->GetBinCenter(tmp4_2->GetMinimumBin()-1));
	  }
	  tmp4_2->SetFillColor(0);
	  tmp4_2->DrawClone("HISTLSAME");

	  h_data->SetMarkerColor(kBlack);
	  h_data->SetLineColor(kBlack);
	  h_data->SetMarkerStyle(20);
	  h_data->SetMarkerSize (0.7);
	  h_data->Draw("EPX0SAME");

	  if (ilepton==1) {
	    leg->AddEntry(h_data,"Z(#rightarrow ee) DATA","p");
	    leg->AddEntry(h_data_b,"Z(#rightarrow ee)+b DATA","p");
	    leg->AddEntry(h_mc1,"Z(#rightarrow ee) MC","l");
	    leg->AddEntry(h_mcg,"Z(#rightarrow ee) MadGraph","l");
	    leg->AddEntry(h_mcg1,"Z(#rightarrow ee) Sherpa","l");
	    leg->AddEntry(h_mcg2,"Z(#rightarrow ee) Powheg","l");
	  }
	  if (ilepton==2){
	    leg->AddEntry(h_data,"Z(#rightarrow #mu#mu) DATA","p");
	    leg->AddEntry(h_data_b,"Z(#rightarrow #mu#mu)+b DATA","p");
	    leg->AddEntry(h_mc1,"Z(#rightarrow #mu#mu) MC","l");
	    leg->AddEntry(h_mcg,"Z(#rightarrow #mu#mu) MadGraph","l");
	    leg->AddEntry(h_mcg1,"Z(#rightarrow #mu#mu) Sherpa","l");
	    leg->AddEntry(h_mcg2,"Z(#rightarrow #mu#mu) Powheg","l");
	  }
	}

	if (isratio==1) {
	  if (ilepton==1) {
	    leg->AddEntry(h_data_b,"Z(#rightarrow ee) DATA","p");
	    leg->AddEntry(h_mc1b_b,"Z(#rightarrow ee) MC","l");
	    leg->AddEntry(h_mcg_b,"Z(#rightarrow ee) MadGraph","l");
	    leg->AddEntry(h_mcg1_b,"Z(#rightarrow ee) Sherpa","l");
	    leg->AddEntry(h_mcg2_b,"Z(#rightarrow ee) Powheg","l");
	  }
	  if (ilepton==2){
	    leg->AddEntry(h_data_b,"Z(#rightarrow #mu#mu) DATA","p");
	    leg->AddEntry(h_mc1b_b,"Z(#rightarrow #mu#mu) MC","l");
	    leg->AddEntry(h_mcg_b,"Z(#rightarrow #mu#mu) MadGraph","l");
	    leg->AddEntry(h_mcg1_b,"Z(#rightarrow #mu#mu) Sherpa","l");
	    leg->AddEntry(h_mcg2_b,"Z(#rightarrow #mu#mu) Powheg","l");
	  }
	}

	leg->Draw();

	c1->cd();

 	TLatex *latexLabel = CMSPrel(Lumi2012/1000.,"",0.15,0.94);
	latexLabel->Draw("same");

	TPad *pad2 = new TPad("pad2","pad2",0,0.29,1,0.4);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.001);
	pad2->Draw();
	pad2->cd();

	TH1F *h_M = h_data_b->Clone("h_M");
	h_M->SetTitle("");
	h_M->SetStats(0);
	h_M->GetXaxis()->SetTitleOffset(0.9);
	h_M->GetXaxis()->SetTitleSize(0.14);
	h_M->GetXaxis()->SetLabelFont(42);
	h_M->GetXaxis()->SetLabelSize(0.12);
	h_M->GetXaxis()->SetTitleFont(42);
	h_M->GetXaxis()->SetTickLength(0.1);
	h_M->GetYaxis()->SetTitle("Data / Theory");
	h_M->GetYaxis()->SetNdivisions(013);
	h_M->GetYaxis()->SetTitleSize(0.17);
	h_M->GetYaxis()->SetLabelSize(0.17);
	h_M->GetYaxis()->SetRangeUser(-0.2, 2.2);
	h_M->GetYaxis()->SetTitleOffset(0.21);
	h_M->GetYaxis()->SetTickLength(0.02);
	h_M->Divide(h_mcg_b);
	h_M->SetMarkerStyle(24);
	h_M->Draw("EPX0");
	if (isratio==0) {
	  TH1F *h_M2= h_data->Clone("h_M2");
	  h_M2->GetXaxis()->SetTitleOffset(0.9);
	  h_M2->GetXaxis()->SetTitleSize(0.14);
	  h_M2->GetXaxis()->SetLabelFont(42);
	  h_M2->GetXaxis()->SetLabelSize(0.12);
	  h_M2->GetXaxis()->SetTitleFont(42);
	  h_M2->GetXaxis()->SetTickLength(0.1);
	  h_M2->GetYaxis()->SetTitle("Data / Theory");
	  h_M2->GetYaxis()->SetNdivisions(013);
	  h_M2->GetYaxis()->SetTitleSize(0.17);
	  h_M2->GetYaxis()->SetLabelSize(0.17);
	  h_M2->GetYaxis()->SetRangeUser(-0.2, 2.2);
	  h_M2->GetYaxis()->SetTitleOffset(0.21);
	  h_M2->GetYaxis()->SetTickLength(0.02);
	  h_M2->Divide(h_mcg);
	  h_M2->SetMarkerStyle(20);
	  h_M2->Draw("EPX0SAME");
	}

	TLatex *t2 = new TLatex();
	t2->SetTextSize(0.2);
	t2->SetTextFont(42);
	t2->SetLineWidth(2);
	t2->SetNDC();
	t2->DrawLatex(0.15,0.7,"MadGraph");

	TLine *OLine2 = new TLine(h_M->GetXaxis()->GetXmin(),1.,h_M->GetXaxis()->GetXmax(),1.);
	OLine2->SetLineColor(kGreen+2);
	OLine2->SetLineWidth(2);
	OLine2->Draw();

	c1->cd();

	TPad *pad3 = new TPad("pad3","pad3",0,0.18,1,0.29);
	pad3->SetTopMargin(0);
	pad3->SetBottomMargin(0.001);
	pad3->Draw();
	pad3->cd();

	TH1F *h_S = h_data_b->Clone("h_S");
	h_S->SetTitle("");
	h_S->SetStats(0);
	h_S->GetXaxis()->SetTitleOffset(0.9);
	h_S->GetXaxis()->SetTitleSize(0.14);
	h_S->GetXaxis()->SetLabelFont(42);
	h_S->GetXaxis()->SetLabelSize(0.12);
	h_S->GetXaxis()->SetTitleFont(42);
	h_S->GetXaxis()->SetTickLength(0.1);
	h_S->GetYaxis()->SetTitle("Data / Theory");
	h_S->GetYaxis()->SetNdivisions(013);
	h_S->GetYaxis()->SetTitleSize(0.17);
	h_S->GetYaxis()->SetLabelSize(0.17);
	h_S->GetYaxis()->SetRangeUser(-0.2, 2.2);
	h_S->GetYaxis()->SetTitleOffset(0.21);
	h_S->GetYaxis()->SetTickLength(0.02);
	h_S->Divide(h_mcg1_b);
	h_S->SetMarkerStyle(24);
	h_S->Draw("EPX0");
	if (isratio==0) {
	  TH1F *h_S2= h_data->Clone("h_S2");
	  h_S2->GetXaxis()->SetTitleOffset(0.9);
	  h_S2->GetXaxis()->SetTitleSize(0.14);
	  h_S2->GetXaxis()->SetLabelFont(42);
	  h_S2->GetXaxis()->SetLabelSize(0.12);
	  h_S2->GetXaxis()->SetTitleFont(42);
	  h_S2->GetXaxis()->SetTickLength(0.1);
	  h_S2->GetYaxis()->SetTitle("Data / Theory");
	  h_S2->GetYaxis()->SetNdivisions(013);
	  h_S2->GetYaxis()->SetTitleSize(0.17);
	  h_S2->GetYaxis()->SetLabelSize(0.17);
	  h_S2->GetYaxis()->SetRangeUser(-0.2, 2.2);
	  h_S2->GetYaxis()->SetTitleOffset(0.21);
	  h_S2->GetYaxis()->SetTickLength(0.02);
	  h_S2->Divide(h_mcg1);
	  h_S2->SetMarkerStyle(20);
	  h_S2->Draw("EPX0SAME");
	}

	TLatex *t3 = new TLatex();
	t3->SetTextSize(0.2);
	t3->SetTextFont(42);
	t3->SetLineWidth(2);
	t3->SetNDC();
	t3->DrawLatex(0.15,0.7,"Sherpa");

	TLine *OLine3 = new TLine(h_S->GetXaxis()->GetXmin(),1.,h_S->GetXaxis()->GetXmax(),1.);
	OLine3->SetLineColor(kMagenta-6);
	OLine3->SetLineWidth(2);
	OLine3->Draw();

	c1->cd();

	TPad *pad4 = new TPad("pad4","pad4",0,0.0,1,0.18);
	pad4->SetTopMargin(0);
	pad4->SetBottomMargin(0.3);
	pad4->Draw();
	pad4->cd();

	TH1F *h_P = h_data_b->Clone("h_P");
	h_P->SetTitle("");
	h_P->SetStats(0);
	h_P->GetXaxis()->SetTitleOffset(0.9);
	h_P->GetXaxis()->SetTitleSize(0.14);
	h_P->GetXaxis()->SetLabelFont(42);
	h_P->GetXaxis()->SetLabelSize(0.12);
	h_P->GetXaxis()->SetTitleFont(42);
	h_P->GetXaxis()->SetTickLength(0.1);
	h_P->GetYaxis()->SetTitle("Data / Theory");
	h_P->GetYaxis()->SetNdivisions(013);
	h_P->GetYaxis()->SetTitleSize(0.11);
	h_P->GetYaxis()->SetLabelSize(0.11);
	h_P->GetYaxis()->SetRangeUser(-0.2, 2.2);
	h_P->GetYaxis()->SetTitleOffset(0.32);
	h_P->GetYaxis()->SetTickLength(0.02);
	h_P->Divide(h_mcg2_b);
	h_P->SetMarkerStyle(24);
	h_P->Draw("EPX0");
	if (isratio==0) {
	  TH1F *h_P2= h_data->Clone("h_P2");
	  h_P2->GetXaxis()->SetTitleOffset(0.9);
	  h_P2->GetXaxis()->SetTitleSize(0.14);
	  h_P2->GetXaxis()->SetLabelFont(42);
	  h_P2->GetXaxis()->SetLabelSize(0.12);
	  h_P2->GetXaxis()->SetTitleFont(42);
	  h_P2->GetXaxis()->SetTickLength(0.1);
	  h_P2->GetYaxis()->SetTitle("Data / Theory");
	  h_P2->GetYaxis()->SetNdivisions(013);
	  h_P2->GetYaxis()->SetTitleSize(0.11);
	  h_P2->GetYaxis()->SetLabelSize(0.11);
	  h_P2->GetYaxis()->SetRangeUser(-0.2, 2.2);
	  h_P2->GetYaxis()->SetTitleOffset(0.32);
	  h_P2->GetYaxis()->SetTickLength(0.02);
	  h_P2->Divide(h_mcg2);
	  h_P2->SetMarkerStyle(20);
	  h_P2->Draw("EPX0SAME");
	}

	TLatex *t4 = new TLatex();
	t4->SetTextSize(0.13);
	t4->SetTextFont(42);
	t4->SetLineWidth(2);
	t4->SetNDC();
	t4->DrawLatex(0.15,0.8,"Powheg");

	TLine *OLine4 = new TLine(h_P->GetXaxis()->GetXmin(),1.,h_P->GetXaxis()->GetXmax(),1.);
	OLine4->SetLineColor(kBlue-4);
	OLine4->SetLineWidth(2);
	OLine4->Draw();

	c1->cd();

	if (isratio==1) {
	  h_mc1b_b->GetYaxis()->SetRangeUser(-0.5, 10);
	}

	if (title_b=="w_first_jet_pt_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dp_{T} [pb]");
	  h_P->GetXaxis()->SetTitle("leading jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp_{T} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(-0.5, 20);
	  }
	} else if (title_b=="w_first_jet_eta_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#eta [pb]");
	  h_P->GetXaxis()->SetTitle("leading jet #eta");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_first_bjet_pt") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dp^{b}_{T} [pb]");
	  h_P->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{b}_{T} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(-0.5, 10);
	  }
	} else if (title_b=="w_first_bjet_eta") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#eta^{b} [pb]");
	  h_P->GetXaxis()->SetTitle("leading b-jet #eta");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta^{b} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_pt_Z_ee_b"||title_b =="w_pt_Z_mm_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dp^{Z}_{T} [pb]");
	  h_P->GetXaxis()->SetTitle("Z boson p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{Z}_{T} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(-0.5, 20);
	  }
	} else if (title_b=="w_Ht_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dH_{T} [pb]");
	  h_P->GetXaxis()->SetTitle("H_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dH_{T} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_delta_phi_ee_b" || title_b=="w_delta_phi_mm_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#Delta#phi_{bZ} [pb]");
	  h_P->GetXaxis()->SetTitle("#Delta#phi(bZ) [rad]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta#phi_{Zb} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	}

	if (plot) {
	  if (unfold==0) {
	    if (isratio==0) {
	      if (ilepton==1) {
	        gSystem->mkdir((path + "/electrons/" + version + "/xsecs/").c_str(), kTRUE);
	        c1->SaveAs((path + "/electrons/" + version + "/xsecs/" + title_b + "_xsecs.pdf").c_str());
	        TFile f((path + "/electrons/" + version + "/xsecs/" + title_b + "_xsecs.root").c_str(),"RECREATE");
	        h_data_raw->Write((title+"_raw").c_str());
                h_data_b_raw->Write((title_b+"_raw").c_str());
                h_data->Write(title.c_str());
                h_data_b->Write(title_b.c_str());
                f.Close();
	      }
	      if (ilepton==2) {
	        gSystem->mkdir((path + "/muons/" + version + "/xsecs/").c_str(), kTRUE);
	        c1->SaveAs((path + "/muons/" + version + "/xsecs/" + title_b + "_xsecs.pdf").c_str());
	        TFile f((path + "/muons/" + version + "/xsecs/" + title_b + "_xsecs.root").c_str(),"RECREATE");
	        h_data_raw->Write((title+"_raw").c_str());
                h_data_b_raw->Write((title_b+"_raw").c_str());
                h_data->Write(title.c_str());
                h_data_b->Write(title_b.c_str());
                f.Close();
	      }
	    }
	    if (isratio==1) {
	      if (ilepton==1) {
	        gSystem->mkdir((path + "/electrons/" + version + "/ratios/").c_str(), kTRUE);
	        c1->SaveAs((path + "/electrons/" + version + "/ratios/" + title_b + "_ratio.pdf").c_str());
	      }
	      if (ilepton==2) {
	        gSystem->mkdir((path + "/muons/" + version + "/ratios/").c_str(), kTRUE);
	        c1->SaveAs((path + "/muons/" + version + "/ratios/" + title_b + "_ratio.pdf").c_str());
	      }
	    }
	  }
	  if (unfold==1) {
	    if (isratio==0) {
	      if (ilepton==1) {
	        gSystem->mkdir((path + "/electrons/" + version + "/xsecs_unfolding/").c_str(), kTRUE);
	        c1->SaveAs((path + "/electrons/" + version + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.pdf").c_str());
	      }
	      if (ilepton==2) {
	        gSystem->mkdir((path + "/muons/" + version + "/xsecs_unfolding/").c_str(), kTRUE);
	        c1->SaveAs((path + "/muons/" + version + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.pdf").c_str());
	      }
	    }
	    if (isratio==1) {
	      if (ilepton==1) {
	        gSystem->mkdir((path + "/electrons/" + version + "/ratios_unfolding/").c_str(), kTRUE);
	        c1->SaveAs((path + "/electrons/" + version + "/ratios_unfolding/" + title_b + "_ratio_unfolding.pdf").c_str());
	      }
	      if (ilepton==2) {
	        gSystem->mkdir((path + "/muons/" + version + "/ratios_unfolding/").c_str(), kTRUE);
	        c1->SaveAs((path + "/muons/" + version + "/ratios_unfolding/" + title_b + "_ratio_unfolding.pdf").c_str());
	      }
	    }
	  }
	}
}
