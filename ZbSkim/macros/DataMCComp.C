#include "LumiLabel.C"
#include "LumiInfo_v11.h"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";

TH1F* h_data_fit = 0;
TH1F* h_mc_fit0 = 0;
TH1F* h_mc_fit1 = 0;
TH1F* h_mc_fit2 = 0;

void fcn(int& npar, double* gin, double& fun, double* par, int iflag) {
  double chisq = 0.0;
  for (int i=1; i<=h_data_fit->GetNbinsX(); i++) {
    double xn = h_data_fit->GetBinContent(i);
    double xd = h_data_fit->GetBinError(i)**2;
    if (npar>0) {
      xn = xn - par[0]*h_mc_fit0->GetBinContent(i);
      xd = xd + (par[0]*h_mc_fit0->GetBinError(i))**2;
    }
    if (npar>1) {
      xn = xn - par[1]*h_mc_fit1->GetBinContent(i);
      xd = xd + (par[1]*h_mc_fit1->GetBinError(i))**2;
    }
    if (npar>2) {
      xn = xn - par[2]*h_mc_fit2->GetBinContent(i);
      xd = xd + (par[2]*h_mc_fit2->GetBinError(i))**2;
    }
    if (xd!=0) chisq = chisq + (xn*xn)/xd;
  }
  fun = chisq;
}

void DataMCComp(int irun=0, string& title="", int plot=0, int ilepton=1, int doBkg=0, int doFit=0) {

//int useFitResults=0; // use MC predictions for c_t
int useFitResults=1;  // use fit results for c_t

//int useEleMuo = 0; // use MC or fit results for c_t
int useEleMuo = 1; // use e-mu fit results for c_t

string subdir="0";
string postfix="";
if (irun==1) {             // irun==1 => JEC Up
  string subdir="1";
  string postfix="Up";
}
if (irun==2) {             // irun==2 => JEC Down
  string subdir="2";
  string postfix="Down";  
}
if (irun==3) {             // irun==3 => PU Up
  string subdir="3";
  string postfix="Pup";  
}
if (irun==4) {             // irun==4 => PU Down
  string subdir="4";
  string postfix="Pum";  
}
if (irun==5) {             // irun==5 => top bkg
  string subdir="5";
  string postfix="";  
}
if (irun==6) {             // irun==6 => b purity
  string subdir="6";
  string postfix="";  
}
if (irun==7) {             // irun==7 => unfolding
  string subdir="7";
  string postfix="";  
}
if (irun==8) {             // irun==8 => unfolding with Sherpa
  string subdir="8";
  string postfix="";
}
if (irun==9) {             // irun==9 => unfolding with Powheg
  string subdir="9";
  string postfix="";
}
if (irun==10) {            // irun==10 => bkg
  string subdir="10";
  string postfix="";
}

	/* top background */

	double c1_t=1.0;
	double ec1_t=0.0;
	double c2_t=1.0;
	double ec2_t=0.0;

	if (doFit==1) {
	  if (title=="w_MET") { useFitResults=0; useEleMuo=0; }
	  if (title=="w_MET_b") { useFitResults=0; useEleMuo=0; }
	}

	ifstream in4, in5, in6, in7;
	if (ilepton==1) {
	  if (useFitResults) {
	    in4.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_MET_doFit" + ".dat").c_str());
	    in5.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_MET_b_doFit" + ".dat").c_str());
	    if (useEleMuo) {
	      in6.open((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub/" + "w_mass_ee_wide_doFit" + ".dat").c_str());
	      in7.open((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub/" + "w_mass_ee_b_wide_doFit" + ".dat").c_str());
	    }
	  }
	}
	if (ilepton==2) {
	  if (useFitResults) {
	    in4.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_MET_doFit" + ".dat").c_str());
	    in5.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_MET_b_doFit" + ".dat").c_str());
	    if (useEleMuo) {
	      in6.open((path + "/muons/" + version + "/" + subdir + "/ttbar_sub/" + "w_mass_ee_wide_doFit" + ".dat").c_str());
	      in7.open((path + "/muons/" + version + "/" + subdir + "/ttbar_sub/" + "w_mass_ee_b_wide_doFit" + ".dat").c_str());
	    }
	  }
	}
	if (ilepton==3) {
	  if (useFitResults) {
	    in4.open((path + "/electrons+muons/" + version + "/" + subdir + "/distributions/" + "w_MET_doFit" + ".dat").c_str());
	    in5.open((path + "/electrons+muons/" + version + "/" + subdir + "/distributions/" + "w_MET_b_doFit" + ".dat").c_str());
	  }
	}
	if (useFitResults) {
	  in4 >> c1_t >> ec1_t;
	  in5 >> c2_t >> ec2_t;
	  in4.close();
	  in5.close();
	  if (useEleMuo) {
	    in6 >> c1_t >> ec1_t;
	    in7 >> c2_t >> ec2_t;
	    in6.close();
	    in7.close();
	  }
	}

	double Lumi2012;

	if (ilepton==1) Lumi2012 = Lumi2012_ele;
	if (ilepton==2) Lumi2012 = Lumi2012_muon;
	if (ilepton==3) Lumi2012 = Lumi2012_ele_muon;

	double norm1 = ((Lumi2012 * Xsec_dy) / Ngen_dy);
	double norm2 = ((Lumi2012 * Xsec_tt) / Ngen_tt);
	if (useEleMuo && ilepton!=3) norm2 = (Lumi2012 / Lumi2012_ele_muon);
	double norm3 = ((Lumi2012 * Xsec_zz) / Ngen_zz);
	double norm4 = ((Lumi2012 * Xsec_wz) / Ngen_wz);
	double norm5 = ((Lumi2012 * Xsec_qcd) / Ngen_qcd);
	double norm6 = ((Lumi2012 * Xsec_ww) / Ngen_ww);
	double norm7 = ((Lumi2012 * Xsec_wj) / Ngen_wj);

	double enorm1 = ((Lumi2012 * eXsec_dy) / Ngen_dy);
	double enorm2 = ((Lumi2012 * eXsec_tt) / Ngen_tt);
	double enorm3 = ((Lumi2012 * eXsec_zz) / Ngen_zz);
	double enorm4 = ((Lumi2012 * eXsec_wz) / Ngen_wz);
	double enorm5 = ((Lumi2012 * eXsec_qcd) / Ngen_qcd);
	double enorm6 = ((Lumi2012 * eXsec_ww) / Ngen_ww);
	double enorm7 = ((Lumi2012 * eXsec_wj) / Ngen_wj);

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
	if (ilepton==3) data = TFile::Open((path + "/" + version + "/" + "MuEG_2012_merge.root").c_str());

	TFile *mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL.root").c_str());
	TFile *mc2 = TFile::Open((path + "/" + version + "/" + "TTbar.root").c_str());
	TFile *mc3 = TFile::Open((path + "/" + version + "/" + "ZZ.root").c_str());
	TFile *mc4 = TFile::Open((path + "/" + version + "/" + "WZ.root").c_str());
//	TFile *mc5 = TFile::Open((path + "/" + version + "/" + "QCD.root").c_str());
	TFile *mc6 = TFile::Open((path + "/" + version + "/" + "WW.root").c_str());
	TFile *mc7 = TFile::Open((path + "/" + version + "/" + "Wj.root").c_str());

	if (ilepton==1) data->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) data->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) data->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_data = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc1->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc1->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc1->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc1b = (TH1F*)gDirectory->Get(("b"+title.substr(1)).c_str());
	TH1F* h_mc1c = (TH1F*)gDirectory->Get(("c"+title.substr(1)).c_str());
	TH1F* h_mc1t = (TH1F*)gDirectory->Get(("t"+title.substr(1)).c_str());

	if (ilepton==1) mc2->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc2->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc2->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());

	if (useEleMuo) {
	  if (ilepton==1) {
	    mc2 = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub/" + title + ".root").c_str());
	    h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	  }
	  if (ilepton==2) {
	    mc2 = TFile::Open((path + "/muons/" + version + "/" + subdir + "/ttbar_sub/" + title + ".root").c_str());
	    h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	  }
	}

	if (ilepton==1) mc3->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc3->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc3->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc4->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc4->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc4->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());

//	if (ilepton==1) mc5->cd(("demoEle"+postfix).c_str());
//	if (ilepton==2) mc5->cd(("demoMuo"+postfix).c_str());
//	if (ilepton==3) mc5->cd(("demoEleMuo"+postfix).c_str());
//	TH1F* h_mc5 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc6->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc6->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc6->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc7->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc7->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc7->cd(("demoEleMuo"+postfix).c_str());
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
	if (h_mc1t) {
	  h_mc1t -> Sumw2();
	  h_mc1t->SetLineColor(kBlack);
	  h_mc1t->SetFillColor(kOrange-4);
	  //h_mc1t->SetFillStyle(3004);
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

	if (irun==10) {
	  norm1 = norm1 + enorm1;
	  norm2 = norm2 + enorm2;
	  norm3 = norm3 + enorm3;
	  norm4 = norm4 + enorm4;
	  norm5 = norm5 + enorm5;
	  norm6 = norm6 + enorm6;
	  norm7 = norm7 + enorm7;
	}

	h_mc1->Scale(norm1);
	if (h_mc1b) h_mc1b->Scale(norm1);
	if (h_mc1c) h_mc1c->Scale(norm1);
	if (h_mc1t) h_mc1t->Scale(norm1);
	h_mc2->Scale(norm2);
	h_mc3->Scale(norm3);
	h_mc4->Scale(norm4);
//	h_mc5->Scale(norm5);
	h_mc6->Scale(norm6);
	h_mc7->Scale(norm7);

	if (useFitResults) {
	  h_mc2->Scale(1./norm2);
	  if (irun==10) norm2 = norm2 - enorm2;
	  if (title.find("_b")==string::npos) {
	    h_mc2->Scale(norm2*c1_t);
	    if (irun==5) h_mc2->Scale((c1_t+ec1_t)/c1_t);
	  } else {
	    h_mc2->Scale(norm2*c2_t);
	    if (irun==5) h_mc2->Scale((c2_t+ec2_t)/c2_t);
	  }
	}

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
	if (h_mc1t) h_mc1->Add(h_mc1t, -1.);
	for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
	  float e = h_mc1->GetBinError(i)**2;
	  if (h_mc1b) e = e - h_mc1b->GetBinError(i)**2;
	  if (h_mc1c) e = e - h_mc1c->GetBinError(i)**2;
	  if (h_mc1t) e = e - h_mc1t->GetBinError(i)**2;
	  h_mc1->SetBinError(i, TMath::Sqrt(e));
	}

	TVirtualFitter* fitter;
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
	    h_data_fit->SetBinError(i, TMath::Sqrt(e));
	  }
	  for (int i=0; i<=h_data_fit->GetNbinsX()+1; i++) {
	    bool skip = false;
	    if (title=="w_MET") {
	      if (h_data_fit->GetXaxis()->GetBinCenter(i) < 125) {
		skip = true;
	      }
	    }
	    if (title=="w_MET_sign") {
	      if (h_data_fit->GetXaxis()->GetBinCenter(i) < 50) {
		skip = true;
	      }
	    }
	    if (title=="w_MET_b") {
	      if (h_data_fit->GetXaxis()->GetBinCenter(i) < 90) {
		skip = true;
	      }
	    }
	    if (title=="w_MET_sign_b") {
	      if (h_data_fit->GetXaxis()->GetBinCenter(i) < 30) {
		skip = true;
	      }
	    }
	    if (skip) {
	      h_data->SetBinContent(i, 0);
	      h_data->SetBinError(i, 0);
	      h_data_fit->SetBinContent(i, 0);
	      h_data_fit->SetBinError(i, 0);
	      h_mc1->SetBinContent(i, 0);
	      h_mc1->SetBinError(i, 0);
	      if (h_mc1b) h_mc1b->SetBinContent(i, 0);
	      if (h_mc1b) h_mc1b->SetBinError(i, 0);
	      if (h_mc1c) h_mc1c->SetBinContent(i, 0);
	      if (h_mc1c) h_mc1c->SetBinError(i, 0);
	      if (h_mc1t) h_mc1t->SetBinContent(i, 0);
	      if (h_mc1t) h_mc1t->SetBinError(i, 0);
	      h_mc2->SetBinContent(i, 0);
	      h_mc2->SetBinError(i, 0);
	      h_mc3->SetBinContent(i, 0);
	      h_mc3->SetBinError(i, 0);
	      h_mc4->SetBinContent(i, 0);
	      h_mc4->SetBinError(i, 0);
//	      h_mc5->SetBinContent(i, 0);
//	      h_mc5->SetBinError(i, 0);
	      h_mc6->SetBinContent(i, 0);
	      h_mc6->SetBinError(i, 0);
	      h_mc7->SetBinContent(i, 0);
	      h_mc7->SetBinError(i, 0);
	    }
	  }
	  fitter = TVirtualFitter::Fitter(0, 1);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(t)", 1.0, 0.1, 0.0, 100.0);
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  h_mc_fit0->Scale(fitter->GetParameter(0));
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
	  fitter = TVirtualFitter::Fitter(0, 2);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(Z+jets)", 1.0, 0.1, 0.0, 100.0);
	  fitter->SetParameter(1, "c(t)", 1.0, 0.1, 0.0, 100.0);
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  if (h_mc1b) h_mc_fit0->Add(h_mc1b, -1.);
	  if (h_mc1c) h_mc_fit0->Add(h_mc1c, -1.);
	  h_mc_fit0->Scale(fitter->GetParameter(0));
	  if (h_mc1b) h_mc1b->Scale(fitter->GetParameter(0));
	  if (h_mc1c) h_mc1c->Scale(fitter->GetParameter(0));
	  h_mc_fit1->Scale(fitter->GetParameter(1));
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
	    bool skip = false;
	    if (title=="w_SVTX_mass") {
	      if (h_data_fit->GetXaxis()->GetBinCenter(i) < 0.2) {
	        skip = true;
	      }
	    }
	    if (skip) {
	      h_data->SetBinContent(i, 0);
	      h_data->SetBinError(i, 0);
	      h_data_fit->SetBinContent(i, 0);
	      h_data_fit->SetBinError(i, 0);
	      h_mc1->SetBinContent(i, 0);
	      h_mc1->SetBinError(i, 0);
	      if (h_mc1b) h_mc1b->SetBinContent(i, 0);
	      if (h_mc1b) h_mc1b->SetBinError(i, 0);
	      if (h_mc1c) h_mc1c->SetBinContent(i, 0);
	      if (h_mc1c) h_mc1c->SetBinError(i, 0);
	      if (h_mc1t) h_mc1t->SetBinContent(i, 0);
	      if (h_mc1t) h_mc1t->SetBinError(i, 0);
	      h_mc2->SetBinContent(i, 0);
	      h_mc2->SetBinError(i, 0);
	      h_mc3->SetBinContent(i, 0);
	      h_mc3->SetBinError(i, 0);
	      h_mc4->SetBinContent(i, 0);
	      h_mc4->SetBinError(i, 0);
//	      h_mc5->SetBinContent(i, 0);
//	      h_mc5->SetBinError(i, 0);
	      h_mc6->SetBinContent(i, 0);
	      h_mc6->SetBinError(i, 0);
	      h_mc7->SetBinContent(i, 0);
	      h_mc7->SetBinError(i, 0);
	    }
	  }
	  fitter = TVirtualFitter::Fitter(0, 3);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(uds)", 1.0, 0.1, 0.0, 100.0);
	  fitter->SetParameter(1, "c(b)", 1.0, 0.1, 0.0, 100.0);
	  fitter->SetParameter(2, "c(c)", 1.0, 0.1, 0.0, 100.0);
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  h_mc_fit0->Scale(fitter->GetParameter(0));
	  h_mc_fit1->Scale(fitter->GetParameter(1));
	  h_mc_fit2->Scale(fitter->GetParameter(2));
	}

	TH1F *ht = h_mc1->Clone("ht");
	ht->Reset();
	if (h_mc1t) ht->Add(h_mc1t);
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
	if (h_mc1t) hs->Add(h_mc1t);
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
	if (ilepton==3) leg->AddEntry(h_data,"Z(#rightarrow e#mu)+jets","p");

	if (h_mc1c && h_mc1b) {
	  leg->AddEntry(h_mc1,"Z+uds-jets","f");
	} else {
	  leg->AddEntry(h_mc1,"Z+jets","f");
	}
	if (h_mc1c) leg->AddEntry(h_mc1c,"Z+c-jets","f");
	if (h_mc1b) leg->AddEntry(h_mc1b,"Z+b-jets","f");
	if (h_mc1t) leg->AddEntry(h_mc1t,"#tau^{+}#tau^{-}+jets","f");
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
	} else if (title=="w_BJP"||title=="w_JBP") {
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
	    sprintf(buff, "c_{t} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	  }
	  if (doFit==2) {
	    sprintf(buff, "c_{Z+jets} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	    sprintf(buff, "c_{t} = %5.3f #pm %5.3f", fitter->GetParameter(1), fitter->GetParError(1));
	    fitLabel->DrawLatex(0.68, 0.43, buff);
	  }
	  if (doFit==3) {
	    sprintf(buff, "c_{uds} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    fitLabel->DrawLatex(0.38, 0.48, buff);
	    sprintf(buff, "c_{b}   = %5.3f #pm %5.3f", fitter->GetParameter(1), fitter->GetParError(1));
	    fitLabel->DrawLatex(0.38, 0.43, buff);
	    sprintf(buff, "c_{c}   = %5.3f #pm %5.3f", fitter->GetParameter(2), fitter->GetParError(2));
	    fitLabel->DrawLatex(0.38, 0.38, buff);
	    float f_uds = 100*h_mc_fit0->Integral()/(h_mc_fit0->Integral()+h_mc_fit1->Integral()+h_mc_fit2->Integral());
	    float ef_uds = f_uds*(fitter->GetParError(0)/fitter->GetParameter(0));
	    sprintf(buff, "f_{uds} = %4.1f #pm %3.1f %%", f_uds, ef_uds);
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	    float f_b = 100*h_mc_fit1->Integral()/(h_mc_fit0->Integral()+h_mc_fit1->Integral()+h_mc_fit2->Integral());
	    float ef_b = f_b*(fitter->GetParError(1)/fitter->GetParameter(1));
	    sprintf(buff, "f_{b}   = %4.1f #pm %3.1f %%", f_b, ef_b);
	    fitLabel->DrawLatex(0.68, 0.43, buff);
	    float f_c = 100*h_mc_fit2->Integral()/(h_mc_fit0->Integral()+h_mc_fit1->Integral()+h_mc_fit2->Integral());
	    float ef_c = f_c*(fitter->GetParError(2)/fitter->GetParameter(2));
	    sprintf(buff, "f_{c}   = %4.1f #pm %3.1f %%", f_c, ef_c);
	    fitLabel->DrawLatex(0.68, 0.38, buff);
	  }
	}

	if (plot) {
	  if (doBkg) title = title + "_doBkg";
	  if (doFit) title = title + "_doFit";
	  ofstream out;
	  if (ilepton==1) {
	    gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electrons/" + version + "/" + subdir + "/distributions/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + title + ".dat").c_str());
	  }
	  if (ilepton==2) {
	    gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/muons/" + version + "/" + subdir + "/distributions/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + title + ".dat").c_str());
	  }
	  if (ilepton==3) {
	    gSystem->mkdir((path + "/electrons+muons/" + version + "/" + subdir + "/distributions/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electrons+muons/" + version + "/" + subdir + "/distributions/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/electrons+muons/" + version + "/" + subdir + "/distributions/" + title + ".dat").c_str());
	  }
	  if (doFit==1) {
	    out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
	    out.close();
	  }
	  if (doFit==2) {
	    out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
	    out << fitter->GetParameter(1) << " " << fitter->GetParError(1) << endl;
	    out.close();
	  }
	  if (doFit==3) {
	    out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
	    out << fitter->GetParameter(1) << " " << fitter->GetParError(1) << endl;
	    out << fitter->GetParameter(2) << " " << fitter->GetParError(2) << endl;
	    out.close();
	  }
	}
}

