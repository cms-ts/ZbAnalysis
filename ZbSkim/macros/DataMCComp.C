#include "DataMCComp.h"
#include "LumiInfo_v14.h"

#include "CMS_lumi.C"
#include "CMS_process.C"

#include "fixrange.C"
#include "rebin.C"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";
//string path = "/gpfs/cms/users/lalicata/work/test/data/";

TH1F* h_data_fit = 0;
TH1F* h_mc_fit0 = 0;
TH1F* h_mc_fit1 = 0;
TH1F* h_mc_fit2 = 0;

void fcn(int& npar, double* gin, double& fun, double* par, int iflag) {
  double chisq = 0.0;
  if (iflag) {};
  if (gin) {};
  for (int i=1; i<=h_data_fit->GetNbinsX(); i++) {
    double xn = h_data_fit->GetBinContent(i);
    double xd = TMath::Power(h_data_fit->GetBinError(i),2);
    if (npar>0) {
      xn = xn - par[0]*h_mc_fit0->GetBinContent(i);
      xd = xd + TMath::Power(par[0]*h_mc_fit0->GetBinError(i),2);
    }
    if (npar>1) {
      xn = xn - par[1]*h_mc_fit1->GetBinContent(i);
      xd = xd + TMath::Power(par[1]*h_mc_fit1->GetBinError(i),2);
    }
    if (npar>2) {
      xn = xn - par[2]*h_mc_fit2->GetBinContent(i);
      xd = xd + TMath::Power(par[2]*h_mc_fit2->GetBinError(i),2);
    }
    if (xd!=0) chisq = chisq + (xn*xn)/xd;
  }
  fun = chisq;

}

void DataMCComp(int irun=0, string title="", int plot=0, int ilepton=1, int doBkg=0, int doFit=0, int numB=0) {

//int useFitResults=0; // use MC predictions for c_t
int useFitResults=1;  // use fit results for c_t

//int useEleMuo = 0; // use MC or fit results for c_t
int useEleMuo = 1; // use e-mu fit results for c_t

int useDY = 0; // use MadGraph DY
//int useDY = 1; // use Sherpa DY
//int useDY = 2; // use Powheg DY
//int useDY = 3; // use MG_aMC@NLO+P8

int useWeights = 0; // do not use weights for numB=0
//int useWeights = 1; // use weights for numB=0

string bSel=" ";
string subdir="0";
string postfix="";
string dirbSel="";

bool bbBkg = false;
bool bbSig = false;

if (irun==1) {             // irun==1 => JEC Up
  subdir="1";
  postfix="Up";
}
if (irun==2) {             // irun==2 => JEC Down
  subdir="2";
  postfix="Down";
}
if (irun==3) {             // irun==3 => PU Up
  subdir="3";
  postfix="Pup";
}
if (irun==4) {             // irun==4 => PU Down
  subdir="4";
  postfix="Pum";
}
if (irun==5) {             // irun==5 => top bkg
  subdir="5";
  postfix="";
}
if (irun==6) {             // irun==6 => b purity
  subdir="6";
  postfix="";
}
if (irun==7) {             // irun==7 => unfolding
  subdir="7";
  postfix="";
}
if (irun==8) {             // irun==8 => unfolding with Sherpa
  subdir="8";
  postfix="";
}
if (irun==9) {             // irun==9 => unfolding with Powheg
  subdir="9";
  postfix="";
}
if (irun==10) {            // irun==10 => bkg systematics
  subdir="10";
  postfix="";
}
if (irun==11) {            // irun==11 => JER Up
  subdir="11";
  postfix="JerUp";
}
if (irun==12) {            // irun==12 => JER Down
  subdir="12";
  postfix="JerDown";
}
if (irun==13) {            // irun==13 => bkg statistics
  subdir="13";
  postfix="";
}
if (irun==14) {            // irun==14 => unfolding with MadGraph 4FS
  subdir="14";
  postfix="";
}
if (irun==15) {            // irun==15 => unfolding with data weight
  subdir="15";
  postfix="";
}
if (irun==16) {            // irun==16 => unfolding with MadGraph aMC@NLO
  subdir="16";
  postfix="";
}
if (irun==17) {            // irun==17 => templates from data
  subdir="17";
  postfix="";
}
if (irun==18) {            // irun==18 => templates from Sherpa
  subdir="18";
  postfix="";
}
if (irun==19) {            // irun==19 => templates from MadGraph aMC@NLO
  subdir="19";
  postfix="";
}
if (numB==1) {
  postfix = postfix + "1b";
  dirbSel = "_1b";
#ifdef PAPER
  bSel = "Z + (= 1) b jet";
#else
  bSel = "Z + (= 1) b-jet";
#endif
}
if (numB==2) {
  postfix = postfix + "2b";
  dirbSel = "_2b";
#ifdef PAPER
  bSel = "Z + (#geq 2) b jet";
#else
  bSel = "Z + (#geq 2) b-jet";
#endif
}

#ifdef PAPER
bSel = "";
#endif

if (numB==1) bbBkg = true;

	if (irun==19) useDY = 3;

        if (doFit==4) bbSig = true;

	/* top background */

	double c1_t=1.0;
	double ec1_t=0.0;
	double c2_t=1.0;
	double ec2_t=0.0;

	if (doFit==1 || doFit==2) {
	  if (title=="w_MET") { useFitResults=0; useEleMuo=0; }
	  if (title=="w_MET_b") { useFitResults=0; useEleMuo=0; }
	  if (title=="w_JBP") { useFitResults=0; useEleMuo=0; }
	  if (title=="w_BJP") { useFitResults=0; useEleMuo=0; }
	  if (title=="w_SVTX_mass") { useFitResults=0; useEleMuo=0; }
	}

	ifstream in4, in5, in6, in7;
	if (ilepton==1) {
	  if (useFitResults) {
	    in4.open((path + "/electrons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_doFit" + ".dat").c_str());
	    in5.open((path + "/electrons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_b_doFit" + ".dat").c_str());
	    if (useEleMuo) {
	      in6.open((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + "w_mass_ee_wide_doFit" + ".dat").c_str());
	      in7.open((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + "w_mass_ee_b_wide_doFit" + ".dat").c_str());
	    }
	  }
	}
	if (ilepton==2) {
	  if (useFitResults) {
	    in4.open((path + "/muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_doFit" + ".dat").c_str());
	    in5.open((path + "/muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_b_doFit" + ".dat").c_str());
	    if (useEleMuo) {
	      in6.open((path + "/muons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + "w_mass_mm_wide_doFit" + ".dat").c_str());
	      in7.open((path + "/muons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + "w_mass_mm_b_wide_doFit" + ".dat").c_str());
	    }
	  }
	}
	if (ilepton==3) {
	  if (useFitResults) {
	    in4.open((path + "/electrons+muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_doFit" + ".dat").c_str());
	    in5.open((path + "/electrons+muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_b_doFit" + ".dat").c_str());
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

        double fScal=1.0;
        double efScal=0.0;

        if (bbBkg) {
          ifstream in8;
          if (ilepton==1) {
            in8.open((path + "/electrons/" + version + "/" + subdir + "/distributions_2b" + "/" + "w_SVTX_mass_doFit" + ".dat").c_str());
          }
          if (ilepton==2) {
            in8.open((path + "/muons/" + version + "/" + subdir + "/distributions_2b" + "/" + "w_SVTX_mass_doFit" + ".dat").c_str());
          }
          in8 >> fScal >> efScal;
        }

	double Lumi2012=0;

	if (ilepton==1) Lumi2012 = Lumi2012_ele;
	if (ilepton==2) Lumi2012 = Lumi2012_muon;
	if (ilepton==3) Lumi2012 = Lumi2012_ele_muon;

	double norm1 = ((Lumi2012 * Xsec_dy) / Ngen_dy);
	double norm1_amc = ((Lumi2012 * Xsec_dy_amc) / Ngen_dy_amc);
	double norm1_1 = ((Lumi2012 * Xsec_dy_1) / Ngen_dy_1);
	double norm1_2=0;
	if (ilepton==1) norm1_2 = ((Lumi2012 * Xsec_dy_2) / Ngen_dy_2_ee);
	if (ilepton==2) norm1_2 = ((Lumi2012 * Xsec_dy_2) / Ngen_dy_2_mm);
	double norm2 = ((Lumi2012 * Xsec_tt) / Ngen_tt);
	if (useEleMuo && ilepton!=3) norm2 = (Lumi2012 / Lumi2012_ele_muon);
	double norm3 = ((Lumi2012 * Xsec_zz) / Ngen_zz);
	double norm4 = ((Lumi2012 * Xsec_wz) / Ngen_wz);
	double norm5 = ((Lumi2012 * Xsec_qcd) / Ngen_qcd);
	double norm6 = ((Lumi2012 * Xsec_ww) / Ngen_ww);
	double norm7 = ((Lumi2012 * Xsec_wj) / Ngen_wj);
	double norm8 = ((Lumi2012 * Xsec_tS) / Ngen_tS);
	double norm9 = ((Lumi2012 * Xsec_tT) / Ngen_tT);
	double norm10 = ((Lumi2012 * Xsec_tW) / Ngen_tW);
	double norm11 = ((Lumi2012 * Xsec_tSb) / Ngen_tSb);
	double norm12 = ((Lumi2012 * Xsec_tTb) / Ngen_tTb);
	double norm13 = ((Lumi2012 * Xsec_tWb) / Ngen_tWb);

	double enorm1 = ((Lumi2012 * eXsec_dy) / Ngen_dy);
	double enorm1_amc = ((Lumi2012 * eXsec_dy_amc) / Ngen_dy_amc);
	double enorm1_1 = ((Lumi2012 * eXsec_dy_1) / Ngen_dy_1);
	double enorm1_2=0;
	if (ilepton==1) enorm1_2 = ((Lumi2012 * eXsec_dy_2) / Ngen_dy_2_ee);
	if (ilepton==2) enorm1_2 = ((Lumi2012 * eXsec_dy_2) / Ngen_dy_2_mm);
	double enorm2 = ((Lumi2012 * eXsec_tt) / Ngen_tt);
	if (useEleMuo && ilepton!=3) enorm2 = 0;
	double enorm3 = ((Lumi2012 * eXsec_zz) / Ngen_zz);
	double enorm4 = ((Lumi2012 * eXsec_wz) / Ngen_wz);
	double enorm5 = ((Lumi2012 * eXsec_qcd) / Ngen_qcd);
	double enorm6 = ((Lumi2012 * eXsec_ww) / Ngen_ww);
	double enorm7 = ((Lumi2012 * eXsec_wj) / Ngen_wj);
	double enorm8 = ((Lumi2012 * eXsec_tS) / Ngen_tS);
	double enorm9 = ((Lumi2012 * eXsec_tT) / Ngen_tT);
	double enorm10 = ((Lumi2012 * eXsec_tW) / Ngen_tW);
	double enorm11 = ((Lumi2012 * eXsec_tSb) / Ngen_tSb);
	double enorm12 = ((Lumi2012 * eXsec_tTb) / Ngen_tTb);
	double enorm13 = ((Lumi2012 * eXsec_tWb) / Ngen_tWb);

	if (title.empty()) title = "w_jetmultiplicity";

	if (ilepton==1) {
	  if (title.find("muon")!=string::npos) return;
	  if (title.find("mm")!=string::npos) return;
	}
	if (ilepton==2) {
	  if (title.find("ele")!=string::npos) return;
	  if (title.find("ee")!=string::npos) return;
	}

	if (useDY!=0 && ilepton==3) return;

	TFile* data=0;
	if (ilepton==1) data = TFile::Open((path + "/" + version + "/" + "DoubleElectron_2012_merge.root").c_str());
	if (ilepton==2) data = TFile::Open((path + "/" + version + "/" + "DoubleMu_2012_merge.root").c_str());
	if (ilepton==3) data = TFile::Open((path + "/" + version + "/" + "MuEG_2012_merge.root").c_str());

	TFile* mc1=0;
	if (useDY==0) {
	  mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL.root").c_str());
	  if (useWeights && numB==0) {
	    mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_weights.root").c_str());
	  }
	}
	if (useDY==1) {
	  norm1 = norm1_1;
	  enorm1 = enorm1_1;
	  mc1 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa.root").c_str());
	}
	if (useDY==2) {
	  norm1 = norm1_2;
	  enorm1 = enorm1_2;
	  if (ilepton==1) mc1 = TFile::Open((path + "/" + version + "/" + "DYToEE_powheg_gen.root").c_str());
	  if (ilepton==2) mc1 = TFile::Open((path + "/" + version + "/" + "DYToMuMu_powheg_gen.root").c_str());
	}
	if (useDY==3) {
	  norm1 = norm1_amc;
          enorm1 = enorm1_amc;
	  mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_aMC.root").c_str());
	  if (useWeights && numB==0) {
	    mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_aMC_weights.root").c_str());
	  }
	}

	TFile* mc2 = TFile::Open((path + "/" + version + "/" + "TTbar.root").c_str());
	TFile* mc3 = TFile::Open((path + "/" + version + "/" + "ZZ.root").c_str());
	TFile* mc4 = TFile::Open((path + "/" + version + "/" + "WZ.root").c_str());
//	TFile* mc5 = TFile::Open((path + "/" + version + "/" + "QCD.root").c_str());
	TFile* mc6 = TFile::Open((path + "/" + version + "/" + "WW.root").c_str());
	TFile* mc7 = TFile::Open((path + "/" + version + "/" + "Wj.root").c_str());
	TFile* mc8 = TFile::Open((path + "/" + version + "/" + "T_s.root").c_str());
	TFile* mc9 = TFile::Open((path + "/" + version + "/" + "T_t.root").c_str());
	TFile* mc10 = TFile::Open((path + "/" + version + "/" + "T_tW.root").c_str());
	TFile* mc11 = TFile::Open((path + "/" + version + "/" + "TBar_s.root").c_str());
	TFile* mc12 = TFile::Open((path + "/" + version + "/" + "TBar_t.root").c_str());
	TFile* mc13 = TFile::Open((path + "/" + version + "/" + "TBar_tW.root").c_str());

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
        TH1F* h_mc1bb = 0;
        if (bbBkg) h_mc1bb = (TH1F*)gDirectory->Get(("bbBkg"+title.substr(1)).c_str());
        if (bbSig) h_mc1bb = (TH1F*)gDirectory->Get(("bbSig"+title.substr(1)).c_str());

	if (!h_mc1bb) bbBkg = 0;

	if (ilepton==1) mc2->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc2->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc2->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());

	if (useEleMuo) {
	  if (ilepton==1) {
	    mc2 = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + title + ".root").c_str());
	    h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	  }
	  if (ilepton==2) {
	    mc2 = TFile::Open((path + "/muons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + title + ".root").c_str());
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

	if (ilepton==1) mc8->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc8->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc8->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_mc8 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc9->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc9->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc9->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_mc9 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc10->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc10->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc10->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_mc10 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc11->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc11->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc11->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_mc11 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc12->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc12->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc12->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_mc12 = (TH1F*)gDirectory->Get(title.c_str());

	if (ilepton==1) mc13->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc13->cd(("demoMuo"+postfix).c_str());
	if (ilepton==3) mc13->cd(("demoEleMuo"+postfix).c_str());
	TH1F* h_mc13 = (TH1F*)gDirectory->Get(title.c_str());

	h_data = fixrange(h_data, numB);
	h_mc1 = fixrange(h_mc1, numB);
	if (h_mc1b) h_mc1b = fixrange(h_mc1b, numB);
	if (h_mc1c) h_mc1c = fixrange(h_mc1c, numB);
	if (h_mc1t) h_mc1t = fixrange(h_mc1t, numB);
	if (bbBkg || bbSig) h_mc1bb = fixrange(h_mc1bb, numB);
	h_mc2 = fixrange(h_mc2, numB);
	h_mc3 = fixrange(h_mc3, numB);
	h_mc4 = fixrange(h_mc4, numB);
//	h_mc5 = fixrange(h_mc5, numB);
	h_mc6 = fixrange(h_mc6, numB);
	h_mc7 = fixrange(h_mc7, numB);
	h_mc8 = fixrange(h_mc8, numB);
	h_mc9 = fixrange(h_mc9, numB);
	h_mc10 = fixrange(h_mc10, numB);
	h_mc11 = fixrange(h_mc11, numB);
	h_mc12 = fixrange(h_mc12, numB);
	h_mc13 = fixrange(h_mc13, numB);

	h_data = rebin(h_data, numB);
	h_mc1 = rebin(h_mc1, numB);
	if (h_mc1b) h_mc1b = rebin(h_mc1b, numB);
	if (h_mc1c) h_mc1c = rebin(h_mc1c, numB);
	if (h_mc1t) h_mc1t = rebin(h_mc1t, numB);
	if (bbBkg || bbSig) h_mc1bb = rebin(h_mc1bb, numB);
	h_mc2 = rebin(h_mc2, numB);
	h_mc3 = rebin(h_mc3, numB);
	h_mc4 = rebin(h_mc4, numB);
//	h_mc5 = rebin(h_mc5, numB);
	h_mc6 = rebin(h_mc6, numB);
	h_mc7 = rebin(h_mc7, numB);
	h_mc8 = rebin(h_mc8, numB);
	h_mc9 = rebin(h_mc9, numB);
	h_mc10 = rebin(h_mc10, numB);
	h_mc11 = rebin(h_mc11, numB);
	h_mc12 = rebin(h_mc12, numB);
	h_mc13 = rebin(h_mc13, numB);

	h_mc1 -> SetLineColor(kBlack);
	h_mc1 -> SetFillColor(kYellow-4);
	//h_mc1 -> SetFillStyle(3004);

	if (h_mc1b) {
	  h_mc1b->SetLineColor(kBlack);
	  h_mc1b->SetFillColor(kYellow-4);
	  h_mc1b->SetFillStyle(3254);
	}
	if (h_mc1c) {
	  h_mc1c->SetLineColor(kBlack);
	  h_mc1c->SetFillColor(kOrange);
	  h_mc1c->SetFillStyle(3245);
	}
	if (h_mc1t) {
	  h_mc1t->SetLineColor(kBlack);
	  h_mc1t->SetFillColor(kOrange-4);
	  //h_mc1t->SetFillStyle(3004);
	}
        if (bbBkg || bbSig) {
          h_mc1bb->SetLineColor(kBlack);
          h_mc1bb->SetFillColor(kYellow+2);
          h_mc1bb->SetFillStyle(3254);
	}

	h_mc2 -> SetLineColor(kBlack);
	h_mc2 -> SetFillColor(kBlue);
	//h_mc2 -> SetFillStyle(3004);

	h_mc3 -> SetLineColor(kBlack);
	h_mc3 -> SetFillColor(kGray+2);
	//h_mc3 -> SetFillStyle(3004);

	h_mc4 -> SetLineColor(kBlack);
	h_mc4 -> SetFillColor(kGray+3);
	//h_mc4 -> SetFillStyle(3004);

//	h_mc5 -> SetLineColor(kBlack);
//	h_mc5 -> SetFillColor(kGray+3);
//	//h_mc5 -> SetFillStyle(3004);

	h_mc6 -> SetLineColor(kBlack);
	h_mc6 -> SetFillColor(kRed+2);
	//h_mc6 -> SetFillStyle(3004);

	h_mc7 -> SetLineColor(kBlack);
	h_mc7 -> SetFillColor(kGray);
	//h_mc7 -> SetFillStyle(3004);

	if (irun==10) {
	  norm1 = norm1 + 0.1*enorm1;
	  norm2 = norm2 + 0.1*enorm2;
	  norm3 = norm3 + 0.1*enorm3;
	  norm4 = norm4 + 0.1*enorm4;
	  norm5 = norm5 + 0.1*enorm5;
	  norm6 = norm6 + 0.1*enorm6;
	  norm7 = norm7 + 0.1*enorm7;
	  norm8 = norm8 + 0.1*enorm8;
	  norm9 = norm9 + 0.1*enorm9;
	  norm10 = norm10 + 0.1*enorm10;
	  norm11 = norm11 + 0.1*enorm11;
	  norm12 = norm12 + 0.1*enorm12;
	  norm13 = norm13 + 0.1*enorm13;
	}

	h_mc1->Scale(norm1);
	if (h_mc1b) h_mc1b->Scale(norm1);
	if (h_mc1c) h_mc1c->Scale(norm1);
	if (h_mc1t) h_mc1t->Scale(norm1);
	if (bbBkg || bbSig) h_mc1bb->Scale(norm1);

	h_mc2->Scale(norm2);
	h_mc3->Scale(norm3);
	h_mc4->Scale(norm4);
//	h_mc5->Scale(norm5);
	h_mc6->Scale(norm6);
	h_mc7->Scale(norm7);
	h_mc8->Scale(norm8);
	h_mc9->Scale(norm9);
	h_mc10->Scale(norm10);
	h_mc11->Scale(norm11);
	h_mc12->Scale(norm12);
	h_mc13->Scale(norm13);

        TH1F* h_mcO = (TH1F*)h_mc8->Clone("h_mcO");
	h_mcO->Reset();

	TH1F* h_mcD = (TH1F*)h_mc3->Clone("h_mcD");
        h_mcD->Reset();

        h_mcO->Add(h_mc13);
        h_mcO->Add(h_mc12);
        h_mcO->Add(h_mc11);
        h_mcO->Add(h_mc10);
        h_mcO->Add(h_mc9);
        h_mcO->Add(h_mc8);
        h_mcO->Add(h_mc7);
        if (h_mc1t) h_mcO->Add(h_mc1t);

        h_mcD->Add(h_mc6);
        h_mcD->Add(h_mc4);
        h_mcD->Add(h_mc3);

        h_mcO -> SetLineColor(kBlack);
        h_mcO -> SetFillColor(kOrange+1);
        //h_mcO -> SetFillStyle(3004);

	h_mcD -> SetLineColor(kBlack);
        h_mcD -> SetFillColor(kBlue-2);
#ifdef PAPER
	h_mcD -> SetFillColor(kGreen+2);
#endif
        //h_mcD -> SetFillStyle(3004);

        if (bbBkg || bbSig) {
	  for (int i=0; i<=h_mc1b->GetNbinsX()+1; i++) {
	    double c = h_mc1b->GetBinContent(i);
	    c = c - h_mc1bb->GetBinContent(i);
	    if (c<0) c = 0.0;
	    h_mc1b->SetBinContent(i,c);
	    double e = TMath::Power(h_mc1b->GetBinError(i),2);
	    e = e - TMath::Power(h_mc1bb->GetBinError(i),2);
	    if (e<0) e = 0.0;
	    h_mc1b->SetBinError(i,TMath::Sqrt(e));
	  }
	}

	if (useFitResults) {
	  int use_c1_t=1;
	  if (title.find("JBP")==string::npos) use_c1_t=0;
	  if (title.find("BJP")==string::npos) use_c1_t=0;
	  if (title.find("SVTX")==string::npos) use_c1_t=0;
	  if (title.find("secondvtx")==string::npos) use_c1_t=0;
	  if (title.find("dxy")==string::npos) use_c1_t=0;
	  if (title.find("flightd")==string::npos) use_c1_t=0;
	  if (title.find("SV")==string::npos) use_c1_t=0;
	  if (title.find("_b")!=string::npos) use_c1_t=0;
	  if (title.find("_2b")!=string::npos) use_c1_t=0;
	  if (title.find("_bb")!=string::npos) use_c1_t=0;
	  if (title.find("eeb_min")==string::npos) use_c1_t=0;
	  if (title.find("mmb_min")==string::npos) use_c1_t=0;
	  if (title.find("emb_min")==string::npos) use_c1_t=0;
	  if (title.find("eeb_max")==string::npos) use_c1_t=0;
	  if (title.find("mmb_max")==string::npos) use_c1_t=0;
	  if (title.find("emb_max")==string::npos) use_c1_t=0;
	  if (title.find("eebb_mass")==string::npos) use_c1_t=0;
	  if (title.find("mmbb_mass")==string::npos) use_c1_t=0;
	  if (title.find("embb_mass")==string::npos) use_c1_t=0;
	  if (title.find("_eeb")==string::npos) use_c1_t=0;
	  if (title.find("_mmb")==string::npos) use_c1_t=0;
	  if (title.find("_emb")==string::npos) use_c1_t=0;
	  if (title.find("_delta_j")==string::npos) use_c1_t=0;
	  if (use_c1_t) {
	    if (irun==5) {
	      h_mc2->Scale(c1_t+0.1*ec1_t);
	    } else {
	      h_mc2->Scale(c1_t);
	    }
	  } else {
	    if (irun==5) {
	      h_mc2->Scale(c2_t+0.1*ec2_t);
	    } else {
	      h_mc2->Scale(c2_t);
	    }
	  }
	}

	if (irun==13) {
	  for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
	    h_mc1->SetBinError(i, 1.1*h_mc1->GetBinError(i));
	    if (h_mc1b) h_mc1b->SetBinError(i, 1.1*h_mc1b->GetBinError(i));
	    if (h_mc1c) h_mc1c->SetBinError(i, 1.1*h_mc1c->GetBinError(i));
            if (bbBkg || bbSig) h_mc1bb->SetBinError(i, 1.1*h_mc1bb->GetBinError(i));
	    h_mc2->SetBinError(i, 1.1*h_mc2->GetBinError(i));
	    h_mcD->SetBinError(i, 1.1*h_mcD->GetBinError(i));
	    h_mcO->SetBinError(i, 1.1*h_mcO->GetBinError(i));
	  }
	}

	for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
	  double c = h_mc1->GetBinContent(i);
	  if (h_mc1b) c = c - h_mc1b->GetBinContent(i);
	  if (h_mc1c) c = c - h_mc1c->GetBinContent(i);
	  if (h_mc1t) c = c - h_mc1t->GetBinContent(i);
	  if (bbBkg || bbSig) c = c - h_mc1bb->GetBinContent(i);
	  if (c<0) c = 0.0;
	  h_mc1->SetBinContent(i,c);
	  double e = TMath::Power(h_mc1->GetBinError(i),2);
	  if (h_mc1b) e = e - TMath::Power(h_mc1b->GetBinError(i),2);
	  if (h_mc1c) e = e - TMath::Power(h_mc1c->GetBinError(i),2);
	  if (h_mc1t) e = e - TMath::Power(h_mc1t->GetBinError(i),2);
          if (bbBkg || bbSig) e = e - TMath::Power(h_mc1bb->GetBinError(i),2);
	  if (e<0) e = 0.0;
	  h_mc1->SetBinError(i,TMath::Sqrt(e));
	}

        if (bbBkg) h_mc1bb->Scale(fScal);

	if (doBkg) {
	  h_data->Add(h_mcO, -1.);
	  h_data->Add(h_mcD, -1.);
	  h_data->Add(h_mc2, -1.);
          if (bbBkg) h_data->Add(h_mc1bb, -1.);
	}

	if (doBkg) {
	  if (bbSig) {
            h_data->Add(h_mc1, -1.);
	    if (h_mc1b) h_data->Add(h_mc1b, -1.);
	    if (h_mc1c) h_data->Add(h_mc1c, -1.);
	    if (h_mc1t) h_data->Add(h_mc1t, -1.);
	  }
	}

	if (doFit==0 && irun==18) {
	  double xval=0;
	  //double xvalb=0;
	  double xvalc=0;

	  if (h_mc1)  xval = h_mc1->Integral(0,h_mc1->GetNbinsX()+1);
	  //if (h_mc1b) xvalb = h_mc1b->Integral(0,h_mc1b->GetNbinsX()+1);
          if (h_mc1c) xvalc = h_mc1c->Integral(0,h_mc1c->GetNbinsX()+1);

	  mc1 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa.root").c_str());
          if (ilepton==1) mc1->cd(("demoEle"+postfix).c_str());
          if (ilepton==2) mc1->cd(("demoMuo"+postfix).c_str());
          if (ilepton==3) mc1->cd(("demoEleMuo"+postfix).c_str());

	  h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	  h_mc1 = fixrange(h_mc1, numB);

	  TH1F* h_mc1b_tmp = h_mc1b;

	  h_mc1b = (TH1F*)gDirectory->Get(("b"+title.substr(1)).c_str());
	  if (h_mc1b) h_mc1b = fixrange(h_mc1b, numB);

	  h_mc1c = (TH1F*)gDirectory->Get(("c"+title.substr(1)).c_str());
	  if (h_mc1c) h_mc1c = fixrange(h_mc1c, numB);

          if (bbBkg) {
            h_mc1bb = (TH1F*)gDirectory->Get(("bbBkg"+title.substr(1)).c_str());
	    h_mc1bb = fixrange(h_mc1bb, numB);
            if (h_mc1b) {
	      for (int i=0; i<=h_mc1b->GetNbinsX()+1; i++) {
	        double c = h_mc1b->GetBinContent(i);
	        c = c - h_mc1bb->GetBinContent(i);
	        if (c<0) c = 0.0;
	        h_mc1b->SetBinContent(i,c);
	        double e = TMath::Power(h_mc1b->GetBinError(i),2);
	        e = e - TMath::Power(h_mc1bb->GetBinError(i),2);
	        if (e<0) e = 0.0;
	        h_mc1b->SetBinError(i,TMath::Sqrt(e));
	      }
	    }
          }

          if (bbSig) {
            h_mc1bb = (TH1F*)gDirectory->Get(("bbSig"+title.substr(1)).c_str());
	    h_mc1bb = fixrange(h_mc1bb, numB);
          }

	  for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
	    double c = h_mc1->GetBinContent(i);
	    if (h_mc1b) c = c - h_mc1b->GetBinContent(i);
	    if (h_mc1c) c = c - h_mc1c->GetBinContent(i);
	    if (h_mc1t) c = c - h_mc1t->GetBinContent(i);
	    if (bbBkg || bbSig) c = c - h_mc1bb->GetBinContent(i);
	    if (c<0) c = 0.0;
	    h_mc1->SetBinContent(i,c);
	    double e = TMath::Power(h_mc1->GetBinError(i),2);
	    if (h_mc1b) e = e - TMath::Power(h_mc1b->GetBinError(i),2);
	    if (h_mc1c) e = e - TMath::Power(h_mc1c->GetBinError(i),2);
	    if (h_mc1t) e = e - TMath::Power(h_mc1t->GetBinError(i),2);
            if (bbBkg || bbSig) e = e - TMath::Power(h_mc1bb->GetBinError(i),2);
	    if (e<0) e = 0.0;
	    h_mc1->SetBinError(i,TMath::Sqrt(e));
	  }

	  if (h_mc1->Integral(0,h_mc1->GetNbinsX()+1)!=0) {
	    xval = xval / h_mc1->Integral(0,h_mc1->GetNbinsX()+1);
	    h_mc1->Scale(xval);
	  }

	  //if (h_mc1b) {
	  //  xvalb = xvalb / h_mc1b->Integral(0,h_mc1b->GetNbinsX()+1);
	  //  h_mc1b->Scale(xvalb);
	  //}

	  h_mc1b = h_mc1b_tmp;

	  if (h_mc1c) {
	    if (h_mc1c->Integral(0,h_mc1c->GetNbinsX()+1)!=0) {
	      xvalc = xvalc / h_mc1c->Integral(0,h_mc1c->GetNbinsX()+1);
	      h_mc1c->Scale(xvalc);
	    }
	  }

	  h_mc1 -> SetLineColor(kBlack);
	  h_mc1 -> SetFillColor(kYellow-4);
	  //h_mc1 -> SetFillStyle(3004);

	  if (h_mc1b) {
	    h_mc1b->SetLineColor(kBlack);
	    h_mc1b->SetFillColor(kYellow-4);
	    h_mc1b->SetFillStyle(3254);
	  }
	  if (h_mc1c) {
	    h_mc1c->SetLineColor(kBlack);
	    h_mc1c->SetFillColor(kOrange);
	    h_mc1c->SetFillStyle(3245);
	  }
	  if (h_mc1bb) {
	    h_mc1bb->SetLineColor(kBlack);
	    h_mc1bb->SetFillColor(kYellow+2);
	    h_mc1bb->SetFillStyle(3254);
	  }
	}

	if (doFit==0 && irun==17) {

	  double xval=0;
	  double xvalc=0;

          string tmp = title;

	  if (title.find("_bjet_")!=string::npos) {
	    tmp.erase(tmp.find("_bjet_")+1, 1);
	  } else if (title.find("_b\0")!=string::npos) {
	    tmp.erase(tmp.find("_b\0"), 2);
	  }

	  if (h_mc1)  xval = h_mc1->Integral(0,h_mc1->GetNbinsX()+1);
          if (h_mc1c) xvalc = h_mc1c->Integral(0,h_mc1c->GetNbinsX()+1);

	  TFile* data1=0;
          if (ilepton==1) data1 = TFile::Open((path + "/" + version + "/" + "DoubleElectron_2012_merge.root").c_str());
          if (ilepton==2) data1 = TFile::Open((path + "/" + version + "/" + "DoubleMu_2012_merge.root").c_str());
          if (ilepton==3) data1 = TFile::Open((path + "/" + version + "/" + "MuEG_2012_merge.root").c_str());

	  if (ilepton==1) data1->cd(("demoEleAntib"+postfix).c_str());
          if (ilepton==2) data1->cd(("demoMuoAntib"+postfix).c_str());
          if (ilepton==3) data1->cd(("demoEleMuoAntib"+postfix).c_str());

	  h_mc1 = (TH1F*)gDirectory->Get(tmp.c_str());
          h_mc1->SetLineColor(kBlack);
          h_mc1->SetFillColor(kYellow-4);

	  xval = xval / h_mc1->Integral(0,h_mc1->GetNbinsX()+1);
	  h_mc1->Scale(xval);

	  xvalc = xvalc / h_mc1c->Integral(0,h_mc1c->GetNbinsX()+1);
	  h_mc1c->Scale(xvalc);
   	}

	TVirtualFitter::SetDefaultFitter("Minuit2");
	TVirtualFitter* fitter=0;
	if (doFit==1) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    h_data_fit->Add(h_mcO, -1.);
	    h_data_fit->Add(h_mcD, -1.);
            if (bbBkg) h_data_fit->Add(h_mc1bb, -1.);
	  }
	  h_data_fit->Add(h_mc1, -1.);
	  if (h_mc1b) h_data_fit->Add(h_mc1b, -1.);
	  if (h_mc1c) h_data_fit->Add(h_mc1c, -1.);
	  h_mc_fit0 = h_mc2;
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
	      if (bbBkg || bbSig) h_mc1bb->SetBinContent(i, 0);
              if (bbBkg || bbSig) h_mc1bb->SetBinError(i, 0);
	      if (h_mc1c) h_mc1c->SetBinContent(i, 0);
	      if (h_mc1c) h_mc1c->SetBinError(i, 0);
	      if (h_mc1t) h_mc1t->SetBinContent(i, 0);
	      if (h_mc1t) h_mc1t->SetBinError(i, 0);
	      h_mc2->SetBinContent(i, 0);
	      h_mc2->SetBinError(i, 0);
	      h_mcD->SetBinContent(i, 0);
	      h_mcD->SetBinError(i, 0);
	      h_mcO->SetBinContent(i, 0);
	      h_mcO->SetBinError(i, 0);
	    }
	  }
	  fitter = TVirtualFitter::Fitter(0, 1);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(t)", 1.00, 0.01, 0.00, 100.00);
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  h_mc_fit0->Scale(fitter->GetParameter(0));
	}

	if (doFit==2) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    h_data_fit->Add(h_mcD, -1.);
	    h_data_fit->Add(h_mcO, -1.);
	    if (bbBkg) h_data_fit->Add(h_mc1bb, -1.);
	  }
	  h_mc_fit0 = h_mc1;
	  if (h_mc1b) h_mc_fit0->Add(h_mc1b, 1.);
	  if (h_mc1c) h_mc_fit0->Add(h_mc1c, 1.);
	  h_mc_fit1 = h_mc2;
	  fitter = TVirtualFitter::Fitter(0, 2);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(Z+jets)", 1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(1, "c(t)", 1.00, 0.01, 0.00, 100.00);
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  if (h_mc1b) h_mc_fit0->Add(h_mc1b, -1.);
	  if (h_mc1c) h_mc_fit0->Add(h_mc1c, -1.);
	  h_mc_fit0->Scale(fitter->GetParameter(0));
	  if (h_mc1b) h_mc1b->Scale(fitter->GetParameter(0));
	  if (h_mc1c) h_mc1c->Scale(fitter->GetParameter(0));
	  h_mc_fit1->Scale(fitter->GetParameter(1));
	}

	/* fit the top and "others" with doFit==2 instead of top+Zjets */
	/*
	if (doFit==2) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    //h_data_fit->Add(h_mcO, -1.);
	    h_data_fit->Add(h_mc7, -1.);
	    h_data_fit->Add(h_mc6, -1.);
//	    h_data_fit->Add(h_mc5, -1.);
	    h_data_fit->Add(h_mc4, -1.);
	    h_data_fit->Add(h_mc3, -1.);
	    h_data_fit->Add(h_mc1t, -1.);
	    if (bbBkg) h_data_fit->Add(h_mc1bb, -1.);
	  }
	  h_mc_fit0 = h_mcO;
	  h_mc_fit1 = h_mc2;
	  fitter = TVirtualFitter::Fitter(0, 2);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(Others)", 1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(1, "c(t)", 1.00, 0.01, 0.00, 100.00);
	  fitter->ExecuteCommand("MIGRAD", arglist, 0);
	  h_mc_fit0->Scale(fitter->GetParameter(0));
	  h_mc_fit1->Scale(fitter->GetParameter(1));
	}
	*/

	if (doFit==3) {
	  h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
	  if (!doBkg) {
	    h_data_fit->Add(h_mcD, -1.);
	    h_data_fit->Add(h_mcO, -1.);
	    h_data_fit->Add(h_mc2, -1.);
	    if (bbBkg) h_data_fit->Add(h_mc1bb, -1.);
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
	      if (bbBkg || bbSig) h_mc1bb->SetBinContent(i, 0);
	      if (bbBkg || bbSig) h_mc1bb->SetBinError(i, 0);
	      h_mc2->SetBinContent(i, 0);
	      h_mc2->SetBinError(i, 0);
	      h_mcD->SetBinContent(i, 0);
	      h_mcD->SetBinError(i, 0);
	      h_mcO->SetBinContent(i, 0);
	      h_mcO->SetBinError(i, 0);
	    }
	  }
	  fitter = TVirtualFitter::Fitter(0, 3);
	  fitter->SetFCN(fcn);
	  double arglist[1] = {-1.0};
	  //fitter->ExecuteCommand("SET PRINT", arglist, 1);
	  fitter->SetParameter(0, "c(uds)", 1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(1, "c(b)",   1.00, 0.01, 0.00, 100.00);
	  fitter->SetParameter(2, "c(c)",   1.00, 0.01, 0.00, 100.00);
	  /*
	  if (irun==99) {
	    fitter->FixParameter(0);
	    fitter->FixParameter(2);
	  }
	  */
          fitter->ExecuteCommand("MIGRAD", arglist, 0);
          h_mc_fit0->Scale(fitter->GetParameter(0));
	  h_mc_fit1->Scale(fitter->GetParameter(1));
	  h_mc_fit2->Scale(fitter->GetParameter(2));
	}

        if (doFit==4) {
          h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
          if (!doBkg) {
            h_data_fit->Add(h_mcO, -1.);
            h_data_fit->Add(h_mcD, -1.);
            h_data_fit->Add(h_mc2, -1.);
            h_data_fit->Add(h_mc1, -1.);
            if (h_mc1b) h_data_fit->Add(h_mc1b, -1.);
            if (h_mc1c) h_data_fit->Add(h_mc1c, -1.);
          }
          h_mc_fit0 = h_mc1bb;
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
	      if (bbBkg || bbSig) h_mc1bb->SetBinContent(i, 0);
	      if (bbBkg || bbSig) h_mc1bb->SetBinError(i, 0);
	      h_mc2->SetBinContent(i, 0);
	      h_mc2->SetBinError(i, 0);
	      h_mcD->SetBinContent(i, 0);
	      h_mcD->SetBinError(i, 0);
	      h_mcO->SetBinContent(i, 0);
	      h_mcO->SetBinError(i, 0);
	    }
	  }
          fitter = TVirtualFitter::Fitter(0, 1);
          fitter->SetFCN(fcn);
          double arglist[1] = {-1.0};
          fitter->ExecuteCommand("SET PRINT", arglist, 1);
          fitter->SetParameter(0, "c(bb)", 1.00, 0.01, 0.00, 100.00);
          fitter->ExecuteCommand("MIGRAD", arglist, 0);
          h_mc_fit0->Scale(fitter->GetParameter(0));
        }

        if (doFit==5) {
          h_data_fit = (TH1F*)h_data->Clone("h_data_fit");
          if (!doBkg) {
            h_data_fit->Add(h_mcD, -1.);
            h_data_fit->Add(h_mcO, -1.);
            h_data_fit->Add(h_mc2, -1.);
          }
          h_mc_fit0 = h_mc1b;

	  h_mc1->Scale(2.355); /* scale uds,c to the svtx fit values before fitting the b */
	  h_mc1c->Scale(1.053);
	  if (h_mc1) h_data_fit->Add(h_mc1, -1.);
          if (h_mc1c) h_data_fit->Add(h_mc1c, -1.);
          fitter = TVirtualFitter::Fitter(0, 1);
          fitter->SetFCN(fcn);
          double arglist[1] = {-1.0};

          fitter->ExecuteCommand("SET PRINT", arglist, 1);
          fitter->SetParameter(0, "c(b)", 1.00, 0.01, 0.00, 100.00);
          fitter->ExecuteCommand("MIGRAD", arglist, 0);
          h_mc_fit0->Scale(fitter->GetParameter(0));
        }

#ifdef PAPER
	h_data->Scale(1., "width");

	h_mcO->Scale(1., "width");
	h_mcD->Scale(1., "width");
	h_mc2->Scale(1., "width");
	if (h_mc1bb) h_mc1bb->Scale(1., "width");
	if (h_mc1b) h_mc1b->Scale(1., "width");
	if (h_mc1c) h_mc1c->Scale(1., "width");
	h_mc1->Scale(1., "width");
#endif

	TH1F *ht = (TH1F*)h_mc1->Clone("ht");
	ht->Reset();
	if (!doBkg) {
	  ht->Add(h_mcO);
	  ht->Add(h_mcD);
	  ht->Add(h_mc2);
          if (bbBkg) ht->Add(h_mc1bb);
	  if (bbSig) {
	    if (h_mc1b) ht->Add(h_mc1b);
	    if (h_mc1c) ht->Add(h_mc1c);
	    ht->Add(h_mc1);
	  }
	}
        if (bbSig) {
	  ht->Add(h_mc1bb);
        } else {
	  if (h_mc1b) ht->Add(h_mc1b);
	  if (h_mc1c) ht->Add(h_mc1c);
	  ht->Add(h_mc1);
	}

	THStack *hs = new THStack("hs","");
	if (!doBkg) {
	  hs->Add(h_mcO);
	  hs->Add(h_mcD);
	  hs->Add(h_mc2);
#ifdef PAPER
	  if (bbSig) {
	    hs->Add(h_mc1);
	    if (h_mc1c) hs->Add(h_mc1c);
	    if (h_mc1b) hs->Add(h_mc1b);
	  }
          if (bbBkg) hs->Add(h_mc1bb);
	}
        if (bbSig) {
	  hs->Add(h_mc1bb);
	} else {
	  hs->Add(h_mc1);
	  if (h_mc1c) hs->Add(h_mc1c);
	  if (h_mc1b) hs->Add(h_mc1b);
	}
#else
          if (bbBkg) hs->Add(h_mc1bb);
	  if (bbSig) {
	    if (h_mc1b) hs->Add(h_mc1b);
	    if (h_mc1c) hs->Add(h_mc1c);
	    hs->Add(h_mc1);
	  }
	}
        if (bbSig) {
	  hs->Add(h_mc1bb);
	} else {
	  if (h_mc1b) hs->Add(h_mc1b);
	  if (h_mc1c) hs->Add(h_mc1c);
	  hs->Add(h_mc1);
	}
#endif

        TCanvas* c1 = 0;
#ifdef PAPER
	c1 = new TCanvas("c","c",10,10,600,600);
#else
	c1 = new TCanvas("c","c",10,10,800,600);
#endif
	c1->cd();

	TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
	pad1->SetBottomMargin(0.001);
	pad1->Draw();
	pad1->cd();
	pad1->SetLogy();
	//if (title.find("MET")!=string::npos)   pad1->SetLogy(0);
	//if (title.find("mass_")!=string::npos) pad1->SetLogy(0);
        if (title.find("y_Z")!=string::npos) pad1->SetLogy(0);
        if (title.find("DR_bb")!=string::npos) pad1->SetLogy(0);
        if (title.find("delta_j")!=string::npos) pad1->SetLogy(0);
        if (title.find("delta_j_n")!=string::npos) pad1->SetLogy(0);
//        if (title.find("DR_jj")!=string::npos) pad1->SetLogy(0);
//        if (title.find("w_eebb_mass")!=string::npos || title.find("w_mmbb_mass")!=string::npos) pad1->SetLogy(0);

	hs->Draw("HIST");
	hs->GetYaxis()->SetTitle("Events");
	hs->GetYaxis()->SetTitleSize(0.05);
 	hs->GetYaxis()->SetLabelSize(0.045);
	hs->GetYaxis()->SetTitleOffset(0.7);
#ifdef PAPER
	hs->GetYaxis()->SetTitleOffset(1.0);
#endif
 	hs->GetXaxis()->SetLabelSize(0.08);
	hs->GetXaxis()->SetTitleOffset(0.7);
	hs->SetMinimum(8);
#ifdef PAPER
	hs->SetMinimum(0.5);
	hs->SetMaximum(1.2*hs->GetMaximum());
	if (title=="w_SVTX_mass") hs->SetMinimum(5);
	if (title.find("w_mass_")!=string::npos) {
	  hs->SetMaximum(1.5*hs->GetMaximum());
	  hs->SetMinimum(8);
	}
#endif

	h_data->Draw("EPX0SAMES");
	h_data->SetMarkerColor(kBlack);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize (1.0);
	h_data->SetStats(0);

	TLegend *leg;
	if (doBkg) {
	  if (h_mc1c && h_mc1b && bbBkg) {
	    leg = new TLegend(0.65, 0.640, 0.91, 0.88);
	  } else if (h_mc1c && h_mc1b && bbSig) {
	    leg = new TLegend(0.65, 0.760, 0.91, 0.88);
	  } else if (h_mc1c && h_mc1b) {
	    leg = new TLegend(0.65, 0.640, 0.91, 0.88);
	  } else if (h_mc1c || h_mc1b) {
	    leg = new TLegend(0.65, 0.740, 0.91, 0.88);
	  } else {
	    leg = new TLegend(0.65, 0.780, 0.91, 0.88);
	  }
	} else {
          if (h_mc1c && h_mc1b && bbBkg) {
            leg = new TLegend(0.65, 0.547, 0.91, 0.88);
          } else if (h_mc1c && h_mc1b && bbSig) {
            leg = new TLegend(0.65, 0.577, 0.91, 0.88);
	  } else if (h_mc1c && h_mc1b) {
#ifdef PAPER
	    leg = new TLegend(0.66, 0.520, 0.94, 0.88);
#else
	    leg = new TLegend(0.65, 0.580, 0.91, 0.88);
#endif
	  } else if (h_mc1c || h_mc1b) {
	    leg = new TLegend(0.65, 0.613, 0.91, 0.88);
	  } else {
	    leg = new TLegend(0.65, 0.647, 0.91, 0.88);
	  }
	}
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (ilepton==1) leg->AddEntry(h_data,"Data","p");
	if (ilepton==2) leg->AddEntry(h_data,"Data","p");
	if (ilepton==3) leg->AddEntry(h_data,"Data","p");

#ifdef PAPER
	if (!bbSig) {
	  if (h_mc1b) leg->AddEntry(h_mc1b,"Z+b jets","f");
	  if (h_mc1c) leg->AddEntry(h_mc1c,"Z+c jets","f");
	  if (h_mc1c && h_mc1b) {
	    leg->AddEntry(h_mc1,"Z+udsg jets","f");
	  } else {
	    leg->AddEntry(h_mc1,"Z+jets","f");
	  }
	}
	if (bbSig) {
	  if (!doBkg) {
	    if (h_mc1b) leg->AddEntry(h_mc1b,"Z+b jets","f");
	    if (h_mc1c) leg->AddEntry(h_mc1c,"Z+c jets","f");
	    if (h_mc1c && h_mc1b) {
	      leg->AddEntry(h_mc1,"Z+udsg jets","f");
	    } else {
	      leg->AddEntry(h_mc1,"Z+jets","f");
	    }
	  }
	  leg->AddEntry(h_mc1bb,"Z+bb jets","f");
	}
#else
	if (!bbSig) {
	  if (h_mc1c && h_mc1b) {
	    leg->AddEntry(h_mc1,"Z+udsg-jets","f");
	  } else {
	    leg->AddEntry(h_mc1,"Z+jets","f");
	  }
	  if (h_mc1c) leg->AddEntry(h_mc1c,"Z+c-jets","f");
	  if (h_mc1b) leg->AddEntry(h_mc1b,"Z+b-jets","f");
	}
	if (bbSig) {
	  if (!doBkg) {
	    if (h_mc1c && h_mc1b) {
	      leg->AddEntry(h_mc1,"Z+udsg-jets","f");
	    } else {
	      leg->AddEntry(h_mc1,"Z+jets","f");
	    }
	    if (h_mc1c) leg->AddEntry(h_mc1c,"Z+c-jets","f");
	    if (h_mc1b) leg->AddEntry(h_mc1b,"Z+b-jets","f");
	  }
	  leg->AddEntry(h_mc1bb,"Z+bb-jets","f");
	}
#endif
	if (!doBkg) {
          if (bbBkg) leg->AddEntry(h_mc1bb,"Z+bb-jets","f");
	  leg->AddEntry(h_mc2,"t#bar{t}","f");
	  leg->AddEntry(h_mcD,"Dibosons","f");
	  leg->AddEntry(h_mcO,"Others","f");
//	  leg->AddEntry(h_mc5,"QCD","f");
	}
	leg->Draw();

	pad1->Update();
	c1->Update();

	c1->cd();

	TH1F *h_ratio = (TH1F*)h_data->Clone("h_ratio");

	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();
	h_ratio->SetTitle("");
	h_ratio->SetStats(0);

	if (title=="w_jetmultiplicity") {
	  h_ratio->GetXaxis ()->SetTitle("jet multiplicity");
	} else if (title=="w_mass_ee"||title=="w_mass_mm"||title=="w_mass_em") {
	  h_ratio->GetXaxis ()->SetTitle("invariant mass [GeV/c^{2}]");
	} else if (title=="w_mass_ee_b"||title=="w_mass_mm_b"||title=="w_mass_em_b") {
	  h_ratio->GetXaxis ()->SetTitle("invariant mass [GeV/c^{2}]");
	} else if (title=="w_first_muon_pt") {
	  h_ratio->GetXaxis ()->SetTitle("muon p_{T} [GeV/c]");
	} else if (title=="w_first_ele_pt") {
	  h_ratio->GetXaxis ()->SetTitle("electron p_{T} [GeV/c]");
	} else if (title=="w_first_jet_pt") {
	  h_ratio->GetXaxis ()->SetTitle("jet p_{T} [GeV/c]");
	} else if (title=="w_secondvtx_N") {
	  h_ratio->GetXaxis ()->SetTitle("CSV discriminator");
	} else if (title=="w_recoVTX" || title=="h_recoVTX") {
	  h_ratio->GetXaxis ()->SetTitle("Number of offline vertices");
	} else if (title=="w_muon_pt") {
	  h_ratio->GetXaxis ()->SetTitle("muon p_{T} [GeV/c]");
	} else if (title=="w_MET"||title=="w_MET_b") {
	  h_ratio->GetXaxis ()->SetTitle("MET [GeV/c]");
	} else if (title=="w_MET_sign"||title=="w_MET_sign_b") {
	  h_ratio->GetXaxis ()->SetTitle("MET Significance");
	} else if (title=="w_Ht"||title=="w_Ht_b") {
	  h_ratio->GetXaxis ()->SetTitle("H_{T} [GeV/c]");
	} else if (title=="w_mass_Zj_ee"||title=="w_mass_Zj_mm"||title=="w_mass_Zj_em") {
	  h_ratio->GetXaxis ()->SetTitle("Zj invariant mass [GeV/c^{2}]");
	} else if (title=="w_mass_Zj_ee_b"||title=="w_mass_Zj_mm_b"||title=="w_mass_Zj_em_b") {
	  h_ratio->GetXaxis ()->SetTitle("Zb invariant mass [GeV/c^{2}]");
	} else if (title=="w_pt_Z_ee"||title=="w_pt_Z_mm"||title=="w_pt_Z_em") {
	  h_ratio->GetXaxis ()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_y_Z_ee_b"||title=="w_y_Z_mm_b"||title=="w_y_Z_em_b") {
          h_ratio->GetXaxis ()->SetTitle("y_{Z}");
        } else if (title=="w_bjetmultiplicity") {
	  h_ratio->GetXaxis ()->SetTitle("jet multiplicity");
        } else if (title=="w_bjetmultiplicity_exc") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("b jet multiplicity");
#else
	  h_ratio->GetXaxis ()->SetTitle("b-jet multiplicity");
#endif
	} else if (title=="w_first_jet_pt_b") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("leading b jet p_{T} [GeV/c]");
#else
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet p_{T} [GeV/c]");
#endif
	} else if (title=="w_first_jet_eta_b") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("leading b jet #eta");
#else
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet #eta");
#endif
	} else if (title=="w_second_jet_pt_b") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("subleading b jet p_{T} [GeV/c]");
#else
	  h_ratio->GetXaxis ()->SetTitle("subleading b-jet p_{T} [GeV/c]");
#endif
	} else if (title=="w_second_jet_eta_b") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("subleading b jet #eta");
#else
	  h_ratio->GetXaxis ()->SetTitle("subleading b-jet #eta");
#endif
	} else if (title=="w_third_jet_pt_b") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("subsubleading b jet p_{T} [GeV/c]");
#else
	  h_ratio->GetXaxis ()->SetTitle("subsubleading b-jet p_{T} [GeV/c]");
#endif
	} else if (title=="w_third_jet_eta_b") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("subsubleading b jet #eta");
#else
	  h_ratio->GetXaxis ()->SetTitle("subsubleading b-jet #eta");
#endif
	} else if (title=="w_mass_ee_b"||title=="w_mm_mass_b"||title=="w_em_mass_b") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("Z mass + (#geq 1 b jet) [GeV/c^{2}]");
#else
	  h_ratio->GetXaxis ()->SetTitle("Z mass + (#geq 1 b-jet) [GeV/c^{2}]");
#endif
	} else if (title=="w_pt_Z_ee_b"||title=="w_pt_Z_mm_b"||title=="w_pt_Z_em_b") {
	  h_ratio->GetXaxis ()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_delta_phi_ee"||title=="w_delta_phi_mm"||title=="w_delta_phi_em") {
	  h_ratio->GetXaxis ()->SetTitle("#Delta#phi(jZ) [rad]");
	} else if (title=="w_delta_phi_ee_b"||title=="w_delta_phi_mm_b"||title=="w_delta_phi_em_b") {
	  h_ratio->GetXaxis ()->SetTitle("#Delta#phi(bZ) [rad]");
	} else if (title=="w_delta_phi_2b") {
          h_ratio->GetXaxis ()->SetTitle("#Delta#phi(bb) [rad]");
        } else if (title=="w_DR_bb") {
          h_ratio->GetXaxis ()->SetTitle("#Delta R(bb) [rad]");
        } else if (title=="w_DR_eeb_min"||title=="w_DR_mmb_min"||title=="w_DR_emb_min") {
          h_ratio->GetXaxis ()->SetTitle("#Delta R(Zb) min [rad]");
        } else if (title=="w_DR_eeb_max"||title=="w_DR_mmb_max"||title=="w_DR_emb_max") {
          h_ratio->GetXaxis ()->SetTitle("#Delta R(Zb) max [rad]");
	} else if (title=="w_A_eeb"||title=="w_A_mmb") {
	  h_ratio->GetXaxis ()->SetTitle("A_{Zbb}");
        } else if (title=="w_bb_mass") {
          h_ratio->GetXaxis ()->SetTitle("m(bb) [GeV/c^{2}]");
        } else if (title=="w_eebb_mass"||title=="w_mmbb_mass"||title=="w_embb_mass") {
          h_ratio->GetXaxis ()->SetTitle("m(Zbb) [GeV/c^{2}]");
        } else if (title=="w_SVTX_mass_jet"||title=="w_SVTX_mass_trk"||title=="w_SVTX_mass") {
	  h_ratio->GetXaxis ()->SetTitle("SV mass [GeV/c^{2}]");
	} else if (title=="w_BJP"||title=="w_BJP0"||title=="w_BJP1"||title=="w_BJP2"||title=="w_BJP_mass"){
	  h_ratio->GetXaxis ()->SetTitle("Jet B Probability Discriminator");
	} else if (title=="w_JBP"||title=="w_JBP_mass"){
	  h_ratio->GetXaxis ()->SetTitle("Jet Probability Discriminator");
	} else if (title=="w_first_bjet_pt") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("leading b jet p_{T} [GeV/c]");
#else
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet p_{T} [GeV/c]");
#endif
	} else if (title=="w_first_bjet_eta") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("leading b jet #eta");
#else
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet #eta");
#endif
	} else if (title=="w_first_bjet_eta_abs") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("leading b jet |#eta|");
#else
	  h_ratio->GetXaxis ()->SetTitle("leading b-jet |#eta|");
#endif
	} else if (title=="w_second_bjet_pt") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("subleading b jet p_{T} [GeV/c]");
#else
	  h_ratio->GetXaxis ()->SetTitle("subleading b-jet p_{T} [GeV/c]");
#endif
	} else if (title=="w_second_bjet_eta") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("subleading b jet #eta");
#else
	  h_ratio->GetXaxis ()->SetTitle("subleading b-jet #eta");
#endif
	} else if (title=="w_second_bjet_eta_abs") {
#ifdef PAPER
	  h_ratio->GetXaxis ()->SetTitle("subleading b jet |#eta|");
#else
	  h_ratio->GetXaxis ()->SetTitle("subleading b-jet |#eta|");
#endif
	} else if (title=="w_flightd") {
	  h_ratio->GetXaxis ()->SetTitle("L_{xy} [cm]");
	} else if (title=="w_flightd_sig") {
	  h_ratio->GetXaxis ()->SetTitle("L_{xy} significance");
	} else if (title=="w_dxy") {
	  h_ratio->GetXaxis ()->SetTitle("d_{xy} [cm]");
	} else if (title=="w_dxy_sig") {
	  h_ratio->GetXaxis ()->SetTitle("d_{xy} significance");
	} else if (title=="w_N_SV") {
	  h_ratio->GetXaxis ()->SetTitle("SV multiplicity");
	} else if (title=="w_SV_NTracks") {
	  h_ratio->GetXaxis ()->SetTitle("Tracks multiplicity at SV");
	}
#ifdef PAPER
	string tmp = h_ratio->GetXaxis()->GetTitle();
	if (tmp.find("/c^{2}")!=string::npos) {
	  tmp.erase(tmp.find("/c^{2}"), 6);
	  h_ratio->GetXaxis ()->SetTitle(tmp.c_str());
	}
	if (tmp.find("/c")!=string::npos) {
	  tmp.erase(tmp.find("/c"), 2);
	  h_ratio->GetXaxis ()->SetTitle(tmp.c_str());
	}
	if (tmp.find("[")!=string::npos) {
	  tmp.replace(tmp.find("["), 1,"(");
	  h_ratio->GetXaxis ()->SetTitle(tmp.c_str());
	}
	if (tmp.find("]")!=string::npos) {
	  tmp.replace(tmp.find("]"), 1, ")");
	  h_ratio->GetXaxis ()->SetTitle(tmp.c_str());
	}
	if (title.find("w_mass_ee")!=string::npos) {
	  h_ratio->GetXaxis()->SetTitle("dielectron invariant mass (GeV)");
	}
	if (title.find("w_mass_mm")!=string::npos) {
	  h_ratio->GetXaxis()->SetTitle("dimuon invariant mass (GeV)");
	}
#endif
	h_ratio->GetXaxis()->SetTitleOffset(0.9);
#ifdef PAPER
	h_ratio->GetXaxis()->SetTitleOffset(1.1);
#endif
 	h_ratio->GetXaxis()->SetTitleSize(0.11);
	h_ratio->GetXaxis()->SetLabelFont(42);
 	h_ratio->GetXaxis()->SetLabelSize(0.10);
	h_ratio->GetXaxis()->SetTitleFont(42);
	h_ratio->GetYaxis()->SetTitle("Data/MC");
	h_ratio->GetYaxis()->SetNdivisions(505);
	h_ratio->GetYaxis()->SetTitleSize(0.11);
	h_ratio->GetYaxis()->SetLabelSize(0.10);
	h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
	h_ratio->GetYaxis()->SetTitleOffset(0.33);
#ifdef PAPER
	h_ratio->GetYaxis()->SetTitleOffset(0.47);
#endif
	h_ratio->Divide(ht);
	h_ratio->SetMarkerStyle(20);
	h_ratio->Draw("E0PX0");

        TH1F *h_ratio_pull =  new TH1F ("h_pull", "h_pull", 50, 0, 6);
        TH1F *h_ratio_pull2 = new TH1F ("h_pull2", "h_pull2", 50, 0, 6);
        TH1F *h_ratio_pull3 = new TH1F ("h_pull3", "h_pull3", 50, 0, 6);

        for (int j=1; j<=h_ratio->GetNbinsX(); j++) {
           h_ratio_pull2->Fill(h_ratio->GetBinContent(j));
        }
        if (doFit==3) {
          for (int i=1; i<=h_ratio->GetNbinsX(); i++) {
               h_ratio_pull->Fill(h_ratio->GetBinContent(i));
          }
	}

        for (int l=1; l<=h_ratio_pull->GetNbinsX(); l++) {
           for (int m=1; m<=h_ratio_pull2->GetNbinsX(); m++) {
       	      double pull = h_ratio_pull2->GetBinContent(m)-h_ratio_pull->GetBinContent(l);
	      h_ratio_pull3->Fill(pull);
	  }
	}

/*
	TFile* dumphistos_file = new TFile("dump.root","RECREATE");
        dumphistos_file->cd();
        //if (title=="w_first_bjet_pt") h_ratio->Write("A");
        if (title=="w_SVTX_mass") h_ratio_pull->Write("A");
        if (title=="w_SVTX_mass") h_ratio_pull2->Write("B");
        if (title=="w_SVTX_mass") h_ratio_pull3->Write("C");
        //if (title=="w_first_bjet_pt_SVTX") h_ratio->Write("B");
	dumphistos_file->Close();
*/

	TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
	OLine->SetLineColor(kRed);
	OLine->SetLineWidth(2);
	OLine->Draw();

#ifdef PAPER
	if (title=="w_bjetmultiplicity_exc") {
	  hs->GetXaxis()->SetRangeUser(0.5, 4.5);
	  hs->SetMinimum(0.01);
	  hs->SetMaximum(100000);
	  hs->GetYaxis()->SetTitle("Events");
	  h_ratio->GetXaxis()->SetRangeUser(0.5, 4.5);
	  h_ratio->GetXaxis()->SetNdivisions(105);
	  h_ratio->GetXaxis()->SetTitle("b jet multiplicity");
	  OLine->SetX2(4.5);
	}
#endif

        writeExtraText = true;
#ifdef PAPER
        writeExtraText = false;
#endif
        lumi_8TeV  = Form("%.1f fb^{-1}", Lumi2012/1000.);
        int iPeriod = 2;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV
        // second parameter drives the position of the CMS logo in the plot
        // iPos=11 : top-left, left-aligned
        // iPos=33 : top-right, right-aligned
        // iPos=22 : center, centered
        // mode generally :
        // iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
        int iPos = 0;
        CMS_lumi(pad1,  iPeriod, iPos);
        c1->cd();

#ifndef PAPER
	if (numB==0 && title.find("_b")==string::npos) {
          if (title=="w_jetmultiplicity") {
            if (ilepton ==1) CMS_process("Z/#gamma*#rightarrow ee selection", 0.44, 0.9);
            if (ilepton ==2) CMS_process("Z/#gamma*#rightarrow #mu#mu selection", 0.44, 0.9);
          }
          if (title=="w_mass_ee" || title=="w_mass_mm") {
            if (ilepton ==1) CMS_process("Z/#gamma*#rightarrow ee selection", 0.135, 0.87);
            if (ilepton ==2) CMS_process("Z/#gamma*#rightarrow #mu#mu selection", 0.135, 0.87);
          }
        }
        if (numB==0 && title.find("_b")!=string::npos) {
          if (title=="w_bjetmultiplicity" || title=="w_bjetmultiplicity_exc") {
            if (ilepton ==1) CMS_process("Z+(#geq1)b-jet selection", 0.34, 0.87);
            if (ilepton ==2) CMS_process("Z+(#geq1)b-jet selection", 0.34, 0.87);
          }
          if (title=="w_mass_ee_b" || title=="w_mass_mm_b") {
            if (ilepton ==1) CMS_process("Z+(#geq1)b-jet selection", 0.135, 0.87);
            if (ilepton ==2) CMS_process("Z+(#geq1)b-jet selection", 0.135, 0.87);
          }
        }
#endif

        TLatex * lab = new TLatex ();
        lab->SetTextSize (0.0275);
        lab->SetTextFont (42);
        lab->DrawLatex (0.63, 0.65, bSel.c_str());

	if (doFit) {
	  TLatex *fitLabel = new TLatex();
	  fitLabel->SetTextSize(0.035);
	  fitLabel->SetTextFont(42);
	  fitLabel->SetLineWidth(2);
	  fitLabel->SetNDC();
	  char buff[100];
	  if (doFit==1) {
	    sprintf(buff, "c_{t} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	  }
	  if (doFit==2) {
	    //sprintf(buff, "c_{Z+jets} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    sprintf(buff, "c_{Others} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
	    fitLabel->DrawLatex(0.68, 0.48, buff);
	    sprintf(buff, "c_{t} = %5.3f #pm %5.3f", fitter->GetParameter(1), fitter->GetParError(1));
	    fitLabel->DrawLatex(0.68, 0.43, buff);
	  }
	  if (doFit==3) {
	    sprintf(buff, "c_{uds} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
#ifdef PAPER
	    fitLabel->DrawLatex(0.15, 0.48, buff);
#else
	    fitLabel->DrawLatex(0.19, 0.48, buff);
#endif
	    sprintf(buff, "c_{b}   = %5.3f #pm %5.3f", fitter->GetParameter(1), fitter->GetParError(1));
#ifdef PAPER
	    fitLabel->DrawLatex(0.15, 0.43, buff);
#else
	    fitLabel->DrawLatex(0.19, 0.43, buff);
#endif
	    sprintf(buff, "c_{c}   = %5.3f #pm %5.3f", fitter->GetParameter(2), fitter->GetParError(2));
#ifdef PAPER
	    fitLabel->DrawLatex(0.15, 0.38, buff);
#else
	    fitLabel->DrawLatex(0.19, 0.38, buff);
#endif
	    float f_uds = 100*h_mc_fit0->Integral(0,h_mc_fit0->GetNbinsX()+1)/(h_mc_fit0->Integral(0,h_mc_fit0->GetNbinsX()+1)+h_mc_fit1->Integral(0,h_mc_fit1->GetNbinsX()+1)+h_mc_fit2->Integral(0,h_mc_fit2->GetNbinsX()+1));
	    float ef_uds = f_uds*(fitter->GetParError(0)/fitter->GetParameter(0));
	    sprintf(buff, "f_{uds} = %4.1f #pm %3.1f %%", f_uds, ef_uds);
	    fitLabel->DrawLatex(0.46, 0.48, buff);
	    float f_b = 100*h_mc_fit1->Integral(0,h_mc_fit1->GetNbinsX()+1)/(h_mc_fit0->Integral(0,h_mc_fit0->GetNbinsX()+1)+h_mc_fit1->Integral(0,h_mc_fit1->GetNbinsX()+1)+h_mc_fit2->Integral(0,h_mc_fit2->GetNbinsX()+1));
	    float ef_b = f_b*(fitter->GetParError(1)/fitter->GetParameter(1));
	    sprintf(buff, "f_{b}   = %4.1f #pm %3.1f %%", f_b, ef_b);
	    fitLabel->DrawLatex(0.46, 0.43, buff);
	    float f_c = 100*h_mc_fit2->Integral(0,h_mc_fit2->GetNbinsX()+1)/(h_mc_fit0->Integral(0,h_mc_fit0->GetNbinsX()+1)+h_mc_fit1->Integral(0,h_mc_fit1->GetNbinsX()+1)+h_mc_fit2->Integral(0,h_mc_fit2->GetNbinsX()+1));
	    float ef_c = f_c*(fitter->GetParError(2)/fitter->GetParameter(2));
	    sprintf(buff, "f_{c}   = %4.1f #pm %3.1f %%", f_c, ef_c);
	    fitLabel->DrawLatex(0.46, 0.38, buff);
	  }
          if (doFit==4) {
            sprintf(buff, "c_{bb} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
#ifdef PAPER
            fitLabel->DrawLatex(0.20, 0.48, buff);
#else
            fitLabel->DrawLatex(0.68, 0.48, buff);
#endif
          }
          if (doFit==5) {
            sprintf(buff, "c_{b} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
            fitLabel->DrawLatex(0.68, 0.48, buff);
          }
	}

//	if (useWeights) subdir = subdir + "_weights";
//	if (useDY==1) subdir = subdir + "_sherpa";
//	if (useDY==2) version = version + "_powheg";
//	if (useDY==3) version = version + "_aMC@NLO";

#ifdef PAPER
	path = path + "/paper/";
#endif

	if (plot) {
	  if (doBkg) title = title + "_doBkg";
	  if (doFit) title = title + "_doFit";
	  ofstream out,out2;
	  if (ilepton==1) {
	    gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/distributions" + dirbSel + "/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electrons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/electrons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + title + ".dat").c_str());
	  }
	  if (ilepton==2) {
	    gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/").c_str(), kTRUE);
	    c1->SaveAs((path + "/muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + title + ".dat").c_str());
	  }
	  if (ilepton==3) {
	    gSystem->mkdir((path + "/electrons+muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electrons+muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + title + ".pdf").c_str());
	    if (doFit) out.open((path + "/electrons+muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + title + ".dat").c_str());
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
          if (doFit==4) {
            out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
            out.close();
          }

	}
}

