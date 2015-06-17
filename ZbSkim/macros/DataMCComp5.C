#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v14.h"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";
//string path = "/gpfs/cms/users/lalicata/work/test/data";

TH1F* h_data = 0;
TH1F* h_data_fit = 0;

void fcn(int& npar, double* gin, double& fun, double* par, int iflag) {
  double chisq = 0.0;
  if (iflag) {};
  if (gin) {};
  for (int i=1; i<=h_data->GetNbinsX(); i++) {
    double xn = h_data->GetBinContent(i);
    double xd = TMath::Power(h_data->GetBinError(i),2);
    if (npar>0) {
      xn = xn - par[0]*h_data_fit->GetBinContent(i);
      xd = xd + TMath::Power(par[0]*h_data_fit->GetBinError(i),2);
    }
    if (xd!=0) chisq = chisq + (xn*xn)/xd;
  }
  fun = chisq;
}

void DataMCComp5(int irun=0, string title="", int plot=0, int ilepton=1, int doFit=0, int numB=0) {

//int useFitResults=0; // use MC predictions for c_t
int useFitResults=1;  // use fit results for c_t

//int useEleMuo = 0; // use MC or fit results for c_t
int useEleMuo = 1; // use e-mu fit results for c_t

string bSel="Z + (#geq 1) b-jet";
string subdir="0";
string postfix="";
string dirbSel="";

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
  bSel = "Z + (= 1) b-jet";
}
if (numB==2) {
  postfix = postfix + "2b";
  dirbSel = "_2b";
  bSel = "Z + (#geq 2) b-jet";
}

      /* top background */

      double c1_t=1.0;
      double ec1_t=0.0;
      double c2_t=1.0;
      double ec2_t=0.0;

      ifstream in4, in5;
      if (ilepton==1) {
        if (useFitResults) {
          in4.open((path + "/electrons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_doFit" + ".dat").c_str());
          in5.open((path + "/electrons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_b_doFit" + ".dat").c_str());
        }
      }
      if (ilepton==2) {
        if (useFitResults) {
          in4.open((path + "/muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_doFit" + ".dat").c_str());
          in5.open((path + "/muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_b_doFit" + ".dat").c_str());
        }
      }
      if (useFitResults) {
        in4 >> c1_t >> ec1_t;
        in5 >> c2_t >> ec2_t;
        in4.close();
        in5.close();
      }

      double Lumi2012=0;

      if (ilepton==1) Lumi2012 = Lumi2012_ele;
      if (ilepton==2) Lumi2012 = Lumi2012_muon;

      double norm1 = ((Lumi2012 * Xsec_dy) / Ngen_dy);
      double norm2 = ((Lumi2012 * Xsec_tt) / Ngen_tt);
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

      double norm1_fit = ((Lumi2012_ele_muon * Xsec_dy) / Ngen_dy);
      double norm2_fit = ((Lumi2012_ele_muon * Xsec_tt) / Ngen_tt);
      double norm3_fit = ((Lumi2012_ele_muon * Xsec_zz) / Ngen_zz);
      double norm4_fit = ((Lumi2012_ele_muon * Xsec_wz) / Ngen_wz);
      double norm5_fit = ((Lumi2012_ele_muon * Xsec_qcd) / Ngen_qcd);
      double norm6_fit = ((Lumi2012_ele_muon * Xsec_ww) / Ngen_ww);
      double norm7_fit = ((Lumi2012_ele_muon * Xsec_wj) / Ngen_wj);
      double norm8_fit = ((Lumi2012 * Xsec_tS) / Ngen_tS);
      double norm9_fit = ((Lumi2012 * Xsec_tT) / Ngen_tT);
      double norm10_fit = ((Lumi2012 * Xsec_tW) / Ngen_tW);
      double norm11_fit = ((Lumi2012 * Xsec_tSb) / Ngen_tSb);
      double norm12_fit = ((Lumi2012 * Xsec_tTb) / Ngen_tTb);
      double norm13_fit = ((Lumi2012 * Xsec_tWb) / Ngen_tWb);

      double enorm1_fit = ((Lumi2012_ele_muon * eXsec_dy) / Ngen_dy);
      double enorm2_fit = ((Lumi2012_ele_muon * eXsec_tt) / Ngen_tt);
      if (useEleMuo && ilepton!=3) enorm2 = 0;
      double enorm3_fit = ((Lumi2012_ele_muon * eXsec_zz) / Ngen_zz);
      double enorm4_fit = ((Lumi2012_ele_muon * eXsec_wz) / Ngen_wz);
      double enorm5_fit = ((Lumi2012_ele_muon * eXsec_qcd) / Ngen_qcd);
      double enorm6_fit = ((Lumi2012_ele_muon * eXsec_ww) / Ngen_ww);
      double enorm7_fit = ((Lumi2012_ele_muon * eXsec_wj) / Ngen_wj);
      double enorm8_fit = ((Lumi2012 * eXsec_tS) / Ngen_tS);
      double enorm9_fit = ((Lumi2012 * eXsec_tT) / Ngen_tT);
      double enorm10_fit = ((Lumi2012 * eXsec_tW) / Ngen_tW);
      double enorm11_fit = ((Lumi2012 * eXsec_tSb) / Ngen_tSb);
      double enorm12_fit = ((Lumi2012 * eXsec_tTb) / Ngen_tTb);
      double enorm13_fit = ((Lumi2012 * eXsec_tWb) / Ngen_tWb);

      if (title.empty()) title = "w_jetmultiplicity";

      if (ilepton==1) {
        if (title.find("muon")!=string::npos) return;
	if (title.find("mm")!=string::npos) return;
      }
      if (ilepton==2) {
	if (title.find("ele")!=string::npos) return;
	if (title.find("ee")!=string::npos) return;
      }

      TFile *data=0;
      if (ilepton==1) data = TFile::Open((path + "/" + version + "/" + "DoubleElectron_2012_merge.root").c_str());
      if (ilepton==2) data = TFile::Open((path + "/" + version + "/" + "DoubleMu_2012_merge.root").c_str());

      TFile *data_fit=0;
      data_fit = TFile::Open((path + "/" + version + "/" + "MuEG_2012_merge.root").c_str());

      TFile *mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL.root").c_str());
      TFile *mc2 = TFile::Open((path + "/" + version + "/" + "TTbar.root").c_str());
      TFile *mc3 = TFile::Open((path + "/" + version + "/" + "ZZ.root").c_str());
      TFile *mc4 = TFile::Open((path + "/" + version + "/" + "WZ.root").c_str());
//    TFile *mc5 = TFile::Open((path + "/" + version + "/" + "QCD.root").c_str());
      TFile *mc6 = TFile::Open((path + "/" + version + "/" + "WW.root").c_str());
      TFile *mc7 = TFile::Open((path + "/" + version + "/" + "Wj.root").c_str());
      TFile *mc8 = TFile::Open((path + "/" + version + "/" + "T_s.root").c_str());
      TFile *mc9 = TFile::Open((path + "/" + version + "/" + "T_t.root").c_str());
      TFile *mc10 = TFile::Open((path + "/" + version + "/" + "T_tW.root").c_str());
      TFile *mc11 = TFile::Open((path + "/" + version + "/" + "TBar_s.root").c_str());
      TFile *mc12 = TFile::Open((path + "/" + version + "/" + "TBar_t.root").c_str());
      TFile *mc13 = TFile::Open((path + "/" + version + "/" + "TBar_tW.root").c_str());

      string title_fit = title;

      if (title_fit.find("ee")!=string::npos) {
        title_fit.replace(title_fit.find("ee")+1, 1, "m");
      }
      if (title_fit.find("mm")!=string::npos) {
        title_fit.replace(title_fit.find("mm"), 1, "e");
      }

      if (ilepton==1) data->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) data->cd(("demoMuo"+postfix).c_str());
      h_data = (TH1F*)gDirectory->Get(title.c_str());
      data_fit->cd(("demoEleMuo"+postfix).c_str());
      h_data_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc1->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc1->cd(("demoMuo"+postfix).c_str());
      TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
      mc1->cd(("demoEleMuo"+postfix).c_str());
      TH1F* h_mc1_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc2->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc2->cd(("demoMuo"+postfix).c_str());
      TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
      mc2->cd(("demoEleMuo"+postfix).c_str());
      TH1F* h_mc2_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc3->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc3->cd(("demoMuo"+postfix).c_str());
      TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());
      mc3->cd(("demoEleMuo"+postfix).c_str());
      TH1F* h_mc3_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc4->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc4->cd(("demoMuo"+postfix).c_str());
      TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());
      mc4->cd(("demoEleMuo"+postfix).c_str());
      TH1F* h_mc4_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

//    if (ilepton==1) mc5->cd(("demoEle"+postfix).c_str());
//    if (ilepton==2) mc5->cd(("demoMuo"+postfix).c_str());
//    TH1F* h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
//    mc5->cd(("demoEleMuo"+postfix).c_str());
//    TH1F* h_mc5_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc6->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc6->cd(("demoMuo"+postfix).c_str());
      TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());
      mc6->cd(("demoEleMuo"+postfix).c_str());
      TH1F* h_mc6_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc7->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc7->cd(("demoMuo"+postfix).c_str());
      TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());
      mc7->cd(("demoEleMuo"+postfix).c_str());
      TH1F* h_mc7_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc8->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc8->cd(("demoMuo"+postfix).c_str());
      TH1F* h_mc8 = (TH1F*)gDirectory->Get(title.c_str());
      mc8->cd(("demoEleMuo"+postfix).c_str());
      TH1F* h_mc8_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc9->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc9->cd(("demoMuo"+postfix).c_str());
      TH1F* h_mc9 = (TH1F*)gDirectory->Get(title.c_str());
      mc9->cd(("demoEleMuo"+postfix).c_str());
      TH1F* h_mc9_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc10->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc10->cd(("demoMuo"+postfix).c_str());
      TH1F* h_mc10 = (TH1F*)gDirectory->Get(title.c_str());
      mc10->cd(("demoEleMuo"+postfix).c_str());
      TH1F* h_mc10_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc11->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc11->cd(("demoMuo"+postfix).c_str());
      TH1F* h_mc11 = (TH1F*)gDirectory->Get(title.c_str());
      mc11->cd(("demoEleMuo"+postfix).c_str());
      TH1F* h_mc11_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc12->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc12->cd(("demoMuo"+postfix).c_str());
      TH1F* h_mc12 = (TH1F*)gDirectory->Get(title.c_str());
      mc12->cd(("demoEleMuo"+postfix).c_str());
      TH1F* h_mc12_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      if (ilepton==1) mc13->cd(("demoEle"+postfix).c_str());
      if (ilepton==2) mc13->cd(("demoMuo"+postfix).c_str());
      TH1F* h_mc13 = (TH1F*)gDirectory->Get(title.c_str());
      mc13->cd(("demoEleMuo"+postfix).c_str());
      TH1F* h_mc13_fit = (TH1F*)gDirectory->Get(title_fit.c_str());

      h_data->Sumw2();
      h_data_fit->Sumw2();

      h_mc1->Sumw2();
      h_mc2->Sumw2();
      h_mc3->Sumw2();
      h_mc4->Sumw2();
//    h_mc5->Sumw2();
      h_mc6->Sumw2();
      h_mc7->Sumw2();
      h_mc8->Sumw2();
      h_mc9->Sumw2();
      h_mc10->Sumw2();
      h_mc11->Sumw2();
      h_mc12->Sumw2();
      h_mc13->Sumw2();

      h_mc1_fit->Sumw2();
      h_mc2_fit->Sumw2();
      h_mc3_fit->Sumw2();
      h_mc4_fit->Sumw2();
//    h_mc5_fit->Sumw2();
      h_mc6_fit->Sumw2();
      h_mc7_fit->Sumw2();
      h_mc8_fit->Sumw2();
      h_mc9_fit->Sumw2();
      h_mc10_fit->Sumw2();
      h_mc11_fit->Sumw2();
      h_mc12_fit->Sumw2();
      h_mc13_fit->Sumw2();

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
      h_mc2->Scale(norm2);
      h_mc3->Scale(norm3);
      h_mc4->Scale(norm4);
//    h_mc5->Scale(norm5);
      h_mc6->Scale(norm6);
      h_mc7->Scale(norm7);
      h_mc8->Scale(norm8);
      h_mc9->Scale(norm9);
      h_mc10->Scale(norm10);
      h_mc11->Scale(norm11);
      h_mc12->Scale(norm12);
      h_mc13->Scale(norm13);

      if (irun==10) {
        norm1_fit = norm1_fit + 0.1*enorm1_fit;
        norm2_fit = norm2_fit + 0.1*enorm2_fit;
        norm3_fit = norm3_fit + 0.1*enorm3_fit;
        norm4_fit = norm4_fit + 0.1*enorm4_fit;
        norm5_fit = norm5_fit + 0.1*enorm5_fit;
        norm6_fit = norm6_fit + 0.1*enorm6_fit;
        norm7_fit = norm7_fit + 0.1*enorm7_fit;
        norm8_fit = norm8_fit + 0.1*enorm8_fit;
        norm9_fit = norm9_fit + 0.1*enorm9_fit;
        norm10_fit = norm10_fit + 0.1*enorm10_fit;
        norm11_fit = norm11_fit + 0.1*enorm11_fit;
        norm12_fit = norm12_fit + 0.1*enorm12_fit;
        norm13_fit = norm13_fit + 0.1*enorm13_fit;
      }

      h_mc1_fit->Scale(norm1_fit);
      h_mc2_fit->Scale(norm2_fit);
      h_mc3_fit->Scale(norm3_fit);
      h_mc4_fit->Scale(norm4_fit);
//    h_mc5_fit->Scale(norm5_fit);
      h_mc6_fit->Scale(norm6_fit);
      h_mc7_fit->Scale(norm7_fit);
      h_mc8_fit->Scale(norm8_fit);
      h_mc9_fit->Scale(norm9_fit);
      h_mc10_fit->Scale(norm10_fit);
      h_mc11_fit->Scale(norm11_fit);
      h_mc12_fit->Scale(norm12_fit);
      h_mc13_fit->Scale(norm13_fit);

      if (title.find("_b")==string::npos) {
        if (irun==5) {
          h_mc2->Scale(c1_t+0.1*ec1_t);
          h_mc2_fit->Scale(c1_t+0.1*ec1_t);
        } else {
          h_mc2->Scale(c1_t);
          h_mc2_fit->Scale(c1_t);
        }
      } else {
        if (irun==5) {
          h_mc2->Scale(c2_t+0.1*ec2_t);
          h_mc2_fit->Scale(c2_t+0.1*ec2_t);
        } else {
          h_mc2->Scale(c2_t);
          h_mc2_fit->Scale(c2_t);
        }
      }

      if (irun==13) {
        for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
          h_mc1->SetBinError(i, 1.1*h_mc1->GetBinError(i));
          h_mc2->SetBinError(i, 1.1*h_mc2->GetBinError(i));
          h_mc3->SetBinError(i, 1.1*h_mc3->GetBinError(i));
          h_mc4->SetBinError(i, 1.1*h_mc4->GetBinError(i));
//        h_mc5->SetBinError(i, 1.1*h_mc5->GetBinError(i));
          h_mc6->SetBinError(i, 1.1*h_mc6->GetBinError(i));
          h_mc7->SetBinError(i, 1.1*h_mc7->GetBinError(i));
          h_mc8->SetBinError(i, 1.1*h_mc8->GetBinError(i));
          h_mc9->SetBinError(i, 1.1*h_mc9->GetBinError(i));
          h_mc10->SetBinError(i, 1.1*h_mc10->GetBinError(i));
          h_mc11->SetBinError(i, 1.1*h_mc11->GetBinError(i));
          h_mc12->SetBinError(i, 1.1*h_mc12->GetBinError(i));
          h_mc13->SetBinError(i, 1.1*h_mc13->GetBinError(i));
        }
      }

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

      h_mcD->Add(h_mc6);
      h_mcD->Add(h_mc4);
      h_mcD->Add(h_mc3);

      TH1F* h_mcO_fit = (TH1F*)h_mc8->Clone("h_mcO_fit");
      h_mcO_fit->Reset();

      TH1F* h_mcD_fit = (TH1F*)h_mc3->Clone("h_mcD_fit");
      h_mcD_fit->Reset();

      h_mcO_fit->Add(h_mc13_fit);
      h_mcO_fit->Add(h_mc12_fit);
      h_mcO_fit->Add(h_mc11_fit);
      h_mcO_fit->Add(h_mc10_fit);
      h_mcO_fit->Add(h_mc9_fit);
      h_mcO_fit->Add(h_mc8_fit);
      h_mcO_fit->Add(h_mc7_fit);

      h_mcD_fit->Add(h_mc6_fit);
      h_mcD_fit->Add(h_mc4_fit);
      h_mcD_fit->Add(h_mc3_fit);

      h_data->Add(h_mcO, -1.);
      h_data->Add(h_mcD, -1.);
      h_data->Add(h_mc1, -1.);

      h_data_fit->Add(h_mcO_fit, -1.);
      h_data_fit->Add(h_mcD_fit, -1.);
      h_data_fit->Add(h_mc1_fit, -1.);

      for (int i=0; i<=h_data->GetNbinsX()+1; i++) {
        if (h_data->GetBinContent(i) < 0) {
          h_data->SetBinContent(i, 0);
          h_data->SetBinError(i, 0);
        }
        if (h_data_fit->GetBinContent(i) < 0) {
          h_data_fit->SetBinContent(i, 0);
          h_data_fit->SetBinError(i, 0);
        }
      }

      TH1F* h_data_fit_raw = (TH1F*)h_data_fit->Clone();

      h_data_fit->Scale(Lumi2012/Lumi2012_ele_muon);

      TVirtualFitter::SetDefaultFitter("Minuit2");
      TVirtualFitter* fitter=0;
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
	fitter = TVirtualFitter::Fitter(0, 1);
	fitter->SetFCN(fcn);
	double arglist[1] = {-1.0};
	fitter->ExecuteCommand("SET PRINT", arglist, 1);
	fitter->SetParameter(0, "c(t)", 0.50, 0.01, 0.00, 1.00);
	fitter->ExecuteCommand("MIGRAD", arglist, 0);
	h_data_fit->Scale(fitter->GetParameter(0));
      }

      TCanvas* c1=0;
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
        h_data->SetStats(0);

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

	TH1F* h_ratio = (TH1F*)h_data->Clone("h_ratio");
	h_ratio->Divide(h_data_fit);

	h_ratio->SetTitle("");
        //h_ratio->SetStats(0);

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

        h_ratio->Draw("E0PX");

        TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
        OLine->SetLineColor(kGreen);
        OLine->SetLineWidth(1);
        OLine->Draw();

        c1->cd();

	TLegend *leg;
	leg = new TLegend(0.5, 0.8, 0.65, 0.9);

	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (ilepton==1) leg->AddEntry(h_data,"Z(#rightarrow ee)+jets","p");
	if (ilepton==2) leg->AddEntry(h_data,"Z(#rightarrow #mu#mu)+jets","p");

	leg->AddEntry(h_data_fit,"Z(#rightarrow e#mu)+jets","l");
	leg->AddEntry(h_mc2,"t#bar{t}","l");

	leg->Draw();

        TLatex *latexLabel = CMSPrel(Lumi2012/1000.,"",0.15,0.94);
        latexLabel->Draw("same");

        TLatex * lab = new TLatex ();
  	lab->SetTextSize (0.0275);
  	lab->SetTextFont (42);
  	lab->DrawLatex (0.63, 0.65, bSel.c_str());

        TLatex *fitLabel = new TLatex();
        fitLabel->SetTextSize(0.0275);
        fitLabel->SetTextFont(42);
        fitLabel->SetLineWidth(2);
        fitLabel->SetNDC();
        char buff[100];
        sprintf(buff, "c_{t} = %5.3f #pm %5.3f", fitter->GetParameter(0), fitter->GetParError(0));
        fitLabel->DrawLatex(0.68, 0.68, buff);
        if (ilepton==1) sprintf(buff, "I_{ee} = %5.1f", h_data->Integral(0,h_data->GetNbinsX()+1));
        if (ilepton==2) sprintf(buff, "I_{#mu#mu} = %5.1f", h_data->Integral(0,h_data->GetNbinsX()+1));
        fitLabel->DrawLatex(0.68, 0.63, buff);
        sprintf(buff, "I_{e#mu} = %5.1f #pm %5.1f", h_data_fit->Integral(0,h_data_fit->GetNbinsX()+1), h_data_fit->Integral(0,h_data_fit->GetNbinsX()+1)*fitter->GetParError(0)/fitter->GetParameter(0));
        fitLabel->DrawLatex(0.68, 0.58, buff);
        sprintf(buff, "I_{t#bar{t}}  = %5.1f", h_mc2->Integral(0,h_mc2->GetNbinsX()+1));
        fitLabel->DrawLatex(0.68, 0.53, buff);
      }

      if (plot) {
        if (doFit) title = title + "_doFit";
	ofstream out;
        if (ilepton==1) {
          gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/").c_str(), kTRUE);
          if (c1) c1->SaveAs((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + title + ".pdf").c_str());
	  if (doFit==0) {
	    TFile f((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + title + ".root").c_str(),"RECREATE");
            h_data_fit_raw->Write(title.c_str());
            f.Close();
	  }
	  if (doFit) out.open((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + title + ".dat").c_str());
        }
        if (ilepton==2) {
          gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/").c_str(), kTRUE);
          if (c1) c1->SaveAs((path + "/muons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + title + ".pdf").c_str());
	  if (doFit==0) {
	    TFile f((path + "/muons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + title + ".root").c_str(),"RECREATE");
            h_data_fit_raw->Write(title.c_str());
            f.Close();
	  }
	  if (doFit) out.open((path + "/muons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + title + ".dat").c_str());
        }
	if (doFit==1) {
	  out << fitter->GetParameter(0) << " " << fitter->GetParError(0) << endl;
	  out.close();
	}
      }
}
