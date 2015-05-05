#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v14.h"

#include "fixrange.C"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";
//string path = "/gpfs/cms/users/lalicata/work/test/data/";

void DataMCComp2(int irun=0, string title="", int plot=0, int ilepton=1, int isratio=1, int unfold=0, int numB=0, int bb=0) {

//int useBinnedEfficiency=0; // use average efficiencies
int useBinnedEfficiency=1; // use bin-by-bin efficiencies

//int useFitResults=0; // use MC predictions for c_b, c_c, c_uds, c_t
int useFitResults=1;  // use fit results for c_b, c_c, c_uds, c_t

//int useEleMuo = 0; // use MC or fit results for c_t
int useEleMuo = 1; // use e-mu fit results for c_t

int drawInclusive = 1; // do plot the "inclusive" histogram

string subdir="0";
string postfix="";
string dirbSel="";
string bSel="";
string genPostfix="";

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
if (irun==55) {            // irun==55 => template syst
  subdir="55";
  postfix="";
}
if (irun==66) {            // irun==66 => unfolding with data weight
  subdir="66";
  postfix="";
}
if (irun==77) {            // irun==77 => unfolding with MadGraph 4FS
  subdir="77";
  postfix="";
}
if (irun==88) {            // irun==88 => deltaR
  subdir="88";
  postfix="DR";
}
if (irun==99) {            // irun==99 => pur
  subdir="99";
  postfix="";
}
if (numB==1) {
  postfix = postfix + "1b";
  dirbSel="_1b";
  bSel="Z + (= 1) b-jet";
  drawInclusive = 0;
  genPostfix = "1b";
}
if (numB==2) {
  postfix = postfix + "2b";
  dirbSel="_2b";
  bSel="Z + (#geq 2) b-jet";
  drawInclusive = 0;
  genPostfix = "2b";
}

if (bb==1 && numB==1) bbBkg = true;
if (bb==1 && numB==2) bbSig = true;

	if (gROOT->GetVersionInt() >= 53401) {
	  //gROOT->GetColor(kRed)->SetAlpha(0.5);
	  gROOT->GetColor(kRed)->SetAlpha(0.0);
	  gROOT->GetColor(kGreen+2)->SetAlpha(0.5);
	  gROOT->GetColor(kMagenta-6)->SetAlpha(0.5);
	  gROOT->GetColor(kBlue-4)->SetAlpha(0.5);
	}

	/* efficiency */

	double e_Z=1.0;
	double ee_Z=0.0;
	double e_Zb=1.0;
	double ee_Zb=0.0;

	/* purity */

	double c_b=1.0;
	double ec_b=0.0;
	double c_c=1.0;
	double ec_c=0.0;
	double c_uds=1.0;
	double ec_uds=0.0;
        double fScal=1.0;
        double efScal=0.0;

	/* top background */

	double c1_t=1.0;
	double ec1_t=0.0;
	double c2_t=1.0;
	double ec2_t=0.0;

	ifstream in1, in2, in3, in4, in5, in6, in7, in8;
	if (ilepton==1) {
	  in1.open((path + "/electrons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + "w_first_jet_eta" + "_efficiency.dat").c_str());
	  in2.open((path + "/electrons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + "w_first_bjet_eta" + "_efficiency.dat").c_str());
          if (useFitResults) {
	    in3.open((path + "/electrons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_SVTX_mass_doFit" + ".dat").c_str());
	    in4.open((path + "/electrons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_doFit" + ".dat").c_str());
	    in5.open((path + "/electrons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_b_doFit" + ".dat").c_str());
	    in8.open((path + "/electrons/" + version + "/" + subdir + "/distributions_2b" + "/" + "w_SVTX_mass_doFit" + ".dat").c_str());
	    if (useEleMuo) {
	      in6.open((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + "w_mass_ee_wide_doFit" + ".dat").c_str());
	      in7.open((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + "w_mass_ee_b_wide_doFit" + ".dat").c_str());
	    }
	  }
	}
	if (ilepton==2) {
	  in1.open((path + "/muons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + "w_first_jet_eta" + "_efficiency.dat").c_str());
	  in2.open((path + "/muons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + "w_first_bjet_eta" + "_efficiency.dat").c_str());
          if (useFitResults) {
	    in3.open((path + "/muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_SVTX_mass_doFit" + ".dat").c_str());
	    in4.open((path + "/muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_doFit" + ".dat").c_str());
	    in5.open((path + "/muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_MET_b_doFit" + ".dat").c_str());
	    in8.open((path + "/muons/" + version + "/" + subdir + "/distributions_2b" + "/" + "w_SVTX_mass_doFit" + ".dat").c_str());
	    if (useEleMuo) {
	      in6.open((path + "/muons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + "w_mass_mm_wide_doFit" + ".dat").c_str());
	      in7.open((path + "/muons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + "w_mass_mm_b_wide_doFit" + ".dat").c_str());
	    }
	  }
	}
	in1 >> e_Z >> ee_Z;
	in2 >> e_Zb >> ee_Zb;
	in1.close();
	in2.close();
	if (useFitResults) {
	  if (numB!=2) {
            in3 >> c_uds >> ec_uds;
            in3 >> c_b >> ec_b;
            in3 >> c_c >> ec_c;
          }
          if (numB==2) {
            in3 >> c_b >> ec_b;
          }
	  in3.close();
	  in4 >> c1_t >> ec1_t;
	  in5 >> c2_t >> ec2_t;
	  in4.close();
	  in5.close();
	  in8 >> fScal >> efScal;
          in8.close();
	  if (useEleMuo) {
	    in6 >> c1_t >> ec1_t;
	    in7 >> c2_t >> ec2_t;
	    in6.close();
	    in7.close();
	  }
	}

	double Lumi2012=0;

	if (ilepton==1) Lumi2012 = Lumi2012_ele;
	if (ilepton==2) Lumi2012 = Lumi2012_muon;

	double norm1 = ((Lumi2012 * Xsec_dy) / Ngen_dy);
	double norm1_1 = ((Lumi2012 * Xsec_dy_1) / Ngen_dy_1);
	double norm1_2=0;
	if (ilepton==1) norm1_2 = ((Lumi2012 * Xsec_dy_2) / Ngen_dy_2_ee);
	if (ilepton==2) norm1_2 = ((Lumi2012 * Xsec_dy_2) / Ngen_dy_2_mm);
	double norm2 = ((Lumi2012 * Xsec_tt) / Ngen_tt);
	if (useEleMuo) norm2 = (Lumi2012 / Lumi2012_ele_muon);
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

	TFile *data=0;
	if (ilepton==1) data = TFile::Open((path + "/" + version + "/" + "DoubleElectron_2012_merge.root").c_str());
	if (ilepton==2) data = TFile::Open((path + "/" + version + "/" + "DoubleMu_2012_merge.root").c_str());

	TFile *mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL.root").c_str());
	TFile *mcg = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	TFile *mcg1 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());
	TFile *mcg2=0;
	if (ilepton==1) mcg2 = TFile::Open((path + "/" + version + "/" + "DYToEE_powheg_gen.root").c_str());
	if (ilepton==2) mcg2 = TFile::Open((path + "/" + version + "/" + "DYToMuMu_powheg_gen.root").c_str());
	TFile *mc2 = TFile::Open((path + "/" + version + "/" + "TTbar.root").c_str());
	TFile *mc3 = TFile::Open((path + "/" + version + "/" + "ZZ.root").c_str());
	TFile *mc4 = TFile::Open((path + "/" + version + "/" + "WZ.root").c_str());
//	TFile *mc5 = TFile::Open((path + "/" + version + "/" + "QCD.root").c_str());
	TFile *mc6 = TFile::Open((path + "/" + version + "/" + "WW.root").c_str());
	TFile *mc7 = TFile::Open((path + "/" + version + "/" + "Wj.root").c_str());
        TFile *mc8 = TFile::Open((path + "/" + version + "/" + "T_s.root").c_str());
        TFile *mc9 = TFile::Open((path + "/" + version + "/" + "T_t.root").c_str());
        TFile *mc10 = TFile::Open((path + "/" + version + "/" + "T_tW.root").c_str());
        TFile *mc11 = TFile::Open((path + "/" + version + "/" + "TBar_s.root").c_str());
        TFile *mc12 = TFile::Open((path + "/" + version + "/" + "TBar_t.root").c_str());
        TFile *mc13 = TFile::Open((path + "/" + version + "/" + "TBar_tW.root").c_str());

	string title_b = title;

        if (drawInclusive) {

	  if (title.find("_bjet_")!=string::npos) {
	    title.erase(title.find("_bjet_")+1, 1);
	  } else {
	    title_b = title + "_b";
          }

	  if (title.find("_single_")!=string::npos) {
	    if (title.find("_jet_")!=string::npos) {
	      title.replace(title.find("_single_"), 8, "_first_");
	    } else {
	      title.erase(title.find("_single_")+1, 7);
	    }
	  }
        }        

        if (!drawInclusive) {
          if (title.find("_abs")!=string::npos) {
            if (title.find("_bjet_")==string::npos) {
              title_b = title;
              title_b = title_b.replace(title_b.find("_abs"), 4, "_b_abs");
            }
          }
        }

	if (ilepton==1) data->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) data->cd(("demoMuo"+postfix).c_str());
	TH1F* h_data=0;
	TH1F* h_data_b=0;

	if (unfold==0) {
	  h_data = (TH1F*)gDirectory->Get(title.c_str());
	  h_data_b = (TH1F*)gDirectory->Get(title_b.c_str());
	}

	if (unfold==1) {
          if (ilepton==1) {
	    TFile f((path + "/electrons/" + version + "/" + subdir + "/unfolding" + dirbSel + "/" + title + "_unfolding.root").c_str());
	    TFile f_b((path + "/electrons/" + version + "/" + subdir + "/unfolding" + dirbSel + "/" + title_b + "_unfolding.root").c_str());
	    h_data = (TH1F*)f.Get(title.c_str())->Clone();
	    h_data_b = (TH1F*)f_b.Get(title_b.c_str())->Clone();
	    h_data->SetDirectory(0);
	    h_data_b->SetDirectory(0);
	    f.Close();
	    f_b.Close();
          }
          if (ilepton==2) {
	    TFile f((path + "/muons/" + version + "/" + subdir + "/unfolding" + dirbSel + "/" + title + "_unfolding.root").c_str());
	    TFile f_b((path + "/muons/" + version + "/" + subdir + "/unfolding" + dirbSel + "/" + title_b + "_unfolding.root").c_str());
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

	if (ilepton==1) mc1->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc1->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc1t = (TH1F*)gDirectory->Get(("t"+title.substr(1)).c_str());
	TH1F* h_mc1_b = (TH1F*)gDirectory->Get(title_b.c_str());
	TH1F* h_mc1b_b = (TH1F*)gDirectory->Get(("b"+title_b.substr(1)).c_str());
	TH1F* h_mc1c_b = (TH1F*)gDirectory->Get(("c"+title_b.substr(1)).c_str());
	TH1F* h_mc1t_b = (TH1F*)gDirectory->Get(("t"+title_b.substr(1)).c_str());
        TH1F* h_mc1bb = 0;
        if (bbSig) h_mc1bb = (TH1F*)gDirectory->Get(("bbSig"+title_b.substr(1)).c_str());
        if (bbBkg) h_mc1bb = (TH1F*)gDirectory->Get(("bbBkg"+title_b.substr(1)).c_str());

	if (ilepton==1) mcg->cd(("demoEleGen"+genPostfix).c_str());
	if (ilepton==2) mcg->cd(("demoMuoGen"+genPostfix).c_str());
	TH1F* h_mcg = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mcg_b = (TH1F*)gDirectory->Get(title_b.c_str());

        bool cdmcg1 = false;
      
	if (ilepton==1) cdmcg1 = mcg1->cd(("demoEleGen"+genPostfix).c_str());
	if (ilepton==2) cdmcg1 = mcg1->cd(("demoMuoGen"+genPostfix).c_str());
	TH1F* h_mcg1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mcg1_b = (TH1F*)gDirectory->Get(title_b.c_str());
        
        if (!h_mcg1 || !cdmcg1) h_mcg1 = (TH1F*)h_mcg->Clone();
	if (!h_mcg1_b || !cdmcg1) h_mcg1_b = (TH1F*)h_mcg_b->Clone();

	if (ilepton==1) mcg2->cd(("demoEleGen"+genPostfix).c_str());
	if (ilepton==2) mcg2->cd(("demoMuoGen"+genPostfix).c_str());
	TH1F* h_mcg2 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mcg2_b = (TH1F*)gDirectory->Get(title_b.c_str());
 
	if (!h_mcg2) h_mcg2 = (TH1F*)h_mcg->Clone();
	if (!h_mcg2_b) h_mcg2_b = (TH1F*)h_mcg_b->Clone();

	if (ilepton==1) mc2->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc2->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc2_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (useEleMuo) {
	  if (ilepton==1) {
	    mc2 = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + title + ".root").c_str());
	    h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	    mc2 = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + title_b + ".root").c_str());
	    h_mc2_b = (TH1F*)gDirectory->Get(title_b.c_str());
	  }
	  if (ilepton==2) {
	    mc2 = TFile::Open((path + "/muons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + title + ".root").c_str());
	    h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
	    mc2 = TFile::Open((path + "/muons/" + version + "/" + subdir + "/ttbar_sub" + dirbSel + "/" + title_b + ".root").c_str());
	    h_mc2_b = (TH1F*)gDirectory->Get(title_b.c_str());
	  }
	}

	if (ilepton==1) mc3->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc3->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc3_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc4->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc4->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc4_b = (TH1F*)gDirectory->Get(title_b.c_str());

//	if (ilepton==1) mc5->cd(("demoEle"+postfix).c_str());
//	if (ilepton==2) mc5->cd(("demoMuo"+postfix).c_str());
//	TH1F* h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
//	TH1F* h_mc5_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc6->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc6->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc6_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mc7->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc7->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc7_b = (TH1F*)gDirectory->Get(title_b.c_str());

        if (ilepton==1) mc8->cd(("demoEle"+postfix).c_str());
        if (ilepton==2) mc8->cd(("demoMuo"+postfix).c_str());
        TH1F* h_mc8 = (TH1F*)gDirectory->Get(title.c_str());
        TH1F* h_mc8_b = (TH1F*)gDirectory->Get(title_b.c_str());

        if (ilepton==1) mc9->cd(("demoEle"+postfix).c_str());
        if (ilepton==2) mc9->cd(("demoMuo"+postfix).c_str());
        TH1F* h_mc9 = (TH1F*)gDirectory->Get(title.c_str());
        TH1F* h_mc9_b = (TH1F*)gDirectory->Get(title_b.c_str());

        if (ilepton==1) mc10->cd(("demoEle"+postfix).c_str());
        if (ilepton==2) mc10->cd(("demoMuo"+postfix).c_str());
        TH1F* h_mc10 = (TH1F*)gDirectory->Get(title.c_str());
        TH1F* h_mc10_b = (TH1F*)gDirectory->Get(title_b.c_str());

        if (ilepton==1) mc11->cd(("demoEle"+postfix).c_str());
        if (ilepton==2) mc11->cd(("demoMuo"+postfix).c_str());
        TH1F* h_mc11 = (TH1F*)gDirectory->Get(title.c_str());
        TH1F* h_mc11_b = (TH1F*)gDirectory->Get(title_b.c_str());

        if (ilepton==1) mc12->cd(("demoEle"+postfix).c_str());
        if (ilepton==2) mc12->cd(("demoMuo"+postfix).c_str());
        TH1F* h_mc12 = (TH1F*)gDirectory->Get(title.c_str());
        TH1F* h_mc12_b = (TH1F*)gDirectory->Get(title_b.c_str());

        if (ilepton==1) mc13->cd(("demoEle"+postfix).c_str());
        if (ilepton==2) mc13->cd(("demoMuo"+postfix).c_str());
        TH1F* h_mc13 = (TH1F*)gDirectory->Get(title.c_str());
        TH1F* h_mc13_b = (TH1F*)gDirectory->Get(title_b.c_str());

        
        if (!drawInclusive) {
          h_data = (TH1F*)h_data->Clone();
          h_mc1 = (TH1F*)h_mc1->Clone();
          h_mc1t = (TH1F*)h_mc1t->Clone();
          h_mcg = (TH1F*)h_mcg->Clone();
          h_mcg1 = (TH1F*)h_mcg1->Clone();
          h_mcg2 = (TH1F*)h_mcg2->Clone();
          h_mc2 = (TH1F*)h_mc2->Clone();
          h_mc3 = (TH1F*)h_mc3->Clone();
          h_mc4 = (TH1F*)h_mc4->Clone();
//        h_mc5 = (TH1F*)h_mc5->Clone();
          h_mc6 = (TH1F*)h_mc6->Clone();
          h_mc7 = (TH1F*)h_mc7->Clone();
          h_mc8 = (TH1F*)h_mc8->Clone();
          h_mc9 = (TH1F*)h_mc9->Clone();
          h_mc10 = (TH1F*)h_mc10->Clone();
          h_mc11 = (TH1F*)h_mc11->Clone();
          h_mc12 = (TH1F*)h_mc12->Clone();
          h_mc13 = (TH1F*)h_mc13->Clone();
        }

	if (unfold==0) {
	  h_data->Sumw2();
	  h_data_b->Sumw2();
	}

	h_mc1->Sumw2();
	if (h_mc1t) h_mc1t->Sumw2();
	h_mcg->Sumw2();
	h_mcg1->Sumw2();
	h_mcg2->Sumw2();
	h_mc2->Sumw2();
	h_mc3->Sumw2();
	h_mc4->Sumw2();
//	h_mc5->Sumw2();
	h_mc6->Sumw2();
	h_mc7->Sumw2();
	h_mc8->Sumw2();
	h_mc9->Sumw2();
	h_mc10->Sumw2();
	h_mc11->Sumw2();
	h_mc12->Sumw2();
	h_mc13->Sumw2();

	h_mc1_b->Sumw2();
	if (h_mc1b_b) h_mc1b_b->Sumw2();
	if (h_mc1c_b) h_mc1c_b->Sumw2();
	if (h_mc1t_b) h_mc1t_b->Sumw2();
        if (bbBkg)  h_mc1bb->Sumw2();
        if (bbSig)  h_mc1bb->Sumw2();
	h_mcg_b->Sumw2();
	h_mcg1_b->Sumw2();
	h_mcg2_b->Sumw2();
	h_mc2_b->Sumw2();
	h_mc3_b->Sumw2();
	h_mc4_b->Sumw2();
//	h_mc5_b->Sumw2();
	h_mc6_b->Sumw2();
	h_mc7_b->Sumw2();
	h_mc8_b->Sumw2();
	h_mc9_b->Sumw2();
	h_mc10_b->Sumw2();
	h_mc11_b->Sumw2();
	h_mc12_b->Sumw2();
	h_mc13_b->Sumw2();

	if (irun==10) {
	  norm1 = norm1 + 0.1*enorm1;
	  norm1_1 = norm1_1 + 0.1*enorm1_1;
	  norm1_2 = norm1_2 + 0.1*enorm1_2;
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
	if (h_mc1t) h_mc1t->Scale(norm1);
	h_mcg->Scale(norm1);
	h_mcg1->Scale(norm1_1);
	h_mcg2->Scale(norm1_2);
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

	h_mc1_b->Scale(norm1);
	if (h_mc1b_b) h_mc1b_b->Scale(norm1);
	if (h_mc1c_b) h_mc1c_b->Scale(norm1);
	if (h_mc1t_b) h_mc1t_b->Scale(norm1);
        if (bbSig)  h_mc1bb->Scale(norm1);
        if (bbBkg)  h_mc1bb->Scale(norm1);
	h_mcg_b->Scale(norm1);
	h_mcg1_b->Scale(norm1_1);
	h_mcg2_b->Scale(norm1_2);
	h_mc2_b->Scale(norm2);
	h_mc3_b->Scale(norm3);
	h_mc4_b->Scale(norm4);
//	h_mc5_b->Scale(norm5);
	h_mc6_b->Scale(norm6);
	h_mc7_b->Scale(norm7);
	h_mc8_b->Scale(norm8);
	h_mc9_b->Scale(norm9);
	h_mc10_b->Scale(norm10);
	h_mc11_b->Scale(norm11);
	h_mc12_b->Scale(norm12);
	h_mc13_b->Scale(norm13);

        TH1F* h_mcO = (TH1F*)h_mc8->Clone("h_mcO");
        TH1F* h_mcO_b = (TH1F*)h_mc8_b->Clone("h_mcO_b");

        h_mcO->Add(h_mc13);
        h_mcO->Add(h_mc12);
        h_mcO->Add(h_mc11);
        h_mcO->Add(h_mc10);
        h_mcO->Add(h_mc9);
        h_mcO->Add(h_mc8);
        
        h_mcO_b->Add(h_mc13_b);
        h_mcO_b->Add(h_mc12_b);
        h_mcO_b->Add(h_mc11_b);
        h_mcO_b->Add(h_mc10_b);
        h_mcO_b->Add(h_mc9_b);
        h_mcO_b->Add(h_mc8_b);

	if (useFitResults) {
	  if (irun==5) {
	    h_mc2->Scale(c1_t+0.1*ec1_t);
	    h_mc2_b->Scale(c2_t+0.1*ec2_t);
	  } else {
	    h_mc2->Scale(c1_t);
	    h_mc2_b->Scale(c2_t);
	  }
	}

	if (irun==13) {
	  for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
	    h_mc1->SetBinError(i, 1.1*h_mc1->GetBinError(i));
	    if (h_mc1t) h_mc1t->SetBinError(i, 1.1*h_mc1t->GetBinError(i));
	    h_mc1_b->SetBinError(i, 1.1*h_mc1_b->GetBinError(i));
	    if (h_mc1b_b) h_mc1b_b->SetBinError(i, 1.1*h_mc1b_b->GetBinError(i));
	    if (h_mc1c_b) h_mc1c_b->SetBinError(i, 1.1*h_mc1c_b->GetBinError(i));
	    if (h_mc1t_b) h_mc1t_b->SetBinError(i, 1.1*h_mc1t_b->GetBinError(i));
            if (bbSig)  h_mc1bb->SetBinError(i, 1.1*h_mc1bb->GetBinError(i));
	    if (bbBkg)  h_mc1bb->SetBinError(i, 1.1*h_mc1bb->GetBinError(i));
	    h_mc2->SetBinError(i, 1.1*h_mc2->GetBinError(i));
	    h_mc2_b->SetBinError(i, 1.1*h_mc2_b->GetBinError(i));
	    h_mc3->SetBinError(i, 1.1*h_mc3->GetBinError(i));
	    h_mc3_b->SetBinError(i, 1.1*h_mc3_b->GetBinError(i));
	    h_mc4->SetBinError(i, 1.1*h_mc4->GetBinError(i));
	    h_mc4_b->SetBinError(i, 1.1*h_mc4_b->GetBinError(i));
//	    h_mc5->SetBinError(i, 1.1*h_mc5->GetBinError(i));
//	    h_mc5_b->SetBinError(i, 1.1*h_mc5_b->GetBinError(i));
	    h_mc6->SetBinError(i, 1.1*h_mc6->GetBinError(i));
	    h_mc6_b->SetBinError(i, 1.1*h_mc6_b->GetBinError(i));
	    h_mc7->SetBinError(i, 1.1*h_mc7->GetBinError(i));
	    h_mc7_b->SetBinError(i, 1.1*h_mc7_b->GetBinError(i));
            h_mcO->SetBinError(i, 1.1*h_mcO->GetBinError(i));
            h_mcO_b->SetBinError(i, 1.1*h_mcO_b->GetBinError(i));
	  }
	}

	if (unfold==0) {
	  h_data->Add(h_mcO, -1.);
	  h_data->Add(h_mc7, -1.);
	  h_data->Add(h_mc6, -1.);
//	  h_data->Add(h_mc5, -1.);
	  h_data->Add(h_mc4, -1.);
	  h_data->Add(h_mc3, -1.);
	  h_data->Add(h_mc2, -1.);
	  h_data->Add(h_mc1t, -1.);

          h_data_b->Add(h_mcO_b, -1.);
	  h_data_b->Add(h_mc7_b, -1.);
	  h_data_b->Add(h_mc6_b, -1.);
//	  h_data_b->Add(h_mc5_b, -1.);
	  h_data_b->Add(h_mc4_b, -1.);
	  h_data_b->Add(h_mc3_b, -1.);
	  h_data_b->Add(h_mc2_b, -1.);
	  h_data_b->Add(h_mc1t_b, -1.);
	}

	TH1F *h_mc1uds_b = (TH1F*)h_mc1_b->Clone("h_mc1uds_b");
	if (h_mc1b_b) h_mc1uds_b->Add(h_mc1b_b, -1);
	if (h_mc1c_b) h_mc1uds_b->Add(h_mc1c_b, -1);
	if (h_mc1t_b) h_mc1uds_b->Add(h_mc1t_b, -1);
	//if (bbBkg)  h_mc1uds_b->Add(h_mc1bb, -1);
        for (int i=0; i<=h_mc1uds_b->GetNbinsX()+1; i++) {
	  float e = TMath::Power(h_mc1uds_b->GetBinError(i),2);
	  if (h_mc1b_b) e = e - TMath::Power(h_mc1b_b->GetBinError(i),2);
	  if (h_mc1c_b) e = e - TMath::Power(h_mc1c_b->GetBinError(i),2);
	  if (h_mc1t_b) e = e - TMath::Power(h_mc1t_b->GetBinError(i),2);
	  h_mc1uds_b->SetBinError(i, TMath::Sqrt(e));
	}

	if (irun==99) {
	  float xval=0;
	  //float xvalb=0;
	  float xvalc=0;

	  if (h_mc1_b)  xval = h_mc1_b->Integral(0,h_mc1_b->GetNbinsX()+1);
	  //if (h_mc1b_b) xvalb = h_mc1b_b->Integral(0,h_mc1b_b->GetNbinsX()+1);
          if (h_mc1c_b) xvalc = h_mc1c_b->Integral(0,h_mc1c_b->GetNbinsX()+1);
	  
	  mc1 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());
	  //mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_aMC.root").c_str());
          if (ilepton==1) mc1->cd(("demoEle"+postfix).c_str());
          if (ilepton==2) mc1->cd(("demoMuo"+postfix).c_str());
          //if (ilepton==3) mc1->cd(("demoEleMuo"+postfix).c_str());
	   
	  h_mc1_b = (TH1F*)gDirectory->Get(title_b.c_str());
	  h_mc1_b->Sumw2();

	  TH1F* h_mc1b_b_tmp = h_mc1b_b;

	  h_mc1b_b = (TH1F*)gDirectory->Get(("b"+title_b.substr(1)).c_str());
	  h_mc1b_b->Sumw2(); 

	  h_mc1c_b = (TH1F*)gDirectory->Get(("c"+title_b.substr(1)).c_str());
          h_mc1c_b->Sumw2(); 

          if (bbBkg) {
            h_mc1bb = (TH1F*)gDirectory->Get(("bbBkg"+title_b.substr(1)).c_str());
            h_mc1bb->Sumw2();
            h_mc1b_b->Add(h_mc1bb, -1.);
          }

	  if (h_mc1b_b) h_mc1_b->Add(h_mc1b_b, -1.);
	  if (h_mc1c_b) h_mc1_b->Add(h_mc1c_b, -1.);
	  if (h_mc1t_b) h_mc1_b->Add(h_mc1t_b, -1.);
//        if (bbBkg || bbSig) h_mc1->Add(h_mc1bb, -1.);
          if (bbBkg) h_mc1->Add(h_mc1bb, -1.);
	  for (int i=0; i<=h_mc1->GetNbinsX()+1; i++) {
	    float e = TMath::Power(h_mc1->GetBinError(i),2);
	    if (h_mc1b_b) e = e - TMath::Power(h_mc1b_b->GetBinError(i),2);
//          if (bbBkg || bbSig) e = e - TMath::Power(h_mc1bb->GetBinError(i),2);
            if (bbBkg) e = e - TMath::Power(h_mc1bb->GetBinError(i),2);
	    if (h_mc1c_b) e = e - TMath::Power(h_mc1c_b->GetBinError(i),2);
	    if (h_mc1t_b) e = e - TMath::Power(h_mc1t_b->GetBinError(i),2);
	    h_mc1_b->SetBinError(i, TMath::Sqrt(e));
	  }

	  xval = xval / h_mc1_b->Integral(0,h_mc1_b->GetNbinsX()+1);
	  h_mc1_b->Scale(xval);

	  //xvalb = xvalb / h_mc1b->Integral(0,h_mc1b->GetNbinsX()+1);
	  //h_mc1b->Scale(xvalb);

	  h_mc1b_b = h_mc1b_b_tmp;

	  xvalc = xvalc / h_mc1c_b->Integral(0,h_mc1c_b->GetNbinsX()+1);
	  h_mc1c_b->Scale(xvalc);

	  h_mc1_b -> SetLineColor(kBlack);

	}

        if (bbBkg) {
          h_mc1b_b->Add(h_mc1bb, -1);
	  if (irun==6) {
            h_mc1bb->Scale(fScal+0.1*efScal);
	  } else {
            h_mc1bb->Scale(fScal);
	  }
        }

	if (h_mc1uds_b) {
	  if (irun==6) {
	    h_mc1uds_b->Scale(c_uds+0.1*ec_uds);
	  } else {
	    h_mc1uds_b->Scale(c_uds);
	  }
	}
	if (h_mc1b_b) {
	  if (irun==6) {
	    h_mc1b_b->Scale(c_b+0.1*ec_b);
	  } else {
	    h_mc1b_b->Scale(c_b);
	  }
	}
	if (h_mc1c_b) {
	  if (irun==6) {
	    h_mc1c_b->Scale(c_c+0.1*ec_c);
	  } else {
	    h_mc1c_b->Scale(c_c);
	  }
        }

	if (unfold==0) {
	  h_data_b->Add(h_mc1c_b, -1.);
	  h_data_b->Add(h_mc1uds_b, -1.);
          if (bbBkg) h_data_b->Add(h_mc1bb, -1.);
          if (bbSig) {
            h_data_b->Add(h_mc1b_b, -1.);
            if (irun==6) {
              h_data_b->Add(h_mc1bb, 1.);
            } else {
              h_data_b->Add(h_mc1bb, 1.);
            }
          }
	}

	TH1F *h_data_raw=0;
	TH1F *h_data_b_raw=0;
        TH1F *h_data_raw2=0;
        TH1F *h_data_b_raw2=0;
	if (unfold==0) {
	  h_data_raw = (TH1F*)h_data->Clone();
	  h_data_b_raw = (TH1F*)h_data_b->Clone();
          h_data_raw2 = (TH1F*)h_data->Clone();
          h_data_b_raw2 = (TH1F*)h_data_b->Clone();
	}

	if (useBinnedEfficiency==0) {
	  if (unfold==0) {
	    h_data->Scale(1./e_Z);
	    h_data_b->Scale(1./e_Zb);
	  }
	  h_mc1->Scale(1./e_Z);
	  h_mc1b_b->Scale(1./e_Zb);
	}

	if (useBinnedEfficiency==1) {
          if (ilepton==1) {
            TFile f((path + "/electrons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + string(h_data->GetName()) + "_efficiency.root").c_str());
            TFile f_b((path + "/electrons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + string(h_data_b->GetName()) + "_efficiency.root").c_str());
	    TH1F* h = (TH1F*)f.Get(h_data->GetName())->Clone();
	    TH1F* h_b = (TH1F*)f_b.Get(h_data_b->GetName())->Clone();
	    h->SetDirectory(0);
	    h_b->SetDirectory(0);
	    f.Close();
	    f_b.Close();
	    if (unfold==0) {
	      h_data->Divide(h);
	      h_data_b->Divide(h_b);
              h_data_raw2->Divide(h);
              h_data_b_raw2->Divide(h_b);
	    }
	    h_mc1->Divide(h);
	    h_mc1b_b->Divide(h_b);
          }
	  if (ilepton==2) {
	    TFile f((path + "/muons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + string(h_data->GetName()) + "_efficiency.root").c_str());
	    TFile f_b((path + "/muons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + string(h_data_b->GetName()) + "_efficiency.root").c_str());
            TH1F* h = (TH1F*)f.Get(h_data->GetName())->Clone();
	    TH1F* h_b = (TH1F*)f_b.Get(h_data_b->GetName())->Clone();
	    h->SetDirectory(0);
	    h_b->SetDirectory(0);
	    f.Close();
	    f_b.Close();
	    if (unfold==0) {
	      h_data->Divide(h);
	      h_data_b->Divide(h_b);
              h_data_raw2->Divide(h);
              h_data_b_raw2->Divide(h_b);
	    }
	    h_mc1->Divide(h);
	    h_mc1b_b->Divide(h_b);
          }
	}

        if (unfold==0) {
          h_data_raw2->Scale(1./Lumi2012, "width");
          h_data_b_raw2->Scale(1./Lumi2012, "width");
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

        if (unfold==0) { 
          h_data_raw2 = fixrange(h_data_raw2);
          h_data_b_raw2 = fixrange(h_data_b_raw2);
        }       
	h_data = fixrange(h_data);
	h_data_b = fixrange(h_data_b);
	h_mc1 = fixrange(h_mc1);
 	h_mc1b_b = fixrange(h_mc1b_b);
        if (bbBkg) h_mc1bb = fixrange(h_mc1bb);
        if (bbSig) h_mc1bb = fixrange(h_mc1bb);
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
	  TH1F* tmp1 = (TH1F*)h_mc1b_b->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp1->GetMinimum()==0) tmp1->GetXaxis()->SetRangeUser(0, tmp1->GetBinCenter(tmp1->GetMinimumBin()-1));
	  }
	  tmp1->SetFillColor(0);
	  tmp1->DrawClone("HISTLSAME");

	  h_mcg_b->Draw("E5SAME");
	  TH1F* tmp2 = (TH1F*)h_mcg_b->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp2->GetMinimum()==0) tmp2->GetXaxis()->SetRangeUser(0, tmp2->GetBinCenter(tmp2->GetMinimumBin()-1));
	  }
	  tmp2->SetFillColor(0);
	  tmp2->DrawClone("HISTLSAME");

          h_mcg1_b->Draw("E5SAME");
          TH1F* tmp2_1 = (TH1F*)h_mcg1_b->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp2_1->GetMinimum()==0) tmp2_1->GetXaxis()->SetRangeUser(0, tmp2_1->GetBinCenter(tmp2_1->GetMinimumBin()-1));
	  }
	  tmp2_1->SetFillColor(0);
	  tmp2_1->DrawClone("HISTLSAME");

	  h_mcg2_b->Draw("E5SAME");
	  TH1F* tmp2_2 = (TH1F*)h_mcg2_b->Clone();
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
	  if (drawInclusive) h_mc1->Draw("E5SAME");
	  TH1F* tmp3 = (TH1F*)h_mc1->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp3->GetMinimum()==0) tmp3->GetXaxis()->SetRangeUser(0, tmp3->GetBinCenter(tmp3->GetMinimumBin()-1));
	  }
	  tmp3->SetFillColor(0);
	  if (drawInclusive) tmp3->DrawClone("HISTLSAME");

	  h_mcg->SetLineColor(kGreen+2);
	  h_mcg->SetLineWidth(2);
	  h_mcg->SetMarkerColor(kGreen+2);
	  h_mcg->SetFillColor(kGreen+2);
	  if (drawInclusive) h_mcg->Draw("E5SAME");
	  TH1F* tmp4 = (TH1F*)h_mcg->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp4->GetMinimum()==0) tmp4->GetXaxis()->SetRangeUser(0, tmp4->GetBinCenter(tmp4->GetMinimumBin()-1));
	  }
	  tmp4->SetFillColor(0);
	  if (drawInclusive) tmp4->DrawClone("HISTLSAME");

	  h_mcg1->SetLineColor(kMagenta-6);
	  h_mcg1->SetLineWidth(2);
	  h_mcg1->SetMarkerColor(kMagenta-6);
	  h_mcg1->SetFillColor(kMagenta-6);
	  if (drawInclusive) h_mcg1->Draw("E5SAME");
	  TH1F* tmp4_1 = (TH1F*)h_mcg1->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp4_1->GetMinimum()==0) tmp4_1->GetXaxis()->SetRangeUser(0, tmp4_1->GetBinCenter(tmp4_1->GetMinimumBin()-1));
	  }
	  tmp4_1->SetFillColor(0);
	  if (drawInclusive) tmp4_1->DrawClone("HISTLSAME");

	  h_mcg2->SetLineColor(kBlue-4);
	  h_mcg2->SetLineWidth(2);
	  h_mcg2->SetMarkerColor(kBlue-4);
	  h_mcg2->SetFillColor(kBlue-4);
	  if (drawInclusive) h_mcg2->Draw("E5SAME");
	  TH1F* tmp4_2 = (TH1F*)h_mcg2->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp4_2->GetMinimum()==0) tmp4_2->GetXaxis()->SetRangeUser(0, tmp4_2->GetBinCenter(tmp4_2->GetMinimumBin()-1));
	  }
	  tmp4_2->SetFillColor(0);
	  if (drawInclusive) tmp4_2->DrawClone("HISTLSAME");

	  h_data->SetMarkerColor(kBlack);
	  h_data->SetLineColor(kBlack);
	  h_data->SetMarkerStyle(20);
	  h_data->SetMarkerSize (0.7);
	  if (drawInclusive) h_data->Draw("EPX0SAME");

	  if (ilepton==1) {
	    if (drawInclusive) leg->AddEntry(h_data,"Z(#rightarrow ee) DATA","p");
	    leg->AddEntry(h_data_b,"Z(#rightarrow ee)+b DATA","p");
	    //leg->AddEntry(h_mc1,"Z(#rightarrow ee) MC","l");
	    leg->AddEntry(h_mcg,"Z(#rightarrow ee) MadGraph","l");
	    leg->AddEntry(h_mcg1,"Z(#rightarrow ee) Sherpa","l");
	    leg->AddEntry(h_mcg2,"Z(#rightarrow ee) Powheg","l");
	  }
	  if (ilepton==2){
	    if (drawInclusive) leg->AddEntry(h_data,"Z(#rightarrow #mu#mu) DATA","p");
	    leg->AddEntry(h_data_b,"Z(#rightarrow #mu#mu)+b DATA","p");
	    //leg->AddEntry(h_mc1,"Z(#rightarrow #mu#mu) MC","l");
	    leg->AddEntry(h_mcg,"Z(#rightarrow #mu#mu) MadGraph","l");
	    leg->AddEntry(h_mcg1,"Z(#rightarrow #mu#mu) Sherpa","l");
	    leg->AddEntry(h_mcg2,"Z(#rightarrow #mu#mu) Powheg","l");
	  }
	}

	if (isratio==1) {
	  if (ilepton==1) {
	    leg->AddEntry(h_data_b,"Z(#rightarrow ee) DATA","p");
	    //leg->AddEntry(h_mc1b_b,"Z(#rightarrow ee) MC","l");
	    leg->AddEntry(h_mcg_b,"Z(#rightarrow ee) MadGraph","l");
	    leg->AddEntry(h_mcg1_b,"Z(#rightarrow ee) Sherpa","l");
	    leg->AddEntry(h_mcg2_b,"Z(#rightarrow ee) Powheg","l");
	  }
	  if (ilepton==2){
	    leg->AddEntry(h_data_b,"Z(#rightarrow #mu#mu) DATA","p");
	    //leg->AddEntry(h_mc1b_b,"Z(#rightarrow #mu#mu) MC","l");
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

	TH1F *h_M = (TH1F*)h_data_b->Clone();
	h_M->Divide(h_mcg_b);
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

	h_M->SetMarkerStyle(24);
	h_M->Draw("EPX0");

	if (isratio==0) {
	  TH1F *h_M2= (TH1F*)h_data->Clone();
	  h_M2->Divide(h_mcg);
	  TGraphErrors *g_M2 = new TGraphErrors(h_M2);

	  float dx = 0.1*(g_M2->GetXaxis()->GetXmax()-g_M2->GetXaxis()->GetXmin())/g_M2->GetN();
	  for (int i=0; i<g_M2->GetN(); i++) {
	    g_M2->SetPoint(i, g_M2->GetX()[i]-dx, g_M2->GetY()[i]);
	    g_M2->SetPointError(i, 0, g_M2->GetEY()[i]);
	  }

	  g_M2->SetMarkerStyle(20);
	  if (drawInclusive) g_M2->Draw("EP0SAME");
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

	TH1F *h_S = (TH1F*)h_data_b->Clone();
	h_S->Divide(h_mcg1_b);

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

	h_S->SetMarkerStyle(24);
	h_S->Draw("EPX0");

	if (isratio==0) {
	  TH1F *h_S2= (TH1F*)h_data->Clone();
	  h_S2->Divide(h_mcg1);

	  TGraphErrors *g_S2 = new TGraphErrors(h_S2);

	  float dx = 0.1*(g_S2->GetXaxis()->GetXmax()-g_S2->GetXaxis()->GetXmin())/g_S2->GetN();
	  for (int i=0; i<g_S2->GetN(); i++) {
	    g_S2->SetPoint(i, g_S2->GetX()[i]-dx, g_S2->GetY()[i]);
	    g_S2->SetPointError(i, 0, g_S2->GetEY()[i]);
	  }

	  g_S2->SetMarkerStyle(20);
	  if (drawInclusive) g_S2->Draw("EP0SAME");
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

	TH1F *h_P = (TH1F*)h_data_b->Clone();
	h_P->Divide(h_mcg2_b);

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

	h_P->SetMarkerStyle(24);
	h_P->Draw("EPX0");

	if (isratio==0) {
	  TH1F *h_P2= (TH1F*)h_data->Clone();
	  h_P2->Divide(h_mcg2);

	  TGraphErrors *g_P2 = new TGraphErrors(h_P2);

	  float dx = 0.1*(g_P2->GetXaxis()->GetXmax()-g_P2->GetXaxis()->GetXmin())/g_P2->GetN();
	  for (int i=0; i<g_P2->GetN(); i++) {
	    g_P2->SetPoint(i, g_P2->GetX()[i]-dx, g_P2->GetY()[i]);
	    g_P2->SetPointError(i, 0, g_P2->GetEY()[i]);
	  }

	  g_P2->SetMarkerStyle(20);
	  if (drawInclusive) g_P2->Draw("EP0SAME");
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
	} else if (title_b=="w_mass_Zj_ee_b" || title_b=="w_mass_Zj_mm_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dM_{Zj} [pb]");
	  h_P->GetXaxis()->SetTitle("M(Zj) [rad]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dM_{Zj} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	}

	if (plot) {
	  if (unfold==0) {
	    if (isratio==0) {
	      ofstream out;
	      if (ilepton==1) {
	        gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/xsecs" + dirbSel + "/").c_str(), kTRUE);
	        c1->SaveAs((path + "/electrons/" + version + "/" + subdir + "/xsecs" + dirbSel + "/" + title_b + "_xsecs.pdf").c_str());
	        out.open((path + "/electrons/" + version + "/" + subdir + "/xsecs" + dirbSel + "/" + title + "_xsecs.dat").c_str());
	        TFile f((path + "/electrons/" + version + "/" + subdir + "/xsecs" + dirbSel + "/" + title_b + "_xsecs.root").c_str(),"RECREATE");
	        h_data_raw->Write((title+"_raw").c_str());
                h_data_b_raw->Write((title_b+"_raw").c_str());
//                h_data_raw2->Write((title+"_raw2").c_str());
//                h_data_b_raw2->Write((title_b+"_raw2").c_str());
                h_data->Write(title.c_str());
                h_data_b->Write(title_b.c_str());
                h_mc1->Write((title+"_MC").c_str());
                h_mc1b_b->Write((title_b+"_MC").c_str());
                f.Close();
	      }
	      if (ilepton==2) {
	        gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/xsecs" + dirbSel + "/").c_str(), kTRUE);
	        c1->SaveAs((path + "/muons/" + version + "/" + subdir + "/xsecs" + dirbSel + "/" + title_b + "_xsecs.pdf").c_str());
	        out.open((path + "/muons/" + version + "/" + subdir + "/xsecs" + dirbSel + "/" + title + "_xsecs.dat").c_str());
	        TFile f((path + "/muons/" + version + "/" + subdir + "/xsecs" + dirbSel + "/" + title_b + "_xsecs.root").c_str(),"RECREATE");
	        h_data_raw->Write((title+"_raw").c_str());
                h_data_b_raw->Write((title_b+"_raw").c_str());
//                h_data_raw2->Write((title+"_raw2").c_str());
//                h_data_b_raw2->Write((title_b+"_raw2").c_str());
                h_data->Write(title.c_str());
                h_data_b->Write(title_b.c_str());
                h_mc1->Write((title+"_MC").c_str());
                h_mc1b_b->Write((title_b+"_MC").c_str());
                f.Close();
	      }
	      out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	      out << h_data->Integral(0, h_data->GetNbinsX()+1, "width") << endl;
	      out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	      out << h_data_b->Integral(0, h_data_b->GetNbinsX()+1, "width") << endl;
	      out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	      out << h_mc1->Integral(0, h_mc1->GetNbinsX()+1, "width") << endl;
	      out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	      out << h_mc1b_b->Integral(0, h_mc1b_b->GetNbinsX()+1, "width") << endl;
	      out.close();
	    }
	    if (isratio==1) {
	      if (ilepton==1) {
	        gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/ratios" + dirbSel + "/").c_str(), kTRUE);
	        c1->SaveAs((path + "/electrons/" + version + "/" + subdir + "/ratios" + dirbSel + "/" + title_b + "_ratio.pdf").c_str());
	      }
	      if (ilepton==2) {
	        gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/ratios" + dirbSel + "/").c_str(), kTRUE);
	        c1->SaveAs((path + "/muons/" + version + "/" + subdir + "/ratios" + dirbSel + "/" + title_b + "_ratio.pdf").c_str());
	      }
	    }
	  }
	  if (unfold==1) {
	    if (isratio==0) {
	      ofstream out;
	      if (ilepton==1) {
	        gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/").c_str(), kTRUE);
	        c1->SaveAs((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.pdf").c_str());
	        out.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + title + "_xsecs_unfolding.dat").c_str());
//                TFile f((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.root").c_str(),"RECREATE");
//                h_data->Write(title.c_str());
//                h_data_b->Write(title_b.c_str());
//                h_mc1->Write((title+"_MC").c_str());
//                h_mc1b_b->Write((title_b+"_MC").c_str());
//                f.Close();
	      }
	      if (ilepton==2) {
	        gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/").c_str(), kTRUE);
	        c1->SaveAs((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.pdf").c_str());
	        out.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + title + "_xsecs_unfolding.dat").c_str());
//                TFile f((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.root").c_str(),"RECREATE");
//                h_data->Write(title.c_str());
//                h_data_b->Write(title_b.c_str());
//                h_mc1->Write((title+"_MC").c_str());
//                h_mc1b_b->Write((title_b+"_MC").c_str());
//                f.Close();
	      }
	      out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	      out << h_data->Integral(0, h_data->GetNbinsX()+1, "width") << endl;
	      out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	      out << h_data_b->Integral(0, h_data_b->GetNbinsX()+1, "width") << endl;
	      out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	      out << h_mc1->Integral(0, h_mc1->GetNbinsX()+1, "width") << endl;
	      out << std::fixed << std::setw( 11 ) << std::setprecision( 4 );
	      out << h_mc1b_b->Integral(0, h_mc1b_b->GetNbinsX()+1, "width") << endl;
	      out.close();
	    }
	    if (isratio==1) {
	      if (ilepton==1) {
	        gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/ratios_unfolding" + dirbSel + "/").c_str(), kTRUE);
	        c1->SaveAs((path + "/electrons/" + version + "/" + subdir + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.pdf").c_str());
	      }
	      if (ilepton==2) {
	        gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/ratios_unfolding" + dirbSel + "/").c_str(), kTRUE);
	        c1->SaveAs((path + "/muons/" + version + "/" + subdir + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.pdf").c_str());
	      }
	    }
	  }
	}
}

