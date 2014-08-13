#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v13.h"

#include "fixrange.C"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";
//string path = "/gpfs/cms/users/lalicata/work/test/data/";

TH1F* read(string subdir, string title, int ilepton, TFile* infile=0, string dirbSel="") {
  TH1F* hist;
  TFile* file = infile;
  string title_tmp = title;
  string postfix = "";
  if (dirbSel=="_1b") postfix = "1b";
  if (dirbSel=="_2b") postfix = "2b";
  if (ilepton==1) {
    if (title=="w_pt_Z") title_tmp="w_pt_Z_ee";
    if (title=="w_pt_Z_b") title_tmp="w_pt_Z_ee_b";
    if (title=="w_delta_phi") title_tmp="w_delta_phi_ee";
    if (title=="w_delta_phi_b") title_tmp="w_delta_phi_ee_b";
    if (title=="w_mass_Zj") title_tmp="w_mass_Zj_ee";
    if (title=="w_mass_Zj_b") title_tmp="w_mass_Zj_ee_b";
    if (title=="w_Zbb_mass") title_tmp="w_eebb_mass";
    if (title=="w_DR_Zb_min") title_tmp="w_DR_eeb_min";
    if (title=="w_DR_Zb_max") title_tmp="w_DR_eeb_max";
    if (title=="w_A_Zb") title_tmp="w_A_eeb";
    if (file) {
      file->cd(("demoEleGen"+postfix).c_str());
    } else {
      file = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/unfolding" + dirbSel + "/" + title_tmp + "_unfolding.root").c_str());
    }
  }
  if (ilepton==2) {
    if (title=="w_pt_Z") title_tmp="w_pt_Z_mm";
    if (title=="w_pt_Z_b") title_tmp="w_pt_Z_mm_b";
    if (title=="w_delta_phi") title_tmp="w_delta_phi_mm";
    if (title=="w_delta_phi_b") title_tmp="w_delta_phi_mm_b";
    if (title=="w_mass_Zj") title_tmp="w_mass_Zj_mm";
    if (title=="w_mass_Zj_b") title_tmp="w_mass_Zj_mm_b";
    if (title=="w_Zbb_mass") title_tmp="w_mmbb_mass";
    if (title=="w_DR_Zb_min") title_tmp="w_DR_mmb_min";
    if (title=="w_DR_Zb_max") title_tmp="w_DR_mmb_max";
    if (title=="w_A_Zb") title_tmp="w_A_mmb";

    if (file) {
      file->cd(("demoMuoGen"+postfix).c_str());
    } else {
      file = TFile::Open((path + "/muons/" + version + "/" + subdir + "/unfolding" + dirbSel + "/" + title_tmp + "_unfolding.root").c_str());
    }
  }
  hist = (TH1F*)gDirectory->Get(title_tmp.c_str())->Clone();
  hist->SetDirectory(0);
  if (!infile) file->Close();
  return hist;
}

double calc(int iflag, double cont1, double cont2, double stat1, double stat2, double stat_bkg1, double stat_bkg2, double syst_eff1, double syst_eff2, double syst_jec1, double syst_jec2, double syst_jer1, double syst_jer2, double syst_pu1, double syst_pu2, double syst_bkg1, double syst_bkg2, double stat_top1, double stat_top2, double stat_bfit1, double stat_bfit2, double syst_btag1, double syst_btag2, double stat_unfold1, double stat_unfold2, double syst_unfold1, double syst_unfold2, double syst_lumi1, double syst_lumi2) {
  double val = 0.0;

  if (iflag == 0) {
    if (cont1*cont2 != 0) {
      val = (cont1/(TMath::Power(stat1,2)+TMath::Power(stat_bkg1,2)+TMath::Power(syst_eff1,2)+TMath::Power(stat_top1,2)+TMath::Power(stat_bfit1,2)+TMath::Power(stat_unfold1,2))+cont2/(TMath::Power(stat2,2)+TMath::Power(stat_bkg2,2)+TMath::Power(syst_eff2,2)+TMath::Power(stat_top2,2)+TMath::Power(stat_bfit2,2)+TMath::Power(stat_unfold2,2)))/(1./(TMath::Power(stat1,2)+TMath::Power(stat_bkg1,2)+TMath::Power(syst_eff1,2)+TMath::Power(stat_top1,2)+TMath::Power(stat_bfit1,2)+TMath::Power(stat_unfold1,2))+1./(TMath::Power(stat2,2)+TMath::Power(stat_bkg2,2)+TMath::Power(syst_eff2,2)+TMath::Power(stat_top2,2)+TMath::Power(stat_bfit2,2)+TMath::Power(stat_unfold2,2)));
    }
  }

  if (iflag == 1) {
    if (cont1*cont2 != 0) {
      val = TMath::Sqrt(TMath::Power((syst_jec1+syst_jec2)/2.,2)+TMath::Power((syst_jer1+syst_jer2)/2.,2)+TMath::Power((syst_pu1+syst_pu2)/2.,2)+TMath::Power((syst_bkg1+syst_bkg2)/2.,2)+TMath::Power((syst_btag1+syst_btag2)/2.,2)+TMath::Power((syst_unfold1+syst_unfold2)/2.,2)+TMath::Power((syst_lumi1+syst_lumi2)/2.,2));
      val = TMath::Sqrt(TMath::Power(val,2)+1./(1./(TMath::Power(stat1,2)+TMath::Power(stat_bkg1,2)+TMath::Power(syst_eff1,2)+TMath::Power(stat_top1,2)+TMath::Power(stat_bfit1,2)+TMath::Power(stat_unfold1,2))+1./(TMath::Power(stat2,2)+TMath::Power(stat_bkg2,2)+TMath::Power(syst_eff2,2)+TMath::Power(stat_top2,2)+TMath::Power(stat_bfit2,2)+TMath::Power(stat_unfold2,2))));
    }
  }

  return val;
}

void DataMCComp8(string title="", int plot=0, int isratio=1, int numB=0) {

int useSherpa=0;
//int useSherpa=1; // use Sherpa MC prediction

//int useNewPowheg=0;
int useNewPowheg=1; // use new Powheg MC prediction

int drawInclusive = 1; // do plot the "inclusive" histogram

string subdir="0";
string postfix="";
string dirbSel="";
string bSel="";

if (numB==1) {
  postfix = postfix + "1b";
  dirbSel="_1b";
  bSel="Z + (= 1) b-jet";
  drawInclusive = 0;
}
if (numB==2) {
  postfix = postfix + "2b";
  dirbSel="_2b";
  bSel="Z + (#geq 2) b-jet";
  drawInclusive = 0;
}

	if (gROOT->GetVersionInt() >= 53401) {
	  gROOT->GetColor(kRed)->SetAlpha(0.5);
	  gROOT->GetColor(kGreen+2)->SetAlpha(0.5);
	  gROOT->GetColor(kMagenta-6)->SetAlpha(0.5);
	  if (!useSherpa) gROOT->GetColor(kMagenta-6)->SetAlpha(0.0);
	  gROOT->GetColor(kBlue-4)->SetAlpha(0.5);
	  gROOT->GetColor(kOrange+7)->SetAlpha(0.5);
	}

	double Lumi2012=0;

	Lumi2012 = (Lumi2012_ele+Lumi2012_muon)/2.;

	double norm1 = ((Lumi2012 * Xsec_dy) / Ngen_dy);
	double norm1_1 = ((Lumi2012 * Xsec_dy_1) / Ngen_dy_1);
	double norm1_2 = ((Lumi2012 * Xsec_dy_2) / ((Ngen_dy_2_ee+Ngen_dy_2_mm)/2.));
	double norm1_3 = ((Lumi2012 * Xsec_dy_3) / Ngen_dy_3);

	if (title.empty()) title = "w_jetmultiplicity";

	TFile *mcg = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	TFile *mcg1 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());
	TFile *mcg2[2];
	mcg2[0] = TFile::Open((path + "/" + version + "/" + "DYToEE_powheg_gen.root").c_str());
	mcg2[1] = TFile::Open((path + "/" + version + "/" + "DYToMuMu_powheg_gen.root").c_str());
	if (useNewPowheg) {
          Ngen_dy_2_ee = 100*10000;
          Ngen_dy_2_mm = 100*10000;
	  Xsec_dy_2 = 333.866;
          norm1_2 = ((Lumi2012 * Xsec_dy_2) / ((Ngen_dy_2_ee+Ngen_dy_2_mm)/2.));
	  mcg2[0] = TFile::Open(("/gpfs/cms/users/candelis/work/ZbSkim/powheg/data/" + version + "/" + "powheg_ele.root").c_str());
	  mcg2[1] = TFile::Open(("/gpfs/cms/users/candelis/work/ZbSkim/powheg/data/" + version + "/" + "powheg_muo.root").c_str());
	}
	TFile *mcg3 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL2_gen.root").c_str());

	string title_b = title;
        
        if (numB==0) {
	  if (title.find("_bjet_")!=string::npos) {
	    title.erase(title.find("_bjet_")+1, 1);
	  } else {
	    title_b = title + "_b";
	  }
        }

	TH1F* w_data[2];
	TH1F* w_data_b[2];
	for (int i=0; i<2; i++) {
	  w_data[i] = read(subdir, title, i+1, 0, dirbSel);
	  w_data_b[i] = read(subdir, title_b, i+1, 0, dirbSel);
	}

	TH1F* w_mcg[2];
	TH1F* w_mcg_b[2];
	TH1F* w_mcg1[2];
	TH1F* w_mcg1_b[2];
	TH1F* w_mcg2[2];
	TH1F* w_mcg2_b[2];
	TH1F* w_mcg3[2];
	TH1F* w_mcg3_b[2];
	for (int i=0; i<2; i++) {
	  w_mcg[i] = read(subdir, title, i+1, mcg, dirbSel);
	  w_mcg_b[i] = read(subdir, title_b, i+1, mcg, dirbSel);
	  w_mcg1[i] = read(subdir, title, i+1, mcg1, dirbSel);
	  w_mcg1_b[i] = read(subdir, title_b, i+1, mcg1, dirbSel);
	  w_mcg2[i] = read(subdir, title, i+1, mcg2[i], dirbSel);
	  w_mcg2_b[i] = read(subdir, title_b, i+1, mcg2[i], dirbSel);
	  w_mcg3[i] = read(subdir, title, i+1, mcg3, dirbSel);
	  w_mcg3_b[i] = read(subdir, title_b, i+1, mcg3, dirbSel);
	}

	TH1F* h_mcg = (TH1F*)w_mcg[0]->Clone();
	TH1F* h_mcg_b = (TH1F*)w_mcg_b[0]->Clone();
	TH1F* h_mcg1 = (TH1F*)w_mcg1[0]->Clone();
	TH1F* h_mcg1_b = (TH1F*)w_mcg1_b[0]->Clone();
	TH1F* h_mcg2 = (TH1F*)w_mcg2[0]->Clone();
	TH1F* h_mcg2_b = (TH1F*)w_mcg2_b[0]->Clone();
	TH1F* h_mcg3 = (TH1F*)w_mcg3[0]->Clone();
	TH1F* h_mcg3_b = (TH1F*)w_mcg3_b[0]->Clone();

	h_mcg->Sumw2();
	h_mcg1->Sumw2();
	h_mcg2->Sumw2();
	h_mcg3->Sumw2();

	h_mcg_b->Sumw2();
	h_mcg1_b->Sumw2();
	h_mcg2_b->Sumw2();
	h_mcg3_b->Sumw2();

	for (int i=0;i<=h_mcg->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (w_mcg[0]->GetBinContent(i)*w_mcg[1]->GetBinContent(i) != 0) {
	    val = (w_mcg[0]->GetBinContent(i)/TMath::Power(w_mcg[0]->GetBinError(i),2)+w_mcg[1]->GetBinContent(i)/TMath::Power(w_mcg[1]->GetBinError(i),2))/(1./TMath::Power(w_mcg[0]->GetBinError(i),2)+1./TMath::Power(w_mcg[1]->GetBinError(i),2));
	    h_mcg->SetBinContent(i, val);
	    val = TMath::Sqrt(1./(1./TMath::Power(w_mcg[0]->GetBinError(i),2)+1./TMath::Power(w_mcg[1]->GetBinError(i),2)));
	    h_mcg->SetBinError(i, val);
	  }
	  if (w_mcg_b[0]->GetBinContent(i)*w_mcg_b[1]->GetBinContent(i) != 0) {
	    val = (w_mcg_b[0]->GetBinContent(i)/TMath::Power(w_mcg_b[0]->GetBinError(i),2)+w_mcg_b[1]->GetBinContent(i)/TMath::Power(w_mcg_b[1]->GetBinError(i),2))/(1./TMath::Power(w_mcg_b[0]->GetBinError(i),2)+1./TMath::Power(w_mcg_b[1]->GetBinError(i),2));
	    h_mcg_b->SetBinContent(i, val);
	    val = TMath::Sqrt(1./(1./TMath::Power(w_mcg_b[0]->GetBinError(i),2)+1./TMath::Power(w_mcg_b[1]->GetBinError(i),2)));
	    h_mcg_b->SetBinError(i, val);
	  }
	  if (w_mcg1[0]->GetBinContent(i)*w_mcg1[1]->GetBinContent(i) != 0) {
	    val = (w_mcg1[0]->GetBinContent(i)/TMath::Power(w_mcg1[0]->GetBinError(i),2)+w_mcg1[1]->GetBinContent(i)/TMath::Power(w_mcg1[1]->GetBinError(i),2))/(1./TMath::Power(w_mcg1[0]->GetBinError(i),2)+1./TMath::Power(w_mcg1[1]->GetBinError(i),2));
	    h_mcg1->SetBinContent(i, val);
	    val = TMath::Sqrt(1./(1./TMath::Power(w_mcg1[0]->GetBinError(i),2)+1./TMath::Power(w_mcg1[1]->GetBinError(i),2)));
	    h_mcg1->SetBinError(i, val);
	  }
	  if (w_mcg1_b[0]->GetBinContent(i)*w_mcg1_b[1]->GetBinContent(i) != 0) {
	    val = (w_mcg1_b[0]->GetBinContent(i)/TMath::Power(w_mcg1_b[0]->GetBinError(i),2)+w_mcg1_b[1]->GetBinContent(i)/TMath::Power(w_mcg1_b[1]->GetBinError(i),2))/(1./TMath::Power(w_mcg1_b[0]->GetBinError(i),2)+1./TMath::Power(w_mcg1_b[1]->GetBinError(i),2));
	    h_mcg1_b->SetBinContent(i, val);
	    val = TMath::Sqrt(1./(1./TMath::Power(w_mcg1_b[0]->GetBinError(i),2)+1./TMath::Power(w_mcg1_b[1]->GetBinError(i),2)));
	    h_mcg1_b->SetBinError(i, val);
	  }
	  if (w_mcg2[0]->GetBinContent(i)*w_mcg2[1]->GetBinContent(i) != 0) {
	    val = (w_mcg2[0]->GetBinContent(i)/TMath::Power(w_mcg2[0]->GetBinError(i),2)+w_mcg2[1]->GetBinContent(i)/TMath::Power(w_mcg2[1]->GetBinError(i),2))/(1./TMath::Power(w_mcg2[0]->GetBinError(i),2)+1./TMath::Power(w_mcg2[1]->GetBinError(i),2));
	    h_mcg2->SetBinContent(i, val);
	    val = TMath::Sqrt(1./(1./TMath::Power(w_mcg2[0]->GetBinError(i),2)+1./TMath::Power(w_mcg2[1]->GetBinError(i),2)));
	    h_mcg2->SetBinError(i, val);
	  }
	  if (w_mcg2_b[0]->GetBinContent(i)*w_mcg2_b[1]->GetBinContent(i) != 0) {
	    val = (w_mcg2_b[0]->GetBinContent(i)/TMath::Power(w_mcg2_b[0]->GetBinError(i),2)+w_mcg2_b[1]->GetBinContent(i)/TMath::Power(w_mcg2_b[1]->GetBinError(i),2))/(1./TMath::Power(w_mcg2_b[0]->GetBinError(i),2)+1./TMath::Power(w_mcg2_b[1]->GetBinError(i),2));
	    h_mcg2_b->SetBinContent(i, val);
	    val = TMath::Sqrt(1./(1./TMath::Power(w_mcg2_b[0]->GetBinError(i),2)+1./TMath::Power(w_mcg2_b[1]->GetBinError(i),2)));
	    h_mcg2_b->SetBinError(i, val);
	  }
	  if (w_mcg3[0]->GetBinContent(i)*w_mcg3[1]->GetBinContent(i) != 0) {
	    val = (w_mcg3[0]->GetBinContent(i)/TMath::Power(w_mcg3[0]->GetBinError(i),2)+w_mcg3[1]->GetBinContent(i)/TMath::Power(w_mcg3[1]->GetBinError(i),2))/(1./TMath::Power(w_mcg3[0]->GetBinError(i),2)+1./TMath::Power(w_mcg3[1]->GetBinError(i),2));
	    h_mcg3->SetBinContent(i, val);
	    val = TMath::Sqrt(1./(1./TMath::Power(w_mcg3[0]->GetBinError(i),2)+1./TMath::Power(w_mcg3[1]->GetBinError(i),2)));
	    h_mcg3->SetBinError(i, val);
	  }
	  if (w_mcg3_b[0]->GetBinContent(i)*w_mcg3_b[1]->GetBinContent(i) != 0) {
	    val = (w_mcg3_b[0]->GetBinContent(i)/TMath::Power(w_mcg3_b[0]->GetBinError(i),2)+w_mcg3_b[1]->GetBinContent(i)/TMath::Power(w_mcg3_b[1]->GetBinError(i),2))/(1./TMath::Power(w_mcg3_b[0]->GetBinError(i),2)+1./TMath::Power(w_mcg3_b[1]->GetBinError(i),2));
	    h_mcg3_b->SetBinContent(i, val);
	    val = TMath::Sqrt(1./(1./TMath::Power(w_mcg3_b[0]->GetBinError(i),2)+1./TMath::Power(w_mcg3_b[1]->GetBinError(i),2)));
	    h_mcg3_b->SetBinError(i, val);
	  }
	}

	for (int i=0; i<2; i++) {
	  w_data[i] = fixrange(w_data[i]);
	  w_data_b[i] = fixrange(w_data_b[i]);
	}

	h_mcg = fixrange(h_mcg);
	h_mcg_b = fixrange(h_mcg_b);
	h_mcg1 = fixrange(h_mcg1);
	h_mcg1_b = fixrange(h_mcg1_b);
	h_mcg2 = fixrange(h_mcg2);
	h_mcg3 = fixrange(h_mcg3);
	h_mcg2_b = fixrange(h_mcg2_b);
	h_mcg3_b = fixrange(h_mcg3_b);

	h_mcg->Scale(norm1);
	h_mcg1->Scale(norm1_1);
	h_mcg2->Scale(norm1_2);
	h_mcg3->Scale(norm1_3);

	h_mcg_b->Scale(norm1);
	h_mcg1_b->Scale(norm1_1);
	h_mcg2_b->Scale(norm1_2);
	h_mcg3_b->Scale(norm1_3);

	for (int i=0; i<2; i++) {
	  w_data[i]->Scale(1./Lumi2012, "width");
	  w_data_b[i]->Scale(1./Lumi2012, "width");
	  if (isratio==1) {
	    w_data_b[i]->Divide(w_data[i]);
	    w_data_b[i]->Scale(100.);
	  }
	}

	h_mcg->Scale(1./Lumi2012, "width");
	h_mcg1->Scale(1./Lumi2012, "width");
	h_mcg2->Scale(1./Lumi2012, "width");
	h_mcg3->Scale(1./Lumi2012, "width");
	h_mcg_b->Scale(1./Lumi2012, "width");
	h_mcg1_b->Scale(1./Lumi2012, "width");
	h_mcg2_b->Scale(1./Lumi2012, "width");
	h_mcg3_b->Scale(1./Lumi2012, "width");
	if (isratio==1) {
	  h_mcg_b->Divide(h_mcg);
	  h_mcg1_b->Divide(h_mcg1);
	  h_mcg2_b->Divide(h_mcg2);
	  h_mcg3_b->Divide(h_mcg);
	  h_mcg_b->Scale(100.);
	  h_mcg1_b->Scale(100.);
	  h_mcg2_b->Scale(100.);
	  h_mcg3_b->Scale(100.);
	}

	TH1F* w_stat_bkg[2];
	TH1F* w_stat_b_bkg[2];

	TH1F* w_syst_eff[2];
	TH1F* w_syst_b_eff[2];

	TH1F* w_syst_jec[2];
	TH1F* w_syst_b_jec[2];

	TH1F* w_syst_jer[2];
	TH1F* w_syst_b_jer[2];

	TH1F* w_syst_pu[2];
	TH1F* w_syst_b_pu[2];

	TH1F* w_syst_bkg[2];
	TH1F* w_syst_b_bkg[2];

	TH1F* w_stat_top[2];
	TH1F* w_stat_b_top[2];

	TH1F* w_stat_bfit[2];
	TH1F* w_stat_b_bfit[2];

	TH1F* w_syst_btag[2];
	TH1F* w_syst_b_btag[2];

	TH1F* w_stat_unfold[2];
	TH1F* w_stat_b_unfold[2];

	TH1F* w_syst_unfold[2];
	TH1F* w_syst_b_unfold[2];

	TH1F* w_syst_lumi[2];
	TH1F* w_syst_b_lumi[2];

	TH1F* w_stat_tot[2];
	TH1F* w_stat_b_tot[2];

	TH1F* w_syst_tot[2];
	TH1F* w_syst_b_tot[2];

	for (int i=0; i<2; i++) {

	  w_stat_bkg[i] = (TH1F*)w_data[0]->Clone();
	  w_stat_b_bkg[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_eff[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_eff[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_jec[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_jec[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_jer[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_jer[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_pu[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_pu[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_bkg[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_bkg[i] = (TH1F*)w_data_b[0]->Clone();

	  w_stat_top[i] = (TH1F*)w_data[0]->Clone();
	  w_stat_b_top[i] = (TH1F*)w_data_b[0]->Clone();

	  w_stat_bfit[i] = (TH1F*)w_data[0]->Clone();
	  w_stat_b_bfit[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_btag[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_btag[i] = (TH1F*)w_data_b[0]->Clone();

	  w_stat_unfold[i] = (TH1F*)w_data[0]->Clone();
	  w_stat_b_unfold[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_unfold[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_unfold[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_lumi[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_lumi[i] = (TH1F*)w_data_b[0]->Clone();

	  w_stat_tot[i] = (TH1F*)w_data[0]->Clone();
	  w_stat_b_tot[i] = (TH1F*)w_data_b[0]->Clone();

	  w_syst_tot[i] = (TH1F*)w_data[0]->Clone();
	  w_syst_b_tot[i] = (TH1F*)w_data_b[0]->Clone();

	  ifstream in;
	  string title_b_tmp = title_b;
	  if (i==0) {
	    if (title_b=="w_pt_Z") title_b_tmp="w_pt_Z_ee";
	    if (title_b=="w_pt_Z_b") title_b_tmp="w_pt_Z_ee_b";
	    if (title_b=="w_delta_phi") title_b_tmp="w_delta_phi_ee";
	    if (title_b=="w_delta_phi_b") title_b_tmp="w_delta_phi_ee_b";
	    if (title_b=="w_mass_Zj") title_b_tmp="w_mass_Zj_ee";
	    if (title_b=="w_mass_Zj_b") title_b_tmp="w_mass_Zj_ee_b";
            if (title_b=="w_Zbb_mass") title_b_tmp="w_eebb_mass";
            if (title_b=="w_DR_Zb_min") title_b_tmp="w_DR_eeb_min";
            if (title_b=="w_DR_Zb_max") title_b_tmp="w_DR_eeb_max";
            if (title_b=="w_A_Zb") title_b_tmp="w_A_eeb";
	    if (isratio==0) in.open((path + "/electrons/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b_tmp + "_xsecs_unfolding.dat").c_str());
	    if (isratio==1) in.open((path + "/electrons/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b_tmp + "_ratio_unfolding.dat").c_str());
	  }
	  if (i==1) {
	    if (title_b=="w_pt_Z") title_b_tmp="w_pt_Z_mm";
	    if (title_b=="w_pt_Z_b") title_b_tmp="w_pt_Z_mm_b";
	    if (title_b=="w_delta_phi") title_b_tmp="w_delta_phi_mm";
	    if (title_b=="w_delta_phi_b") title_b_tmp="w_delta_phi_mm_b";
	    if (title_b=="w_mass_Zj") title_b_tmp="w_mass_Zj_mm";
	    if (title_b=="w_mass_Zj_b") title_b_tmp="w_mass_Zj_mm_b";
            if (title_b=="w_Zbb_mass") title_b_tmp="w_mmbb_mass";
            if (title_b=="w_DR_Zb_min") title_b_tmp="w_DR_mmb_min";
            if (title_b=="w_DR_Zb_max") title_b_tmp="w_DR_mmb_max";
            if (title_b=="w_A_Zb") title_b_tmp="w_A_mmb";
            if (isratio==0) in.open((path + "/muons/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b_tmp + "_xsecs_unfolding.dat").c_str());
	    if (isratio==1) in.open((path + "/muons/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b_tmp + "_ratio_unfolding.dat").c_str());
	  }

	  string tmp;

	  if (isratio==0) {
	    getline(in, tmp);
	    getline(in, tmp);
	    getline(in, tmp);
	    for (int j=0; j<w_data[0]->GetNbinsX()+2; j++) {
	      in >> tmp;
	      double val = 0.0;
	      in >> val; w_data[i]->SetBinContent(j, val); in >> tmp;
	      in >> val; w_data[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_stat_bkg[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_eff[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_jec[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_jer[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_pu[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_bkg[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_stat_top[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_stat_bfit[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_btag[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_stat_unfold[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_unfold[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_lumi[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_stat_tot[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_tot[i]->SetBinError(j, val); in >> tmp;
	      in >> val; in >> tmp; in >> val;
	      in.ignore();
	    }
	  }

	  getline(in, tmp);
	  getline(in, tmp);
	  getline(in, tmp);
	  for (int j=0; j<w_data_b[0]->GetNbinsX()+2; j++) {
	    in >> tmp;
	    double val = 0.0;
	    in >> val; w_data_b[i]->SetBinContent(j, val); in >> tmp;
	    in >> val; w_data_b[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_b_bkg[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_eff[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_jec[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_jer[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_pu[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_bkg[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_b_top[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_b_bfit[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_btag[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_b_unfold[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_unfold[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_lumi[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_b_tot[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_tot[i]->SetBinError(j, val); in >> tmp;
	    in >> val; in >> tmp; in >> val;
	    in.ignore();
	  }

	  in.close();

	}

	TH1F* h_data = (TH1F*)w_data[0]->Clone();
	TH1F* h_data_b = (TH1F*)w_data_b[0]->Clone();
	TH1F* h_data_stat = (TH1F*)w_data[0]->Clone();
	TH1F* h_data_b_stat = (TH1F*)w_data_b[0]->Clone();
	TH1F* h_data_syst = (TH1F*)w_data[0]->Clone();
	TH1F* h_data_b_syst = (TH1F*)w_data[0]->Clone();
	TH1F* h_data_tot = (TH1F*)w_data[0]->Clone();
	TH1F* h_data_b_tot = (TH1F*)w_data[0]->Clone();

	TH1F* stat_bkg = (TH1F*)w_data[0]->Clone();
	TH1F* stat_b_bkg = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_eff = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_eff = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_jec = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_jec = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_jer = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_jer = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_pu = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_pu = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_bkg = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_bkg = (TH1F*)w_data_b[0]->Clone();

	TH1F* stat_top = (TH1F*)w_data[0]->Clone();
	TH1F* stat_b_top = (TH1F*)w_data_b[0]->Clone();

	TH1F* stat_bfit = (TH1F*)w_data[0]->Clone();
	TH1F* stat_b_bfit = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_btag = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_btag = (TH1F*)w_data_b[0]->Clone();

	TH1F* stat_unfold = (TH1F*)w_data[0]->Clone();
	TH1F* stat_b_unfold = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_unfold = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_unfold = (TH1F*)w_data_b[0]->Clone();

	TH1F* syst_lumi = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_lumi = (TH1F*)w_data_b[0]->Clone();

	for (int i=0;i<=h_data_stat->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = calc(0, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  h_data->SetBinContent(i, val);
	  h_data_stat->SetBinContent(i, val);
	  h_data_syst->SetBinContent(i, val);
	  h_data_tot->SetBinContent(i, val);
	  double ref = 0.0;
	  ref = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			1.1*w_data[0]->GetBinError(i), 1.1*w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  h_data->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			1.1*w_stat_bkg[0]->GetBinError(i), 1.1*w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_bkg->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			1.1*w_syst_eff[0]->GetBinError(i), 1.1*w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_eff->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			1.1*w_syst_jer[0]->GetBinError(i), 1.1*w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_jer->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			1.1*w_syst_jec[0]->GetBinError(i), 1.1*w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_jec->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			1.1*w_syst_pu[0]->GetBinError(i), 1.1*w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_pu->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			1.1*w_syst_bkg[0]->GetBinError(i), 1.1*w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_bkg->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			1.1*w_stat_top[0]->GetBinError(i), 1.1*w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_top->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			1.1*w_stat_bfit[0]->GetBinError(i), 1.1*w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_bfit->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			1.1*w_syst_btag[0]->GetBinError(i), 1.1*w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_btag->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			1.1*w_stat_unfold[0]->GetBinError(i), 1.1*w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_unfold->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			1.1*w_syst_unfold[0]->GetBinError(i), 1.1*w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_unfold->SetBinError(i, val);
	  val = calc(1, w_data[0]->GetBinContent(i), w_data[1]->GetBinContent(i),
			w_data[0]->GetBinError(i), w_data[1]->GetBinError(i),
			w_stat_bkg[0]->GetBinError(i), w_stat_bkg[1]->GetBinError(i),
			w_syst_eff[0]->GetBinError(i), w_syst_eff[1]->GetBinError(i),
			w_syst_jer[0]->GetBinError(i), w_syst_jer[1]->GetBinError(i),
			w_syst_jec[0]->GetBinError(i), w_syst_jec[1]->GetBinError(i),
			w_syst_pu[0]->GetBinError(i), w_syst_pu[1]->GetBinError(i),
			w_syst_bkg[0]->GetBinError(i), w_syst_bkg[1]->GetBinError(i),
			w_stat_top[0]->GetBinError(i), w_stat_top[1]->GetBinError(i),
			w_stat_bfit[0]->GetBinError(i), w_stat_bfit[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			1.1*w_syst_lumi[0]->GetBinError(i), 1.1*w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_lumi->SetBinError(i, val);

	  val = TMath::Sqrt(TMath::Power(h_data->GetBinError(i),2)+TMath::Power(stat_top->GetBinError(i),2)+TMath::Power(stat_bfit->GetBinError(i),2));
	  h_data_stat->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Power(stat_bkg->GetBinError(i),2)+TMath::Power(syst_eff->GetBinError(i),2)+TMath::Power(syst_jec->GetBinError(i),2)+TMath::Power(syst_jer->GetBinError(i),2)+TMath::Power(syst_pu->GetBinError(i),2)+TMath::Power(syst_bkg->GetBinError(i),2)+TMath::Power(syst_btag->GetBinError(i),2)+TMath::Power(stat_unfold->GetBinError(i),2)+TMath::Power(syst_unfold->GetBinError(i),2)+TMath::Power(syst_lumi->GetBinError(i),2));
	  h_data_syst->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Power(h_data_stat->GetBinError(i),2)+TMath::Power(h_data_syst->GetBinError(i),2));
	  h_data_tot->SetBinError(i, val);
	}

	for (int i=0;i<=h_data_b_stat->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = calc(0, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  h_data_b->SetBinContent(i, val);
	  h_data_b_stat->SetBinContent(i, val);
	  h_data_b_syst->SetBinContent(i, val);
	  h_data_b_tot->SetBinContent(i, val);
	  double ref = 0.0;
	  ref = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			1.1*w_data_b[0]->GetBinError(i), 1.1*w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  h_data_b->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			1.1*w_stat_b_bkg[0]->GetBinError(i), 1.1*w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_b_bkg->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			1.1*w_syst_b_eff[0]->GetBinError(i), 1.1*w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_eff->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			1.1*w_syst_b_jer[0]->GetBinError(i), 1.1*w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_jer->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			1.1*w_syst_b_jec[0]->GetBinError(i), 1.1*w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_jec->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			1.1*w_syst_b_pu[0]->GetBinError(i), 1.1*w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_pu->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			1.1*w_syst_b_bkg[0]->GetBinError(i), 1.1*w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_bkg->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			1.1*w_stat_b_top[0]->GetBinError(i), 1.1*w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_b_top->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			1.1*w_stat_b_bfit[0]->GetBinError(i), 1.1*w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_b_bfit->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			1.1*w_syst_b_btag[0]->GetBinError(i), 1.1*w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_btag->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			1.1*w_stat_b_unfold[0]->GetBinError(i), 1.1*w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  stat_b_unfold->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			1.1*w_syst_b_unfold[0]->GetBinError(i), 1.1*w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_unfold->SetBinError(i, val);
	  val = calc(1, w_data_b[0]->GetBinContent(i), w_data_b[1]->GetBinContent(i),
			w_data_b[0]->GetBinError(i), w_data_b[1]->GetBinError(i),
			w_stat_b_bkg[0]->GetBinError(i), w_stat_b_bkg[1]->GetBinError(i),
			w_syst_b_eff[0]->GetBinError(i), w_syst_b_eff[1]->GetBinError(i),
			w_syst_b_jer[0]->GetBinError(i), w_syst_b_jer[1]->GetBinError(i),
			w_syst_b_jec[0]->GetBinError(i), w_syst_b_jec[1]->GetBinError(i),
			w_syst_b_pu[0]->GetBinError(i), w_syst_b_pu[1]->GetBinError(i),
			w_syst_b_bkg[0]->GetBinError(i), w_syst_b_bkg[1]->GetBinError(i),
			w_stat_b_top[0]->GetBinError(i), w_stat_b_top[1]->GetBinError(i),
			w_stat_b_bfit[0]->GetBinError(i), w_stat_b_bfit[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			1.1*w_syst_b_lumi[0]->GetBinError(i), 1.1*w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_lumi->SetBinError(i, val);

	  val = TMath::Sqrt(TMath::Power(h_data_b->GetBinError(i),2)+TMath::Power(stat_b_top->GetBinError(i),2)+TMath::Power(stat_b_bfit->GetBinError(i),2));
	  h_data_b_stat->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Power(stat_b_bkg->GetBinError(i),2)+TMath::Power(syst_b_eff->GetBinError(i),2)+TMath::Power(syst_b_jec->GetBinError(i),2)+TMath::Power(syst_b_jer->GetBinError(i),2)+TMath::Power(syst_b_pu->GetBinError(i),2)+TMath::Power(syst_b_bkg->GetBinError(i),2)+TMath::Power(syst_b_btag->GetBinError(i),2)+TMath::Power(stat_b_unfold->GetBinError(i),2)+TMath::Power(syst_b_unfold->GetBinError(i),2)+TMath::Power(syst_b_lumi->GetBinError(i),2));
	  h_data_b_syst->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Power(h_data_b_stat->GetBinError(i),2)+TMath::Power(h_data_b_syst->GetBinError(i),2));
	  h_data_b_tot->SetBinError(i, val);
	}

	TCanvas* c1 = new TCanvas("c", "c", 10, 10, 800, 600);
	c1->cd();

	h_mcg_b->SetTitle("");
	if (isratio==1) {
	  h_mcg_b->GetYaxis()->SetTitle("#sigma(Z+b) / #sigma(Z+j) [%]");
	} else {
	  h_mcg_b->GetYaxis()->SetTitle("#sigma [pb]");
	}

	TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
	pad1->SetTopMargin(0.115);
	pad1->SetBottomMargin(0.0001);
	pad1->Draw();
	pad1->cd();

	h_mcg_b->SetLineColor(kGreen+2);
	h_mcg_b->SetLineWidth(2);
	h_mcg_b->SetFillColor(kGreen+2);
	h_mcg_b->SetMarkerColor(kGreen+2);
	h_mcg_b->SetStats(0);
	if (isratio==1) {
	  h_mcg_b->Draw("E5");
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

	h_mcg3_b->SetLineColor(kOrange+7);
	h_mcg3_b->SetLineWidth(2);
	h_mcg3_b->SetFillColor(kOrange+7);
	h_mcg3_b->SetMarkerColor(kOrange+7);
	if (isratio==1) {
	  h_mcg3_b->Draw("E5SAME");
	}

	if (isratio==1) {
	  h_data_b_tot->GetYaxis()->SetTitle("#sigma_{Z+b-jets}/#sigma_{Z+jets} [%]");
	}
	h_data_b_tot->GetYaxis()->SetTitleOffset(1.2);
	h_data_b_tot->GetXaxis()->SetTitleOffset(1.3);
	h_data_b_tot->SetMarkerColor(kRed+1);
	h_data_b_tot->SetLineColor(kRed+1);
	h_data_b_tot->SetMarkerStyle(24);
	h_data_b_tot->SetMarkerSize(0.7);
	h_data_b_tot->SetStats(0);
	h_data_b_stat->GetYaxis()->SetTitleOffset(1.2);
	h_data_b_stat->GetXaxis()->SetTitleOffset(1.3);
	h_data_b_stat->SetMarkerColor(kBlack);
	h_data_b_stat->SetLineColor(kBlack);
	h_data_b_stat->SetMarkerStyle(24);
	h_data_b_stat->SetMarkerSize(0.7);
	h_data_b_stat->SetStats(0);
	if (isratio==1) {
	  h_data_b_tot->Draw("E1PX0SAME");
	  h_data_b_stat->Draw("E1PX0SAME");
	}

	TLegend *leg = new TLegend(0.62, 0.580, 0.88, 0.88);
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (isratio==0) {
	  pad1->SetLogy();

	  h_mcg_b->SetMaximum(4*h_data_tot->GetMaximum());
	  h_mcg_b->SetMinimum(TMath::Max(0.000002,0.25*h_data_b_tot->GetBinContent(h_data_b_tot->GetMinimumBin())));

	  h_mcg_b->Draw("E5");
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

	  h_mcg3_b->Draw("E5SAME");
	  TH1F* tmp2_3 = (TH1F*)h_mcg3_b->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp2_3->GetMinimum()==0) tmp2_3->GetXaxis()->SetRangeUser(0, tmp2_3->GetBinCenter(tmp2_3->GetMinimumBin()-1));
	  }
	  tmp2_3->SetFillColor(0);
	  tmp2_3->DrawClone("HISTLSAME");

	  h_data_b_tot->Draw("E1PX0SAME");
	  h_data_b_stat->Draw("E1PX0SAME");

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

	  h_mcg3->SetLineColor(kOrange+7);
	  h_mcg3->SetLineWidth(2);
	  h_mcg3->SetMarkerColor(kOrange+7);
	  h_mcg3->SetFillColor(kOrange+7);

	  h_data_tot->SetMarkerColor(kRed+1);
	  h_data_tot->SetLineColor(kRed+1);
	  h_data_tot->SetMarkerStyle(20);
	  h_data_tot->SetMarkerSize (0.7);
	  h_data_stat->SetLineColor(kBlack);
	  h_data_stat->SetMarkerColor(kBlack);
	  h_data_stat->SetMarkerStyle(20);
	  h_data_stat->SetMarkerSize (0.7);
	  if (drawInclusive) h_data_tot->Draw("E1PX0SAME");
	  if (drawInclusive) h_data_stat->Draw("E1PX0SAME");

	  if (drawInclusive) leg->AddEntry(h_data_stat,"Z(#rightarrow ll) DATA","p");
	  leg->AddEntry(h_data_b_stat,"Z(#rightarrow ll)+b DATA","p");
	  leg->AddEntry(h_mcg,"Z(#rightarrow ll) MadGraph 5FS","l");
	  leg->AddEntry(h_mcg3,"Z(#rightarrow ll) MadGraph 4FS","l");
	  if (useSherpa) leg->AddEntry(h_mcg1,"Z(#rightarrow ll) Sherpa","l");
	  leg->AddEntry(h_mcg2,"Z(#rightarrow ll) Powheg","l");
	}

	if (isratio==1) {
	  leg->AddEntry(h_data_b_stat,"Z(#rightarrow ll) DATA","p");
	  leg->AddEntry(h_mcg_b,"Z(#rightarrow ll) MadGraph 5FS","l");
	  leg->AddEntry(h_mcg3_b,"Z(#rightarrow ll) MadGraph 4FS","l");
	  if (useSherpa) leg->AddEntry(h_mcg1_b,"Z(#rightarrow ll) Sherpa","l");
	  leg->AddEntry(h_mcg2_b,"Z(#rightarrow ll) Powheg","l");
	}

	leg->Draw();

	c1->cd();
       
 	TLatex *latexLabel = CMSFinal (Lumi2012/1000., "Z/#gamma*#rightarrow ll selection" , 0, 0.135, 0.51);
 	
        if (isratio==1) {
          latexLabel = CMSFinal (Lumi2012/1000., "Z/#gamma*#rightarrow ll selection", 0, 0.135, 0.85);
        }
        if (isratio==0) {
          if (title_b=="w_Ht_b" || title_b=="w_first_bjet_pt" || title_b=="w_pt_Z_b" || title_b=="w_DR_bb" || title_b=="w_bb_mass" || title_b=="w_Zbb_mass"|| title_b=="w_DR_Zb_min"|| title_b=="w_DR_Zb_max"|| title_b=="w_A_Zb") {
            latexLabel = CMSFinal (Lumi2012/1000., "Z/#gamma*#rightarrow ll selection", 0, 0.135, 0.51);
          }
          if (title_b=="w_delta_phi_b" || title_b=="w_delta_phi_2b") {
            latexLabel = CMSFinal (Lumi2012/1000., "Z/#gamma*#rightarrow ll selection", 0, 0.68, 0.51);
          }
          if (title_b=="w_first_bjet_eta") {
            latexLabel = CMSFinal (Lumi2012/1000., "Z/#gamma*#rightarrow ll selection", 0, 0.68, 0.51);
          }
        }
	latexLabel->Draw("same");

	TPad *pad2 = new TPad("pad2","pad2",0,0.29,1,0.4);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.001);
	pad2->Draw();
	pad2->cd();

	TH1F *h_M_tot = (TH1F*)h_mcg_b->Clone();
	TH1F *h_M_stat = (TH1F*)h_mcg_b->Clone();

	h_M_tot->Divide(h_data_b_tot);
	h_M_stat->Divide(h_data_b_stat);

	h_M_tot->SetTitle("");
	h_M_tot->SetStats(0);
	h_M_tot->GetXaxis()->SetTitleOffset(0.9);
	h_M_tot->GetXaxis()->SetTitleSize(0.14);
	h_M_tot->GetXaxis()->SetLabelFont(42);
	h_M_tot->GetXaxis()->SetLabelSize(0.12);
	h_M_tot->GetXaxis()->SetTitleFont(42);
	h_M_tot->GetXaxis()->SetTickLength(0.1);
	h_M_tot->GetYaxis()->SetTitle("Theory / Data");
	h_M_tot->GetYaxis()->SetNdivisions(013);
	h_M_tot->GetYaxis()->SetTitleSize(0.17);
	h_M_tot->GetYaxis()->SetLabelSize(0.17);
	h_M_tot->GetYaxis()->SetRangeUser(-0.2, 2.2);
	h_M_tot->GetYaxis()->SetTitleOffset(0.21);
	h_M_tot->GetYaxis()->SetTickLength(0.02);

	h_M_tot->SetMarkerColor(kRed+1);
	h_M_tot->SetLineColor(kRed+1);
	h_M_tot->SetLineWidth(1);
	h_M_tot->SetMarkerSize(0.7);
	h_M_stat->GetXaxis()->SetTitleOffset(0.7);
	h_M_stat->SetMarkerColor(kBlack);
	h_M_stat->SetLineColor(kBlack);
	h_M_stat->SetLineWidth(1);
	h_M_stat->SetMarkerSize(0.7);

	h_M_tot->SetMarkerStyle(24);
	h_M_tot->Draw("E1PX0");
	h_M_stat->SetMarkerStyle(24);
	h_M_stat->Draw("E1PX0SAME");

	if (isratio==0) {
	  TH1F *h_M2_tot= (TH1F*)h_mcg->Clone();
	  TH1F *h_M2_stat= (TH1F*)h_mcg->Clone();

	  h_M2_tot->Divide(h_data_tot);
	  h_M2_stat->Divide(h_data_stat);

	  TGraphErrors *g_M2_tot = new TGraphErrors(h_M2_tot);
	  TGraphErrors *g_M2_stat = new TGraphErrors(h_M2_stat);

	  float dx = 0.1*(g_M2_tot->GetXaxis()->GetXmax()-g_M2_tot->GetXaxis()->GetXmin())/g_M2_tot->GetN();
	  for (int i=0; i<g_M2_tot->GetN(); i++) {
	    g_M2_stat->SetPoint(i, g_M2_stat->GetX()[i]-dx, g_M2_stat->GetY()[i]);
	    g_M2_stat->SetPointError(i, 0, g_M2_stat->GetEY()[i]);
	    g_M2_tot->SetPoint(i, g_M2_tot->GetX()[i]-dx, g_M2_tot->GetY()[i]);
	    g_M2_tot->SetPointError(i, 0, g_M2_tot->GetEY()[i]);
	  }

	  g_M2_tot->SetMarkerColor(kRed+1);
	  g_M2_tot->SetLineColor(kRed+1);
	  g_M2_tot->SetLineWidth(1);
	  g_M2_tot->SetMarkerSize(0.7);
	  g_M2_stat->GetXaxis()->SetTitleOffset(0.7);
	  g_M2_stat->SetMarkerColor(kBlack);
	  g_M2_stat->SetLineColor(kBlack);
	  g_M2_stat->SetLineWidth(1);
	  g_M2_stat->SetMarkerSize(0.7);

	  g_M2_tot->SetMarkerStyle(20);
	  if (drawInclusive) g_M2_tot->Draw("E1PX0SAME");
	  g_M2_stat->SetMarkerStyle(20);
	  if (drawInclusive) g_M2_stat->Draw("E1PX0SAME");
	}

	TLatex *t2 = new TLatex();
	t2->SetTextSize(0.2);
	t2->SetTextFont(42);
	t2->SetLineWidth(2);
	t2->SetNDC();
	if (useSherpa) {
	  t2->DrawLatex(0.15,0.7,"MadGraph 5FS / MadGraph 4FS");
	} else {
	  t2->DrawLatex(0.15,0.7,"MadGraph 5FS");
	}

	TLine *OLine2 = new TLine(h_M_tot->GetXaxis()->GetXmin(),1.,h_M_tot->GetXaxis()->GetXmax(),1.);
	OLine2->SetLineColor(kGreen+2);
	OLine2->SetLineWidth(2);
	OLine2->Draw();

	c1->cd();

	TPad *pad3 = new TPad("pad3","pad3",0,0.18,1,0.29);
	pad3->SetTopMargin(0);
	pad3->SetBottomMargin(0.001);
	pad3->Draw();
	pad3->cd();

	TH1F *h_S_tot = (TH1F*)h_mcg1_b->Clone();
	TH1F *h_S_stat = (TH1F*)h_mcg1_b->Clone();

	h_S_tot->Divide(h_data_b_tot);
	h_S_stat->Divide(h_data_b_stat);

	h_S_tot->SetTitle("");
	h_S_tot->SetStats(0);
	h_S_tot->GetXaxis()->SetTitleOffset(0.9);
	h_S_tot->GetXaxis()->SetTitleSize(0.14);
	h_S_tot->GetXaxis()->SetLabelFont(42);
	h_S_tot->GetXaxis()->SetLabelSize(0.12);
	h_S_tot->GetXaxis()->SetTitleFont(42);
	h_S_tot->GetXaxis()->SetTickLength(0.1);
	h_S_tot->GetYaxis()->SetTitle("Theory / Data");
	h_S_tot->GetYaxis()->SetNdivisions(013);
	h_S_tot->GetYaxis()->SetTitleSize(0.17);
	h_S_tot->GetYaxis()->SetLabelSize(0.17);
	h_S_tot->GetYaxis()->SetRangeUser(-0.2, 2.2);
	h_S_tot->GetYaxis()->SetTitleOffset(0.21);
	h_S_tot->GetYaxis()->SetTickLength(0.02);

	h_S_tot->SetMarkerColor(kRed+1);
	h_S_tot->SetLineColor(kRed+1);
	h_S_tot->SetLineWidth(1);
	h_S_tot->SetMarkerSize(0.7);
	h_S_stat->GetXaxis()->SetTitleOffset(0.7);
	h_S_stat->SetMarkerColor(kBlack);
	h_S_stat->SetLineColor(kBlack);
	h_S_stat->SetLineWidth(1);
	h_S_stat->SetMarkerSize(0.7);

	h_S_tot->SetMarkerStyle(24);
	if (useSherpa) {
	  h_S_tot->Draw("E1PX0");
	} else {
	  for (int i=0;i<=h_S_tot->GetNbinsX()+1;i++) {
	    h_S_tot->SetBinContent(i, -0.5);
	  }
	  h_S_tot->Draw("E1PX0");
	}
	h_S_stat->SetMarkerStyle(24);
	if (useSherpa) h_S_stat->Draw("E1PX0SAME");

	if (isratio==0) {
	  TH1F *h_S2_tot= (TH1F*)h_mcg1->Clone();
	  TH1F *h_S2_stat= (TH1F*)h_mcg1->Clone();

	  h_S2_tot->Divide(h_data_tot);
	  h_S2_stat->Divide(h_data_stat);

	  TGraphErrors *g_S2_tot = new TGraphErrors(h_S2_tot);
	  TGraphErrors *g_S2_stat = new TGraphErrors(h_S2_stat);

	  float dx = 0.1*(g_S2_tot->GetXaxis()->GetXmax()-g_S2_tot->GetXaxis()->GetXmin())/g_S2_tot->GetN();
	  for (int i=0; i<g_S2_tot->GetN(); i++) {
	    g_S2_stat->SetPoint(i, g_S2_stat->GetX()[i]-dx, g_S2_stat->GetY()[i]);
	    g_S2_stat->SetPointError(i, 0, g_S2_stat->GetEY()[i]);
	    g_S2_tot->SetPoint(i, g_S2_tot->GetX()[i]-dx, g_S2_tot->GetY()[i]);
	    g_S2_tot->SetPointError(i, 0, g_S2_tot->GetEY()[i]);
	  }

	  g_S2_tot->SetMarkerColor(kRed+1);
	  g_S2_tot->SetLineColor(kRed+1);
	  g_S2_tot->SetLineWidth(1);
	  g_S2_tot->SetMarkerSize(0.7);
	  g_S2_stat->GetXaxis()->SetTitleOffset(0.7);
	  g_S2_stat->SetMarkerColor(kBlack);
	  g_S2_stat->SetLineColor(kBlack);
	  g_S2_stat->SetLineWidth(1);
	  g_S2_stat->SetMarkerSize(0.7);

	  g_S2_tot->SetMarkerStyle(20);
	  if (useSherpa) g_S2_tot->Draw("E1PX0SAME");
	  g_S2_stat->SetMarkerStyle(20);
	  if (useSherpa) g_S2_stat->Draw("E1PX0SAME");
	}

	TLatex *t3 = new TLatex();
	t3->SetTextSize(0.2);
	t3->SetTextFont(42);
	t3->SetLineWidth(2);
	t3->SetNDC();
	if (useSherpa) {
	  t3->DrawLatex(0.15,0.7,"Sherpa");
	} else {
	  t3->DrawLatex(0.15,0.7,"MadGraph 4FS");
	}

	if (useSherpa) {
	  TLine *OLine3 = new TLine(h_S_tot->GetXaxis()->GetXmin(),1.,h_S_tot->GetXaxis()->GetXmax(),1.);
	  OLine3->SetLineColor(kMagenta-6);
	  OLine3->SetLineWidth(2);
	  OLine3->Draw();
	}

	c1->cd();

	TPad *pad4 = new TPad("pad4","pad4",0,0.0,1,0.18);
	pad4->SetTopMargin(0);
	pad4->SetBottomMargin(0.3);
	pad4->Draw();
	pad4->cd();

	TH1F *h_P_tot = (TH1F*)h_mcg2_b->Clone();
	TH1F *h_P_stat = (TH1F*)h_mcg2_b->Clone();

	h_P_tot->Divide(h_data_b_tot);
	h_P_stat->Divide(h_data_b_stat);

	h_P_tot->SetTitle("");
	h_P_tot->SetStats(0);
	h_P_tot->GetXaxis()->SetTitleOffset(0.9);
	h_P_tot->GetXaxis()->SetTitleSize(0.14);
	h_P_tot->GetXaxis()->SetLabelFont(42);
	h_P_tot->GetXaxis()->SetLabelSize(0.12);
	h_P_tot->GetXaxis()->SetTitleFont(42);
	h_P_tot->GetXaxis()->SetTickLength(0.1);
	h_P_tot->GetYaxis()->SetTitle("Theory / Data");
	h_P_tot->GetYaxis()->SetNdivisions(013);
	h_P_tot->GetYaxis()->SetTitleSize(0.11);
	h_P_tot->GetYaxis()->SetLabelSize(0.11);
	h_P_tot->GetYaxis()->SetRangeUser(-0.2, 2.2);
	h_P_tot->GetYaxis()->SetTitleOffset(0.32);
	h_P_tot->GetYaxis()->SetTickLength(0.02);

	h_P_tot->SetMarkerColor(kRed+1);
	h_P_tot->SetLineColor(kRed+1);
	h_P_tot->SetLineWidth(1);
	h_P_tot->SetMarkerSize(0.7);
	h_P_stat->GetXaxis()->SetTitleOffset(0.7);
	h_P_stat->SetMarkerColor(kBlack);
	h_P_stat->SetLineColor(kBlack);
	h_P_stat->SetLineWidth(1);
	h_P_stat->SetMarkerSize(0.7);

	h_P_tot->SetMarkerStyle(24);
	h_P_tot->Draw("E1PX0");
	h_P_stat->SetMarkerStyle(24);
	h_P_stat->Draw("E1PX0SAME");

	if (isratio==0) {
	  TH1F *h_P2_tot= (TH1F*)h_mcg2->Clone();
	  TH1F *h_P2_stat= (TH1F*)h_mcg2->Clone();

	  h_P2_tot->Divide(h_data_tot);
	  h_P2_stat->Divide(h_data_stat);

	  TGraphErrors *g_P2_tot = new TGraphErrors(h_P2_tot);
	  TGraphErrors *g_P2_stat = new TGraphErrors(h_P2_stat);

	  float dx = 0.1*(g_P2_stat->GetXaxis()->GetXmax()-g_P2_stat->GetXaxis()->GetXmin())/g_P2_stat->GetN();
	  for (int i=0; i<g_P2_stat->GetN(); i++) {
	    g_P2_stat->SetPoint(i, g_P2_stat->GetX()[i]-dx, g_P2_stat->GetY()[i]);
	    g_P2_stat->SetPointError(i, 0, g_P2_stat->GetEY()[i]);
	    g_P2_tot->SetPoint(i, g_P2_tot->GetX()[i]-dx, g_P2_tot->GetY()[i]);
	    g_P2_tot->SetPointError(i, 0, g_P2_tot->GetEY()[i]);
	  }

	  g_P2_tot->SetMarkerColor(kRed+1);
	  g_P2_tot->SetLineColor(kRed+1);
	  g_P2_tot->SetLineWidth(1);
	  g_P2_tot->SetMarkerSize(0.7);
	  g_P2_stat->GetXaxis()->SetTitleOffset(0.7);
	  g_P2_stat->SetMarkerColor(kBlack);
	  g_P2_stat->SetLineColor(kBlack);
	  g_P2_stat->SetLineWidth(1);
	  g_P2_stat->SetMarkerSize(0.7);

	  g_P2_tot->SetMarkerStyle(20);
	  if (drawInclusive) g_P2_tot->Draw("E1PX0SAME");
	  g_P2_stat->SetMarkerStyle(20);
	  if (drawInclusive) g_P2_stat->Draw("E1PX0SAME");
	}

	TLatex *t4 = new TLatex();
	t4->SetTextSize(0.13);
	t4->SetTextFont(42);
	t4->SetLineWidth(2);
	t4->SetNDC();
	t4->DrawLatex(0.15,0.8,"Powheg");

	TLine *OLine4 = new TLine(h_P_tot->GetXaxis()->GetXmin(),1.,h_P_tot->GetXaxis()->GetXmax(),1.);
	OLine4->SetLineColor(kBlue-4);
	OLine4->SetLineWidth(2);
	OLine4->Draw();

	if (useSherpa) {
	  pad2->cd();
	} else {
	  pad3->cd();
	}

	TH1F *h_M3_tot = (TH1F*)h_mcg3_b->Clone();
	TH1F *h_M3_stat = (TH1F*)h_mcg3_b->Clone();

	h_M3_tot->Divide(h_data_b_tot);
	h_M3_stat->Divide(h_data_b_stat);

	TGraphErrors *g_M3_tot = new TGraphErrors(h_M3_tot);
	TGraphErrors *g_M3_stat = new TGraphErrors(h_M3_stat);

	float dx = 0.1*(g_M3_tot->GetXaxis()->GetXmax()-g_M3_tot->GetXaxis()->GetXmin())/g_M3_tot->GetN();
	for (int i=0; i<g_M3_tot->GetN(); i++) {
	  g_M3_stat->SetPoint(i, g_M3_stat->GetX()[i]+dx, g_M3_stat->GetY()[i]);
	  g_M3_stat->SetPointError(i, 0, g_M3_stat->GetEY()[i]);
	  g_M3_tot->SetPoint(i, g_M3_tot->GetX()[i]+dx, g_M3_tot->GetY()[i]);
	  g_M3_tot->SetPointError(i, 0, g_M3_tot->GetEY()[i]);
	}

	g_M3_tot->SetMarkerColor(kRed+1);
	g_M3_tot->SetLineColor(kRed+1);
	g_M3_tot->SetLineWidth(1);
	g_M3_tot->SetMarkerSize(0.7);
	g_M3_stat->GetXaxis()->SetTitleOffset(0.7);
	g_M3_stat->SetMarkerColor(kBlack);
	g_M3_stat->SetLineColor(kBlack);
	g_M3_stat->SetLineWidth(1);
	g_M3_stat->SetMarkerSize(0.7);

	if (useSherpa) {
	  g_M3_stat->SetMarkerStyle(25);
	  g_M3_tot->SetMarkerStyle(25);
	} else {
	  g_M3_stat->SetMarkerStyle(24);
	  g_M3_stat->SetMarkerColor(kBlack);
	  g_M3_tot->SetMarkerStyle(24);
	  g_M3_tot->SetMarkerColor(kBlack);
	}
	if (useSherpa) {
	  g_M3_tot->SetMarkerStyle(25);
	} else {
	  g_M3_tot->SetMarkerStyle(24);
	}
	g_M3_tot->Draw("E1P");
	g_M3_stat->Draw("E1PSAME");

	TLine *OLine5 = new TLine(h_P_tot->GetXaxis()->GetXmin(),0.93,h_P_tot->GetXaxis()->GetXmax(),0.93);
	OLine5->SetLineColor(kOrange+7);
	OLine5->SetLineWidth(2);
	OLine5->Draw();

	c1->cd();

	if (isratio==1) {
	  h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 10);
	}

	if (title_b=="w_first_jet_pt_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp_{T} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("leading jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 20);
	  }
	} else if (title_b=="w_first_jet_eta_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#eta [pb]");
	  h_P_tot->GetXaxis()->SetTitle("leading jet #eta");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_first_bjet_pt") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp^{b}_{T} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{b}_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 10);
	  }
	} else if (title_b=="w_first_bjet_eta") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#eta^{b} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("leading b-jet #eta");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta^{b} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_pt_Z_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp^{Z}_{T} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("Z boson p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{Z}_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 20);
	  }
	} else if (title_b=="w_Ht_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dH_{T} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("H_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dH_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_delta_phi_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#Delta#phi_{bZ} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta#phi(bZ) [rad]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta#phi_{Zb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	}

	if (plot) {
	  ofstream out, out1, out2;
	  if (isratio==0) {
	    gSystem->mkdir((path + "/combined/" + version + "/xsecs_unfolding" + dirbSel + "/").c_str(), kTRUE);
	    c1->SaveAs((path + "/combined/" + version + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.pdf").c_str());
	    out.open((path + "/combined/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.dat").c_str());
	    out1.open((path + "/combined/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.txt").c_str());
	    out2.open((path + "/combined/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.tex").c_str());
	  }
	  if (isratio==1) {
	    gSystem->mkdir((path + "/combined/" + version + "/ratios_unfolding" + dirbSel + "/").c_str(), kTRUE);
	    c1->SaveAs((path + "/combined/" + version + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.pdf").c_str());
	    out.open((path + "/combined/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.dat").c_str());
	    out1.open((path + "/combined/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.txt").c_str());
	    out2.open((path + "/combined/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.tex").c_str());
	  }
	  if (isratio==0) {
	    out << h_data->GetName();
	    out << endl;
	    out << std::setw(25) << "data";
	    out << std::setw(12) << "bkg";
	    out << std::setw(12) << "eff";
	    out << std::setw(12) << "jec";
	    out << std::setw(12) << "jer";
	    out << std::setw(12) << "pu";
	    out << std::setw(12) << "bkg";
	    out << std::setw(12) << "ttbar";
	    out << std::setw(12) << "bfit";
	    out << std::setw(12) << "btag";
	    out << std::setw(12) << "unfold";
	    out << std::setw(12) << "unfold";
	    out << std::setw(12) << "lumi";
	    out << std::setw(12) << "total";
	    out << std::setw(12) << "total";
	    out << std::setw(12) << "total";
	    out << endl;
	    out << std::setw(25) << "stat";
	    out << std::setw(12) << "stat";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "stat";
	    out << std::setw(12) << "stat";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "stat";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "stat";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "error";
	    out << std::setw(8) << "%";
	    out << endl;
	    for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	      out << std::fixed;
	      out << std::setw(2) << i;
	      out << " ";
	      out << std::setprecision(6);
	      out << std::setw(10) << h_data->GetBinContent(i);
	      out << " +- ";
	      out << std::setw(8) << h_data->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << stat_bkg->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_eff->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_jec->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_jer->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_pu->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_bkg->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << stat_top->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << stat_bfit->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_btag->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << stat_unfold->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_unfold->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << syst_lumi->GetBinError(i);
	      out << " => ";
	      out << std::setw(8) << h_data_stat->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << h_data_syst->GetBinError(i);
	      out << " => ";
	      out << std::setw(8) << h_data_tot->GetBinError(i);
	      out << " => ";
	      out << std::setprecision(1);
	      out << std::setw(4) << 100.*(h_data_stat->GetBinContent(i)==0 ? 0 : h_data_tot->GetBinError(i)/h_data_stat->GetBinContent(i));
	      out << endl;
	    }
	  }
	  out << h_data_b->GetName();
	  out << endl;
	  out << std::setw(25) << "data";
	  out << std::setw(12) << "bkg";
	  out << std::setw(12) << "eff";
	  out << std::setw(12) << "jec";
	  out << std::setw(12) << "jer";
	  out << std::setw(12) << "pu";
	  out << std::setw(12) << "bkg";
	  out << std::setw(12) << "ttbar";
	  out << std::setw(12) << "bfit";
	  out << std::setw(12) << "btag";
	  out << std::setw(12) << "unfold";
	  out << std::setw(12) << "unfold";
	  out << std::setw(12) << "lumi";
	  out << std::setw(12) << "total";
	  out << std::setw(12) << "total";
	  out << std::setw(12) << "total";
	  out << endl;
	  out << std::setw(25) << "stat";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "stat";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "error";
	  out << std::setw(8) << "%";
	  out << endl;
	  for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	    out << std::fixed;
	    out << std::setw(2) << i;
	    out << " ";
	    out << std::setprecision(6);
	    out << std::setw(10) << h_data_b->GetBinContent(i);
	    out << " +- ";
	    out << std::setw(8) << h_data_b->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_b_bkg->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_eff->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_jec->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_jer->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_pu->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_bkg->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_b_top->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_b_bfit->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_btag->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << stat_b_unfold->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_unfold->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << syst_b_lumi->GetBinError(i);
	    out << " => ";
	    out << std::setw(8) << h_data_b_stat->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << h_data_b_syst->GetBinError(i);
	    out << " => ";
	    out << std::setw(8) << h_data_b_tot->GetBinError(i);
	    out << " => ";
	    out << std::setprecision(1);
	    out << std::setw(4) << 100.*(h_data_b_stat->GetBinContent(i)==0 ? 0 : h_data_b_tot->GetBinError(i)/h_data_b_stat->GetBinContent(i));
	    out << endl;
	  }
	  out.close();
	  if (isratio==0) {
	    out1 << h_data->GetName() << " - RELATIVE ERRORS";
	    out1 << endl;
	    out1 << std::setw(7) << "data";
	    out1 << std::setw(8) << "bkg";
	    out1 << std::setw(8) << "eff";
	    out1 << std::setw(8) << "jec";
	    out1 << std::setw(8) << "jer";
	    out1 << std::setw(8) << "pu";
	    out1 << std::setw(8) << "bkg";
	    out1 << std::setw(8) << "ttbar";
	    out1 << std::setw(8) << "bfit";
	    out1 << std::setw(8) << "btag";
	    out1 << std::setw(8) << "unfold";
	    out1 << std::setw(8) << "unfold";
	    out1 << std::setw(8) << "lumi";
	    out1 << std::setw(8) << "total";
	    out1 << std::setw(8) << "total";
	    out1 << std::setw(8) << "total";
	    out1 << endl;
	    out1 << std::setw(7) << "stat";
	    out1 << std::setw(8) << "stat";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "stat";
	    out1 << std::setw(8) << "stat";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "stat";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "stat";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "error";
	    out1 << endl;
	    for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	      double val = 100.*(h_data->GetBinContent(i)==0 ? 0 : 1./h_data->GetBinContent(i));
	      out1 << std::fixed;
	      out1 << std::setw(2) << i;
	      out1 << " ";
	      out1 << std::setprecision(1);
	      out1 << std::setw(4) << h_data->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << stat_bkg->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_eff->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_jec->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_jer->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_pu->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_bkg->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << stat_top->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << stat_bfit->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_btag->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << stat_unfold->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_unfold->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << syst_lumi->GetBinError(i)*val;
	      out1 << " => ";
	      out1 << std::setw(4) << h_data_stat->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << h_data_syst->GetBinError(i)*val;
	      out1 << " => ";
	      out1 << std::setw(4) << h_data_tot->GetBinError(i)*val;
	      out1 << endl;
	    }
	  }
	  out1 << h_data_b->GetName() << " - RELATIVE ERRORS";
	  out1 << endl;
	  out1 << std::setw(7) << "data";
	  out1 << std::setw(8) << "bkg";
	  out1 << std::setw(8) << "eff";
	  out1 << std::setw(8) << "jec";
	  out1 << std::setw(8) << "jer";
	  out1 << std::setw(8) << "pu";
	  out1 << std::setw(8) << "bkg";
	  out1 << std::setw(8) << "ttbar";
	  out1 << std::setw(8) << "bfit";
	  out1 << std::setw(8) << "btag";
	  out1 << std::setw(8) << "unfold";
	  out1 << std::setw(8) << "unfold";
	  out1 << std::setw(8) << "lumi";
	  out1 << std::setw(8) << "total";
	  out1 << std::setw(8) << "total";
	  out1 << std::setw(8) << "total";
	  out1 << endl;
	  out1 << std::setw(7) << "stat";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "error";
	  out1 << endl;
	  for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	    double val = 100.*(h_data_b->GetBinContent(i)==0 ? 0 : 1./h_data_b->GetBinContent(i));
	    out1 << std::fixed;
	    out1 << std::setw(2) << i;
	    out1 << " ";
	    out1 << std::setprecision(1);
	    out1 << std::setw(4) << h_data_b->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_b_bkg->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_eff->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_jec->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_jer->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_pu->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_bkg->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_b_top->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_b_bfit->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_btag->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << stat_b_unfold->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_unfold->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << syst_b_lumi->GetBinError(i)*val;
	    out1 << " => ";
	    out1 << std::setw(4) << h_data_b_stat->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(4) << h_data_b_syst->GetBinError(i)*val;
	    out1 << " => ";
	    out1 << std::setw(4) << h_data_b_tot->GetBinError(i)*val;
	    out1 << endl;
	  }
	  out1.close();
	  if (isratio==0) {
	    //out2 << h_data->GetName() << " - RELATIVE ERRORS";
	    //out2 << endl;
	    out2 << std::setw(7) << "\\textbf{data} &"  ;
	    out2 << std::setw(8) << "\\textbf{bkg} &"   ;
	    out2 << std::setw(8) << "\\textbf{eff} &"   ;
	    out2 << std::setw(8) << "\\textbf{jec} &"   ;
	    out2 << std::setw(8) << "\\textbf{jer} &"   ;
	    out2 << std::setw(8) << "\\textbf{pu} &"    ;
	    out2 << std::setw(8) << "\\textbf{bkg} &"   ;
	    out2 << std::setw(8) << "\\textbf{ttbar} &" ;
	    out2 << std::setw(8) << "\\textbf{bfit} &"  ;
	    out2 << std::setw(8) << "\\textbf{btag} &"  ;
	    out2 << std::setw(8) << "\\textbf{unfold} &";
	    out2 << std::setw(8) << "\\textbf{unfold} &";
	    out2 << std::setw(8) << "\\textbf{lumi} &"  ;
	    out2 << std::setw(8) << "\\textbf{total} &" ;
	    out2 << std::setw(8) << "\\textbf{total} &" ;
	    out2 << std::setw(8) << "\\textbf{total} &" ;
	    out2 << endl;
	    out2 << std::setw(7) << "\\textbf{stat} & ";
	    out2 << std::setw(8) << "\\textbf{stat} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{stat} & ";
	    out2 << std::setw(8) << "\\textbf{stat} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{stat} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{stat} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{error & ";
	    out2 << endl;
	    for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	      double val = 100.*(h_data->GetBinContent(i)==0 ? 0 : 1./h_data->GetBinContent(i));
	      out2 << std::fixed;
	      out2 << std::setw(2) << i;
	      out2 << " ";
	      out2 << std::setprecision(1);
	      out2 << std::setw(4) << h_data->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << stat_bkg->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << syst_eff->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << syst_jec->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << syst_jer->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << syst_pu->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << syst_bkg->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << stat_top->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << stat_bfit->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << syst_btag->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << stat_unfold->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << syst_unfold->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << syst_lumi->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << h_data_stat->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << h_data_syst->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << h_data_tot->GetBinError(i)*val;
	      out2 << std::setw(4) << "\\tabularnewline" << "   " << "\\hline";
	      out2 << endl;
	    }
	  }
	  //out2 << h_data_b->GetName() << " - RELATIVE ERRORS";
	  //out2 << endl;
	  out2 << std::setw(7) << "\\textbf{data} &";
	  out2 << std::setw(8) << "\\textbf{bkg} &";
	  out2 << std::setw(8) << "\\textbf{eff} &";
	  out2 << std::setw(8) << "\\textbf{jec} &";
	  out2 << std::setw(8) << "\\textbf{jer} &";
	  out2 << std::setw(8) << "\\textbf{pu} &";
	  out2 << std::setw(8) << "\\textbf{bkg} &";
	  out2 << std::setw(8) << "\\textbf{ttbar} &";
	  out2 << std::setw(8) << "\\textbf{bfit} &";
	  out2 << std::setw(8) << "\\textbf{btag} &";
	  out2 << std::setw(8) << "\\textbf{unfold} &";
	  out2 << std::setw(8) << "\\textbf{unfold} &";
	  out2 << std::setw(8) << "\\textbf{lumi} &";
	  out2 << std::setw(8) << "\\textbf{total} &";
	  out2 << std::setw(8) << "\\textbf{total} &";
	  out2 << std::setw(8) << "\\textbf{total} &";
	  out2 << endl;
	  out2 << std::setw(7) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{error} &";
	  out2 << endl;
	  for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	    double val = 100.*(h_data_b->GetBinContent(i)==0 ? 0 : 1./h_data_b->GetBinContent(i));
	    out2 << std::fixed;
	    out2 << std::setw(2) << i;
	    out2 << " ";
	    out2 << std::setprecision(1);
	    out2 << std::setw(4) << h_data_b->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << stat_b_bkg->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_eff->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_jec->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_jer->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_pu->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_bkg->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << stat_b_top->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << stat_b_bfit->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_btag->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << stat_b_unfold->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_unfold->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << syst_b_lumi->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << h_data_b_stat->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << h_data_b_syst->GetBinError(i)*val;
	    out2 << " & ";
	    out2 << std::setw(4) << h_data_b_tot->GetBinError(i)*val;
	    out2 << std::setw(4) << "\\tabularnewline" << "   " << "\\hline";
	    out2 << endl;
	  }
           
	}
}

