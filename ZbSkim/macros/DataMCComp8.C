#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v14.h"

#include "fixrange.C"
#include "rebin.C"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data";
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

double calc(int iflag, double cont1, double cont2, double stat1, double stat2, double stat_bkg1, double stat_bkg2, double syst_eff1, double syst_eff2, double syst_jec1, double syst_jec2, double syst_jer1, double syst_jer2, double syst_pu1, double syst_pu2, double syst_bkg1, double syst_bkg2, double stat_top1, double stat_top2, double stat_bfit1, double stat_bfit2, double syst_bfit1, double syst_bfit2, double syst_btag1, double syst_btag2, double stat_unfold1, double stat_unfold2, double syst_unfold1, double syst_unfold2, double syst_lumi1, double syst_lumi2) {

  if (cont1*cont2==0) return 0.0;

  double tmp1 = TMath::Power(stat1,2);
  tmp1 = tmp1+TMath::Power(stat_bkg1,2);
  tmp1 = tmp1+TMath::Power(syst_eff1,2);
  tmp1 = tmp1+TMath::Power(stat_top1,2);
  tmp1 = tmp1+TMath::Power(stat_bfit1,2);
  tmp1 = tmp1+TMath::Power(stat_unfold1,2);

  double tmp2 = TMath::Power(stat2,2);
  tmp2 = tmp2+TMath::Power(stat_bkg2,2);
  tmp2 = tmp2+TMath::Power(syst_eff2,2);
  tmp2 = tmp2+TMath::Power(stat_top2,2);
  tmp2 = tmp2+TMath::Power(stat_bfit2,2);
  tmp2 = tmp2+TMath::Power(stat_unfold2,2);

  double val0 = (cont1/tmp1+cont2/tmp2)/(1.0/tmp1+1.0/tmp2);

  double val1 = TMath::Sqrt(1.0/(1.0/tmp1+1.0/tmp2));
  val1 = TMath::Sqrt(TMath::Power(val1,2)+TMath::Power((syst_jec1+syst_jec2)/2.,2));
  val1 = TMath::Sqrt(TMath::Power(val1,2)+TMath::Power((syst_jer1+syst_jer2)/2.,2));
  val1 = TMath::Sqrt(TMath::Power(val1,2)+TMath::Power((syst_pu1+syst_pu2)/2.,2));
  val1 = TMath::Sqrt(TMath::Power(val1,2)+TMath::Power((syst_bkg1+syst_bkg2)/2.,2));
  val1 = TMath::Sqrt(TMath::Power(val1,2)+TMath::Power((syst_bfit1+syst_bfit2)/2.,2));
  val1 = TMath::Sqrt(TMath::Power(val1,2)+TMath::Power((syst_unfold1+syst_unfold2)/2.,2));
  val1 = TMath::Sqrt(TMath::Power(val1,2)+TMath::Power(val0*(TMath::Abs(syst_btag1/cont1)+TMath::Abs(syst_btag2/cont2))/2.0,2));
  val1 = TMath::Sqrt(TMath::Power(val1,2)+TMath::Power(val0*(TMath::Abs(syst_lumi1/cont1)+TMath::Abs(syst_lumi2/cont2))/2.0,2));

  double val = 0.0;
  if (iflag == 0) val = val0;
  if (iflag == 1) val = val1;

  return val;
}

void DataMCComp8(string title="", int plot=0, int isratio=1, int numB=0) {

int useSherpa=0;
//int useSherpa=1; // use Sherpa MC prediction

//int useNewPowheg=0;
int useNewPowheg=1; // use new Powheg MC prediction

//int drawInclusive=0;
int drawInclusive=1; // do plot the "inclusive" histogram

string subdir="0";
string postfix="";
string dirbSel="";
string bSel="";

if (numB==1) {
  postfix = postfix + "1b";
  dirbSel = "_1b";
  bSel = "Z + (= 1) b-jet";
  drawInclusive = 0;
}
if (numB==2) {
  postfix = postfix + "2b";
  dirbSel = "_2b";
  bSel = "Z + (#geq 2) b-jet";
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
	  norm1_2 = ((Lumi2012 * Xsec_dy_4) / ((Ngen_dy_4_ee+Ngen_dy_4_mm)/2.));
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

	h_mcg = fixrange(h_mcg, numB);
	h_mcg_b = fixrange(h_mcg_b, numB);
	h_mcg1 = fixrange(h_mcg1, numB);
	h_mcg1_b = fixrange(h_mcg1_b, numB);
	h_mcg2 = fixrange(h_mcg2, numB);
	h_mcg3 = fixrange(h_mcg3, numB);
	h_mcg2_b = fixrange(h_mcg2_b, numB);
	h_mcg3_b = fixrange(h_mcg3_b, numB);

	h_mcg = rebin(h_mcg, numB);
	h_mcg_b = rebin(h_mcg_b, numB);
	h_mcg1 = rebin(h_mcg1, numB);
	h_mcg1_b = rebin(h_mcg1_b, numB);
	h_mcg2 = rebin(h_mcg2, numB);
	h_mcg3 = rebin(h_mcg3, numB);
	h_mcg2_b = rebin(h_mcg2_b, numB);
	h_mcg3_b = rebin(h_mcg3_b, numB);

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

	TH1F* w_syst_bfit2[2];
	TH1F* w_syst_b_bfit2[2];

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

        double xsec_data[2];
        double xsec_data_b[2];

        double xsec_stat_data[2];
        double xsec_stat_data_b[2];

        double xsec_stat_bkg[2];
        double xsec_stat_b_bkg[2];

        double xsec_syst_eff[2];
        double xsec_syst_b_eff[2];

        double xsec_syst_jec[2];
        double xsec_syst_b_jec[2];

        double xsec_syst_jer[2];
        double xsec_syst_b_jer[2];

        double xsec_syst_pu[2];
        double xsec_syst_b_pu[2];

        double xsec_syst_bkg[2];
        double xsec_syst_b_bkg[2];

        double xsec_stat_top[2];
        double xsec_stat_b_top[2];

        double xsec_stat_bfit[2];
        double xsec_stat_b_bfit[2];

        double xsec_syst_bfit2[2];
        double xsec_syst_b_bfit2[2];

        double xsec_syst_btag[2];
        double xsec_syst_b_btag[2];

        double xsec_stat_unfold[2];
        double xsec_stat_b_unfold[2];

        double xsec_syst_unfold[2];
        double xsec_syst_b_unfold[2];

        double xsec_syst_lumi[2];
        double xsec_syst_b_lumi[2];

        double xsec_stat_tot[2];
        double xsec_stat_b_tot[2];

        double xsec_syst_tot[2];
        double xsec_syst_b_tot[2];

	for (int i=0; i<2; i++) {

	  w_stat_bkg[i] = (TH1F*)w_data[i]->Clone();
	  w_stat_b_bkg[i] = (TH1F*)w_data_b[i]->Clone();

	  w_syst_eff[i] = (TH1F*)w_data[i]->Clone();
	  w_syst_b_eff[i] = (TH1F*)w_data_b[i]->Clone();

	  w_syst_jec[i] = (TH1F*)w_data[i]->Clone();
	  w_syst_b_jec[i] = (TH1F*)w_data_b[i]->Clone();

	  w_syst_jer[i] = (TH1F*)w_data[i]->Clone();
	  w_syst_b_jer[i] = (TH1F*)w_data_b[i]->Clone();

	  w_syst_pu[i] = (TH1F*)w_data[i]->Clone();
	  w_syst_b_pu[i] = (TH1F*)w_data_b[i]->Clone();

	  w_syst_bkg[i] = (TH1F*)w_data[i]->Clone();
	  w_syst_b_bkg[i] = (TH1F*)w_data_b[i]->Clone();

	  w_stat_top[i] = (TH1F*)w_data[i]->Clone();
	  w_stat_b_top[i] = (TH1F*)w_data_b[i]->Clone();

	  w_stat_bfit[i] = (TH1F*)w_data[i]->Clone();
	  w_stat_b_bfit[i] = (TH1F*)w_data_b[i]->Clone();

	  w_syst_bfit2[i] = (TH1F*)w_data[i]->Clone();
	  w_syst_b_bfit2[i] = (TH1F*)w_data_b[i]->Clone();

	  w_syst_btag[i] = (TH1F*)w_data[i]->Clone();
	  w_syst_b_btag[i] = (TH1F*)w_data_b[i]->Clone();

	  w_stat_unfold[i] = (TH1F*)w_data[i]->Clone();
	  w_stat_b_unfold[i] = (TH1F*)w_data_b[i]->Clone();

	  w_syst_unfold[i] = (TH1F*)w_data[i]->Clone();
	  w_syst_b_unfold[i] = (TH1F*)w_data_b[i]->Clone();

	  w_syst_lumi[i] = (TH1F*)w_data[i]->Clone();
	  w_syst_b_lumi[i] = (TH1F*)w_data_b[i]->Clone();

	  w_stat_tot[i] = (TH1F*)w_data[i]->Clone();
	  w_stat_b_tot[i] = (TH1F*)w_data_b[i]->Clone();

	  w_syst_tot[i] = (TH1F*)w_data[i]->Clone();
	  w_syst_b_tot[i] = (TH1F*)w_data_b[i]->Clone();

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

	  if (isratio==0 && drawInclusive) {
	    getline(in, tmp);
	    getline(in, tmp);
	    getline(in, tmp);
	    for (int j=0; j<=w_data[0]->GetNbinsX()+1; j++) {
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
	      in >> val; w_syst_bfit2[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_btag[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_stat_unfold[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_unfold[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_lumi[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_stat_tot[i]->SetBinError(j, val); in >> tmp;
	      in >> val; w_syst_tot[i]->SetBinError(j, val); in >> tmp;
	      in >> val; in >> tmp; in >> val;
	      in.ignore();
	    }
            in >> tmp;
            in >> xsec_data[i]; in >> tmp;
            in >> xsec_stat_data[i]; in >> tmp;
            in >> xsec_stat_bkg[i]; in >> tmp;
            in >> xsec_syst_eff[i]; in >> tmp;
            in >> xsec_syst_jec[i]; in >> tmp;
            in >> xsec_syst_jer[i]; in >> tmp;
            in >> xsec_syst_pu[i]; in >> tmp;
            in >> xsec_syst_bkg[i]; in >> tmp;
            in >> xsec_stat_top[i]; in >> tmp;
            in >> xsec_stat_bfit[i]; in >> tmp;
            in >> xsec_syst_bfit2[i]; in >> tmp;
            in >> xsec_syst_btag[i]; in >> tmp;
            in >> xsec_stat_unfold[i]; in >> tmp;
            in >> xsec_syst_unfold[i]; in >> tmp;
            in >> xsec_syst_lumi[i]; in >> tmp;
            in >> xsec_stat_tot[i]; in >> tmp;
            in >> xsec_syst_tot[i]; in >> tmp;
            in >> tmp; in >> tmp; in >> tmp;
            in.ignore();
	  }

	  getline(in, tmp);
	  getline(in, tmp);
	  getline(in, tmp);
	  for (int j=0; j<=w_data_b[0]->GetNbinsX()+1; j++) {
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
	    in >> val; w_syst_b_bfit2[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_btag[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_b_unfold[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_unfold[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_lumi[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_stat_b_tot[i]->SetBinError(j, val); in >> tmp;
	    in >> val; w_syst_b_tot[i]->SetBinError(j, val); in >> tmp;
	    in >> val; in >> tmp; in >> val;
	    in.ignore();
	  }
          in >> tmp;
          in >> xsec_data_b[i]; in >> tmp;
          in >> xsec_stat_data_b[i]; in >> tmp;
          in >> xsec_stat_b_bkg[i]; in >> tmp;
          in >> xsec_syst_b_eff[i]; in >> tmp;
          in >> xsec_syst_b_jec[i]; in >> tmp;
          in >> xsec_syst_b_jer[i]; in >> tmp;
          in >> xsec_syst_b_pu[i]; in >> tmp;
          in >> xsec_syst_b_bkg[i]; in >> tmp;
          in >> xsec_stat_b_top[i]; in >> tmp;
          in >> xsec_stat_b_bfit[i]; in >> tmp;
          in >> xsec_syst_b_bfit2[i]; in >> tmp;
          in >> xsec_syst_b_btag[i]; in >> tmp;
          in >> xsec_stat_b_unfold[i]; in >> tmp;
          in >> xsec_syst_b_unfold[i]; in >> tmp;
          in >> xsec_syst_b_lumi[i]; in >> tmp;
          in >> xsec_stat_b_tot[i]; in >> tmp;
          in >> xsec_syst_b_tot[i]; in >> tmp;
          in >> tmp; in >> tmp; in >> tmp;
          in.ignore();

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

	TH1F* syst_bfit2 = (TH1F*)w_data[0]->Clone();
	TH1F* syst_b_bfit2 = (TH1F*)w_data_b[0]->Clone();

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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			1.1*w_syst_bfit2[0]->GetBinError(i), 1.1*w_syst_bfit2[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			w_syst_lumi[0]->GetBinError(i), w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_bfit2->SetBinError(i, val);
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
	                w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
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
			w_syst_bfit2[0]->GetBinError(i), w_syst_bfit2[1]->GetBinError(i),
			w_syst_btag[0]->GetBinError(i), w_syst_btag[1]->GetBinError(i),
			w_stat_unfold[0]->GetBinError(i), w_stat_unfold[1]->GetBinError(i),
			w_syst_unfold[0]->GetBinError(i), w_syst_unfold[1]->GetBinError(i),
			1.1*w_syst_lumi[0]->GetBinError(i), 1.1*w_syst_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_lumi->SetBinError(i, val);

	  val = 0.0;
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(h_data->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_top->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_bfit->GetBinError(i),2));
	  h_data_stat->SetBinError(i, val);
	  val = 0.0;
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_bkg->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_eff->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_jec->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_jer->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_pu->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_bkg->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_bfit2->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_btag->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_unfold->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_unfold->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_lumi->GetBinError(i),2));
	  h_data_syst->SetBinError(i, val);
	  val = 0.0;
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(h_data_stat->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(h_data_syst->GetBinError(i),2));
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			1.1*w_syst_b_bfit2[0]->GetBinError(i), 1.1*w_syst_b_bfit2[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			w_syst_b_lumi[0]->GetBinError(i), w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_bfit2->SetBinError(i, val);
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
	                w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
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
			w_syst_b_bfit2[0]->GetBinError(i), w_syst_b_bfit2[1]->GetBinError(i),
			w_syst_b_btag[0]->GetBinError(i), w_syst_b_btag[1]->GetBinError(i),
			w_stat_b_unfold[0]->GetBinError(i), w_stat_b_unfold[1]->GetBinError(i),
			w_syst_b_unfold[0]->GetBinError(i), w_syst_b_unfold[1]->GetBinError(i),
			1.1*w_syst_b_lumi[0]->GetBinError(i), 1.1*w_syst_b_lumi[1]->GetBinError(i));
	  val = TMath::Sqrt((TMath::Power(val,2)-TMath::Power(ref,2))/(TMath::Power(1.1,2)-1));
	  syst_b_lumi->SetBinError(i, val);

	  val = 0.0;
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(h_data_b->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_b_top->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_b_bfit->GetBinError(i),2));
	  h_data_b_stat->SetBinError(i, val);
	  val = 0.0;
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_b_bkg->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_eff->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_jec->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_jer->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_pu->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_bkg->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_bfit2->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_btag->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_b_unfold->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_unfold->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_lumi->GetBinError(i),2));
	  h_data_b_syst->SetBinError(i, val);
	  val = 0.0;
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(h_data_b_stat->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(h_data_b_syst->GetBinError(i),2));
	  h_data_b_tot->SetBinError(i, val);
	}

        double xsec_tot_data;
        double xsec_tot_data_b;

        double xsec_tot_stat_data;
        double xsec_tot_stat_data_b;

        double xsec_tot_stat_bkg;
        double xsec_tot_stat_b_bkg;

        double xsec_tot_syst_eff;
        double xsec_tot_syst_b_eff;

        double xsec_tot_syst_jec;
        double xsec_tot_syst_b_jec;

        double xsec_tot_syst_jer;
        double xsec_tot_syst_b_jer;

        double xsec_tot_syst_pu;
        double xsec_tot_syst_b_pu;

        double xsec_tot_syst_bkg;
        double xsec_tot_syst_b_bkg;

        double xsec_tot_stat_top;
        double xsec_tot_stat_b_top;

        double xsec_tot_stat_bfit;
        double xsec_tot_stat_b_bfit;

        double xsec_tot_syst_bfit2;
        double xsec_tot_syst_b_bfit2;

        double xsec_tot_syst_btag;
        double xsec_tot_syst_b_btag;

        double xsec_tot_stat_unfold;
        double xsec_tot_stat_b_unfold;

        double xsec_tot_syst_unfold;
        double xsec_tot_syst_b_unfold;

        double xsec_tot_syst_lumi;
        double xsec_tot_syst_b_lumi;

        double xsec_tot_stat_tot;
        double xsec_tot_stat_b_tot;

        double xsec_tot_syst_tot;
        double xsec_tot_syst_b_tot;

        double xsec_tot_data_tot;
        double xsec_tot_data_b_tot;

        double xval = 0.0;
        xval = calc(0, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xsec_tot_data = xval;

        double xref = 0.0;
        xref = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = calc(1, xsec_data[0], xsec_data[1],
                       1.1*xsec_stat_data[0], 1.1*xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_stat_data = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       1.1*xsec_stat_bkg[0], 1.1*xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_stat_bkg = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       1.1*xsec_syst_eff[0], 1.1*xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_eff = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       1.1*xsec_syst_jer[0], 1.1*xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_jer = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       1.1*xsec_syst_jec[0], 1.1*xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_jec = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       1.1*xsec_syst_pu[0], 1.1*xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_pu = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       1.1*xsec_syst_bkg[0], 1.1*xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_bkg = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       1.1*xsec_stat_top[0], 1.1*xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_stat_top = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       1.1*xsec_stat_bfit[0], 1.1*xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_stat_bfit = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       1.1*xsec_syst_bfit2[0], 1.1*xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_bfit2 = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       1.1*xsec_syst_btag[0], 1.1*xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_btag = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       1.1*xsec_stat_unfold[0], 1.1*xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_stat_unfold = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       1.1*xsec_syst_unfold[0], 1.1*xsec_syst_unfold[1],
                       xsec_syst_lumi[0], xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_unfold = xval;
        xval = calc(1, xsec_data[0], xsec_data[1],
                       xsec_stat_data[0], xsec_stat_data[1],
                       xsec_stat_bkg[0], xsec_stat_bkg[1],
                       xsec_syst_eff[0], xsec_syst_eff[1],
                       xsec_syst_jer[0], xsec_syst_jer[1],
                       xsec_syst_jec[0], xsec_syst_jec[1],
                       xsec_syst_pu[0], xsec_syst_pu[1],
                       xsec_syst_bkg[0], xsec_syst_bkg[1],
                       xsec_stat_top[0], xsec_stat_top[1],
                       xsec_stat_bfit[0], xsec_stat_bfit[1],
                       xsec_syst_bfit2[0], xsec_syst_bfit2[1],
                       xsec_syst_btag[0], xsec_syst_btag[1],
                       xsec_stat_unfold[0], xsec_stat_unfold[1],
                       xsec_syst_unfold[0], xsec_syst_unfold[1],
                       1.1*xsec_syst_lumi[0], 1.1*xsec_syst_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_lumi = xval;

	xval = 0.0;
        xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_stat_data,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_stat_top,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_stat_bfit,2));
        xsec_tot_stat_tot = xval;
	xval = 0.0;
        xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_stat_bkg,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_eff,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_jec,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_jer,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_pu,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_bkg,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_bfit2,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_btag,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_stat_unfold,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_unfold,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_lumi,2));
        xsec_tot_syst_tot = xval;
	xval = 0.0;
        xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_stat_tot,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_tot,2));
        xsec_tot_data_tot = xval;

        xval = 0.0;
        xval = calc(0, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xsec_tot_data_b = xval;

        xref = 0.0;
        xref = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       1.1*xsec_stat_data_b[0], 1.1*xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_stat_data_b = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       1.1*xsec_stat_b_bkg[0], 1.1*xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_stat_b_bkg = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       1.1*xsec_syst_b_eff[0], 1.1*xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_b_eff = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       1.1*xsec_syst_b_jer[0], 1.1*xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_b_jer = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       1.1*xsec_syst_b_jec[0], 1.1*xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_b_jec = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       1.1*xsec_syst_b_pu[0], 1.1*xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_b_pu = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       1.1*xsec_syst_b_bkg[0], 1.1*xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_b_bkg = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       1.1*xsec_stat_b_top[0], 1.1*xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_stat_b_top = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       1.1*xsec_stat_b_bfit[0], 1.1*xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_stat_b_bfit = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       1.1*xsec_syst_b_bfit2[0], 1.1*xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_b_bfit2 = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       1.1*xsec_syst_b_btag[0], 1.1*xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_b_btag = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       1.1*xsec_stat_b_unfold[0], 1.1*xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_stat_b_unfold = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       1.1*xsec_syst_b_unfold[0], 1.1*xsec_syst_b_unfold[1],
                       xsec_syst_b_lumi[0], xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_b_unfold = xval;
        xval = calc(1, xsec_data_b[0], xsec_data_b[1],
                       xsec_stat_data_b[0], xsec_stat_data_b[1],
                       xsec_stat_b_bkg[0], xsec_stat_b_bkg[1],
                       xsec_syst_b_eff[0], xsec_syst_b_eff[1],
                       xsec_syst_b_jer[0], xsec_syst_b_jer[1],
                       xsec_syst_b_jec[0], xsec_syst_b_jec[1],
                       xsec_syst_b_pu[0], xsec_syst_b_pu[1],
                       xsec_syst_b_bkg[0], xsec_syst_b_bkg[1],
                       xsec_stat_b_top[0], xsec_stat_b_top[1],
                       xsec_stat_b_bfit[0], xsec_stat_b_bfit[1],
                       xsec_syst_b_bfit2[0], xsec_syst_b_bfit2[1],
                       xsec_syst_b_btag[0], xsec_syst_b_btag[1],
                       xsec_stat_b_unfold[0], xsec_stat_b_unfold[1],
                       xsec_syst_b_unfold[0], xsec_syst_b_unfold[1],
                       1.1*xsec_syst_b_lumi[0], 1.1*xsec_syst_b_lumi[1]);
        xval = TMath::Sqrt((TMath::Power(xval,2)-TMath::Power(xref,2))/(TMath::Power(1.1,2)-1));
        xsec_tot_syst_b_lumi = xval;

	xval = 0.0;
        xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_stat_data_b,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_stat_b_top,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_stat_b_bfit,2));
        xsec_tot_stat_b_tot = xval;
	xval = 0.0;
        xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_stat_b_bkg,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_b_eff,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_b_jec,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_b_jer,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_b_pu,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_b_bkg,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_b_bfit2,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_b_btag,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_stat_b_unfold,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_b_unfold,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_b_lumi,2));
        xsec_tot_syst_b_tot = xval;
	xval = 0.0;
        xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_stat_b_tot,2));
	xval = TMath::Sqrt(TMath::Power(xval,2)+TMath::Power(xsec_tot_syst_b_tot,2));
        xsec_tot_data_b_tot = xval;

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
	  TH1F* tmp2 = (TH1F*)h_mcg_b->Clone();
	  tmp2->SetFillColor(0);
	  tmp2->Draw("HISTLSAME");
	}

	h_mcg1_b->SetLineColor(kMagenta-6);
	h_mcg1_b->SetLineWidth(2);
	h_mcg1_b->SetFillColor(kMagenta-6);
	h_mcg1_b->SetMarkerColor(kMagenta-6);
	if (isratio==1) {
	  h_mcg1_b->Draw("E5SAME");
	  TH1F* tmp2_1 = (TH1F*)h_mcg1_b->Clone();
	  tmp2_1->SetFillColor(0);
	  tmp2_1->Draw("HISTLSAME");
	}

	h_mcg2_b->SetLineColor(kBlue-4);
	h_mcg2_b->SetLineWidth(2);
	h_mcg2_b->SetFillColor(kBlue-4);
	h_mcg2_b->SetMarkerColor(kBlue-4);
	if (isratio==1) {
	  h_mcg2_b->Draw("E5SAME");
	  TH1F* tmp2_2 = (TH1F*)h_mcg2_b->Clone();
	  tmp2_2->SetFillColor(0);
	  tmp2_2->Draw("HISTLSAME");
	}

	h_mcg3_b->SetLineColor(kOrange+7);
	h_mcg3_b->SetLineWidth(2);
	h_mcg3_b->SetFillColor(kOrange+7);
	h_mcg3_b->SetMarkerColor(kOrange+7);
	if (isratio==1) {
	  h_mcg3_b->Draw("E5SAME");
	  TH1F* tmp2_3 = (TH1F*)h_mcg3_b->Clone();
	  tmp2_3->SetFillColor(0);
	  tmp2_3->Draw("HISTLSAME");
	}

	if (isratio==1) {
	  h_data_b_tot->GetYaxis()->SetTitle("#sigma_{Z+b-jets}/#sigma_{Z+jets} [%]");
	}
	h_data_b_tot->GetYaxis()->SetTitleOffset(1.2);
	h_data_b_tot->GetXaxis()->SetTitleOffset(1.3);
	h_data_b_tot->SetMarkerColor(kRed+1);
	h_data_b_tot->SetLineColor(kRed+1);
	//h_data_b_tot->SetMarkerSize(0.7);
	h_data_b_tot->SetStats(0);
	if (isratio==0) {
	  h_data_b_tot->SetMarkerStyle(24);
	  h_data_b_tot->SetMarkerSize(0.7);
	}
	if (isratio==1) {
	  h_data_b_tot->SetMarkerStyle(22);
	  h_data_b_tot->SetMarkerSize(0.9);
	}
	h_data_b_stat->GetYaxis()->SetTitleOffset(1.2);
	h_data_b_stat->GetXaxis()->SetTitleOffset(1.3);
	h_data_b_stat->SetMarkerColor(kBlack);
	h_data_b_stat->SetLineColor(kBlack);
	//h_data_b_stat->SetMarkerSize(0.7);
	h_data_b_stat->SetStats(0);
	if (isratio==0) {
	  h_data_b_stat->SetMarkerStyle(24);
	  h_data_b_stat->SetMarkerSize(0.7);
	}
	if (isratio==1) {
	  h_data_b_stat->SetMarkerStyle(22);
	  h_data_b_stat->SetMarkerSize(0.9);
	}
	if (isratio==1) {
	  h_data_b_tot->Draw("E1PX0SAME");
	  h_data_b_stat->Draw("E1PX0SAME");
	}

	TLegend *leg = NULL;
	if (isratio==0) {
	  leg = new TLegend(0.64, 0.590, 0.88, 0.88);
	}
	if (isratio) {
	  leg = new TLegend(0.52, 0.510, 0.90, 0.88);
	}
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (isratio==0) {
	  pad1->SetLogy();

	  if (title=="w_delta_phi" || title_b=="w_first_bjet_eta" || title_b=="w_first_bjet_eta_abs") {
	    h_mcg_b->SetMaximum(18*h_data_tot->GetMaximum());
	  } else {
	    h_mcg_b->SetMaximum(4*h_data_tot->GetMaximum());
	  }
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

	  if (drawInclusive) leg->AddEntry(h_data_stat,"Z(#rightarrow ll)+j DATA","p");
	  leg->AddEntry(h_data_b_stat,"Z(#rightarrow ll)+b DATA","p");
	  leg->AddEntry(h_mcg,"Z(#rightarrow ll) MadGraph 5FS","l");
	  leg->AddEntry(h_mcg3,"Z(#rightarrow ll) MadGraph 4FS","l");
	  if (useSherpa) leg->AddEntry(h_mcg1,"Z(#rightarrow ll) Sherpa","l");
	  leg->AddEntry(h_mcg2,"Z(#rightarrow ll) Powheg","l");
	}

	if (isratio==1) {
	  leg->AddEntry(h_data_b_stat,"[Z(#rightarrow ll)+b] / [Z(#rightarrow ll)+j] DATA","p");
	  leg->AddEntry(h_mcg_b,"[Z(#rightarrow ll)+b] / [Z(#rightarrow ll)+j] MadGraph 5FS","l");
	  leg->AddEntry(h_mcg3_b,"[Z(#rightarrow ll)+b] / [Z(#rightarrow ll)+j] MadGraph 4FS","l");
	  if (useSherpa) leg->AddEntry(h_mcg1_b,"[Z(#rightarrow ll)+b] / [Z(#rightarrow ll)+j] Sherpa","l");
	  leg->AddEntry(h_mcg2_b,"[Z(#rightarrow ll)+b] / [Z(#rightarrow ll)+j] Powheg","l");
	}

	leg->Draw();

	c1->cd();

 	TLatex *latexLabel = 0;

	if (isratio==0) {
	  if (title_b=="w_Ht_b" || title_b=="w_first_bjet_pt" || title_b=="w_pt_Z_b" || title_b=="w_DR_bb" || title_b=="w_bb_mass" || title_b=="w_Zbb_mass"|| title_b=="w_DR_Zb_min"|| title_b=="w_DR_Zb_max"|| title_b=="w_A_Zb") {
	    latexLabel = CMSPrel2 (Lumi2012/1000., "Z/#gamma*#rightarrow ll selection", 0, 0.135, 0.51);
	  }
	  if (title_b=="w_delta_phi_b" || title_b=="w_delta_phi_2b") {
	    latexLabel = CMSPrel2 (Lumi2012/1000., "Z/#gamma*#rightarrow ll selection", 0, 0.68, 0.51);
	  }
	  if (title_b=="w_first_bjet_eta" || title_b=="w_first_bjet_eta_abs") {
	    latexLabel = CMSPrel2 (Lumi2012/1000., "Z/#gamma*#rightarrow ll selection", 0, 0.68, 0.51);
	  }
	}
	if (isratio==1) {
	  latexLabel = CMSPrel2 (Lumi2012/1000., "Z/#gamma*#rightarrow ll selection", 0, 0.135, 0.85);
	}

	if (latexLabel) latexLabel->Draw("same");

	TPad *pad2 = new TPad("pad2","pad2",0,0.29,1,0.4);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.001);
	pad2->Draw();
	pad2->cd();

	TH1F *h_M_tot = (TH1F*)h_mcg_b->Clone();
	TH1F *h_M_stat = (TH1F*)h_mcg_b->Clone();

	h_M_tot->Divide(h_data_b_tot);
	h_M_stat->Divide(h_data_b_stat);

	for (int i=0;i<=h_M_tot->GetNbinsX()+1;i++) {
	  if (h_M_tot->GetBinContent(i)<=0) h_M_tot->SetBinContent(i,999.);
	  if (h_M_stat->GetBinContent(i)<=0) h_M_stat->SetBinContent(i,999.);
	}

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
	//h_M_tot->SetMarkerSize(0.7);
	h_M_stat->GetXaxis()->SetTitleOffset(0.7);
	h_M_stat->SetMarkerColor(kBlack);
	h_M_stat->SetLineColor(kBlack);
	h_M_stat->SetLineWidth(1);
	//h_M_stat->SetMarkerSize(0.7);

	if (isratio==0) {
	  h_M_tot->SetMarkerStyle(24);
	  h_M_tot->SetMarkerSize(0.7);
	}
	if (isratio==1) {
	  h_M_tot->SetMarkerStyle(22);
	  h_M_tot->SetMarkerSize(0.9);
	}
	h_M_tot->Draw("E1PX0");
	h_M_tot->Draw("E0PX0SAME");
	if (isratio==0) {
	  h_M_stat->SetMarkerStyle(24);
	  h_M_stat->SetMarkerSize(0.7);
	}
	if (isratio==1) {
	  h_M_stat->SetMarkerStyle(22);
	  h_M_stat->SetMarkerSize(0.9);
	}
	h_M_stat->Draw("E1PX0SAME");
	h_M_stat->Draw("E0PX0SAME");

	if (isratio==0) {
	  TH1F *h_M2_tot= (TH1F*)h_mcg->Clone();
	  TH1F *h_M2_stat= (TH1F*)h_mcg->Clone();

	  h_M2_tot->Divide(h_data_tot);
	  h_M2_stat->Divide(h_data_stat);

	  for (int i=0;i<=h_M2_tot->GetNbinsX()+1;i++) {
	    if (h_M2_tot->GetBinContent(i)<=0) h_M2_tot->SetBinContent(i,999.);
	    if (h_M2_stat->GetBinContent(i)<=0) h_M2_stat->SetBinContent(i,999.);
	  }

	  TGraphErrors *g_M2_tot = new TGraphErrors(h_M2_tot);
	  TGraphErrors *g_M2_stat = new TGraphErrors(h_M2_stat);

	  float dx = 0.0;
	  if (drawInclusive) dx = 0.1*(g_M2_tot->GetXaxis()->GetXmax()-g_M2_tot->GetXaxis()->GetXmin())/g_M2_tot->GetN();
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
	  if (drawInclusive) g_M2_tot->Draw("E0PX0SAME");
	  g_M2_stat->SetMarkerStyle(20);
	  if (drawInclusive) g_M2_stat->Draw("E1PX0SAME");
	  if (drawInclusive) g_M2_stat->Draw("E0PX0SAME");
	}

	TLatex *t2 = new TLatex();
	t2->SetTextSize(0.2);
	t2->SetTextFont(42);
	t2->SetLineWidth(2);
	t2->SetNDC();
	if (useSherpa) {
	  t2->DrawLatex(0.15,0.7,"MadGraph 5FS / MadGraph 4FS");
	} else {
	  t2->DrawLatex(0.15,0.13,"MadGraph 5FS, normalized to #sigma_{NNLO}");
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

	for (int i=0;i<=h_S_tot->GetNbinsX()+1;i++) {
	  if (h_S_tot->GetBinContent(i)<=0) h_S_tot->SetBinContent(i,999.);
	  if (h_S_stat->GetBinContent(i)<=0) h_S_stat->SetBinContent(i,999.);
	}

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
	//h_S_tot->SetMarkerSize(0.7);
	h_S_stat->GetXaxis()->SetTitleOffset(0.7);
	h_S_stat->SetMarkerColor(kBlack);
	h_S_stat->SetLineColor(kBlack);
	h_S_stat->SetLineWidth(1);
	//h_S_stat->SetMarkerSize(0.7);

	if (isratio==0) {
	  h_S_tot->SetMarkerStyle(24);
	  h_S_tot->SetMarkerSize(0.7);
	}
	if (isratio==1) {
	  h_S_tot->SetMarkerStyle(22);
	  h_S_tot->SetMarkerSize(0.9);
	}
	if (useSherpa) {
	  h_S_tot->Draw("E1PX0");
	  h_S_tot->Draw("E0PX0SAME");
	} else {
	  for (int i=0;i<=h_S_tot->GetNbinsX()+1;i++) {
	    h_S_tot->SetBinContent(i, -999.);
	    h_S_tot->SetBinError(i, 0.);
	  }
	  h_S_tot->Draw("E1PX0");
	  h_S_tot->Draw("E0PX0SAME");
	}
	if (isratio==0) {
	  h_S_stat->SetMarkerStyle(24);
	  h_S_stat->SetMarkerSize(0.7);
	}
	if (isratio==1) {
	  h_S_stat->SetMarkerStyle(22);
	  h_S_stat->SetMarkerSize(0.9);
	}
	if (useSherpa) h_S_stat->Draw("E1PX0SAME");
	if (useSherpa) h_S_stat->Draw("E0PX0SAME");

	if (isratio==0) {
	  TH1F *h_S2_tot= (TH1F*)h_mcg1->Clone();
	  TH1F *h_S2_stat= (TH1F*)h_mcg1->Clone();

	  h_S2_tot->Divide(h_data_tot);
	  h_S2_stat->Divide(h_data_stat);

	  for (int i=0;i<=h_S2_tot->GetNbinsX()+1;i++) {
	    if (h_S2_tot->GetBinContent(i)<=0) h_S2_tot->SetBinContent(i,999.);
	    if (h_S2_stat->GetBinContent(i)<=0) h_S2_stat->SetBinContent(i,999.);
	  }

	  TGraphErrors *g_S2_tot = new TGraphErrors(h_S2_tot);
	  TGraphErrors *g_S2_stat = new TGraphErrors(h_S2_stat);

	  float dx = 0.0;
	  if (drawInclusive) dx = 0.1*(g_S2_tot->GetXaxis()->GetXmax()-g_S2_tot->GetXaxis()->GetXmin())/g_S2_tot->GetN();
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
	  if (useSherpa && drawInclusive) g_S2_tot->Draw("E1PX0SAME");
	  if (useSherpa && drawInclusive) g_S2_tot->Draw("E0PX0SAME");
	  g_S2_stat->SetMarkerStyle(20);
	  if (useSherpa && drawInclusive) g_S2_stat->Draw("E1PX0SAME");
	  if (useSherpa && drawInclusive) g_S2_stat->Draw("E0PX0SAME");
	}

	TLatex *t3 = new TLatex();
	t3->SetTextSize(0.2);
	t3->SetTextFont(42);
	t3->SetLineWidth(2);
	t3->SetNDC();
	if (useSherpa) {
	  t3->DrawLatex(0.15,0.7,"Sherpa");
	} else {
	  t3->DrawLatex(0.15,0.13,"MadGraph 4FS, normalized to #sigma_{NLO}");
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

	for (int i=0;i<=h_P_tot->GetNbinsX()+1;i++) {
	  if (h_P_tot->GetBinContent(i)<=0) h_P_tot->SetBinContent(i,999.);
	  if (h_P_stat->GetBinContent(i)<=0) h_P_stat->SetBinContent(i,999.);
	}

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
	//h_P_tot->SetMarkerSize(0.7);
	h_P_stat->GetXaxis()->SetTitleOffset(0.7);
	h_P_stat->SetMarkerColor(kBlack);
	h_P_stat->SetLineColor(kBlack);
	h_P_stat->SetLineWidth(1);
	//h_P_stat->SetMarkerSize(0.7);

	if (isratio==0) {
	  h_P_tot->SetMarkerStyle(24);
	  h_P_tot->SetMarkerSize(0.7);
	}
	if (isratio==1) {
	  h_P_tot->SetMarkerStyle(22);
	  h_P_tot->SetMarkerSize(0.9);
	}
	h_P_tot->Draw("E1PX0");
	h_P_tot->Draw("E0PX0SAME");
	if (isratio==0) {
	  h_P_stat->SetMarkerStyle(24);
	  h_P_stat->SetMarkerSize(0.7);
	}
	if (isratio==1) {
	  h_P_stat->SetMarkerStyle(22);
	  h_P_stat->SetMarkerSize(0.9);
	}
	h_P_stat->Draw("E1PX0SAME");
	h_P_stat->Draw("E0PX0SAME");

	if (isratio==0) {
	  TH1F *h_P2_tot= (TH1F*)h_mcg2->Clone();
	  TH1F *h_P2_stat= (TH1F*)h_mcg2->Clone();

	  h_P2_tot->Divide(h_data_tot);
	  h_P2_stat->Divide(h_data_stat);

	  for (int i=0;i<=h_P2_tot->GetNbinsX()+1;i++) {
	    if (h_P2_tot->GetBinContent(i)<=0) h_P2_tot->SetBinContent(i,999.);
	    if (h_P2_stat->GetBinContent(i)<=0) h_P2_stat->SetBinContent(i,999.);
	  }

	  TGraphErrors *g_P2_tot = new TGraphErrors(h_P2_tot);
	  TGraphErrors *g_P2_stat = new TGraphErrors(h_P2_stat);

	  float dx = 0.0;
	  if (drawInclusive) dx = 0.1*(g_P2_stat->GetXaxis()->GetXmax()-g_P2_stat->GetXaxis()->GetXmin())/g_P2_stat->GetN();
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
	  if (drawInclusive) g_P2_tot->Draw("E0PX0SAME");
	  g_P2_stat->SetMarkerStyle(20);
	  if (drawInclusive) g_P2_stat->Draw("E1PX0SAME");
	  if (drawInclusive) g_P2_stat->Draw("E0PX0SAME");
	}

	TLatex *t4 = new TLatex();
	t4->SetTextSize(0.13);
	t4->SetTextFont(42);
	t4->SetLineWidth(2);
	t4->SetNDC();
	t4->DrawLatex(0.15,0.40,"Powheg, normalized to #sigma_{NLO}");

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

	for (int i=0;i<=h_M3_tot->GetNbinsX()+1;i++) {
	  if (h_M3_tot->GetBinContent(i)<=0) h_M3_tot->SetBinContent(i,999.);
	  if (h_M3_stat->GetBinContent(i)<=0) h_M3_stat->SetBinContent(i,999.);
	}

	TGraphErrors *g_M3_tot = new TGraphErrors(h_M3_tot);
	TGraphErrors *g_M3_stat = new TGraphErrors(h_M3_stat);

	float dx = 0.0;
	for (int i=0; i<g_M3_tot->GetN(); i++) {
	  g_M3_stat->SetPoint(i, g_M3_stat->GetX()[i]+dx, g_M3_stat->GetY()[i]);
	  g_M3_stat->SetPointError(i, 0, g_M3_stat->GetEY()[i]);
	  g_M3_tot->SetPoint(i, g_M3_tot->GetX()[i]+dx, g_M3_tot->GetY()[i]);
	  g_M3_tot->SetPointError(i, 0, g_M3_tot->GetEY()[i]);
	}

	g_M3_tot->SetMarkerColor(kRed+1);
	g_M3_tot->SetLineColor(kRed+1);
	g_M3_tot->SetLineWidth(1);
	//g_M3_tot->SetMarkerSize(0.7);
	g_M3_stat->GetXaxis()->SetTitleOffset(0.7);
	g_M3_stat->SetMarkerColor(kBlack);
	g_M3_stat->SetLineColor(kBlack);
	g_M3_stat->SetLineWidth(1);
	//g_M3_stat->SetMarkerSize(0.7);

	if (useSherpa) {
	  g_M3_stat->SetMarkerStyle(25);
	  g_M3_tot->SetMarkerStyle(25);
	} else {
	  if (isratio==0) {
	    g_M3_stat->SetMarkerStyle(24);
	    g_M3_stat->SetMarkerSize(0.7);
	  }
	  if (isratio==1) {
	    g_M3_stat->SetMarkerStyle(22);
	    g_M3_stat->SetMarkerSize(0.9);
	  }
	  g_M3_stat->SetMarkerColor(kBlack);
	  if (isratio==0) {
	    g_M3_tot->SetMarkerStyle(24);
	    g_M3_tot->SetMarkerSize(0.7);
	  }
	  if (isratio==1) {
	    g_M3_tot->SetMarkerStyle(22);
	    g_M3_tot->SetMarkerSize(0.9);
	  }
	  g_M3_tot->SetMarkerColor(kBlack);
	}
	if (useSherpa) {
	  g_M3_tot->SetMarkerStyle(25);
	} else {
	  if (isratio==0) {
	    g_M3_tot->SetMarkerStyle(24);
	    g_M3_tot->SetMarkerSize(0.7);
	  }
	  if (isratio==1) {
	    g_M3_tot->SetMarkerStyle(22);
	    g_M3_tot->SetMarkerSize(0.9);
	  }
	}
	g_M3_tot->Draw("E1P");
	g_M3_tot->Draw("E0PSAME");
	g_M3_stat->Draw("E1PSAME");
	g_M3_stat->Draw("E0PSAME");

	TLine *OLine5 = new TLine(h_P_tot->GetXaxis()->GetXmin(),1.,h_P_tot->GetXaxis()->GetXmax(),1.);
	OLine5->SetLineColor(kOrange+7);
	OLine5->SetLineWidth(2);
	OLine5->Draw();

	c1->cd();

	if (isratio==1) {
	  h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 10);
	}

	if (title_b=="w_first_jet_pt_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp_{T} [pb/GeV]");
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
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp^{b}_{T} [pb/GeV]");
	  h_P_tot->GetXaxis()->SetTitle("leading jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{b}_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 10);
	  }
	} else if (title_b=="w_second_bjet_pt") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp^{b}_{T} [pb/GeV]");
	  h_P_tot->GetXaxis()->SetTitle("subleading jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{b}_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 10);
	  }
	} else if (title_b=="w_first_bjet_eta") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#eta^{b} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("leading jet #eta");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta^{b} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_second_bjet_eta") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#eta^{b} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("subleading jet #eta");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta^{b} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_first_bjet_eta_abs") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d|#eta^{b}| [pb]");
	  h_P_tot->GetXaxis()->SetTitle("leading jet |#eta|");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d|#eta^{b}| [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_second_bjet_eta_abs") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d|#eta^{b}| [pb]");
	  h_P_tot->GetXaxis()->SetTitle("subleading jet |#eta|");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d|#eta^{b}| [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_pt_Z_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp^{Z}_{T} [pb/GeV]");
	  h_P_tot->GetXaxis()->SetTitle("Z boson p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{Z}_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 20);
	  }
	} else if (title_b=="w_Ht_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dH_{T} [pb/GeV]");
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
	} else if (title_b=="w_DR_Zb_min") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#Delta R^{min}_{bZ} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta R^{min}(bZ) [rad]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta R^{min}_{Zb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_DR_Zb_max") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#Delta R^{max}_{bZ} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta R^{max}(bZ) [rad]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta R^{max}_{Zb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_delta_phi_2b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#Delta#phi_{bb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta#phi(bb) [rad]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta#phi_{bb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_DR_bb") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d #Delta R_{bb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta R(bb) [rad]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta R_{bb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_Zbb_mass") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d M_{Zbb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("M(Zbb) [GeV]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d M_{Zbb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_bb_mass") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d M_{bb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("M(bb) [GeV]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d M_{bb} [%]");
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
	    TFile f((path + "/combined/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.root").c_str(),"RECREATE");
	    if (drawInclusive) {
	      h_data_stat->Write((title+"_data_stat").c_str());
	      h_data_tot->Write((title+"_data_tot").c_str());
	    }
	    h_data_b_stat->Write((title_b+"_data_stat").c_str());
	    h_data_b_tot->Write((title_b+"_data_tot").c_str());
	    if (drawInclusive) {
	      h_mcg->Write((title+"_mcg").c_str());
	      h_mcg1->Write((title+"_mcg1").c_str());
	      h_mcg2->Write((title+"_mcg2").c_str());
	      h_mcg3->Write((title+"_mcg3").c_str());
	    }
	    h_mcg_b->Write((title_b+"_mcg").c_str());
	    h_mcg1_b->Write((title_b+"_mcg1").c_str());
	    h_mcg2_b->Write((title_b+"_mcg2").c_str());
	    h_mcg3_b->Write((title_b+"_mcg3").c_str());
	    f.Close();
	  }
	  if (isratio==1) {
	    gSystem->mkdir((path + "/combined/" + version + "/ratios_unfolding" + dirbSel + "/").c_str(), kTRUE);
	    c1->SaveAs((path + "/combined/" + version + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.pdf").c_str());
	    out.open((path + "/combined/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.dat").c_str());
	    out1.open((path + "/combined/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.txt").c_str());
	    out2.open((path + "/combined/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.tex").c_str());
	    TFile f((path + "/combined/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.root").c_str(),"RECREATE");
	    h_data_b_stat->Write((title_b+"_data_stat").c_str());
	    h_data_b_tot->Write((title_b+"_data_tot").c_str());
	    h_mcg_b->Write((title_b+"_mcg").c_str());
	    h_mcg1_b->Write((title_b+"_mcg1").c_str());
	    h_mcg2_b->Write((title_b+"_mcg2").c_str());
	    h_mcg3_b->Write((title_b+"_mcg3").c_str());
	    f.Close();
	  }
	  if (isratio==0 && drawInclusive) {
	    out << h_data->GetName();
	    out << endl;
	    out << std::setw(29) << "data";
	    out << std::setw(14) << "bkg";
	    out << std::setw(14) << "eff";
	    out << std::setw(14) << "jec";
	    out << std::setw(14) << "jer";
	    out << std::setw(14) << "pu";
	    out << std::setw(14) << "bkg";
	    out << std::setw(14) << "ttbar";
	    out << std::setw(14) << "bfit";
	    out << std::setw(14) << "btemp";
	    out << std::setw(14) << "btag";
	    out << std::setw(14) << "unfold";
	    out << std::setw(14) << "unfold";
	    out << std::setw(14) << "lumi";
	    out << std::setw(14) << "total";
	    out << std::setw(14) << "total";
	    out << std::setw(14) << "total";
	    out << endl;
	    out << std::setw(29) << "stat";
	    out << std::setw(14) << "stat";
	    out << std::setw(14) << "syst";
	    out << std::setw(14) << "syst";
	    out << std::setw(14) << "syst";
	    out << std::setw(14) << "syst";
	    out << std::setw(14) << "syst";
	    out << std::setw(14) << "stat";
	    out << std::setw(14) << "stat";
	    out << std::setw(14) << "syst";
	    out << std::setw(14) << "syst";
	    out << std::setw(14) << "stat";
	    out << std::setw(14) << "syst";
	    out << std::setw(14) << "syst";
	    out << std::setw(14) << "stat";
	    out << std::setw(14) << "syst";
	    out << std::setw(14) << "error";
	    out << std::setw(8) << "%";
	    out << endl;
	    for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	      out << std::fixed;
	      out << std::setw(2) << i;
	      out << " ";
	      out << std::setprecision(8);
	      out << std::setw(12) << h_data->GetBinContent(i);
	      out << " +- ";
	      out << std::setw(10) << h_data->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << stat_bkg->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << syst_eff->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << syst_jec->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << syst_jer->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << syst_pu->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << syst_bkg->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << stat_top->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << stat_bfit->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << syst_bfit2->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << syst_btag->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << stat_unfold->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << syst_unfold->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << syst_lumi->GetBinError(i);
	      out << " => ";
	      out << std::setw(10) << h_data_stat->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << h_data_syst->GetBinError(i);
	      out << " => ";
	      out << std::setw(10) << h_data_tot->GetBinError(i);
	      out << " => ";
	      out << std::setprecision(1);
	      out << std::setw(4) << TMath::Abs(100.*(h_data_stat->GetBinContent(i)==0 ? 0 : h_data_tot->GetBinError(i)/h_data_stat->GetBinContent(i)));
	      out << endl;
	    }
	    out << "tot";
	    out << std::setprecision(8);
	    out << std::setw(12) << xsec_tot_data;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_stat_data;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_stat_bkg;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_syst_eff;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_syst_jec;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_syst_jer;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_syst_pu;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_syst_bkg;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_stat_top;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_stat_bfit;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_syst_bfit2;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_syst_btag;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_stat_unfold;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_syst_unfold;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_syst_lumi;
	    out << " => ";
	    out << std::setw(10) << xsec_tot_stat_tot;
	    out << " +- ";
	    out << std::setw(10) << xsec_tot_syst_tot;
	    out << " => ";
	    out << std::setw(10) << xsec_tot_data_tot;
	    out << " => ";
	    out << std::setprecision(1);
	    out << std::setw(4) << TMath::Abs(100.*(xsec_tot_data_tot==0 ? 0 : xsec_tot_data_tot/xsec_tot_data));
	    out << endl;
	  }
	  out << h_data_b->GetName();
	  out << endl;
	  out << std::setw(29) << "data";
	  out << std::setw(14) << "bkg";
	  out << std::setw(14) << "eff";
	  out << std::setw(14) << "jec";
	  out << std::setw(14) << "jer";
	  out << std::setw(14) << "pu";
	  out << std::setw(14) << "bkg";
	  out << std::setw(14) << "ttbar";
	  out << std::setw(14) << "bfit";
	  out << std::setw(14) << "btemp";
	  out << std::setw(14) << "btag";
	  out << std::setw(14) << "unfold";
	  out << std::setw(14) << "unfold";
	  out << std::setw(14) << "lumi";
	  out << std::setw(14) << "total";
	  out << std::setw(14) << "total";
	  out << std::setw(14) << "total";
	  out << endl;
	  out << std::setw(29) << "stat";
	  out << std::setw(14) << "stat";
	  out << std::setw(14) << "syst";
	  out << std::setw(14) << "syst";
	  out << std::setw(14) << "syst";
	  out << std::setw(14) << "syst";
	  out << std::setw(14) << "syst";
	  out << std::setw(14) << "stat";
	  out << std::setw(14) << "stat";
	  out << std::setw(14) << "syst";
	  out << std::setw(14) << "syst";
	  out << std::setw(14) << "stat";
	  out << std::setw(14) << "syst";
	  out << std::setw(14) << "syst";
	  out << std::setw(14) << "stat";
	  out << std::setw(14) << "syst";
	  out << std::setw(14) << "error";
	  out << std::setw(8) << "%";
	  out << endl;
	  for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	    out << std::fixed;
	    out << std::setw(2) << i;
	    out << " ";
	    out << std::setprecision(8);
	    out << std::setw(12) << h_data_b->GetBinContent(i);
	    out << " +- ";
	    out << std::setw(10) << h_data_b->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << stat_b_bkg->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << syst_b_eff->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << syst_b_jec->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << syst_b_jer->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << syst_b_pu->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << syst_b_bkg->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << stat_b_top->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << stat_b_bfit->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << syst_b_bfit2->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << syst_b_btag->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << stat_b_unfold->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << syst_b_unfold->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << syst_b_lumi->GetBinError(i);
	    out << " => ";
	    out << std::setw(10) << h_data_b_stat->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << h_data_b_syst->GetBinError(i);
	    out << " => ";
	    out << std::setw(10) << h_data_b_tot->GetBinError(i);
	    out << " => ";
	    out << std::setprecision(1);
	    out << std::setw(4) << TMath::Abs(100.*(h_data_b_stat->GetBinContent(i)==0 ? 0 : h_data_b_tot->GetBinError(i)/h_data_b_stat->GetBinContent(i)));
	    out << endl;
	  }
	  out << "tot";
	  out << std::setprecision(8);
	  out << std::setw(12) << xsec_tot_data_b;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_stat_data_b;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_stat_b_bkg;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_syst_b_eff;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_syst_b_jec;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_syst_b_jer;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_syst_b_pu;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_syst_b_bkg;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_stat_b_top;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_stat_b_bfit;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_syst_b_bfit2;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_syst_b_btag;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_stat_b_unfold;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_syst_b_unfold;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_syst_b_lumi;
	  out << " => ";
	  out << std::setw(10) << xsec_tot_stat_b_tot;
	  out << " +- ";
	  out << std::setw(10) << xsec_tot_syst_b_tot;
	  out << " => ";
	  out << std::setw(10) << xsec_tot_data_b_tot;
	  out << " => ";
	  out << std::setprecision(1);
	  out << std::setw(4) << TMath::Abs(100.*(xsec_tot_data_b==0 ? 0 : xsec_tot_data_b_tot/xsec_tot_data_b));
	  out << endl;
	  out.close();
	  if (isratio==0 && drawInclusive) {
	    out1 << h_data->GetName() << " - RELATIVE ERRORS";
	    out1 << endl;
	    out1 << std::setw(8) << "data";
	    out1 << std::setw(9) << "bkg";
	    out1 << std::setw(9) << "eff";
	    out1 << std::setw(9) << "jec";
	    out1 << std::setw(9) << "jer";
	    out1 << std::setw(9) << "pu";
	    out1 << std::setw(9) << "bkg";
	    out1 << std::setw(9) << "ttbar";
	    out1 << std::setw(9) << "bfit";
	    out1 << std::setw(9) << "btemp";
	    out1 << std::setw(9) << "btag";
	    out1 << std::setw(9) << "unfold";
	    out1 << std::setw(9) << "unfold";
	    out1 << std::setw(9) << "lumi";
	    out1 << std::setw(9) << "total";
	    out1 << std::setw(9) << "total";
	    out1 << std::setw(9) << "total";
	    out1 << endl;
	    out1 << std::setw(8) << "stat";
	    out1 << std::setw(9) << "stat";
	    out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "stat";
	    out1 << std::setw(9) << "stat";
	    out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "stat";
	    out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "stat";
	    out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "error";
	    out1 << endl;
	    for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	      double val = TMath::Abs(100.*(h_data->GetBinContent(i)==0 ? 0 : 1./h_data->GetBinContent(i)));
	      out1 << std::fixed;
	      out1 << std::setw(2) << i;
	      out1 << " ";
	      out1 << std::setprecision(1);
	      out1 << std::setw(5) << h_data->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << stat_bkg->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_eff->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_jec->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_jer->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_pu->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_bkg->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << stat_top->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << stat_bfit->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_bfit2->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_btag->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << stat_unfold->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_unfold->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_lumi->GetBinError(i)*val;
	      out1 << " => ";
	      out1 << std::setw(5) << h_data_stat->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << h_data_syst->GetBinError(i)*val;
	      out1 << " => ";
	      out1 << std::setw(5) << h_data_tot->GetBinError(i)*val;
	      out1 << endl;
	    }
	    xval = TMath::Abs(100.*(xsec_tot_data==0 ? 0 : 1./xsec_tot_data));
	    out1 << "tot";
	    out1 << std::setprecision(1);
	    out1 << std::setw(5) << xsec_tot_stat_data*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_stat_bkg*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_syst_eff*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_syst_jec*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_syst_jer*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_syst_pu*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_syst_bkg*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_stat_top*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_stat_bfit*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_syst_bfit2*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_syst_btag*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_stat_unfold*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_syst_unfold*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_syst_lumi*xval;
	    out1 << " => ";
	    out1 << std::setw(5) << xsec_tot_stat_tot*xval;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_tot_syst_tot*xval;
	    out1 << " => ";
	    out1 << std::setw(5) << xsec_tot_data_tot*xval;
	    out1 << endl;
	  }
	  out1 << h_data_b->GetName() << " - RELATIVE ERRORS";
	  out1 << endl;
	  out1 << std::setw(8) << "data";
	  out1 << std::setw(9) << "bkg";
	  out1 << std::setw(9) << "eff";
	  out1 << std::setw(9) << "jec";
	  out1 << std::setw(9) << "jer";
	  out1 << std::setw(9) << "pu";
	  out1 << std::setw(9) << "bkg";
	  out1 << std::setw(9) << "ttbar";
	  out1 << std::setw(9) << "bfit";
	  out1 << std::setw(9) << "btemp";
	  out1 << std::setw(9) << "btag";
	  out1 << std::setw(9) << "unfold";
	  out1 << std::setw(9) << "unfold";
	  out1 << std::setw(9) << "lumi";
	  out1 << std::setw(9) << "total";
	  out1 << std::setw(9) << "total";
	  out1 << std::setw(9) << "total";
	  out1 << endl;
	  out1 << std::setw(8) << "stat";
	  out1 << std::setw(9) << "stat";
	  out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "stat";
	  out1 << std::setw(9) << "stat";
	  out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "stat";
	  out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "stat";
	  out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "error";
	  out1 << endl;
	  for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	    double val = TMath::Abs(100.*(h_data_b->GetBinContent(i)==0 ? 0 : 1./h_data_b->GetBinContent(i)));
	    out1 << std::fixed;
	    out1 << std::setw(2) << i;
	    out1 << " ";
	    out1 << std::setprecision(1);
	    out1 << std::setw(5) << h_data_b->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << stat_b_bkg->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << syst_b_eff->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << syst_b_jec->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << syst_b_jer->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << syst_b_pu->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << syst_b_bkg->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << stat_b_top->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << stat_b_bfit->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << syst_b_bfit2->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << syst_b_btag->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << stat_b_unfold->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << syst_b_unfold->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << syst_b_lumi->GetBinError(i)*val;
	    out1 << " => ";
	    out1 << std::setw(5) << h_data_b_stat->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << h_data_b_syst->GetBinError(i)*val;
	    out1 << " => ";
	    out1 << std::setw(5) << h_data_b_tot->GetBinError(i)*val;
	    out1 << endl;
	  }
	  xval = TMath::Abs(100.*(xsec_tot_data_b==0 ? 0 : 1./xsec_tot_data_b));
	  out1 << "tot";
	  out1 << std::setprecision(1);
	  out1 << std::setw(5) << xsec_tot_stat_data_b*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_stat_b_bkg*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_syst_b_eff*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_syst_b_jec*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_syst_b_jer*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_syst_b_pu*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_syst_b_bkg*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_stat_b_top*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_stat_b_bfit*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_syst_b_bfit2*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_syst_b_btag*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_stat_b_unfold*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_syst_b_unfold*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_syst_b_lumi*xval;
	  out1 << " => ";
	  out1 << std::setw(5) << xsec_tot_stat_b_tot*xval;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_tot_syst_b_tot*xval;
	  out1 << " => ";
	  out1 << std::setw(5) << xsec_tot_data_b_tot*xval;
	  out1 << endl;
	  out1.close();
	  /*
	  if (isratio==0 && drawInclusive) {
	    out2 << std::setw(7) << "\\textbf{data} & ";
	    out2 << std::setw(8) << "\\textbf{bkg} & ";
	    out2 << std::setw(8) << "\\textbf{eff} & ";
	    out2 << std::setw(8) << "\\textbf{jec} & ";
	    out2 << std::setw(8) << "\\textbf{jer} & ";
	    out2 << std::setw(8) << "\\textbf{pu} & ";
	    out2 << std::setw(8) << "\\textbf{bkg} & ";
	    out2 << std::setw(8) << "\\textbf{ttbar} & ";
	    out2 << std::setw(8) << "\\textbf{bfit} & ";
	    out2 << std::setw(8) << "\\textbf{btemp} & ";
	    out2 << std::setw(8) << "\\textbf{btag} & ";
	    out2 << std::setw(8) << "\\textbf{unfold} &";
	    out2 << std::setw(8) << "\\textbf{unfold} &";
	    out2 << std::setw(8) << "\\textbf{lumi} & ";
	    out2 << std::setw(8) << "\\textbf{total} & ";
	    out2 << std::setw(8) << "\\textbf{total} & ";
	    out2 << std::setw(8) << "\\textbf{total} ";
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
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{stat} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{stat} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{error} \\\\";
	    out2 << endl;
	    out2 << "\\tabularnewline" << " " << "\\hline";
	    out2 << endl;
	    for (int i=1;i<=h_data->GetNbinsX();i++) {
	      double val = TMath::Abs(100.*(h_data->GetBinContent(i)==0 ? 0 : 1./h_data->GetBinContent(i)));
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
	      out2 << std::setw(4) << syst_bfit2->GetBinError(i)*val;
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
	      out2 << " \\tabularnewline" << " " << "\\hline";
	      out2 << endl;
	    }
	  }
	  */
	  out2 << std::setw(7) << "\\textbf{data} & ";
	  out2 << std::setw(8) << "\\textbf{bkg} & ";
	  out2 << std::setw(8) << "\\textbf{eff} & ";
	  out2 << std::setw(8) << "\\textbf{jec} & ";
	  out2 << std::setw(8) << "\\textbf{jer} & ";
	  out2 << std::setw(8) << "\\textbf{pu} & ";
	  out2 << std::setw(8) << "\\textbf{bkg} & ";
	  out2 << std::setw(8) << "\\textbf{ttbar} & ";
	  out2 << std::setw(8) << "\\textbf{bfit} & ";
	  out2 << std::setw(8) << "\\textbf{btemp} & ";
	  out2 << std::setw(8) << "\\textbf{btag} & ";
	  out2 << std::setw(8) << "\\textbf{unfold} & ";
	  out2 << std::setw(8) << "\\textbf{unfold} & ";
	  out2 << std::setw(8) << "\\textbf{lumi} & ";
	  out2 << std::setw(8) << "\\textbf{total} & ";
	  out2 << std::setw(8) << "\\textbf{total} & ";
	  out2 << std::setw(8) << "\\textbf{total} ";
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
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{stat} & ";
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{stat} & ";
	  out2 << std::setw(8) << "\\textbf{syst} & ";
	  out2 << std::setw(8) << "\\textbf{error} \\\\";
	  out2 << endl;
	  out2 << "\\tabularnewline" << " " << "\\hline";
	  out2 << endl;
	  for (int i=1;i<=h_data_b->GetNbinsX();i++) {
	    double val = TMath::Abs(100.*(h_data_b->GetBinContent(i)==0 ? 0 : 1./h_data_b->GetBinContent(i)));
	    out2 << std::fixed;
	    out2 << std::setw(2) << i;
	    out2 << " & ";
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
	    out2 << std::setw(4) << syst_b_bfit2->GetBinError(i)*val;
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
	    out2 << " \\tabularnewline" << " " << "\\hline";
	    out2 << endl;
	  }

	}
}

