#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v14.h"

#include "fixrange.C"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";
//string path = "/gpfs/cms/users/lalicata/work/test/data/";

TH1F* read(string subdir, string title, int ilepton, string dirbSel) {
  TH1F* hist;
  TFile* file=0;
  if (ilepton==1) {
    file = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/unfolding" + dirbSel + "/" + title + "_unfolding.root").c_str());
  }
  if (ilepton==2) {
    file = TFile::Open((path + "/muons/" + version + "/" + subdir + "/unfolding" + dirbSel + "/" + title + "_unfolding.root").c_str());
  }
  hist = (TH1F*)gDirectory->Get(title.c_str())->Clone();
  hist->SetDirectory(0);
  file->Close();
  return hist;
}

void DataMCComp7(string title="", int plot=0, int ilepton=1, int isratio=1, int numB=0) {

//int useSysBfit2=0;
//int useSysBfit2=1; // include Bfit2 systematics, templates from data
//int useSysBfit2=2; // include Bfit2 systematics, templates from Sherpa
int useSysBfit2=3; // include Bfit2 systematics, templates from MadGraph aMC@NLO

int useSysRMS=0;
//int useSysRMS=1; // include xsec RMS systematics

//int useSysUnfold=0;
//int useSysUnfold=1; // use Sherpa for unfolding systematics
//int useSysUnfold=2; // use Powheg for unfolding systematics
//int useSysUnfold=3; // use MadGraph 4FS for unfolding systematics
//int useSysUnfold=4; // use data weighted MC for unfolding systematics
int useSysUnfold=5; // use MadGraph aMC@NLO for unfolding systematics

int useMC=0;
//int useMC=1; // use MC prediction

int useSherpa=0;
//int useSherpa=1; // use Sherpa MC prediction

//int useNewPowheg=0;
int useNewPowheg=1; // use new Powheg MC prediction

//int drawInclusive=0;
int drawInclusive=1; // do plot the "inclusive" histogram

double ele_eff_sys=0.015;
double muo_eff_sys=0.020;
double btag_sys=0.030;
double lumi_sys=0.026;

string bSel=" ";
string dirbSel="";
string subdir="0";
string postfix="";

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
	  if (!useMC) gROOT->GetColor(kRed)->SetAlpha(0.0);
	  gROOT->GetColor(kGreen+2)->SetAlpha(0.5);
	  gROOT->GetColor(kMagenta-6)->SetAlpha(0.5);
	  if (!useSherpa) gROOT->GetColor(kMagenta-6)->SetAlpha(0.0);
	  gROOT->GetColor(kBlue-4)->SetAlpha(0.5);
	  gROOT->GetColor(kOrange+7)->SetAlpha(0.5);
	}

	/* purity */

	double c_b=1.0;
	double ec_b=0.0;
	double c_c=1.0;
	double ec_c=0.0;
	double c_uds=1.0;
	double ec_uds=0.0;

	ifstream in;
	if (ilepton==1) {
	  in.open((path + "/electrons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_BJP_doFit" + ".dat").c_str());
	}
	if (ilepton==2) {
	  in.open((path + "/muons/" + version + "/" + subdir + "/distributions" + dirbSel + "/" + "w_BJP_doFit" + ".dat").c_str());
	}
	if (numB!=2) {
	  in >> c_uds >> ec_uds;
	  in >> c_b >> ec_b;
	  in >> c_c >> ec_c;
	}
	else {
	  in >> c_b >> ec_b;
	}
	in.close();

	double Lumi2012=0;

	if (ilepton==1) Lumi2012 = Lumi2012_ele;
	if (ilepton==2) Lumi2012 = Lumi2012_muon;

	double norm1 = ((Lumi2012 * Xsec_dy) / Ngen_dy);
	double norm1_1 = ((Lumi2012 * Xsec_dy_1) / Ngen_dy_1);
	double norm1_2=0;
	if (ilepton==1) norm1_2 = ((Lumi2012 * Xsec_dy_2) / Ngen_dy_2_ee);
	if (ilepton==2) norm1_2 = ((Lumi2012 * Xsec_dy_2) / Ngen_dy_2_mm);
	double norm1_3 = ((Lumi2012 * Xsec_dy_3) / Ngen_dy_3);

	if (title.empty()) title = "w_jetmultiplicity";

	if (ilepton==1) {
	  if (title.find("muon")!=string::npos) return;
	  if (title.find("mm")!=string::npos) return;
	}
	if (ilepton==2) {
	  if (title.find("ele")!=string::npos) return;
	  if (title.find("ee")!=string::npos) return;
	}

	TFile *mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL.root").c_str());
	TFile *mcg = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	TFile *mcg1 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());
	TFile *mcg2=0;
	if (ilepton==1) mcg2 = TFile::Open((path + "/" + version + "/" + "DYToEE_powheg_gen.root").c_str());
	if (ilepton==2) mcg2 = TFile::Open((path + "/" + version + "/" + "DYToMuMu_powheg_gen.root").c_str());
	if (useNewPowheg) {
	  Ngen_dy_2_ee = 100*10000;
	  Ngen_dy_2_mm = 100*10000;
	  Xsec_dy_2 = 333.866;
	  norm1_2 = ((Lumi2012 * Xsec_dy_2) / ((Ngen_dy_2_ee+Ngen_dy_2_mm)/2.));
	  if (ilepton==1) mcg2 = TFile::Open(("/gpfs/cms/users/candelis/work/ZbSkim/powheg/data/" + version + "/" + "powheg_ele.root").c_str());
	  if (ilepton==2) mcg2 = TFile::Open(("/gpfs/cms/users/candelis/work/ZbSkim/powheg/data/" + version + "/" + "powheg_muo.root").c_str());
	}
	TFile *mcg3 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL2_gen.root").c_str());

	string title_b = title;

	if (drawInclusive) {
	  if (title.find("_bjet_")!=string::npos) {
	    title.erase(title.find("_bjet_")+1, 1);
	  } else {
	    title_b = title + "_b";
	  }
	}

	TH1F* h_data;
	TH1F* h_data_b;
	h_data = read(subdir, title, ilepton, dirbSel);
	h_data_b = read(subdir, title_b, ilepton, dirbSel);
	if (drawInclusive) {
	  h_data = (TH1F*)h_data->Clone();
	}
	h_data->SetStats(0);
	h_data_b->SetStats(0);

	const int NMAX = 100;
	TH1F* h_data_scan[NMAX];
	TH1F* h_data_b_scan[NMAX];

	for (int i=0;i<NMAX;i++) {
	  h_data_scan[i] = 0;
	  h_data_b_scan[i] = 0;
	  if (i==8 && useSysUnfold!=1) continue;
	  if (i==9 && useSysUnfold!=2) continue;
	  if (i==14 && useSysUnfold!=3) continue;
	  if (i==15 && useSysUnfold!=4) continue;
	  if (i==16 && useSysUnfold!=5) continue;
	  if (i==17 && useSysBfit2!=1) continue;
	  if (i==18 && useSysBfit2!=2) continue;
	  if (i==19 && useSysBfit2!=3) continue;
	  if (i<=19) {
	    stringstream ss;
	    ss << i;
	    h_data_scan[i] = read(ss.str(), title, ilepton, dirbSel);
	    h_data_b_scan[i] = read(ss.str(), title_b, ilepton, dirbSel);
	  }
	}

	if (ilepton==1) mc1->cd(("demoEle"+postfix).c_str());
	if (ilepton==2) mc1->cd(("demoMuo"+postfix).c_str());
	TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mc1b_b = (TH1F*)gDirectory->Get(("b"+title_b.substr(1)).c_str());

	if (ilepton==1) mcg->cd(("demoEleGen"+postfix).c_str());
	if (ilepton==2) mcg->cd(("demoMuoGen"+postfix).c_str());
	TH1F* h_mcg = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mcg_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (ilepton==1) mcg1->cd(("demoEleGen"+postfix).c_str());
	if (ilepton==2) mcg1->cd(("demoMuoGen"+postfix).c_str());
	TH1F* h_mcg1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mcg1_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (!h_mcg1) h_mcg1 = (TH1F*)h_mcg->Clone();
	if (!h_mcg1_b) h_mcg1_b = (TH1F*)h_mcg_b->Clone();

	if (ilepton==1) mcg2->cd(("demoEleGen"+postfix).c_str());
	if (ilepton==2) mcg2->cd(("demoMuoGen"+postfix).c_str());
	TH1F* h_mcg2 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mcg2_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (!h_mcg2) h_mcg2 = (TH1F*)h_mcg->Clone();
	if (!h_mcg2_b) h_mcg2_b = (TH1F*)h_mcg_b->Clone();

	if (ilepton==1) mcg3->cd(("demoEleGen"+postfix).c_str());
	if (ilepton==2) mcg3->cd(("demoMuoGen"+postfix).c_str());
	TH1F* h_mcg3 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_mcg3_b = (TH1F*)gDirectory->Get(title_b.c_str());

	if (!h_mcg3) h_mcg3 = (TH1F*)h_mcg->Clone();
	if (!h_mcg3_b) h_mcg3_b = (TH1F*)h_mcg_b->Clone();

	if (!drawInclusive) {
	  h_mc1 = (TH1F*)h_mc1->Clone();
	  h_mcg = (TH1F*)h_mcg->Clone();
	  h_mcg1 = (TH1F*)h_mcg1->Clone();
	  h_mcg2 = (TH1F*)h_mcg2->Clone();
	  h_mcg3 = (TH1F*)h_mcg3->Clone();
	}

	h_mc1->Sumw2();
	h_mcg->Sumw2();
	h_mcg1->Sumw2();
	h_mcg2->Sumw2();
	h_mcg3->Sumw2();

	h_mc1b_b->Sumw2();
	h_mcg_b->Sumw2();
	h_mcg1_b->Sumw2();
	h_mcg2_b->Sumw2();
	h_mcg3_b->Sumw2();

	h_mc1b_b->Scale(c_b);
	for (int i=0; i<=h_mc1b_b->GetNbinsX()+1; i++) {
	  float e = TMath::Power(h_mc1b_b->GetBinError(i),2);
	  e = e + TMath::Power(h_mc1b_b->GetBinContent(i)*(ec_b/c_b),2);
	  h_mc1b_b->SetBinError(i, TMath::Sqrt(e));
	}

	if (ilepton==1) {
	  TFile f((path + "/electrons/" + version + "/" + subdir +"/efficiency" + dirbSel + "/" + string(h_data->GetName()) + "_efficiency.root").c_str());
	  TFile f_b((path + "/electrons/" + version + "/" + subdir +"/efficiency" + dirbSel + "/" + string(h_data_b->GetName()) + "_efficiency.root").c_str());
	  TH1F* h = (TH1F*)f.Get(h_data->GetName())->Clone();
	  TH1F* h_b = (TH1F*)f_b.Get(h_data_b->GetName())->Clone();
	  h->SetDirectory(0);
	  h_b->SetDirectory(0);
	  f.Close();
	  f_b.Close();
	  h_mc1->Divide(h);
	  h_mc1b_b->Divide(h_b);
	}
	if (ilepton==2) {
	  TFile f((path + "/muons/" + version + "/" + subdir +"/efficiency" + dirbSel + "/" + string(h_data->GetName()) + "_efficiency.root").c_str());
	  TFile f_b((path + "/muons/" + version + "/" + subdir +"/efficiency" + dirbSel + "/" + string(h_data_b->GetName()) + "_efficiency.root").c_str());
	  TH1F* h = (TH1F*)f.Get(h_data->GetName())->Clone();
	  TH1F* h_b = (TH1F*)f_b.Get(h_data_b->GetName())->Clone();
	  h->SetDirectory(0);
	  h_b->SetDirectory(0);
	  f.Close();
	  f_b.Close();
	  h_mc1->Divide(h);
	  h_mc1b_b->Divide(h_b);
	}

	h_data = fixrange(h_data, numB);
	h_data_b = fixrange(h_data_b, numB);
	for (int i=0;i<NMAX;i++) {
	  if (h_data_scan[i]) h_data_scan[i] = fixrange(h_data_scan[i], numB);
	  if (h_data_b_scan[i]) h_data_b_scan[i] = fixrange(h_data_b_scan[i], numB);
	}

	h_mc1 = fixrange(h_mc1, numB);
	h_mc1b_b = fixrange(h_mc1b_b, numB);
	h_mcg = fixrange(h_mcg, numB);
	h_mcg_b = fixrange(h_mcg_b, numB);
	h_mcg1 = fixrange(h_mcg1, numB);
	h_mcg1_b = fixrange(h_mcg1_b, numB);
	h_mcg2 = fixrange(h_mcg2, numB);
	h_mcg3 = fixrange(h_mcg3, numB);
	h_mcg2_b = fixrange(h_mcg2_b, numB);
	h_mcg3_b = fixrange(h_mcg3_b, numB);

	h_mc1->Scale(norm1);
	h_mcg->Scale(norm1);
	h_mcg1->Scale(norm1_1);
	h_mcg2->Scale(norm1_2);
	h_mcg3->Scale(norm1_3);

	h_mc1b_b->Scale(norm1);
	h_mcg_b->Scale(norm1);
	h_mcg1_b->Scale(norm1_1);
	h_mcg2_b->Scale(norm1_2);
	h_mcg3_b->Scale(norm1_3);

	h_data->Scale(1./Lumi2012, "width");
	h_data_b->Scale(1./Lumi2012, "width");
	for (int i=0;i<NMAX;i++) {
	  if (h_data_scan[i]) h_data_scan[i]->Scale(1./Lumi2012, "width");
	  if (h_data_b_scan[i]) h_data_b_scan[i]->Scale(1./Lumi2012, "width");
	}
	h_mc1->Scale(1./Lumi2012, "width");
	h_mc1b_b->Scale(1./Lumi2012, "width");
	if (isratio==1) {
	  h_data_b->Divide(h_data);
	  h_data_b->Scale(100.);
	  for (int i=0;i<NMAX;i++) {
	    if (h_data_scan[i]) h_data_b_scan[i]->Divide(h_data_scan[i]);
	    if (h_data_b_scan[i]) h_data_b_scan[i]->Scale(100.);
	  }
	  h_mc1b_b->Divide(h_mc1);
	  h_mc1b_b->Scale(100.);
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

	TH1F* stat_bkg = (TH1F*)h_data->Clone();
	TH1F* stat_b_bkg = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Sqrt(TMath::Max(0.,TMath::Power(h_data_scan[13]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2))/(TMath::Power(1.1,2)-1));
	  stat_bkg->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Max(0.,TMath::Power(h_data->GetBinError(i),2)-TMath::Power(stat_bkg->GetBinError(i),2)));
	  h_data->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Sqrt(TMath::Max(0.,TMath::Power(h_data_b_scan[13]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))/(TMath::Power(1.1,2)-1));
	  stat_b_bkg->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Max(0.,TMath::Power(h_data_b->GetBinError(i),2)-TMath::Power(stat_b_bkg->GetBinError(i),2)));
	  h_data_b->SetBinError(i, val);
	}

	double xsec_data = 0.0;
	double xsec_stat_data = 0.0;
	double xsec_data_b = 0.0;
	double xsec_stat_data_b = 0.0;
	double xsec_stat_bkg = 0.0;
	double xsec_stat_b_bkg = 0.0;
	xsec_data = h_data->IntegralAndError(0,h_data->GetNbinsX()+1,xsec_stat_data,"width");
	xsec_data_b = h_data_b->IntegralAndError(0,h_data_b->GetNbinsX()+1,xsec_stat_data_b,"width");
	stat_bkg->IntegralAndError(0,h_data->GetNbinsX()+1,xsec_stat_bkg,"width");
	stat_b_bkg->IntegralAndError(0,h_data_b->GetNbinsX()+1,xsec_stat_b_bkg,"width");

	TH1F* syst_eff = (TH1F*)h_data->Clone();
	TH1F* syst_b_eff = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (ilepton==1) val = ele_eff_sys * h_data_scan[0]->GetBinContent(i);
	  if (ilepton==2) val = muo_eff_sys * h_data_scan[0]->GetBinContent(i);
	  if (isratio) val = 0.0;
	  syst_eff->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (ilepton==1) val = ele_eff_sys * h_data_b_scan[0]->GetBinContent(i);
	  if (ilepton==2) val = muo_eff_sys * h_data_b_scan[0]->GetBinContent(i);
	  if (isratio) val = 0.0;
	  syst_b_eff->SetBinError(i, val);
	}
	double xsec_syst_eff = 0.0;
	double xsec_syst_b_eff = 0.0;
	if (ilepton==1) {
	  xsec_syst_eff = ele_eff_sys * h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width");
	  xsec_syst_b_eff = ele_eff_sys * h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width");
	}
	if (ilepton==2) {
	  xsec_syst_eff = muo_eff_sys * h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width");
	  xsec_syst_b_eff = muo_eff_sys * h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width");
	}

	TH1F* syst_jec = (TH1F*)h_data->Clone();
	TH1F* syst_b_jec = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_scan[2]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	  val = TMath::Max(val,TMath::Abs(h_data_scan[1]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i)));
	  syst_jec->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_b_scan[2]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	  val = TMath::Max(val,TMath::Abs(h_data_b_scan[1]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i)));
	  syst_b_jec->SetBinError(i, val);
	}
	double xsec_syst_jec = 0.0;
	double xsec_syst_b_jec = 0.0;
	xsec_syst_jec = TMath::Abs(h_data_scan[2]->Integral(0,h_data_scan[2]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width"));
	xsec_syst_jec = TMath::Max(xsec_syst_jec,TMath::Abs(h_data_scan[1]->Integral(0,h_data_scan[1]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width")));
	xsec_syst_b_jec = TMath::Abs(h_data_b_scan[2]->Integral(0,h_data_b_scan[2]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width"));
	xsec_syst_b_jec = TMath::Max(xsec_syst_b_jec,TMath::Abs(h_data_b_scan[1]->Integral(0,h_data_b_scan[1]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width")));

	TH1F* syst_jer = (TH1F*)h_data->Clone();
	TH1F* syst_b_jer = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_scan[12]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	  val = TMath::Max(val,TMath::Abs(h_data_scan[11]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i)));
	  syst_jer->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_b_scan[12]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	  val = TMath::Max(val,TMath::Abs(h_data_b_scan[11]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i)));
	  syst_b_jer->SetBinError(i, val);
	}
	double xsec_syst_jer = 0.0;
	double xsec_syst_b_jer = 0.0;
	xsec_syst_jer = TMath::Abs(h_data_scan[12]->Integral(0,h_data_scan[12]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width"));
	xsec_syst_jer = TMath::Max(xsec_syst_jer,TMath::Abs(h_data_scan[11]->Integral(0,h_data_scan[11]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width")));
	xsec_syst_b_jer = TMath::Abs(h_data_b_scan[12]->Integral(0,h_data_b_scan[12]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width"));
	xsec_syst_b_jer = TMath::Max(xsec_syst_b_jer,TMath::Abs(h_data_b_scan[11]->Integral(0,h_data_b_scan[11]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width")));

	TH1F* syst_pu = (TH1F*)h_data->Clone();
	TH1F* syst_b_pu = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_scan[3]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	  val = TMath::Max(val,TMath::Abs(h_data_scan[4]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i)));
	  syst_pu->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_b_scan[3]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	  val = TMath::Max(val,TMath::Abs(h_data_b_scan[4]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i)));
	  syst_b_pu->SetBinError(i, val);
	}
	double xsec_syst_pu = 0.0;
	double xsec_syst_b_pu = 0.0;
	xsec_syst_pu = TMath::Abs(h_data_scan[3]->Integral(0,h_data_scan[3]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width"));
	xsec_syst_pu = TMath::Max(xsec_syst_pu,TMath::Abs(h_data_scan[4]->Integral(0,h_data_scan[4]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width")));
	xsec_syst_b_pu = TMath::Abs(h_data_b_scan[3]->Integral(0,h_data_b_scan[3]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width"));
	xsec_syst_b_pu = TMath::Max(xsec_syst_b_pu,TMath::Abs(h_data_b_scan[4]->Integral(0,h_data_b_scan[4]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width")));

	TH1F* syst_bkg = (TH1F*)h_data->Clone();
	TH1F* syst_b_bkg = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_scan[10]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i))/0.1;
	  syst_bkg->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_b_scan[10]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i))/0.1;
	  syst_b_bkg->SetBinError(i, val);
	}
	double xsec_syst_bkg = 0.0;
	double xsec_syst_b_bkg = 0.0;
	xsec_syst_bkg = TMath::Abs(h_data_scan[10]->Integral(0,h_data_scan[10]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width"))/0.1;
	xsec_syst_b_bkg = TMath::Abs(h_data_b_scan[10]->Integral(0,h_data_b_scan[10]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width"))/0.1;

	TH1F* stat_top = (TH1F*)h_data->Clone();
	TH1F* stat_b_top = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_scan[5]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i))/0.1;
	  stat_top->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_b_scan[5]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i))/0.1;
	  stat_b_top->SetBinError(i, val);
	}
	double xsec_stat_top = 0.0;
	double xsec_stat_b_top = 0.0;
	xsec_stat_top = TMath::Abs(h_data_scan[5]->Integral(0,h_data_scan[5]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width"))/0.1;
	xsec_stat_b_top = TMath::Abs(h_data_b_scan[5]->Integral(0,h_data_b_scan[5]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width"))/0.1;

	TH1F* stat_bfit = (TH1F*)h_data->Clone();
	TH1F* stat_b_bfit = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_scan[6]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i))/0.1;
	  stat_bfit->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Abs(h_data_b_scan[6]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i))/0.1;
	  stat_b_bfit->SetBinError(i, val);
	}
	double xsec_stat_bfit = 0.0;
	double xsec_stat_b_bfit = 0.0;
	xsec_stat_bfit = TMath::Abs(h_data_scan[6]->Integral(0,h_data_scan[6]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width"))/0.1;
	xsec_stat_b_bfit = TMath::Abs(h_data_b_scan[6]->Integral(0,h_data_b_scan[6]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width"))/0.1;

	TH1F* syst_bfit2 = (TH1F*)h_data->Clone();
	TH1F* syst_b_bfit2 = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (useSysBfit2==1) {
	    val = TMath::Abs(h_data_scan[17]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Abs(TMath::Power(h_data_scan[17]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2))));
	  }
	  if (useSysBfit2==2) {
	    val = TMath::Abs(h_data_scan[18]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Abs(TMath::Power(h_data_scan[18]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2))));
	  }
	  if (useSysBfit2==3) {
	    val = TMath::Abs(h_data_scan[19]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Abs(TMath::Power(h_data_scan[19]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2))));
	  }
	  syst_bfit2->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (useSysBfit2==1) {
	    val = TMath::Abs(h_data_b_scan[17]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Abs(TMath::Power(h_data_b_scan[17]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))));
	  }
	  if (useSysBfit2==2) {
	    val = TMath::Abs(h_data_b_scan[18]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Abs(TMath::Power(h_data_b_scan[18]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))));
	  }
	  if (useSysBfit2==3) {
	    val = TMath::Abs(h_data_b_scan[19]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Abs(TMath::Power(h_data_b_scan[19]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))));
	  }
	  syst_b_bfit2->SetBinError(i, val);
	}
	double xsec_syst_bfit2 = 0.0;
	double xsec_syst_b_bfit2 = 0.0;
	if (useSysBfit2==1) {
	  xsec_syst_bfit2 = TMath::Abs(h_data_scan[17]->Integral(0,h_data_scan[17]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width"));
	  xsec_syst_b_bfit2 = TMath::Abs(h_data_b_scan[17]->Integral(0,h_data_b_scan[17]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width"));
	}
	if (useSysBfit2==2) {
	  xsec_syst_bfit2 = TMath::Abs(h_data_scan[18]->Integral(0,h_data_scan[18]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width"));
	  xsec_syst_b_bfit2 = TMath::Abs(h_data_b_scan[18]->Integral(0,h_data_b_scan[18]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width"));
	}
	if (useSysBfit2==3) {
	  xsec_syst_bfit2 = TMath::Abs(h_data_scan[19]->Integral(0,h_data_scan[19]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width"));
	  xsec_syst_b_bfit2 = TMath::Abs(h_data_b_scan[19]->Integral(0,h_data_b_scan[19]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width"));
	}

	TH1F* syst_btag = (TH1F*)h_data->Clone();
	TH1F* syst_b_btag = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = 0.0;
	  syst_btag->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (numB!=2) {
	    val = btag_sys * h_data_b_scan[0]->GetBinContent(i);
	  }
	  if (numB==2) {
	    val = 2.0 * btag_sys * h_data_b_scan[0]->GetBinContent(i);
	  }
	  syst_b_btag->SetBinError(i, val);
	}
	double xsec_syst_btag = 0.0;
	double xsec_syst_b_btag = 0.0;
	if (title.find("_b")!=string::npos) {
	  xsec_syst_btag = btag_sys * h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width");
	}
	if (numB!=2) {
	  xsec_syst_b_btag = btag_sys * h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width");
	}
	if (numB==2) {
	  xsec_syst_b_btag = 2.0 * btag_sys * h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width");
	}

	TH1F* stat_unfold = (TH1F*)h_data->Clone();
	TH1F* stat_b_unfold = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Sqrt(TMath::Max(0.,TMath::Power(h_data_scan[7]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2)));
	  stat_unfold->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = TMath::Sqrt(TMath::Max(0.,TMath::Power(h_data_b_scan[7]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2)));
	  stat_b_unfold->SetBinError(i, val);
	}
	double xsec_stat_unfold = 0.0;
	double xsec_stat_b_unfold = 0.0;
	double val7 = 0.0;
	double val0 = 0.0;
	h_data_scan[7]->IntegralAndError(0,h_data_scan[7]->GetNbinsX()+1,val7,"width");
	h_data_scan[0]->IntegralAndError(0,h_data_scan[0]->GetNbinsX()+1,val0,"width");
	xsec_stat_unfold = TMath::Sqrt(TMath::Max(0.,TMath::Power(val7,2)-TMath::Power(val0,2)));
	h_data_b_scan[7]->IntegralAndError(0,h_data_b_scan[7]->GetNbinsX()+1,val7,"width");
	h_data_b_scan[0]->IntegralAndError(0,h_data_b_scan[0]->GetNbinsX()+1,val0,"width");
	xsec_stat_b_unfold = TMath::Sqrt(TMath::Max(0.,TMath::Power(val7,2)-TMath::Power(val0,2)));

	TH1F* syst_unfold = (TH1F*)h_data->Clone();
	TH1F* syst_b_unfold = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (useSysUnfold==1) {
	    val = TMath::Abs(h_data_scan[8]->GetBinContent(i)-h_data_scan[7]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Max(0.,TMath::Power(h_data_scan[8]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2))-TMath::Max(0.,TMath::Power(h_data_scan[7]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2))));
	  }
	  if (useSysUnfold==2) {
	    val = TMath::Abs(h_data_scan[9]->GetBinContent(i)-h_data_scan[7]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Max(0.,TMath::Power(h_data_scan[9]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2))-TMath::Max(0.,TMath::Power(h_data_scan[7]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2))));
	  }
	  if (useSysUnfold==4) {
	    val = TMath::Abs(h_data_scan[15]->GetBinContent(i)-h_data_scan[0]->GetBinContent(i));
	  }
	  if (useSysUnfold==5) {
	    val = TMath::Abs(h_data_scan[16]->GetBinContent(i)-h_data_scan[7]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Max(0.,TMath::Power(h_data_scan[16]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2))-TMath::Max(0.,TMath::Power(h_data_scan[7]->GetBinError(i),2)-TMath::Power(h_data_scan[0]->GetBinError(i),2))));
	  }
	  syst_unfold->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  if (useSysUnfold==1) {
	    val = TMath::Abs(h_data_b_scan[8]->GetBinContent(i)-h_data_b_scan[7]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Max(0.,TMath::Power(h_data_b_scan[8]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))-TMath::Max(0.,TMath::Power(h_data_b_scan[7]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))));
	  }
	  if (useSysUnfold==2) {
	    val = TMath::Abs(h_data_b_scan[9]->GetBinContent(i)-h_data_b_scan[7]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Max(0.,TMath::Power(h_data_b_scan[9]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))-TMath::Max(0.,TMath::Power(h_data_b_scan[7]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))));
	  }
	  if (useSysUnfold==3) {
	    val = TMath::Abs(h_data_b_scan[14]->GetBinContent(i)-h_data_b_scan[7]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Max(0.,TMath::Power(h_data_b_scan[14]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))-TMath::Max(0.,TMath::Power(h_data_b_scan[7]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))));
	  }
	  if (useSysUnfold==4) {
	    val = TMath::Abs(h_data_b_scan[15]->GetBinContent(i)-h_data_b_scan[0]->GetBinContent(i));
	  }
	  if (useSysUnfold==5) {
	    val = TMath::Abs(h_data_b_scan[16]->GetBinContent(i)-h_data_b_scan[7]->GetBinContent(i));
	    val = TMath::Sqrt(TMath::Max(0.,TMath::Power(val,2)-TMath::Max(0.,TMath::Power(h_data_b_scan[16]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))-TMath::Max(0.,TMath::Power(h_data_b_scan[7]->GetBinError(i),2)-TMath::Power(h_data_b_scan[0]->GetBinError(i),2))));
	  }
	  syst_b_unfold->SetBinError(i, val);
	}
	double xsec_syst_unfold = 0.0;
	double xsec_syst_b_unfold = 0.0;
	if (useSysUnfold==1) {
	  xsec_syst_unfold = TMath::Abs(h_data_scan[8]->Integral(0,h_data_scan[8]->GetNbinsX()+1,"width")-h_data_scan[7]->Integral(0,h_data_scan[7]->GetNbinsX()+1,"width"));
	  xsec_syst_b_unfold = TMath::Abs(h_data_b_scan[8]->Integral(0,h_data_b_scan[8]->GetNbinsX()+1,"width")-h_data_b_scan[7]->Integral(0,h_data_b_scan[7]->GetNbinsX()+1,"width"));
	}
	if (useSysUnfold==2) {
	  xsec_syst_unfold = TMath::Abs(h_data_scan[9]->Integral(0,h_data_scan[9]->GetNbinsX()+1,"width")-h_data_scan[7]->Integral(0,h_data_scan[7]->GetNbinsX()+1,"width"));
	  xsec_syst_b_unfold = TMath::Abs(h_data_b_scan[9]->Integral(0,h_data_b_scan[9]->GetNbinsX()+1,"width")-h_data_b_scan[7]->Integral(0,h_data_b_scan[7]->GetNbinsX()+1,"width"));
	}
	if (useSysUnfold==3) {
	  xsec_syst_unfold = 0.0;
	  xsec_syst_b_unfold = TMath::Abs(h_data_b_scan[14]->Integral(0,h_data_b_scan[14]->GetNbinsX()+1,"width")-h_data_b_scan[7]->Integral(0,h_data_b_scan[7]->GetNbinsX()+1,"width"));
	}
	if (useSysUnfold==4) {
	  xsec_syst_unfold = TMath::Abs(h_data_scan[15]->Integral(0,h_data_scan[15]->GetNbinsX()+1,"width")-h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width"));
	  xsec_syst_b_unfold = TMath::Abs(h_data_b_scan[15]->Integral(0,h_data_b_scan[15]->GetNbinsX()+1,"width")-h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width"));
	}
	if (useSysUnfold==5) {
	  xsec_syst_unfold = TMath::Abs(h_data_scan[16]->Integral(0,h_data_scan[16]->GetNbinsX()+1,"width")-h_data_scan[7]->Integral(0,h_data_scan[7]->GetNbinsX()+1,"width"));
	  xsec_syst_b_unfold = TMath::Abs(h_data_b_scan[16]->Integral(0,h_data_b_scan[16]->GetNbinsX()+1,"width")-h_data_b_scan[7]->Integral(0,h_data_b_scan[7]->GetNbinsX()+1,"width"));
	}

	TH1F* syst_lumi = (TH1F*)h_data->Clone();
	TH1F* syst_b_lumi = (TH1F*)h_data_b->Clone();
	for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = lumi_sys * h_data_scan[0]->GetBinContent(i);
	  if (isratio) val = 0.0;
	  syst_lumi->SetBinError(i, val);
	}
	for (int i=0;i<=h_data_b->GetNbinsX()+1;i++) {
	  double val = 0.0;
	  val = lumi_sys * h_data_b_scan[0]->GetBinContent(i);
	  if (isratio) val = 0.0;
	  syst_b_lumi->SetBinError(i, val);
	}
	double xsec_syst_lumi = 0.0;
	double xsec_syst_b_lumi = 0.0;
	xsec_syst_lumi = lumi_sys * h_data_scan[0]->Integral(0,h_data_scan[0]->GetNbinsX()+1,"width");
	xsec_syst_b_lumi = lumi_sys * h_data_b_scan[0]->Integral(0,h_data_b_scan[0]->GetNbinsX()+1,"width");

	float sum1, sum2, sum3, sum4;
	float sum1_b, sum2_b, sum3_b, sum4_b;
	ifstream in1, in2, in3, in4;
	if (ilepton==1) {
	  if (numB==0) {
	    in1.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_first_jet_pt" + "_xsecs_unfolding.dat").c_str());
	    in2.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_first_jet_eta" + "_xsecs_unfolding.dat").c_str());
	    in3.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_pt_Z_ee" + "_xsecs_unfolding.dat").c_str());
	    in4.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_Ht" + "_xsecs_unfolding.dat").c_str());
	  }
	  if (numB!=0) {
	    in1.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_first_bjet_pt" + "_xsecs_unfolding.dat").c_str());
	    in2.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_first_bjet_eta" + "_xsecs_unfolding.dat").c_str());
	    in3.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_pt_Z_ee_b" + "_xsecs_unfolding.dat").c_str());
	    in4.open((path + "/electrons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_Ht_b" + "_xsecs_unfolding.dat").c_str());
	  }
	}
	if (ilepton==2) {
	  if (numB==0) {
	    in1.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_first_jet_pt" + "_xsecs_unfolding.dat").c_str());
	    in2.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_first_jet_eta" + "_xsecs_unfolding.dat").c_str());
	    in3.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_pt_Z_mm" + "_xsecs_unfolding.dat").c_str());
	    in4.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_Ht" + "_xsecs_unfolding.dat").c_str());
	  }
	  if (numB!=0) {
	    in1.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_first_bjet_pt" + "_xsecs_unfolding.dat").c_str());
	    in2.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_first_bjet_eta" + "_xsecs_unfolding.dat").c_str());
	    in3.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_pt_Z_mm_b" + "_xsecs_unfolding.dat").c_str());
	    in4.open((path + "/muons/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + "w_Ht_b" + "_xsecs_unfolding.dat").c_str());
	  }
	}
	in1 >> sum1; in1 >> sum1_b;
	in2 >> sum2; in2 >> sum2_b;
	in3 >> sum3; in3 >> sum3_b;
	in4 >> sum4; in4 >> sum4_b;
	in1.close();
	in2.close();
	in3.close();
	in4.close();

	float tot = (sum1+sum2+sum3+sum4)/4.;
	float tot_b = (sum1_b+sum2_b+sum3_b+sum4_b)/4.;
	float rms = TMath::Sqrt((TMath::Power(sum1-tot,2)+TMath::Power(sum2-tot,2)+TMath::Power(sum3-tot,2)+TMath::Power(sum4-tot,2))/(4-1));
	float rms_b = TMath::Sqrt((TMath::Power(sum1_b-tot_b,2)+TMath::Power(sum2_b-tot_b,2)+TMath::Power(sum3_b-tot_b,2)+TMath::Power(sum4_b-tot_b,2))/(4-1));

	if (isratio==1) {
	  float tmp1 = (tot_b/tot);
	  float tmp2 = (tot_b/tot)*TMath::Sqrt(TMath::Power(rms/tot,2)+TMath::Power(rms_b/tot_b,2));
	  tot_b = tmp1;
	  rms_b = tmp2;
	}

	TH1F* h_data_stat = (TH1F*)h_data->Clone();
	TH1F* h_data_b_stat = (TH1F*)h_data_b->Clone();
	TH1F* h_data_syst = (TH1F*)h_data->Clone();
	TH1F* h_data_b_syst = (TH1F*)h_data_b->Clone();
	TH1F* h_data_tot = (TH1F*)h_data->Clone();
	TH1F* h_data_b_tot = (TH1F*)h_data_b->Clone();

	for (int i=0;i<=h_data_stat->GetNbinsX()+1;i++) {
	  h_data_stat->SetBinError(i, TMath::Sqrt(TMath::Power(h_data_stat->GetBinError(i),2)+TMath::Power(stat_top->GetBinError(i),2)));
	  h_data_stat->SetBinError(i, TMath::Sqrt(TMath::Power(h_data_stat->GetBinError(i),2)+TMath::Power(stat_bfit->GetBinError(i),2)));
	  double val = 0.0;
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_bkg->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_eff->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_jec->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_jer->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_pu->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_bkg->GetBinError(i),2));
	  if (useSysBfit2) val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_bfit2->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_btag->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_unfold->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_unfold->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_lumi->GetBinError(i),2));
	  if (useSysRMS) val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(h_data_stat->GetBinContent(i)*rms/tot,2));
	  h_data_syst->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Power(h_data_stat->GetBinError(i),2)+TMath::Power(h_data_syst->GetBinError(i),2));
	  h_data_tot->SetBinError(i, val);
	}

	for (int i=0;i<=h_data_b_stat->GetNbinsX()+1;i++) {
	  h_data_b_stat->SetBinError(i, TMath::Sqrt(TMath::Power(h_data_b_stat->GetBinError(i),2)+TMath::Power(stat_b_top->GetBinError(i),2)));
	  h_data_b_stat->SetBinError(i, TMath::Sqrt(TMath::Power(h_data_b_stat->GetBinError(i),2)+TMath::Power(stat_b_bfit->GetBinError(i),2)));
	  double val = 0.0;
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_b_bkg->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_eff->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_jec->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_jer->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_pu->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_bkg->GetBinError(i),2));
	  if (useSysBfit2) val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_bfit2->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_btag->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(stat_b_unfold->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_unfold->GetBinError(i),2));
	  val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(syst_b_lumi->GetBinError(i),2));
	  if (useSysRMS) val = TMath::Sqrt(TMath::Power(val,2)+TMath::Power(h_data_b_stat->GetBinContent(i)*rms_b/tot_b,2));
	  h_data_b_syst->SetBinError(i, val);
	  val = TMath::Sqrt(TMath::Power(h_data_b_stat->GetBinError(i),2)+TMath::Power(h_data_b_syst->GetBinError(i),2));
	  h_data_b_tot->SetBinError(i, val);
	}

	double xsec_data_tot_stat = 0.0;
	double xsec_data_b_tot_stat = 0.0;
	double xsec_data_tot_syst = 0;
	double xsec_data_b_tot_syst = 0;
	double xsec_data_tot_tot = 0;
	double xsec_data_b_tot_tot = 0;
	xsec_data_tot_stat = TMath::Sqrt(TMath::Power(xsec_data_tot_stat,2)+TMath::Power(xsec_stat_data,2));
	xsec_data_tot_stat = TMath::Sqrt(TMath::Power(xsec_data_tot_stat,2)+TMath::Power(xsec_stat_top,2));
	xsec_data_tot_stat = TMath::Sqrt(TMath::Power(xsec_data_tot_stat,2)+TMath::Power(xsec_stat_bfit,2));
	xsec_data_b_tot_stat = TMath::Sqrt(TMath::Power(xsec_data_b_tot_stat,2)+TMath::Power(xsec_stat_data_b,2));
	xsec_data_b_tot_stat = TMath::Sqrt(TMath::Power(xsec_data_b_tot_stat,2)+TMath::Power(xsec_stat_b_top,2));
	xsec_data_b_tot_stat = TMath::Sqrt(TMath::Power(xsec_data_b_tot_stat,2)+TMath::Power(xsec_stat_b_bfit,2));
	xsec_data_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_tot_syst,2)+TMath::Power(xsec_stat_bkg,2));
	xsec_data_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_tot_syst,2)+TMath::Power(xsec_syst_eff,2));
	xsec_data_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_tot_syst,2)+TMath::Power(xsec_syst_jec,2));
	xsec_data_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_tot_syst,2)+TMath::Power(xsec_syst_jer,2));
	xsec_data_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_tot_syst,2)+TMath::Power(xsec_syst_pu,2));
	xsec_data_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_tot_syst,2)+TMath::Power(xsec_syst_bkg,2));
	if (useSysBfit2) xsec_data_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_tot_syst,2)+TMath::Power(xsec_syst_bfit2,2));
	xsec_data_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_tot_syst,2)+TMath::Power(xsec_syst_btag,2));
	xsec_data_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_tot_syst,2)+TMath::Power(xsec_stat_unfold,2));
	xsec_data_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_tot_syst,2)+TMath::Power(xsec_syst_unfold,2));
	xsec_data_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_tot_syst,2)+TMath::Power(xsec_syst_lumi,2));
	if (useSysRMS) xsec_data_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_tot_syst,2)+TMath::Power(xsec_data*(tot>0 ? rms/tot : 0),2));
	xsec_data_b_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_b_tot_syst,2)+TMath::Power(xsec_stat_b_bkg,2));
	xsec_data_b_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_b_tot_syst,2)+TMath::Power(xsec_syst_b_eff,2));
	xsec_data_b_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_b_tot_syst,2)+TMath::Power(xsec_syst_b_jec,2));
	xsec_data_b_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_b_tot_syst,2)+TMath::Power(xsec_syst_b_jer,2));
	xsec_data_b_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_b_tot_syst,2)+TMath::Power(xsec_syst_b_pu,2));
	xsec_data_b_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_b_tot_syst,2)+TMath::Power(xsec_syst_b_bkg,2));
	if (useSysBfit2) xsec_data_b_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_b_tot_syst,2)+TMath::Power(xsec_syst_b_bfit2,2));
	xsec_data_b_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_b_tot_syst,2)+TMath::Power(xsec_syst_b_btag,2));
	xsec_data_b_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_b_tot_syst,2)+TMath::Power(xsec_stat_b_unfold,2));
	xsec_data_b_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_b_tot_syst,2)+TMath::Power(xsec_syst_b_unfold,2));
	xsec_data_b_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_b_tot_syst,2)+TMath::Power(xsec_syst_b_lumi,2));
	if (useSysRMS) xsec_data_b_tot_syst = TMath::Sqrt(TMath::Power(xsec_data_b_tot_syst,2)+TMath::Power(xsec_data_b*(tot_b>0 ? rms_b/tot_b : 0),2));
	xsec_data_tot_tot = TMath::Sqrt(TMath::Power(xsec_data_tot_stat,2)+TMath::Power(xsec_data_tot_syst,2));
	xsec_data_b_tot_tot = TMath::Sqrt(TMath::Power(xsec_data_b_tot_stat,2)+TMath::Power(xsec_data_b_tot_syst,2));

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

	  h_mc1b_b->SetMaximum(4*h_data_tot->GetMaximum());
	  h_mc1b_b->SetMinimum(TMath::Max(0.000002,0.25*h_data_b_tot->GetBinContent(h_data_b_tot->GetMinimumBin())));

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

	  h_mcg3_b->Draw("E5SAME");
	  TH1F* tmp2_3 = (TH1F*)h_mcg3_b->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp2_3->GetMinimum()==0) tmp2_3->GetXaxis()->SetRangeUser(0, tmp2_3->GetBinCenter(tmp2_3->GetMinimumBin()-1));
	  }
	  tmp2_3->SetFillColor(0);
	  tmp2_3->DrawClone("HISTLSAME");

	  h_data_b_tot->Draw("E1PX0SAME");
	  h_data_b_stat->Draw("E1PX0SAME");

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

	  if (ilepton==1) {
	    if (drawInclusive) leg->AddEntry(h_data_stat,"Z(#rightarrow ee)+j DATA","p");
	    leg->AddEntry(h_data_b_stat,"Z(#rightarrow ee)+b DATA","p");
	    if (useMC) leg->AddEntry(h_mc1,"Z(#rightarrow ee) MC","l");
	    leg->AddEntry(h_mcg,"Z(#rightarrow ee) MadGraph 5FS","l");
	    leg->AddEntry(h_mcg3,"Z(#rightarrow ee) MadGraph 4FS","l");
	    if (useSherpa) leg->AddEntry(h_mcg1,"Z(#rightarrow ee) Sherpa","l");
	    leg->AddEntry(h_mcg2,"Z(#rightarrow ee) Powheg","l");
	  }
	  if (ilepton==2){
	    if (drawInclusive) leg->AddEntry(h_data_stat,"Z(#rightarrow #mu#mu)+j DATA","p");
	    leg->AddEntry(h_data_b_stat,"Z(#rightarrow #mu#mu)+b DATA","p");
	    if (useMC) leg->AddEntry(h_mc1,"Z(#rightarrow #mu#mu) MC","l");
	    leg->AddEntry(h_mcg,"Z(#rightarrow #mu#mu) MadGraph 5FS","l");
	    leg->AddEntry(h_mcg3,"Z(#rightarrow #mu#mu) MadGraph 4FS","l");
	    if (useSherpa) leg->AddEntry(h_mcg1,"Z(#rightarrow #mu#mu) Sherpa","l");
	    leg->AddEntry(h_mcg2,"Z(#rightarrow #mu#mu) Powheg","l");
	  }
	}

	if (isratio==1) {
	  if (ilepton==1) {
	    leg->AddEntry(h_data_b_stat,"[Z(#rightarrow ee)+b] / [Z(#rightarrow ee)+j] DATA","p");
	    if (useMC) leg->AddEntry(h_mc1b_b,"[Z(#rightarrow ee)+b] / [Z(#rightarrow ee)+j] MC","l");
	    leg->AddEntry(h_mcg_b,"[Z(#rightarrow ee)+b] / [Z(#rightarrow ee)+j] MadGraph 5FS","l");
	    leg->AddEntry(h_mcg3_b,"[Z(#rightarrow ee)+b] / [Z(#rightarrow ee)+j] MadGraph 4FS","l");
	    if (useSherpa) leg->AddEntry(h_mcg1_b,"[Z(#rightarrow ee)+b] / [Z(#rightarrow ee)+j] Sherpa","l");
	    leg->AddEntry(h_mcg2_b,"[Z(#rightarrow ee)+b] / [Z(#rightarrow ee)+j] Powheg","l");
	  }
	  if (ilepton==2){
	    leg->AddEntry(h_data_b_stat,"[Z(#rightarrow #mu#mu)+b] / [Z(#rightarrow #mu#mu)+j] DATA","p");
	    if (useMC) leg->AddEntry(h_mc1b_b,"[Z(#rightarrow #mu#mu)+b] / [Z(#rightarrow #mu#mu)+j] MC","l");
	    leg->AddEntry(h_mcg_b,"[Z(#rightarrow #mu#mu)+b] / [Z(#rightarrow #mu#mu)+j] MadGraph 5FS","l");
	    leg->AddEntry(h_mcg3_b,"[Z(#rightarrow #mu#mu)+b] / [Z(#rightarrow #mu#mu)+j] MadGraph 4FS","l");
	    if (useSherpa) leg->AddEntry(h_mcg1_b,"[Z(#rightarrow #mu#mu)+b] / [Z(#rightarrow #mu#mu)+j] Sherpa","l");
	    leg->AddEntry(h_mcg2_b,"[Z(#rightarrow #mu#mu)+b] / [Z(#rightarrow #mu#mu)+j] Powheg","l");
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
	  t2->DrawLatex(0.15,0.15,"MadGraph 5FS / MadGraph 4FS");
	} else {
	  t2->DrawLatex(0.15,0.13,"MadGraph 5FS, normalized to  #sigma_{NNLO}");
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
	  t3->DrawLatex(0.15,0.15,"Sherpa");
	} else {
	  t3->DrawLatex(0.15,0.13,"MadGraph 4FS, normalized to  #sigma_{LO}");
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
	  if (drawInclusive) dx = 0.1*(g_P2_tot->GetXaxis()->GetXmax()-g_P2_tot->GetXaxis()->GetXmin())/g_P2_tot->GetN();
	  for (int i=0; i<g_P2_tot->GetN(); i++) {
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
	t4->DrawLatex(0.15,0.40,"Powheg, normalized to  #sigma_{NLO}");

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

	TLine *OLine5 = new TLine(h_P_tot->GetXaxis()->GetXmin(),0.93,h_P_tot->GetXaxis()->GetXmax(),0.93);
	OLine5->SetLineColor(kOrange+7);
	OLine5->SetLineWidth(2);
	OLine5->Draw();

	c1->cd();

	if (isratio==1) {
	  h_mc1b_b->GetYaxis()->SetRangeUser(-0.5, 10);
	}

	if (title_b=="w_first_jet_pt_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dp_{T} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("leading jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp_{T} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(-0.5, 20);
	  }
	} else if (title_b=="w_first_jet_eta_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#eta [pb]");
	  h_P_tot->GetXaxis()->SetTitle("leading jet #eta");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_first_bjet_pt") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dp^{b}_{T} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{b}_{T} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(-0.5, 10);
	  }
	} else if (title_b=="w_second_bjet_pt") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dp^{b}_{T} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("subleading b-jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{b}_{T} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(-0.5, 10);
	  }
	} else if (title_b=="w_first_bjet_eta") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#eta^{b} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("leading b-jet #eta");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta^{b} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_second_bjet_eta") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#eta^{b} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("subleading b-jet #eta");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta^{b} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_first_bjet_eta_abs") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d|#eta^{b}| [pb]");
	  h_P_tot->GetXaxis()->SetTitle("leading b-jet |#eta|");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d|#eta^{b}| [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_second_bjet_eta_abs") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d|#eta^{b}| [pb]");
	  h_P_tot->GetXaxis()->SetTitle("subleading b-jet |#eta|");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d|#eta^{b}| [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_pt_Z_ee_b"||title_b =="w_pt_Z_mm_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dp^{Z}_{T} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("Z boson p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{Z}_{T} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(-0.5, 20);
	  }
	} else if (title_b=="w_Ht_b") {
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / dH_{T} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("H_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dH_{T} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_delta_phi_ee_b" || title_b=="w_delta_phi_mm_b") {
	  leg->SetX1(0.12);
	  leg->SetX2(0.38);
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#Delta#phi_{bZ} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta#phi(bZ) [rad]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta#phi_{Zb} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_delta_phi_2b") {
	  leg->SetX1(0.12);
	  leg->SetX2(0.38);
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#Delta#phi_{bb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta#phi(bb) [rad]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta#phi_{Zb} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_DR_eeb_min" || title_b=="w_DR_mmb_min") {
	  leg->SetX1(0.12);
	  leg->SetX2(0.38);
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#Delta R^{min}_{bb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta R^{min}(bb) [rad]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta R^{min}_{Zb} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_DR_eeb_max" || title_b=="w_DR_mmb_max") {
	  leg->SetX1(0.12);
	  leg->SetX2(0.38);
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#Delta R^{max}_{bb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta R^{max}(bb) [rad]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta R^{max}_{Zb} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_DR_bb") {
	  leg->SetX1(0.12);
	  leg->SetX2(0.38);
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d#Delta R{bb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta R(bb) [rad]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta R_{bb} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_eebb_mass" || title_b=="w_mmbb_mass") {
	  leg->SetX1(0.12);
	  leg->SetX2(0.38);
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d M_{Zbb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("M(Zbb) [GeV]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d M_{Zbb} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_eebb_mass" || title_b=="w_mmbb_mass") {
	  leg->SetX1(0.12);
	  leg->SetX2(0.38);
	  h_mc1b_b->GetYaxis()->SetTitle("d#sigma / d M_{bb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("M(bb) [GeV]");
	  if (isratio==1) {
	    h_mc1b_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d M_{bb} [%]");
	    h_mc1b_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	}

	if (plot) {
	  ofstream out, out1;
	  if (isratio==0) {
	    if (ilepton==1) {
	      gSystem->mkdir((path + "/electrons/" + version + "/xsecs_unfolding" + dirbSel + "/").c_str(), kTRUE);
	      c1->SaveAs((path + "/electrons/" + version + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.pdf").c_str());
	      out.open((path + "/electrons/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.dat").c_str());
	      out1.open((path + "/electrons/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.txt").c_str());
	    }
	    if (ilepton==2) {
	      gSystem->mkdir((path + "/muons/" + version + "/xsecs_unfolding" + dirbSel + "/").c_str(), kTRUE);
	      c1->SaveAs((path + "/muons/" + version + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.pdf").c_str());
	      out.open((path + "/muons/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.dat").c_str());
	      out1.open((path + "/muons/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.txt").c_str());
	    }
	  }
	  if (isratio==1) {
	    if (ilepton==1) {
	      gSystem->mkdir((path + "/electrons/" + version + "/ratios_unfolding" + dirbSel + "/").c_str(), kTRUE);
	      c1->SaveAs((path + "/electrons/" + version + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.pdf").c_str());
	      out.open((path + "/electrons/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.dat").c_str());
	      out1.open((path + "/electrons/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.txt").c_str());
	    }
	    if (ilepton==2) {
	      gSystem->mkdir((path + "/muons/" + version + "/ratios_unfolding" + dirbSel + "/").c_str(), kTRUE);
	      c1->SaveAs((path + "/muons/" + version + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.pdf").c_str());
	      out.open((path + "/muons/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.dat").c_str());
	      out1.open((path + "/muons/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.txt").c_str());
	    }
	  }
	  if (isratio==0 && drawInclusive) {
	    out << h_data->GetName();
	    out << std::fixed << std::setprecision(4);
	    out << " : average unfolded total cross section = " << tot << " +- " << rms << " pb (" << 100*(rms/tot) << " %)";
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
	    if (useSysBfit2) out << std::setw(14) << "btemp";
	    out << std::setw(14) << "btag";
	    out << std::setw(14) << "unfold";
	    out << std::setw(14) << "unfold";
	    out << std::setw(14) << "lumi";
	    if (useSysRMS) out << std::setw(14) << "unfold";
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
	    if (useSysBfit2) out << std::setw(14) << "syst";
	    out << std::setw(14) << "syst";
	    out << std::setw(14) << "stat";
	    out << std::setw(14) << "syst";
	    out << std::setw(14) << "syst";
	    if (useSysRMS) out << std::setw(14) << "rms";
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
	      if (useSysBfit2) {
		out << " +- ";
		out << std::setw(10) << syst_bfit2->GetBinError(i);
	      }
	      out << " +- ";
	      out << std::setw(10) << syst_btag->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << stat_unfold->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << syst_unfold->GetBinError(i);
	      out << " +- ";
	      out << std::setw(10) << syst_lumi->GetBinError(i);
	      if (useSysRMS) {
		out << " +- ";
		out << std::setw(10) << h_data->GetBinContent(i)*(rms/tot);
	      }
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
	    out << std::setw(12) << xsec_data;
	    out << " +- ";
	    out << std::setw(10) << xsec_stat_data;
	    out << " +- ";
	    out << std::setw(10) << xsec_stat_bkg;
	    out << " +- ";
	    out << std::setw(10) << xsec_syst_eff;
	    out << " +- ";
	    out << std::setw(10) << xsec_syst_jec;
	    out << " +- ";
	    out << std::setw(10) << xsec_syst_jer;
	    out << " +- ";
	    out << std::setw(10) << xsec_syst_pu;
	    out << " +- ";
	    out << std::setw(10) << xsec_syst_bkg;
	    out << " +- ";
	    out << std::setw(10) << xsec_stat_top;
	    out << " +- ";
	    out << std::setw(10) << xsec_stat_bfit;
	    if (useSysBfit2) {
	      out << " +- ";
	      out << std::setw(10) << xsec_syst_bfit2;
	    }
	    out << " +- ";
	    out << std::setw(10) << xsec_syst_btag;
	    out << " +- ";
	    out << std::setw(10) << xsec_stat_unfold;
	    out << " +- ";
	    out << std::setw(10) << xsec_syst_unfold;
	    out << " +- ";
	    out << std::setw(10) << xsec_syst_lumi;
	    if (useSysRMS) {
	      out << " +- ";
	      out << std::setw(10) << xsec_data*(tot>0 ? rms/tot : 0);
	    }
	    out << " => ";
	    out << std::setw(10) << xsec_data_tot_stat;
	    out << " +- ";
	    out << std::setw(10) << xsec_data_tot_syst;
	    out << " => ";
	    out << std::setw(10) << xsec_data_tot_tot;
	    out << " => ";
	    out << std::setprecision(1);
	    out << std::setw(4) << TMath::Abs(100.*(xsec_data==0 ? 0 : xsec_data_tot_tot/xsec_data));
	    out << endl;
	  }
	  out << h_data_b->GetName();
	  if (isratio==0) {
	    out << std::fixed << std::setprecision(4);
	    out << " : average unfolded total cross section = " << tot_b << " +- " << rms_b << " pb (" << 100*(rms_b/tot_b) << " %)";
	  }
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
	  if (useSysBfit2) out << std::setw(14) << "btemp";
	  out << std::setw(14) << "btag";
	  out << std::setw(14) << "unfold";
	  out << std::setw(14) << "unfold";
	  out << std::setw(14) << "lumi";
	  if (useSysRMS) out << std::setw(14) << "unfold";
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
	  if (useSysBfit2) out << std::setw(14) << "syst";
	  out << std::setw(14) << "syst";
	  out << std::setw(14) << "stat";
	  out << std::setw(14) << "syst";
	  out << std::setw(14) << "syst";
	  if (useSysRMS) out << std::setw(1.) << "rms";
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
	    if (useSysBfit2) {
	      out << " +- ";
	      out << std::setw(10) << syst_b_bfit2->GetBinError(i);
	    }
	    out << " +- ";
	    out << std::setw(10) << syst_b_btag->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << stat_b_unfold->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << syst_b_unfold->GetBinError(i);
	    out << " +- ";
	    out << std::setw(10) << syst_b_lumi->GetBinError(i);
	    if (useSysRMS) {
	      out << " +- ";
	      out << std::setw(10) << h_data_b->GetBinContent(i)*(rms_b/tot_b);
	    }
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
	  out << std::setw(12) << xsec_data_b;
	  out << " +- ";
	  out << std::setw(10) << xsec_stat_data_b;
	  out << " +- ";
	  out << std::setw(10) << xsec_stat_b_bkg;
	  out << " +- ";
	  out << std::setw(10) << xsec_syst_b_eff;
	  out << " +- ";
	  out << std::setw(10) << xsec_syst_b_jec;
	  out << " +- ";
	  out << std::setw(10) << xsec_syst_b_jer;
	  out << " +- ";
	  out << std::setw(10) << xsec_syst_b_pu;
	  out << " +- ";
	  out << std::setw(10) << xsec_syst_b_bkg;
	  out << " +- ";
	  out << std::setw(10) << xsec_stat_b_top;
	  out << " +- ";
	  out << std::setw(10) << xsec_stat_b_bfit;
	  if (useSysBfit2) {
	    out << " +- ";
	    out << std::setw(10) << xsec_syst_b_bfit2;
	  }
	  out << " +- ";
	  out << std::setw(10) << xsec_syst_b_btag;
	  out << " +- ";
	  out << std::setw(10) << xsec_stat_b_unfold;
	  out << " +- ";
	  out << std::setw(10) << xsec_syst_b_unfold;
	  out << " +- ";
	  out << std::setw(10) << xsec_syst_b_lumi;
	  if (useSysRMS) {
	    out << " +- ";
	    out << std::setw(10) << xsec_data_b*(tot_b>0 ? rms_b/tot_b : 0);
	  }
	  out << " => ";
	  out << std::setw(10) << xsec_data_b_tot_stat;
	  out << " +- ";
	  out << std::setw(10) << xsec_data_b_tot_syst;
	  out << " => ";
	  out << std::setw(10) << xsec_data_b_tot_tot;
	  out << " => ";
	  out << std::setprecision(1);
	  out << std::setw(4) << TMath::Abs(100.*(xsec_data_b==0 ? 0 : xsec_data_b_tot_tot/xsec_data_b));
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
	    if (useSysBfit2) out1 << std::setw(9) << "btemp";
	    out1 << std::setw(9) << "btag";
	    out1 << std::setw(9) << "unfold";
	    out1 << std::setw(9) << "unfold";
	    out1 << std::setw(9) << "lumi";
	    if (useSysRMS) out1 << std::setw(9) << "unfold";
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
	    if (useSysBfit2) out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "stat";
	    out1 << std::setw(9) << "syst";
	    out1 << std::setw(9) << "syst";
	    if (useSysRMS) out1 << std::setw(9) << "rms";
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
	      if (useSysBfit2) {
		out1 << " +- ";
		out1 << std::setw(5) << syst_bfit2->GetBinError(i)*val;
	      }
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_btag->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << stat_unfold->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_unfold->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_lumi->GetBinError(i)*val;
	      if (useSysRMS) {
		out1 << " +- ";
		out1 << std::setw(5) << h_data->GetBinContent(i)*(rms/tot)*val;
	      }
	      out1 << " => ";
	      out1 << std::setw(5) << h_data_stat->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(5) << h_data_syst->GetBinError(i)*val;
	      out1 << " => ";
	      out1 << std::setw(5) << h_data_tot->GetBinError(i)*val;
	      out1 << endl;
	    }
	    double val = TMath::Abs(100.*(xsec_data==0 ? 0 : 1./xsec_data));
	    out1 << "tot";
	    out1 << std::setprecision(1);
	    out1 << std::setw(5) << xsec_stat_data*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_stat_bkg*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_syst_eff*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_syst_jec*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_syst_jer*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_syst_pu*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_syst_bkg*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_stat_top*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_stat_bfit*val;
	    if (useSysBfit2) {
	      out1 << " +- ";
	      out1 << std::setw(5) << xsec_syst_bfit2*val;
	    }
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_syst_btag*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_stat_unfold*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_syst_unfold*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_syst_lumi*val;
	    if (useSysRMS) {
	      out1 << " +- ";
	      out1 << std::setw(5) << xsec_data*(tot>0 ? rms/tot : 0)*val;
	    }
	    out1 << " => ";
	    out1 << std::setw(5) << xsec_data_tot_stat*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_data_tot_syst*val;
	    out1 << " => ";
	    out1 << std::setw(5) << xsec_data_tot_tot*val;
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
	  if (useSysBfit2) out1 << std::setw(9) << "btemp";
	  out1 << std::setw(9) << "btag";
	  out1 << std::setw(9) << "unfold";
	  out1 << std::setw(9) << "unfold";
	  out1 << std::setw(9) << "lumi";
	  if (useSysRMS) out1 << std::setw(9) << "unfold";
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
	  if (useSysBfit2) out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "stat";
	  out1 << std::setw(9) << "syst";
	  out1 << std::setw(9) << "syst";
	  if (useSysRMS) out1 << std::setw(9) << "rms";
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
	    if (useSysBfit2) {
	      out1 << " +- ";
	      out1 << std::setw(5) << syst_b_bfit2->GetBinError(i)*val;
	    }
	    out1 << " +- ";
	    out1 << std::setw(5) << syst_b_btag->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << stat_b_unfold->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << syst_b_unfold->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << syst_b_lumi->GetBinError(i)*val;
	    if (useSysRMS) {
	      out1 << " +- ";
	      out1 << std::setw(5) << h_data_b->GetBinContent(i)*(rms_b/tot_b)*val;
	    }
	    out1 << " => ";
	    out1 << std::setw(5) << h_data_b_stat->GetBinError(i)*val;
	    out1 << " +- ";
	    out1 << std::setw(5) << h_data_b_syst->GetBinError(i)*val;
	    out1 << " => ";
	    out1 << std::setw(5) << h_data_b_tot->GetBinError(i)*val;
	    out1 << endl;
	  }
	  double val = TMath::Abs(100.*(xsec_data_b==0 ? 0 : 1./xsec_data_b));
	  out1 << "tot";
	  out1 << std::setprecision(1);
	  out1 << std::setw(5) << xsec_stat_data_b*val;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_stat_b_bkg*val;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_syst_b_eff*val;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_syst_b_jec*val;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_syst_b_jer*val;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_syst_b_pu*val;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_syst_b_bkg*val;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_stat_b_top*val;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_stat_b_bfit*val;
	  if (useSysBfit2) {
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_syst_b_bfit2*val;
	  }
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_syst_b_btag*val;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_stat_b_unfold*val;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_syst_b_unfold*val;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_syst_b_lumi*val;
	  if (useSysRMS) {
	    out1 << " +- ";
	    out1 << std::setw(5) << xsec_data_b*(tot>0 ? rms_b/tot_b : 0)*val;
	  }
	  out1 << " => ";
	  out1 << std::setw(5) << xsec_data_b_tot_stat*val;
	  out1 << " +- ";
	  out1 << std::setw(5) << xsec_data_b_tot_syst*val;
	  out1 << " => ";
	  out1 << std::setw(5) << xsec_data_b_tot_tot*val;
	  out1 << endl;
	  out1.close();
	}
}

