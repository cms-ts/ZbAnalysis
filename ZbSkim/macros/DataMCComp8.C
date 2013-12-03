#include "LumiLabel.C"
#include "LumiInfo_v11.h"

#include "fixrange.C"
#include <iomanip>

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";

TH1F* read(string& path, string& subdir, string& title, int ilepton, TFile* infile=0) {
  TH1F* hist;
  TFile* file = infile;
  string title_tmp = title;
  if (ilepton==1) {
    if (title=="w_pt_Z") title_tmp="w_pt_Z_ee";
    if (title=="w_pt_Z_b") title_tmp="w_pt_Z_ee_b";
    if (title=="w_delta_phi") title_tmp="w_delta_phi_ee";
    if (title=="w_delta_phi_b") title_tmp="w_delta_phi_ee_b";
    if (file) {
      file->cd("demoEleGen");
    } else {
      file = TFile::Open((path + "/electrons/" + version + "/" + subdir +"/unfolding/" + title_tmp + "_unfolding.root").c_str());
    }
  }
  if (ilepton==2) {
    if (title=="w_pt_Z") title_tmp="w_pt_Z_mm";
    if (title=="w_pt_Z_b") title_tmp="w_pt_Z_mm_b";
    if (title=="w_delta_phi") title_tmp="w_delta_phi_mm";
    if (title=="w_delta_phi_b") title_tmp="w_delta_phi_mm_b";
    if (file) {
      file->cd("demoMuoGen");
    } else {
      file = TFile::Open((path + "/muons/" + version + "/" + subdir +"/unfolding/" + title_tmp + "_unfolding.root").c_str());
    }
  }
  hist = (TH1F*)gDirectory->Get(title_tmp.c_str())->Clone();
  hist->SetDirectory(0);
  if (!infile) file->Close();
  return hist;
}

void DataMCComp8(string& title="", int plot=0, int isratio=1) {

int useSherpa=0;
//int useSherpa=1; // use Sherpa MC prediction

double ele_eff_sys=0.015;
double muo_eff_sys=0.020;
double btag_sys=0.030;

string subdir="0";

	if (gROOT->GetVersionInt() >= 53401) {
	  gROOT->GetColor(kRed)->SetAlpha(0.5);
	  gROOT->GetColor(kGreen+2)->SetAlpha(0.5);
	  gROOT->GetColor(kMagenta-6)->SetAlpha(0.5);
	  if (!useSherpa) gROOT->GetColor(kMagenta-6)->SetAlpha(0.0);
	  gROOT->GetColor(kBlue-4)->SetAlpha(0.5);
	  gROOT->GetColor(kOrange+7)->SetAlpha(0.5);
	}

	double Lumi2012;

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
	TFile *mcg3 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL2_gen.root").c_str());

	string title_b = title;

	if (title.find("_bjet_")!=string::npos) {
	  title.erase(title.find("_bjet_")+1, 1);
	} else {
	  title_b = title + "_b";
	}

	TH1F* w_data[2];
	TH1F* w_data_b[2];
	for (int i=0; i<2; i++) {
	  w_data[i] = read(path, subdir, title, i+1);
	  w_data_b[i] = read(path, subdir, title_b, i+1);
	}

	TH1F* h_data = w_data[0]->Clone();
	h_data->Add(w_data[1]);
	h_data->Scale(0.5);
	TH1F* h_data_b = w_data_b[0]->Clone();
	h_data_b->Add(w_data_b[1]);
	h_data_b->Scale(0.5);
	h_data->SetStats(0);
	h_data_b->SetStats(0);

	TH1F* w_mcg[2];
	TH1F* w_mcg_b[2];
	TH1F* w_mcg1[2];
	TH1F* w_mcg1_b[2];
	TH1F* w_mcg2[2];
	TH1F* w_mcg2_b[2];
	TH1F* w_mcg3[2];
	TH1F* w_mcg3_b[2];
	for (int i=0; i<2; i++) {
	  w_mcg[i] = read(path, subdir, title, i+1, mcg);
	  w_mcg_b[i] = read(path, subdir, title_b, i+1, mcg);
	  w_mcg1[i] = read(path, subdir, title, i+1, mcg1);
	  w_mcg1_b[i] = read(path, subdir, title_b, i+1, mcg1);
	  w_mcg2[i] = read(path, subdir, title, i+1, mcg2[i]);
	  w_mcg2_b[i] = read(path, subdir, title_b, i+1, mcg2[i]);
	  w_mcg3[i] = read(path, subdir, title, i+1, mcg3);
	  w_mcg3_b[i] = read(path, subdir, title_b, i+1, mcg3);
	}

	TH1F* h_mcg = w_mcg[0]->Clone();
	h_mcg->Add(w_mcg[1]);
	h_mcg->Scale(0.5);
	TH1F* h_mcg_b = w_mcg_b[0]->Clone();
	h_mcg_b->Add(w_mcg_b[1]);
	h_mcg_b->Scale(0.5);
	TH1F* h_mcg1 = w_mcg1[0]->Clone();
	h_mcg1->Add(w_mcg1[1]);
	h_mcg1->Scale(0.5);
	TH1F* h_mcg1_b = w_mcg1_b[0]->Clone();
	h_mcg1_b->Add(w_mcg1_b[1]);
	h_mcg1_b->Scale(0.5);
	TH1F* h_mcg2 = w_mcg2[0]->Clone();
	h_mcg2->Add(w_mcg2[1]);
	h_mcg2->Scale(0.5);
	TH1F* h_mcg2_b = w_mcg2_b[0]->Clone();
	h_mcg2_b->Add(w_mcg2_b[1]);
	h_mcg2_b->Scale(0.5);
	TH1F* h_mcg3 = w_mcg3[0]->Clone();
	h_mcg3->Add(w_mcg3[1]);
	h_mcg3->Scale(0.5);
	TH1F* h_mcg3_b = w_mcg3_b[0]->Clone();
	h_mcg3_b->Add(w_mcg3_b[1]);
	h_mcg3_b->Scale(0.5);

	h_mcg->Sumw2();
	h_mcg1->Sumw2();
	h_mcg2->Sumw2();
	h_mcg3->Sumw2();

	h_mcg_b->Sumw2();
	h_mcg1_b->Sumw2();
	h_mcg2_b->Sumw2();
	h_mcg3_b->Sumw2();

	h_mcg->Scale(norm1);
	h_mcg1->Scale(norm1_1);
	h_mcg2->Scale(norm1_2);
	h_mcg3->Scale(norm1_3);

	h_mcg_b->Scale(norm1);
	h_mcg1_b->Scale(norm1_1);
	h_mcg2_b->Scale(norm1_2);
	h_mcg3_b->Scale(norm1_3);

	h_data->Scale(1./Lumi2012, "width");
	h_data_b->Scale(1./Lumi2012, "width");
	if (isratio==1) {
	  h_data_b->Divide(h_data);
	  h_data_b->Scale(100.);
	}

	for (int i=0; i<2; i++) {
	  ifstream in;
	  string title_b_tmp = title_b;
	  if (i==0) {
	    if (title_b=="w_pt_Z") title_b_tmp="w_pt_Z_ee";
	    if (title_b=="w_pt_Z_b") title_b_tmp="w_pt_Z_ee_b";
	    if (title_b=="w_delta_phi") title_b_tmp="w_delta_phi_ee";
	    if (title_b=="w_delta_phi_b") title_b_tmp="w_delta_phi_ee_b";
	    in.open((path + "/electrons/" + version + "/" + "/xsecs_unfolding/" + title_b_tmp + "_xsecs_unfolding.dat").c_str());
//cout << path + "/electrons/" + version + "/" + "/xsecs_unfolding/" + title_b_tmp + "_xsecs_unfolding.dat" << endl;
	  }
	  if (i==1) {
	    if (title_b=="w_pt_Z") title_b_tmp="w_pt_Z_mm";
	    if (title_b=="w_pt_Z_b") title_b_tmp="w_pt_Z_mm_b";
	    if (title_b=="w_delta_phi") title_b_tmp="w_delta_phi_mm";
	    if (title_b=="w_delta_phi_b") title_b_tmp="w_delta_phi_mm_b";
	    in.open((path + "/muons/" + version + "/" + "/xsecs_unfolding/" + title_b_tmp + "_xsecs_unfolding.dat").c_str());
//cout << path + "/muons/" + version + "/" + "/xsecs_unfolding/" + title_b_tmp + "_xsecs_unfolding.dat" << endl;
	  }

	  string tmp;
	  const int N=16;
	  float x[2][100][100];
	  float x_b[2][100][100];

	  getline(in, tmp);
	  getline(in, tmp);
	  for (int j=0; j<h_data->GetNbinsX()+2; j++) {
	    for (int k=0; k<N; k++){
	      in >> tmp;
	      in >> x[i][k][j];
//cout << x[i][k][j] << " ";
	    }
	    in.ignore();
//cout << endl;
	  }
	  getline(in, tmp);
	  getline(in, tmp);
	  for (int j=0; j<h_data_b->GetNbinsX()+2; j++) {
	    for (int k=0; k<N; k++){
	      in >> tmp;
	      in >> x_b[i][k][j];
//cout << x_b[i][k][j] << " ";
	    }
	    in.ignore();
//cout << endl;
	  }
	  in.close();
	}

	TH1F* h_data_stat = h_data->Clone();
	TH1F* h_data_b_stat = h_data_b->Clone();
	TH1F* h_data_syst = h_data->Clone();
	TH1F* h_data_b_syst = h_data_b->Clone();
	TH1F* h_data_tot = h_data->Clone();
	TH1F* h_data_b_tot = h_data_b->Clone();

	for (int i=0;i<=h_data_stat->GetNbinsX()+1;i++) {
	  float val = 0.0;
	  if (x[0][1][i]*x[1][1][i] != 0) {
	    val = (x[0][0][i]/(x[0][N-4][i]**2)+x[1][0][i]/(x[1][N-4][i]**2))/(1./(x[0][N-4][i]**2)+1./(x[1][N-4][i]**2));
	    h_data_stat->SetBinContent(i, val);
	    val = sqrt(1./(1./(x[0][N-4][i]**2)+1./(x[1][N-4][i]**2)));
	    h_data_stat->SetBinError(i, val);
	  }
	  val = (x[0][N-3][i]+x[1][N-3][i])/2.;
	  h_data_syst->SetBinError(i, val);
	  val = TMath::Sqrt(h_data_stat->GetBinError(i)**2+h_data_syst->GetBinError(i)**2);
	  h_data_tot->SetBinError(i, val);
	}

	for (int i=0;i<=h_data_b_stat->GetNbinsX()+1;i++) {
	  float val = 0.0;
	  if (x_b[0][1][i]*x_b[1][1][i] != 0) {
	    val = (x_b[0][0][i]/(x_b[0][N-4][i]**2)+x_b[1][0][i]/(x_b[1][N-4][i]**2))/(1./(x_b[0][N-4][i]**2)+1./(x_b[1][N-4][i]**2));
	    h_data_b_stat->SetBinContent(i, val);
	    val = sqrt(1./(1./(x_b[0][N-4][i]**2)+1./(x_b[1][N-4][i]**2)));
	    h_data_b_stat->SetBinError(i, val);
	  }
	  val = (x_b[0][N-3][i]+x_b[1][N-3][i])/2.;
	  h_data_b_syst->SetBinError(i, val);
	  val = TMath::Sqrt(h_data_b_stat->GetBinError(i)**2+h_data_b_syst->GetBinError(i)**2);
	  h_data_b_tot->SetBinError(i, val);
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
	  h_mcg3_b->Divide(h_mcg1);
	  h_mcg_b->Scale(100.);
	  h_mcg1_b->Scale(100.);
	  h_mcg2_b->Scale(100.);
	  h_mcg3_b->Scale(100.);
	}

	h_data = fixrange(h_data);
	h_data_b = fixrange(h_data_b);
	h_data_stat = fixrange(h_data_stat);
	h_data_b_stat = fixrange(h_data_b_stat);
	h_data_tot = fixrange(h_data_tot);
	h_data_b_tot = fixrange(h_data_b_tot);
	h_mcg = fixrange(h_mcg);
	h_mcg_b = fixrange(h_mcg_b);
	h_mcg1 = fixrange(h_mcg1);
	h_mcg1_b = fixrange(h_mcg1_b);
	h_mcg2 = fixrange(h_mcg2);
	h_mcg3 = fixrange(h_mcg3);
	h_mcg2_b = fixrange(h_mcg2_b);
	h_mcg3_b = fixrange(h_mcg3_b);

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
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
	  h_data_b->GetYaxis()->SetTitle("#sigma_{Z+b-jets}/#sigma_{Z+jets} [%]");
	}
	h_data_b->GetYaxis()->SetTitleOffset(1.2);
	h_data_b->GetXaxis()->SetTitleOffset(1.3);
	h_data_b->SetMarkerColor(kBlack);
	h_data_b->SetLineColor(kBlack);
	h_data_b->SetMarkerStyle(24);
	h_data_b->SetMarkerSize(0.7);
	h_data_b->SetStats(0);
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
	  //h_data_b->Draw("EPX0SAME");
	}

	TLegend *leg = new TLegend(0.62, 0.580, 0.88, 0.88);
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	if (isratio==0) {
	  pad1->SetLogy();

	  h_mcg_b->SetMaximum(4*h_data->GetMaximum());
	  h_mcg_b->SetMinimum(TMath::Max(0.000002,0.25*h_data_b->GetBinContent(h_data_b->GetMinimumBin())));

	  h_mcg_b->Draw("E5");
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

	  h_mcg3_b->Draw("E5SAME");
	  TH1F* tmp2_3 = h_mcg3_b->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp2_3->GetMinimum()==0) tmp2_3->GetXaxis()->SetRangeUser(0, tmp2_3->GetBinCenter(tmp2_3->GetMinimumBin()-1));
	  }
	  tmp2_3->SetFillColor(0);
	  tmp2_3->DrawClone("HISTLSAME");

	  h_data_b_tot->Draw("E1PX0SAME");
	  h_data_b_stat->Draw("E1PX0SAME");

	  //h_data_b->Draw("EPX0SAME");

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

	  h_mcg3->SetLineColor(kOrange+7);
	  h_mcg3->SetLineWidth(2);
	  h_mcg3->SetMarkerColor(kOrange+7);
	  h_mcg3->SetFillColor(kOrange+7);

	  h_data_tot->SetMarkerColor(kRed);
	  h_data_tot->SetLineColor(kRed);
	  h_data_tot->SetMarkerStyle(20);
	  h_data_tot->SetMarkerSize (0.7);
	  h_data_stat->SetLineColor(kBlack);
	  h_data_stat->SetMarkerColor(kBlack);
	  h_data_stat->SetMarkerStyle(20);
	  h_data_stat->SetMarkerSize (0.7);
	  h_data_tot->Draw("E1PX0SAME");
	  h_data_stat->Draw("E1PX0SAME");

	  h_data->SetMarkerColor(kBlack);
	  h_data->SetLineColor(kBlack);
	  h_data->SetMarkerStyle(20);
	  h_data->SetMarkerSize (0.7);
	  //h_data->Draw("EPX0SAME");

	  leg->AddEntry(h_data,"Z(#rightarrow ll) DATA","p");
	  leg->AddEntry(h_data_b,"Z(#rightarrow ll)+b DATA","p");
	  leg->AddEntry(h_mcg,"Z(#rightarrow ll) MadGraph 5FS","l");
	  leg->AddEntry(h_mcg3,"Z(#rightarrow ll) MadGraph 4FS","l");
	  if (useSherpa) leg->AddEntry(h_mcg1,"Z(#rightarrow ll) Sherpa","l");
	  leg->AddEntry(h_mcg2,"Z(#rightarrow ll) Powheg","l");
	}

	if (isratio==1) {
	  leg->AddEntry(h_data_b,"Z(#rightarrow ll) DATA","p");
	  leg->AddEntry(h_mcg_b,"Z(#rightarrow ll) MadGraph 5FS","l");
	  leg->AddEntry(h_mcg3_b,"Z(#rightarrow ll) MadGraph 4FS","l");
	  if (useSherpa) leg->AddEntry(h_mcg1_b,"Z(#rightarrow ll) Sherpa","l");
	  leg->AddEntry(h_mcg2_b,"Z(#rightarrow ll) Powheg","l");
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

	TH1F *h_M = h_data_b->Clone();
	TH1F *h_M_tot = h_data_b_tot->Clone();
	TH1F *h_M_stat = h_data_b_stat->Clone();

	h_M->Divide(h_mcg_b);
	h_M_tot->Divide(h_mcg_b);
	h_M_stat->Divide(h_mcg_b);

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
	h_M_tot->SetMarkerStyle(24);
	h_M_tot->Draw("E1PX0SAME");
	h_M_stat->SetMarkerStyle(24);
	h_M_stat->Draw("E1PX0SAME");

	if (isratio==0) {
	  TH1F *h_M2= h_data->Clone();
	  TH1F *h_M2_tot= h_data_tot->Clone();
	  TH1F *h_M2_stat= h_data_stat->Clone();

	  h_M2->Divide(h_mcg);
	  h_M2_tot->Divide(h_mcg);
	  h_M2_stat->Divide(h_mcg);

	  TGraphErrors *g_M2 = new TGraphErrors(h_M2);
	  TGraphErrors *g_M2_tot = new TGraphErrors(h_M2_tot);
	  TGraphErrors *g_M2_stat = new TGraphErrors(h_M2_stat);

	  float dx = 0.1*(g_M2->GetXaxis()->GetXmax()-g_M2->GetXaxis()->GetXmin())/g_M2->GetN();
	  for (int i=0; i<g_M2->GetN(); i++) {
	    g_M2->SetPoint(i, g_M2->GetX()[i]-dx, g_M2->GetY()[i]);
	    g_M2->SetPointError(i, 0, g_M2->GetEY()[i]);
	    g_M2_stat->SetPoint(i, g_M2_stat->GetX()[i]-dx, g_M2_stat->GetY()[i]);
	    g_M2_stat->SetPointError(i, 0, g_M2_stat->GetEY()[i]);
	    g_M2_tot->SetPoint(i, g_M2_tot->GetX()[i]-dx, g_M2_tot->GetY()[i]);
	    g_M2_tot->SetPointError(i, 0, g_M2_tot->GetEY()[i]);
	  }

	  g_M2->SetMarkerStyle(20);
	  g_M2->Draw("EP0SAME");
	  g_M2_tot->SetMarkerStyle(20);
	  g_M2_tot->Draw("E1P0SAME");
	  g_M2_stat->SetMarkerStyle(20);
	  g_M2_stat->Draw("E1P0SAME");
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

	TH1F *h_S = h_data_b->Clone();
	TH1F *h_S_tot = h_data_b_tot->Clone();
	TH1F *h_S_stat = h_data_b_stat->Clone();

	h_S->Divide(h_mcg1_b);
	h_S_tot->Divide(h_mcg1_b);
	h_S_stat->Divide(h_mcg1_b);

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
	if (useSherpa) {
	  h_S->Draw("EPX0");
	} else {
	  for (int i=0;i<=h_S->GetNbinsX()+1;i++) {
	    h_S->SetBinContent(i, -0.5);
	  }
	  h_S->Draw("EPX0");
	}
	h_S_tot->SetMarkerStyle(24);
	if (useSherpa) h_S_tot->Draw("E1PX0SAME");
	h_S_stat->SetMarkerStyle(24);
	if (useSherpa) h_S_stat->Draw("E1PX0SAME");

	if (isratio==0) {
	  TH1F *h_S2= h_data->Clone();
	  TH1F *h_S2_tot= h_data_tot->Clone();
	  TH1F *h_S2_stat= h_data_stat->Clone();

	  h_S2->Divide(h_mcg1);
	  h_S2_tot->Divide(h_mcg1);
	  h_S2_stat->Divide(h_mcg1);

	  TGraphErrors *g_S2 = new TGraphErrors(h_S2);
	  TGraphErrors *g_S2_tot = new TGraphErrors(h_S2_tot);
	  TGraphErrors *g_S2_stat = new TGraphErrors(h_S2_stat);

	  float dx = 0.1*(g_S2->GetXaxis()->GetXmax()-g_S2->GetXaxis()->GetXmin())/g_S2->GetN();
	  for (int i=0; i<g_S2->GetN(); i++) {
	    g_S2->SetPoint(i, g_S2->GetX()[i]-dx, g_S2->GetY()[i]);
	    g_S2->SetPointError(i, 0, g_S2->GetEY()[i]);
	    g_S2_stat->SetPoint(i, g_S2_stat->GetX()[i]-dx, g_S2_stat->GetY()[i]);
	    g_S2_stat->SetPointError(i, 0, g_S2_stat->GetEY()[i]);
	    g_S2_tot->SetPoint(i, g_S2_tot->GetX()[i]-dx, g_S2_tot->GetY()[i]);
	    g_S2_tot->SetPointError(i, 0, g_S2_tot->GetEY()[i]);
	  }

	  g_S2->SetMarkerStyle(20);
	  if (useSherpa) g_S2->Draw("EP0SAME");
	  g_S2_tot->SetMarkerStyle(20);
	  if (useSherpa) g_S2_tot->Draw("E1P0SAME");
	  g_S2_stat->SetMarkerStyle(20);
	  if (useSherpa) g_S2_stat->Draw("E1P0SAME");
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
	  TLine *OLine3 = new TLine(h_S->GetXaxis()->GetXmin(),1.,h_S->GetXaxis()->GetXmax(),1.);
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

	TH1F *h_P = h_data_b->Clone();
	TH1F *h_P_tot = h_data_b_tot->Clone();
	TH1F *h_P_stat = h_data_b_stat->Clone();

	h_P->Divide(h_mcg2_b);
	h_P_tot->Divide(h_mcg2_b);
	h_P_stat->Divide(h_mcg2_b);

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
	h_P_tot->SetMarkerStyle(24);
	h_P_tot->Draw("E1PX0SAME");
	h_P_stat->SetMarkerStyle(24);
	h_P_stat->Draw("E1PX0SAME");

	if (isratio==0) {
	  TH1F *h_P2= h_data->Clone();
	  TH1F *h_P2_tot= h_data_tot->Clone();
	  TH1F *h_P2_stat= h_data_stat->Clone();

	  h_P2->Divide(h_mcg2);
	  h_P2_tot->Divide(h_mcg2);
	  h_P2_stat->Divide(h_mcg2);

	  TGraphErrors *g_P2 = new TGraphErrors(h_P2);
	  TGraphErrors *g_P2_tot = new TGraphErrors(h_P2_tot);
	  TGraphErrors *g_P2_stat = new TGraphErrors(h_P2_stat);

	  float dx = 0.1*(g_P2->GetXaxis()->GetXmax()-g_P2->GetXaxis()->GetXmin())/g_P2->GetN();
	  for (int i=0; i<g_P2->GetN(); i++) {
	    g_P2->SetPoint(i, g_P2->GetX()[i]-dx, g_P2->GetY()[i]);
	    g_P2->SetPointError(i, 0, g_P2->GetEY()[i]);
	    g_P2_stat->SetPoint(i, g_P2_stat->GetX()[i]-dx, g_P2_stat->GetY()[i]);
	    g_P2_stat->SetPointError(i, 0, g_P2_stat->GetEY()[i]);
	    g_P2_tot->SetPoint(i, g_P2_tot->GetX()[i]-dx, g_P2_tot->GetY()[i]);
	    g_P2_tot->SetPointError(i, 0, g_P2_tot->GetEY()[i]);
	  }

	  g_P2->SetMarkerStyle(20);
	  g_P2->Draw("EP0SAME");
	  g_P2_tot->SetMarkerStyle(20);
	  g_P2_tot->Draw("E1P0SAME");
	  g_P2_stat->SetMarkerStyle(20);
	  g_P2_stat->Draw("E1P0SAME");
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

	if (useSherpa) {
	  pad2->cd();
	} else {
	  pad3->cd();
	}

	TH1F *h_M3 = h_data_b->Clone();
	TH1F *h_M3_tot = h_data_b_tot->Clone();
	TH1F *h_M3_stat = h_data_b_stat->Clone();

	h_M3->Divide(h_mcg3_b);
	h_M3_tot->Divide(h_mcg3_b);
	h_M3_stat->Divide(h_mcg3_b);

	TGraphErrors *g_M3 = new TGraphErrors(h_M3);
	TGraphErrors *g_M3_tot = new TGraphErrors(h_M3_tot);
	TGraphErrors *g_M3_stat = new TGraphErrors(h_M3_stat);

	float dx = 0.1*(g_M3->GetXaxis()->GetXmax()-g_M3->GetXaxis()->GetXmin())/g_M3->GetN();
	for (int i=0; i<g_M3->GetN(); i++) {
	  g_M3->SetPoint(i, g_M3->GetX()[i]+dx, g_M3->GetY()[i]);
	  g_M3->SetPointError(i, 0, g_M3->GetEY()[i]);
	  g_M3_stat->SetPoint(i, g_M3_stat->GetX()[i]+dx, g_M3_stat->GetY()[i]);
	  g_M3_stat->SetPointError(i, 0, g_M3_stat->GetEY()[i]);
	  g_M3_tot->SetPoint(i, g_M3_tot->GetX()[i]+dx, g_M3_tot->GetY()[i]);
	  g_M3_tot->SetPointError(i, 0, g_M3_tot->GetEY()[i]);
	}

	if (useSherpa) {
	  g_M3->SetMarkerStyle(25);
	  g_M3_stat->SetMarkerStyle(25);
	  g_M3_tot->SetMarkerStyle(25);
	} else {
	  g_M3->SetMarkerStyle(24);
	  g_M3->SetMarkerColor(kBlack);
	  g_M3_stat->SetMarkerStyle(24);
	  g_M3_stat->SetMarkerColor(kBlack);
	  g_M3_tot->SetMarkerStyle(24);
	  g_M3_tot->SetMarkerColor(kBlack);
	}
	g_M3->Draw("EP0SAME");
	if (useSherpa) {
	  g_M3_tot->SetMarkerStyle(25);
	} else {
	  g_M3_tot->SetMarkerStyle(24);
	}
	g_M3_tot->Draw("EP0SAME");

	TLine *OLine5 = new TLine(h_P->GetXaxis()->GetXmin(),0.93,h_P->GetXaxis()->GetXmax(),0.93);
	OLine5->SetLineColor(kOrange+7);
	OLine5->SetLineWidth(2);
	OLine5->Draw();

	c1->cd();

	if (isratio==1) {
	  h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 10);
	}

	if (title_b=="w_first_jet_pt_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp_{T} [pb]");
	  h_P->GetXaxis()->SetTitle("leading jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 20);
	  }
	} else if (title_b=="w_first_jet_eta_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#eta [pb]");
	  h_P->GetXaxis()->SetTitle("leading jet #eta");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_first_bjet_pt") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp^{b}_{T} [pb]");
	  h_P->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{b}_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 10);
	  }
	} else if (title_b=="w_first_bjet_eta") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#eta^{b} [pb]");
	  h_P->GetXaxis()->SetTitle("leading b-jet #eta");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta^{b} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_pt_Z_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp^{Z}_{T} [pb]");
	  h_P->GetXaxis()->SetTitle("Z boson p_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{Z}_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 20);
	  }
	} else if (title_b=="w_Ht_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dH_{T} [pb]");
	  h_P->GetXaxis()->SetTitle("H_{T} [GeV/c]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dH_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_delta_phi_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#Delta#phi_{bZ} [pb]");
	  h_P->GetXaxis()->SetTitle("#Delta#phi(bZ) [rad]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta#phi_{Zb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	}

	if (plot) {
	  ofstream out;
	  if (isratio==0) {
	    gSystem->mkdir((path + "/combined/" + version + "/xsecs_unfolding/").c_str(), kTRUE);
	    c1->SaveAs((path + "/combined/" + version + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.pdf").c_str());
	    out.open((path + "/combined/" + version + "/" + "/xsecs_unfolding/" + title_b + "_xsecs_unfolding.dat").c_str());
	  }
	  if (isratio==1) {
	    gSystem->mkdir((path + "/combined/" + version + "/ratios_unfolding/").c_str(), kTRUE);
	    c1->SaveAs((path + "/combined/" + version + "/ratios_unfolding/" + title_b + "_ratio_unfolding.pdf").c_str());
	    out.open((path + "/combined/" + version + "/" + "/ratios_unfolding/" + title_b + "_ratio_unfolding.dat").c_str());
	  }
	  out << h_data->GetName();
	  out << endl;
	  out << std::setw(25) << "tot stat";
	  out << std::setw(12) << "tot sys";
	  out << std::setw(12) << "tot error";
	  out << std::setw(8) << "%";
	  out << endl;
	  for (int i=0;i<=h_data_stat->GetNbinsX()+1;i++) {
	    out << std::fixed << std::setw(2);
	    out << i << " ";
	    out << std::fixed << std::setprecision(6) << std::setw(10);
	    out << h_data_stat->GetBinContent(i);
	    out << " +- ";
	    out << std::fixed << std::setprecision(6) << std::setw(8);
	    out << h_data_stat->GetBinError(i);
	    out << " +- ";
	    out << h_data_syst->GetBinError(i);
	    out << " => ";
	    out << h_data_tot->GetBinError(i);
	    out << " => ";
	    out << std::setprecision(1) << std::setw(4);
	    out << 100.*(h_data_stat->GetBinContent(i)==0 ? 0 : h_data_tot->GetBinError(i)/h_data_stat->GetBinContent(i));
	    out << endl;
	  }
	  out << h_data_b->GetName();
	  out << endl;
	  out << std::setw(25) << "tot stat";
	  out << std::setw(12) << "tot sys";
	  out << std::setw(12) << "tot error";
	  out << std::setw(8) << "%";
	  out << endl;
	  for (int i=0;i<=h_data_b_stat->GetNbinsX()+1;i++) {
	    out << std::fixed << std::setw(2);
	    out << i << " ";
	    out << std::fixed << std::setprecision(6) << std::setw(10);
	    out << h_data_b_stat->GetBinContent(i);
	    out << " +- ";
	    out << h_data_b_stat->GetBinError(i);
	    out << " +- ";
	    out << h_data_b_syst->GetBinError(i);
	    out << " => ";
	    out << h_data_b_tot->GetBinError(i);
	    out << " => ";
	    out << std::setprecision(1) << std::setw(4);
	    out << 100.*(h_data_b_stat->GetBinContent(i)==0 ? 0 : h_data_b_tot->GetBinError(i)/h_data_b_stat->GetBinContent(i));
	    out << endl;
	  }
	  out.close();
	}
}

