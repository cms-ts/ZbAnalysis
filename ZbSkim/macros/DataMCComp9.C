#include "DataMCComp.h"
#include "LumiInfo_v14.h"

#include "CMS_lumi.C"
#include "CMS_process.C"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data";
//string path = "/gpfs/cms/users/lalicata/work/test/data/";

void DataMCComp9(string title="", int plot=0, int isratio=1, int numB=0) {

int useSherpa=0;
//int useSherpa=1; // use Sherpa MC prediction

int drawInclusive=0;
//int drawInclusive=1; // do plot the "inclusive" histogram

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

	string title_b = title;

	if (numB==0) {
	  if (title.find("_bjet_")!=string::npos) {
	    title.erase(title.find("_bjet_")+1, 1);
	  } else {
	    title_b = title + "_b";
	  }
	}

	TFile* file = 0;

	if (isratio==0) file = TFile::Open((path + "/combined/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.root").c_str());
	if (isratio==1) file = TFile::Open((path + "/combined/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.root").c_str());

	file->cd();

	TH1F* h_data_tot;
	TH1F* h_data_stat;
	TH1F* h_data_b_tot;
	TH1F* h_data_b_stat;

	string tmp;

	if (drawInclusive) {
	  tmp = title+"_data_tot";
	  h_data_tot = (TH1F*)gDirectory->Get(tmp.c_str())->Clone();
	  h_data_tot->SetDirectory(0);
	  tmp = title+"_data_stat";
	  h_data_stat = (TH1F*)gDirectory->Get(tmp.c_str())->Clone();
	  h_data_stat->SetDirectory(0);
	}
	tmp = title_b+"_data_tot";
	h_data_b_tot = (TH1F*)gDirectory->Get(tmp.c_str())->Clone();
	h_data_b_tot->SetDirectory(0);
	tmp = title_b+"_data_stat";
	h_data_b_stat = (TH1F*)gDirectory->Get(tmp.c_str())->Clone();
	h_data_b_stat->SetDirectory(0);

	TH1F* h_data_syst;

	if (drawInclusive) {
	  h_data_syst = (TH1F*)h_data_tot->Clone("h_data_syst");
	  for (int i=0;i<=h_data_tot->GetNbinsX()+1;i++) {
	    double xval = TMath::Sqrt(TMath::Power(h_data_tot->GetBinError(i),2)-TMath::Power(h_data_stat->GetBinError(i),2));
	    h_data_syst->SetBinError(i, xval);
	  }
	  h_data_syst->SetDirectory(0);
	}

	TH1F* h_data_b_syst = (TH1F*)h_data_b_tot->Clone("h_data_b_syst");
	for (int i=0;i<=h_data_b_tot->GetNbinsX()+1;i++) {
	  double xval = TMath::Sqrt(TMath::Power(h_data_b_tot->GetBinError(i),2)-TMath::Power(h_data_b_stat->GetBinError(i),2));
	  h_data_b_syst->SetBinError(i, xval);
	}
	h_data_b_syst->SetDirectory(0);

	TH1F* h_mcg;
	TH1F* h_mcg1;
	TH1F* h_mcg2;
	TH1F* h_mcg3;
	TH1F* h_mcg_b;
	TH1F* h_mcg1_b;
	TH1F* h_mcg2_b;
	TH1F* h_mcg3_b;

	if (drawInclusive) {
	  tmp = title+"_mcg";
	  h_mcg = (TH1F*)gDirectory->Get(tmp.c_str())->Clone();
	  h_mcg->SetDirectory(0);
	  tmp = title+"_mcg1";
	  h_mcg1 = (TH1F*)gDirectory->Get(tmp.c_str())->Clone();
	  h_mcg1->SetDirectory(0);
	  tmp = title+"_mcg2";
	  h_mcg2 = (TH1F*)gDirectory->Get(tmp.c_str())->Clone();
	  h_mcg2->SetDirectory(0);
	  tmp = title+"_mcg3";
	  h_mcg3 = (TH1F*)gDirectory->Get(tmp.c_str())->Clone();
	  h_mcg3->SetDirectory(0);
	}

	tmp = title_b+"_mcg";
	h_mcg_b = (TH1F*)gDirectory->Get(tmp.c_str())->Clone();
	h_mcg_b->SetDirectory(0);
	tmp = title_b+"_mcg1";
	h_mcg1_b = (TH1F*)gDirectory->Get(tmp.c_str())->Clone();
	h_mcg1_b->SetDirectory(0);
	tmp = title_b+"_mcg2";
	h_mcg2_b = (TH1F*)gDirectory->Get(tmp.c_str())->Clone();
	h_mcg2_b->SetDirectory(0);
	tmp = title_b+"_mcg3";
	h_mcg3_b = (TH1F*)gDirectory->Get(tmp.c_str())->Clone();
	h_mcg3_b->SetDirectory(0);

	file->Close();

	if (!drawInclusive) {
	  h_data_tot = (TH1F*)h_data_b_tot->Clone();
	  h_data_stat = (TH1F*)h_data_b_stat->Clone();
	  h_data_syst = (TH1F*)h_data_b_syst->Clone();
	  h_mcg = (TH1F*)h_mcg_b->Clone();
	  h_mcg1 = (TH1F*)h_mcg1_b->Clone();
	  h_mcg2 = (TH1F*)h_mcg2_b->Clone();
	  h_mcg3 = (TH1F*)h_mcg3_b->Clone();
	}

/*
	if (!drawInclusive) {
	  for (int i=0;i<=h_mcg->GetNbinsX()+1;i++) {
	    double val = 0.;
	    val = TMath::Power(h_mcg->GetBinError(i),2);
	    val = val + TMath::Power(h_mcg->GetBinContent(i)*eXsec_dy/Xsec_dy,2);
	    h_mcg->SetBinError(i, TMath::Sqrt(val));
	    val = TMath::Power(h_mcg_b->GetBinError(i),2);
	    val = val + TMath::Power(h_mcg_b->GetBinContent(i)*eXsec_dy/Xsec_dy,2);
	    h_mcg_b->SetBinError(i, TMath::Sqrt(val));
	    val = TMath::Power(h_mcg1->GetBinError(i),2);
	    val = val + TMath::Power(h_mcg1->GetBinContent(i)*eXsec_dy_1/Xsec_dy_1,2);
	    h_mcg1->SetBinError(i, TMath::Sqrt(val));
	    val = TMath::Power(h_mcg1_b->GetBinError(i),2);
	    val = val + TMath::Power(h_mcg1_b->GetBinContent(i)*eXsec_dy_1/Xsec_dy_1,2);
	    h_mcg1_b->SetBinError(i, TMath::Sqrt(val));
	    val = TMath::Power(h_mcg2->GetBinError(i),2);
	    val = val + TMath::Power(h_mcg2->GetBinContent(i)*eXsec_dy_2/Xsec_dy_2,2);
	    h_mcg2->SetBinError(i, TMath::Sqrt(val));
	    val = TMath::Power(h_mcg2_b->GetBinError(i),2);
	    val = val + TMath::Power(h_mcg2_b->GetBinContent(i)*eXsec_dy_2/Xsec_dy_2,2);
	    h_mcg2_b->SetBinError(i, TMath::Sqrt(val));
	    val = TMath::Power(h_mcg3->GetBinError(i),2);
	    val = val + TMath::Power(h_mcg3->GetBinContent(i)*eXsec_dy_3/Xsec_dy_3,2);
	    h_mcg3->SetBinError(i, TMath::Sqrt(val));
	    val = TMath::Power(h_mcg3_b->GetBinError(i),2);
	    val = val + TMath::Power(h_mcg3_b->GetBinContent(i)*eXsec_dy_3/Xsec_dy_3,2);
	    h_mcg3_b->SetBinError(i, TMath::Sqrt(val));
	  }
	}
*/

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
	  h_mcg_b->Draw("E0");
	  h_mcg_b->Draw("E2SAME");
	}

	h_mcg1_b->SetLineColor(kMagenta-6);
	h_mcg1_b->SetLineWidth(2);
	h_mcg1_b->SetFillColor(kMagenta-6);
	h_mcg1_b->SetMarkerColor(kMagenta-6);
	if (isratio==1) {
	  h_mcg1_b->Draw("E0SAME");
	  h_mcg1_b->Draw("E2SAME");
	}

	h_mcg2_b->SetLineColor(kBlue-4);
	h_mcg2_b->SetLineWidth(2);
	h_mcg2_b->SetFillColor(kBlue-4);
	h_mcg2_b->SetMarkerColor(kBlue-4);
	if (isratio==1) {
	  h_mcg2_b->Draw("E0SAME");
	  h_mcg2_b->Draw("E2SAME");
	}

	h_mcg3_b->SetLineColor(kOrange+7);
	h_mcg3_b->SetLineWidth(2);
	h_mcg3_b->SetFillColor(kOrange+7);
	h_mcg3_b->SetMarkerColor(kOrange+7);
	if (isratio==1) {
	  h_mcg3_b->Draw("E0SAME");
	  h_mcg3_b->Draw("E2SAME");
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
          h_data_b_tot->SetMarkerStyle(26);
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
          h_data_b_stat->SetMarkerStyle(26);
          h_data_b_stat->SetMarkerSize(0.9);
        }
	if (isratio==1) {
	  h_data_b_tot->Draw("E1PX0SAME");
	  h_data_b_stat->Draw("E1PX0SAME");
	}

        TLegend *leg = new TLegend(0.613, 0.590, 0.883, 0.880);
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
	  if (title=="w_DR_Zb_max") {
	    h_mcg_b->SetMinimum(0.00005);
	    h_mcg_b->SetMaximum(10.0);
	  }
	  if (title=="w_DR_Zb_min") {
	    h_mcg_b->SetMinimum(0.0001);
	    h_mcg_b->SetMaximum(10.0);
	  }
	  if (title=="w_DR_bb") {
	    h_mcg_b->SetMinimum(0.0005);
	    h_mcg_b->SetMaximum(2.0);
	  }
	  if (title=="w_Ht_b") {
	    h_mcg_b->SetMinimum(0.00002);
	    h_mcg_b->SetMaximum(0.01);
	    if (numB==1) {
	      h_mcg_b->SetMinimum(0.00005);
	      h_mcg_b->SetMaximum(0.4);
	    }
	  }

	  h_mcg_b->Draw("E0");
	  h_mcg1_b->Draw("E0SAME");
	  h_mcg2_b->Draw("E0SAME");
	  h_mcg3_b->Draw("E0SAME");

	  h_mcg_b->Draw("E2SAME");
	  h_mcg1_b->Draw("E2SAME");
	  h_mcg2_b->Draw("E2SAME");
	  h_mcg3_b->Draw("E2SAME");

	  h_data_b_tot->Draw("E1PX0SAME");
	  h_data_b_stat->Draw("E1PX0SAME");

	  h_mcg->SetLineColor(kGreen+2);
	  h_mcg->SetLineWidth(2);
	  h_mcg->SetMarkerColor(kGreen+2);
	  h_mcg->SetFillColor(kGreen+2);
	  if (drawInclusive) h_mcg->Draw("E2SAME");
	  TH1F* tmp4 = (TH1F*)h_mcg->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp4->GetMinimum()==0) tmp4->GetXaxis()->SetRangeUser(0, tmp4->GetBinCenter(tmp4->GetMinimumBin()-1));
	  }
	  tmp4->SetFillColor(0);
	  if (drawInclusive) tmp4->DrawClone("HISTSAME");

	  h_mcg1->SetLineColor(kMagenta-6);
	  h_mcg1->SetLineWidth(2);
	  h_mcg1->SetMarkerColor(kMagenta-6);
	  h_mcg1->SetFillColor(kMagenta-6);
	  if (drawInclusive) h_mcg1->Draw("E2SAME");
	  TH1F* tmp4_1 = (TH1F*)h_mcg1->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp4_1->GetMinimum()==0) tmp4_1->GetXaxis()->SetRangeUser(0, tmp4_1->GetBinCenter(tmp4_1->GetMinimumBin()-1));
	  }
	  tmp4_1->SetFillColor(0);
	  if (drawInclusive) tmp4_1->DrawClone("HISTSAME");

	  h_mcg2->SetLineColor(kBlue-4);
	  h_mcg2->SetLineWidth(2);
	  h_mcg2->SetMarkerColor(kBlue-4);
	  h_mcg2->SetFillColor(kBlue-4);
	  if (drawInclusive) h_mcg2->Draw("E2SAME");
	  TH1F* tmp4_2 = (TH1F*)h_mcg2->Clone();
	  if (title.find("_pt")!=string::npos || title.find("_Ht")!=string::npos) {
	    if (tmp4_2->GetMinimum()==0) tmp4_2->GetXaxis()->SetRangeUser(0, tmp4_2->GetBinCenter(tmp4_2->GetMinimumBin()-1));
	  }
	  tmp4_2->SetFillColor(0);
	  if (drawInclusive) tmp4_2->DrawClone("HISTSAME");

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
	}

	leg->AddEntry(h_data_b_stat,"DATA","lp");
	leg->AddEntry(h_mcg_b,"MadGraph 5FS + Pythia6","lf");
	leg->AddEntry(h_mcg3_b,"MadGraph 4FS + Pythia6","lf");
	if (useSherpa) leg->AddEntry(h_mcg1_b,"Sherpa","lf");
	leg->AddEntry(h_mcg2_b,"Powheg + Pythia6","lf");

	leg->Draw();

        writeExtraText = true;
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

        if (isratio==0) {
          if (title_b=="w_Ht_b" || title_b=="w_first_bjet_pt" || title_b=="w_second_bjet_pt" || title_b=="w_first_bjet_eta_abs" || title_b=="w_second_bjet_eta_abs" || title_b=="w_pt_Z_b" || title_b=="w_DR_bb" || title_b=="w_bb_mass" || title_b=="w_Zbb_mass"|| title_b=="w_DR_Zb_min"|| title_b=="w_DR_Zb_max"|| title_b=="w_A_Zb" ) {
            //CMS_process("Z/#gamma*#rightarrow ll selection", 0.135, 0.51);
            if (numB==0) CMS_process("Z/#gamma* + at least 1 b jet", 0.135, 0.51);
            if (numB==1) CMS_process("Z/#gamma* + exactly b jet", 0.135, 0.51);
            if (numB==2) CMS_process("Z/#gamma* + at least 2 b jets", 0.135, 0.51);
          }
          if (title_b=="w_delta_phi_b" || title_b=="w_delta_phi_2b" || title_b=="w_mass_Zj_b") {
            //CMS_process("Z/#gamma*#rightarrow ll selection", 0.68, 0.51);
            if (numB==0) CMS_process("Z/#gamma* + at least 1 b jet", 0.68, 0.51);
            if (numB==1) CMS_process("Z/#gamma* + exactly b jet", 0.68, 0.51);
            if (numB==2) CMS_process("Z/#gamma* + at least 2 b jets", 0.68, 0.51);
          }
          if (title_b=="w_first_bjet_eta" || title_b=="w_second_bjet_eta") {
            //CMS_process("Z/#gamma*#rightarrow ll selection", 0.68, 0.51);
            if (numB==0) CMS_process("Z/#gamma* + at least 1 b jet", 0.68, 0.51);
            if (numB==1) CMS_process("Z/#gamma* + exactly b jet", 0.68, 0.51);
            if (numB==2) CMS_process("Z/#gamma* + at least 2 b jets", 0.68, 0.51);
          }
        }
        if (isratio==1) {
          //CMS_process("Z/#gamma*#rightarrow ll selection", 0.135, 0.85);
          CMS_process("Z/#gamma* + at least 1 b jet", 0.135, 0.85);
        }

	TPad *pad2 = new TPad("pad2","pad2",0,0.29,1,0.4);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.001);
	pad2->Draw();
	pad2->cd();

	TH1F *h_M_tot = (TH1F*)h_mcg_b->Clone();
	TH1F *h_M_stat = (TH1F*)h_mcg_b->Clone();

	if (!drawInclusive) {
	  for (int i=0;i<=h_M_tot->GetNbinsX()+1;i++) {
	    h_M_tot->SetBinError(i,0.);
	    h_M_stat->SetBinError(i,0.);
	  }
	}

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
	h_M_tot->GetYaxis()->SetNdivisions(205);
	if (numB==2) h_M_tot->GetYaxis()->SetNdivisions(505);
	h_M_tot->GetYaxis()->SetTitleSize(0.17);
	h_M_tot->GetYaxis()->SetLabelSize(0.17);
	h_M_tot->GetYaxis()->SetRangeUser(0.65, 1.35);
	if (numB==2) h_M_tot->GetYaxis()->SetRangeUser(0.3, 1.7);
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
          h_M_tot->SetMarkerStyle(26);
          h_M_tot->SetMarkerSize(0.9);
        }
	h_M_tot->Draw("E1PX0");
	h_M_tot->Draw("E0PX0SAME");
        if (isratio==0) {
          h_M_stat->SetMarkerStyle(24);
          h_M_stat->SetMarkerSize(0.7);
        }
        if (isratio==1) {
          h_M_stat->SetMarkerStyle(26);
          h_M_stat->SetMarkerSize(0.9);
        }
	h_M_stat->Draw("E1PX0SAME");
	h_M_stat->Draw("E0PX0SAME");

	if (isratio==1) {
	  if (!drawInclusive) {
	    TH1F *h_M = (TH1F*)h_mcg_b->Clone();
	    for (int i=0;i<=h_M->GetNbinsX()+1;i++) {
	      h_M->SetBinError(i, h_M->GetBinContent(i)==0 ? 0 : h_M->GetBinError(i)/h_M->GetBinContent(i));
	      h_M->SetBinContent(i,h_M_tot->GetBinContent(i));
	      h_M_tot->SetBinContent(i,1.);
	      h_M_tot->SetBinError(i,h_data_b_tot->GetBinContent(i)==0 ? 0 : h_data_b_tot->GetBinError(i)/h_data_b_tot->GetBinContent(i));
	      h_M_stat->SetBinContent(i,1.);
	      h_M_stat->SetBinError(i,h_data_b_stat->GetBinContent(i)==0 ? 0 : h_data_b_stat->GetBinError(i)/h_data_b_stat->GetBinContent(i));
	    }
	    h_M->Draw("E2SAME");
	    h_M->Draw("E0SAME");
	  }
	}

	if (isratio==0) {
	  TH1F *h_M2_tot= (TH1F*)h_mcg->Clone();
	  TH1F *h_M2_stat= (TH1F*)h_mcg->Clone();

	  if (!drawInclusive) {
	    for (int i=0;i<=h_M2_tot->GetNbinsX()+1;i++) {
	      h_M2_tot->SetBinError(i,0.);
	      h_M2_stat->SetBinError(i,0.);
	    }
	  }

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

	  if (!drawInclusive) {
	    TH1F *h_M2 = (TH1F*)h_mcg_b->Clone();
	    for (int i=0;i<=h_M2->GetNbinsX()+1;i++) {
	      h_M2->SetBinError(i, h_M2->GetBinContent(i)==0 ? 0 : h_M2->GetBinError(i)/h_M2->GetBinContent(i));
	      h_M2->SetBinContent(i,h_M_tot->GetBinContent(i));
	      h_M_tot->SetBinContent(i,1.);
	      h_M_tot->SetBinError(i,h_data_b_tot->GetBinContent(i)==0 ? 0 : h_data_b_tot->GetBinError(i)/h_data_b_tot->GetBinContent(i));
	      h_M_stat->SetBinContent(i,1.);
	      h_M_stat->SetBinError(i,h_data_b_stat->GetBinContent(i)==0 ? 0 : h_data_b_stat->GetBinError(i)/h_data_b_stat->GetBinContent(i));
	    }
	    if (title=="w_DR_bb") {
	      h_M_tot->SetBinContent(1, -999.);
	      h_M_stat->SetBinContent(1, -999.);
	    }
	    h_M2->SetLineWidth(2);
	    h_M2->Draw("E2SAME");
	    h_M2->Draw("E0SAME");
	  }
	}

	TLatex *t2 = new TLatex();
	t2->SetTextSize(0.2);
	t2->SetTextFont(42);
	t2->SetLineWidth(2);
	t2->SetNDC();
	if (useSherpa) {
	  t2->DrawLatex(0.15,0.7,"MadGraph 5FS + Pythia6, normalized to #sigma_{NNLO} / MadGraph 4FS + Pythia6, normalized to #sigma_{NLO}");
	} else {
	  t2->DrawLatex(0.15,0.13,"MadGraph 5FS + Pythia6, normalized to  #sigma_{NNLO}");
	}

	TLine *OLine2 = new TLine(h_M_tot->GetXaxis()->GetXmin(),1.,h_M_tot->GetXaxis()->GetXmax(),1.);
	if (drawInclusive) {
	  OLine2->SetLineColor(kGreen+2);
	  OLine2->SetLineWidth(2);
	} else {
	  OLine2->SetLineColor(kBlack);
	  OLine2->SetLineWidth(0);
	}
	OLine2->Draw();

	c1->cd();

	TPad *pad3 = new TPad("pad3","pad3",0,0.18,1,0.29);
	pad3->SetTopMargin(0);
	pad3->SetBottomMargin(0.001);
	pad3->Draw();
	pad3->cd();

	TH1F *h_S_tot = (TH1F*)h_mcg1_b->Clone();
	TH1F *h_S_stat = (TH1F*)h_mcg1_b->Clone();

	if (!drawInclusive) {
	  for (int i=0;i<=h_S_tot->GetNbinsX()+1;i++) {
	    h_S_tot->SetBinError(i,0.);
	    h_S_stat->SetBinError(i,0.);
	  }
	}

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
	h_S_tot->GetYaxis()->SetNdivisions(205);
	if (numB==2) h_S_tot->GetYaxis()->SetNdivisions(505);
	h_S_tot->GetYaxis()->SetTitleSize(0.17);
	h_S_tot->GetYaxis()->SetLabelSize(0.17);
	h_S_tot->GetYaxis()->SetRangeUser(0.65, 1.35);
	if (numB==2) h_S_tot->GetYaxis()->SetRangeUser(0.3, 1.7);
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
          h_S_tot->SetMarkerStyle(26);
          h_S_tot->SetMarkerSize(0.9);
        }
	if (useSherpa) {
	  h_S_tot->Draw("E1PX0");
	  h_S_tot->Draw("E0PX0SAME");
	} else {
	  for (int i=0;i<=h_S_tot->GetNbinsX()+1;i++) {
	    h_S_tot->SetBinContent(i,-999.);
	    h_S_tot->SetBinError(i,0.);
	  }
	  h_S_tot->Draw("E1PX0");
	  h_S_tot->Draw("E0PX0SAME");
	}
        if (isratio==0) {
          h_S_stat->SetMarkerStyle(24);
          h_S_stat->SetMarkerSize(0.7);
        }
        if (isratio==1) {
          h_S_stat->SetMarkerStyle(26);
          h_S_stat->SetMarkerSize(0.9);
        }
	if (useSherpa) h_S_stat->Draw("E1PX0SAME");
	if (useSherpa) h_S_stat->Draw("E0PX0SAME");

	if (isratio==1) {
	  if (!drawInclusive) {
	    TH1F *h_S = (TH1F*)h_mcg1_b->Clone();
	    for (int i=0;i<=h_S->GetNbinsX()+1;i++) {
	      h_S->SetBinError(i, h_S->GetBinContent(i)==0 ? 0 : h_S->GetBinError(i)/h_S->GetBinContent(i));
	      h_S->SetBinContent(i,h_S_tot->GetBinContent(i));
	      h_S_tot->SetBinContent(i,1.);
	      h_S_tot->SetBinError(i,h_data_b_tot->GetBinContent(i)==0 ? 0 : h_data_b_tot->GetBinError(i)/h_data_b_tot->GetBinContent(i));
	      h_S_stat->SetBinContent(i,1.);
	      h_S_stat->SetBinError(i,h_data_b_stat->GetBinContent(i)==0 ? 0 : h_data_b_stat->GetBinError(i)/h_data_b_stat->GetBinContent(i));
	    }
	    if (useSherpa) h_S->Draw("E2SAME");
	    if (useSherpa) h_S->Draw("E0SAME");
	  }
	}

	if (isratio==0) {
	  TH1F *h_S2_tot= (TH1F*)h_mcg1->Clone();
	  TH1F *h_S2_stat= (TH1F*)h_mcg1->Clone();

	  if (!drawInclusive) {
	    for (int i=0;i<=h_S2_tot->GetNbinsX()+1;i++) {
	      h_S2_tot->SetBinError(i,0.);
	      h_S2_stat->SetBinError(i,0.);
	    }
	  }

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

	  if (!drawInclusive) {
	    TH1F *h_S2 = (TH1F*)h_mcg1_b->Clone();
	    for (int i=0;i<=h_S2->GetNbinsX()+1;i++) {
	      h_S2->SetBinError(i, h_S2->GetBinContent(i)==0 ? 0 : h_S2->GetBinError(i)/h_S2->GetBinContent(i));
	      h_S2->SetBinContent(i,h_S_tot->GetBinContent(i));
	      h_S_tot->SetBinContent(i,1.);
	      h_S_tot->SetBinError(i,h_data_b_tot->GetBinContent(i)==0 ? 0 : h_data_b_tot->GetBinError(i)/h_data_b_tot->GetBinContent(i));
	      h_S_stat->SetBinContent(i,1.);
	      h_S_stat->SetBinError(i,h_data_b_stat->GetBinContent(i)==0 ? 0 : h_data_b_stat->GetBinError(i)/h_data_b_stat->GetBinContent(i));
	    }
	    if (title=="w_DR_bb") {
	      h_S_tot->SetBinContent(1, -999.);
	      h_S_stat->SetBinContent(1, -999.);
	    }
	    h_S2->SetLineWidth(2);
	    h_S2->Draw("E2SAME");
	    h_S2->Draw("E0SAME");
	  }
	}

	TLatex *t3 = new TLatex();
	t3->SetTextSize(0.2);
	t3->SetTextFont(42);
	t3->SetLineWidth(2);
	t3->SetNDC();
	if (useSherpa) {
	  t3->DrawLatex(0.15,0.7,"Sherpa, normalized to #sigma_{NNLO}");
	} else {
	  t3->DrawLatex(0.15,0.13,"MadGraph 4FS + Pythia6, normalized to  #sigma_{NLO}");
	}

	if (useSherpa) {
	  TLine *OLine3 = new TLine(h_S_tot->GetXaxis()->GetXmin(),1.,h_S_tot->GetXaxis()->GetXmax(),1.);
	  if (drawInclusive) {
	    OLine3->SetLineColor(kMagenta-6);
	    OLine3->SetLineWidth(2);
	  } else {
	    OLine3->SetLineColor(kBlack);
	    OLine3->SetLineWidth(0);
	  }
	  OLine3->Draw();
	}

	c1->cd();

	TPad *pad4 = new TPad("pad4","pad4",0,0.0,1,0.18);
	pad4->SetTopMargin(0);
	pad4->SetBottomMargin(0.35);
	pad4->Draw();
	pad4->cd();

	TH1F *h_P_tot = (TH1F*)h_mcg2_b->Clone();
	TH1F *h_P_stat = (TH1F*)h_mcg2_b->Clone();

	if (!drawInclusive) {
	  for (int i=0;i<=h_P_tot->GetNbinsX()+1;i++) {
	    h_P_tot->SetBinError(i,0.);
	    h_P_stat->SetBinError(i,0.);
	  }
	}

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
	h_P_tot->GetYaxis()->SetNdivisions(205);
	if (numB==2) h_P_tot->GetYaxis()->SetNdivisions(505);
	h_P_tot->GetYaxis()->SetTitleSize(0.11);
	h_P_tot->GetYaxis()->SetLabelSize(0.11);
	h_P_tot->GetYaxis()->SetRangeUser(0.65, 1.35);
	if (numB==2) h_P_tot->GetYaxis()->SetRangeUser(0.3, 1.7);
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
          h_P_tot->SetMarkerStyle(26);
          h_P_tot->SetMarkerSize(0.9);
        }
	h_P_tot->Draw("E1PX0");
	h_P_tot->Draw("E0PX0SAME");
        if (isratio==0) {
          h_P_stat->SetMarkerStyle(24);
          h_P_stat->SetMarkerSize(0.7);
        }
        if (isratio==1) {
          h_P_stat->SetMarkerStyle(26);
          h_P_stat->SetMarkerSize(0.9);
        }
	h_P_stat->Draw("E1PX0SAME");
	h_P_stat->Draw("E0PX0SAME");

	if (isratio==1) {
	  if (!drawInclusive) {
	    TH1F *h_P = (TH1F*)h_mcg2_b->Clone();
	    for (int i=0;i<=h_P->GetNbinsX()+1;i++) {
	      h_P->SetBinError(i, h_P->GetBinContent(i)==0 ? 0 : h_P->GetBinError(i)/h_P->GetBinContent(i));
	      h_P->SetBinContent(i,h_P_tot->GetBinContent(i));
	      h_P_tot->SetBinContent(i,1.);
	      h_P_tot->SetBinError(i,h_data_b_tot->GetBinContent(i)==0 ? 0 : h_data_b_tot->GetBinError(i)/h_data_b_tot->GetBinContent(i));
	      h_P_stat->SetBinContent(i,1.);
	      h_P_stat->SetBinError(i,h_data_b_stat->GetBinContent(i)==0 ? 0 : h_data_b_stat->GetBinError(i)/h_data_b_stat->GetBinContent(i));
	    }
	    h_P->Draw("E2SAME");
	    h_P->Draw("E0SAME");
	  }
	}

	if (isratio==0) {
	  TH1F *h_P2_tot= (TH1F*)h_mcg2->Clone();
	  TH1F *h_P2_stat= (TH1F*)h_mcg2->Clone();

	  if (!drawInclusive) {
	    for (int i=0;i<=h_P2_tot->GetNbinsX()+1;i++) {
	      h_P2_tot->SetBinError(i,0.);
	      h_P2_stat->SetBinError(i,0.);
	    }
	  }

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

	  if (!drawInclusive) {
	    TH1F *h_P2 = (TH1F*)h_mcg2_b->Clone();
	    for (int i=0;i<=h_P2->GetNbinsX()+1;i++) {
	      h_P2->SetBinError(i, h_P2->GetBinContent(i)==0 ? 0 : h_P2->GetBinError(i)/h_P2->GetBinContent(i));
	      h_P2->SetBinContent(i,h_P_tot->GetBinContent(i));
	      h_P_tot->SetBinContent(i,1.);
	      h_P_tot->SetBinError(i,h_data_b_tot->GetBinContent(i)==0 ? 0 : h_data_b_tot->GetBinError(i)/h_data_b_tot->GetBinContent(i));
	      h_P_stat->SetBinContent(i,1.);
	      h_P_stat->SetBinError(i,h_data_b_stat->GetBinContent(i)==0 ? 0 : h_data_b_stat->GetBinError(i)/h_data_b_stat->GetBinContent(i));
	    }
	    if (title=="w_DR_bb") {
	      h_P_tot->SetBinContent(1, -999.);
	      h_P_stat->SetBinContent(1, -999.);
	    }
	    h_P2->SetLineWidth(2);
	    h_P2->Draw("E2SAME");
	    h_P2->Draw("E0SAME");
	  }
	}

	TLatex *t4 = new TLatex();
	t4->SetTextSize(0.13);
	t4->SetTextFont(42);
	t4->SetLineWidth(2);
	t4->SetNDC();
	t4->DrawLatex(0.15,0.43,"Powheg + Pythia6, normalized to #sigma_{NLO}");

	TLine *OLine4 = new TLine(h_P_tot->GetXaxis()->GetXmin(),1.,h_P_tot->GetXaxis()->GetXmax(),1.);
	if (drawInclusive) {
	  OLine4->SetLineColor(kBlue-4);
	  OLine4->SetLineWidth(2);
	} else {
	  OLine4->SetLineColor(kBlack);
	  OLine4->SetLineWidth(0);
	}
	OLine4->Draw();

	if (useSherpa) {
	  pad2->cd();
	} else {
	  pad3->cd();
	}

	if (isratio==1) {
	  if (!drawInclusive) {
	    TH1F *h_M = (TH1F*)h_mcg3_b->Clone();
	    for (int i=0;i<=h_M->GetNbinsX()+1;i++) {
	      h_M->SetBinError(i, h_M->GetBinContent(i)==0 ? 0 : h_M->GetBinError(i)/h_M->GetBinContent(i));
	      h_M->SetBinContent(i,h_M_tot->GetBinContent(i));
	      h_M_tot->SetBinContent(i,1.);
	      h_M_tot->SetBinError(i,h_data_b_tot->GetBinContent(i)==0 ? 0 : h_data_b_tot->GetBinError(i)/h_data_b_tot->GetBinContent(i));
	      h_M_stat->SetBinContent(i,1.);
	      h_M_stat->SetBinError(i,h_data_b_stat->GetBinContent(i)==0 ? 0 : h_data_b_stat->GetBinError(i)/h_data_b_stat->GetBinContent(i));
	    }
	  }
	}

	TH1F *h_M3_tot = (TH1F*)h_mcg3_b->Clone();
	TH1F *h_M3_stat = (TH1F*)h_mcg3_b->Clone();

	if (!drawInclusive) {
	  for (int i=0;i<=h_M3_tot->GetNbinsX()+1;i++) {
	    h_M3_tot->SetBinError(i,0.);
	    h_M3_stat->SetBinError(i,0.);
	  }
	}

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
	  g_M3_stat->SetMarkerSize(0.7);
	  g_M3_tot->SetMarkerStyle(25);
	  g_M3_tot->SetMarkerSize(0.7);
	} else {
          if (isratio==0) {
            g_M3_stat->SetMarkerStyle(24);
            g_M3_stat->SetMarkerSize(0.7);
          }
          if (isratio==1) {
            g_M3_stat->SetMarkerStyle(26);
            g_M3_stat->SetMarkerSize(0.9);
          }
	  g_M3_stat->SetMarkerColor(kBlack);
          if (isratio==0) {
            g_M3_tot->SetMarkerStyle(24);
            g_M3_tot->SetMarkerSize(0.7);
          }
          if (isratio==1) {
            g_M3_tot->SetMarkerStyle(26);
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
            g_M3_tot->SetMarkerStyle(26);
            g_M3_tot->SetMarkerSize(0.9);
          }
	}
	g_M3_tot->Draw("E1P");
	g_M3_tot->Draw("E0PSAME");
	g_M3_stat->Draw("E1PSAME");
	g_M3_stat->Draw("E0PSAME");

	if (!drawInclusive) {
	  TH1F *h_M3 = (TH1F*)h_mcg3_b->Clone();
	  for (int i=0;i<=h_M3->GetNbinsX()+1;i++) {
	    h_M3->SetBinError(i, h_M3->GetBinContent(i)==0 ? 0 : h_M3->GetBinError(i)/h_M3->GetBinContent(i));
	    h_M3->SetBinContent(i, h_M3_tot->GetBinContent(i));
	  }
	  for (int i=0; i<g_M3_tot->GetN(); i++) {
	    g_M3_tot->SetPoint(i,g_M3_tot->GetX()[i],1.);
	    g_M3_tot->SetPointError(i,0.,h_data_b_tot->GetBinContent(i+1)==0 ? 0 : h_data_b_tot->GetBinError(i+1)/h_data_b_tot->GetBinContent(i+1));
	    g_M3_stat->SetPoint(i,g_M3_stat->GetX()[i],1.);
	    g_M3_stat->SetPointError(i,0.,h_data_b_stat->GetBinContent(i+1)==0 ? 0 : h_data_b_stat->GetBinError(i+1)/h_data_b_stat->GetBinContent(i+1));
          }
	  if (title=="w_DR_bb") {
	    g_M3_tot->SetPoint(0, g_M3_tot->GetX()[0], -999.);
	    g_M3_stat->SetPoint(0, g_M3_stat->GetX()[0], -999.);
	  }
	  h_M3->SetLineWidth(2);
	  h_M3->Draw("E2SAME");
	  h_M3->Draw("E0SAME");
	}

	TLine *OLine5 = new TLine(h_P_tot->GetXaxis()->GetXmin(),1.,h_P_tot->GetXaxis()->GetXmax(),1.);
	if (drawInclusive) {
	  OLine5->SetLineColor(kOrange+7);
	  OLine5->SetLineWidth(2);
	} else {
	  OLine5->SetLineColor(kBlack);
	  OLine5->SetLineWidth(0);
	}
	OLine5->Draw();

	c1->cd();

	if (isratio==1) {
	  h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 10);
	}

        h_mcg_b->GetYaxis()->SetTitleSize(0.05);
        h_mcg_b->GetYaxis()->SetTitleOffset(0.6);
        h_P_tot->GetXaxis()->SetTitleSize(0.165);
        h_P_tot->GetXaxis()->SetTitleOffset(0.92);
        h_P_tot->GetXaxis()->SetTickLength(0.07);

        h_mcg_b->GetYaxis()->SetLabelSize(0.036);
        h_M_tot->GetYaxis()->SetLabelSize(0.19);
        h_S_tot->GetYaxis()->SetLabelSize(0.19);
        h_P_tot->GetYaxis()->SetLabelSize(0.12);
        h_P_tot->GetXaxis()->SetLabelOffset(0.02);

        h_M_tot->GetYaxis()->SetTitleSize(0.28);
        h_M_tot->GetYaxis()->SetTitleOffset(0.1);
        h_M_tot->GetYaxis()->SetTitle("/ Data   ");
        h_S_tot->GetYaxis()->SetTitleSize(0.28);
        h_S_tot->GetYaxis()->SetTitleOffset(0.1);
        h_S_tot->GetYaxis()->SetTitle("Theory");
        h_P_tot->GetYaxis()->SetTitle("");

	if (title_b=="w_first_jet_pt_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp_{T} [pb/GeV]");
	  h_P_tot->GetXaxis()->SetTitle("leading jet p_{T} [GeV]");
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
	  h_P_tot->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV]");
	  if (drawInclusive) h_P_tot->GetXaxis()->SetTitle("leading jet p_{T} [GeV]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{b}_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 10);
	  }
	} else if (title_b=="w_second_bjet_pt") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp^{b}_{T} [pb/GeV]");
	  h_P_tot->GetXaxis()->SetTitle("subleading b-jet p_{T} [GeV]");
	  if (drawInclusive) h_P_tot->GetXaxis()->SetTitle("subleading jet p_{T} [GeV]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{b}_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 10);
	  }
	} else if (title_b=="w_first_bjet_eta") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#eta^{b} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("leading b-jet #eta");
	  if (drawInclusive) h_P_tot->GetXaxis()->SetTitle("leading jet #eta");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta^{b} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_second_bjet_eta") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#eta^{b} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("subleading b-jet #eta");
	  if (drawInclusive) h_P_tot->GetXaxis()->SetTitle("subleading jet #eta");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#eta^{b} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 10);
	  }
	} else if (title_b=="w_first_bjet_eta_abs") {
          h_mcg_b->GetYaxis()->SetTitle("d#sigma / d|#eta^{b}| [pb]");
          h_P_tot->GetXaxis()->SetTitle("leading b-jet |#eta|");
          if (drawInclusive) h_P_tot->GetXaxis()->SetTitle("leading jet |#eta|");
          if (isratio==1) {
            h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d|#eta^{b}| [%]");
            h_mcg_b->GetYaxis()->SetRangeUser(0, 10);
          }
	} else if (title_b=="w_second_bjet_eta_abs") {
          h_mcg_b->GetYaxis()->SetTitle("d#sigma / d|#eta^{b}| [pb]");
          h_P_tot->GetXaxis()->SetTitle("subleading b-jet |#eta|");
          if (drawInclusive) h_P_tot->GetXaxis()->SetTitle("subleading jet |#eta|");
          if (isratio==1) {
            h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d|#eta^{b}| [%]");
            h_mcg_b->GetYaxis()->SetRangeUser(0, 10);
          }
        } else if (title_b=="w_pt_Z_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dp^{Z}_{T} [pb/GeV]");
	  h_P_tot->GetXaxis()->SetTitle("Z boson p_{T} [GeV]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dp^{Z}_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(-0.5, 20);
	  }
	} else if (title_b=="w_Ht_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / dH_{T} [pb/GeV]");
	  h_P_tot->GetXaxis()->SetTitle("H_{T} [GeV]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / dH_{T} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_delta_phi_b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#Delta#phi_{Zb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta#phi_{Zb} [rad]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta#phi_{Zb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_DR_Zb_min") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#Delta R^{min}_{Zb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta R^{min}_{Zb} [rad]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta R^{min}_{Zb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_DR_Zb_max") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#Delta R^{max}_{Zb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta R^{max}_{Zb} [rad]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta R^{max}_{Zb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_delta_phi_2b") {
	  h_mcg_b->GetYaxis()->SetTitle("d#sigma / d#Delta#phi_{bb} [pb]");
	  h_P_tot->GetXaxis()->SetTitle("#Delta#phi_{bb} [rad]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta#phi_{bb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_DR_bb") {
          h_mcg_b->GetYaxis()->SetTitle("d#sigma / d #Delta R_{bb} [pb]");
          h_P_tot->GetXaxis()->SetTitle("#Delta R_{bb} [rad]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d#Delta R_{bb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_Zbb_mass") {
          h_mcg_b->GetYaxis()->SetTitle("d#sigma / d M_{Zbb} [pb]");
          h_P_tot->GetXaxis()->SetTitle("M_{Zbb} [GeV]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d M_{Zbb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_bb_mass") {
          h_mcg_b->GetYaxis()->SetTitle("d#sigma / d M_{bb} [pb]");
          h_P_tot->GetXaxis()->SetTitle("M_{bb} [GeV]");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d M_{bb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
	} else if (title_b=="w_A_Zb") {
          h_mcg_b->GetYaxis()->SetTitle("d#sigma / d A_{Zbb} [pb]");
          h_P_tot->GetXaxis()->SetTitle("A_{Zbb}");
	  if (isratio==1) {
	    h_mcg_b->GetYaxis()->SetTitle("d[#sigma(Z+b) / #sigma(Z+j)] / d A_{Zbb} [%]");
	    h_mcg_b->GetYaxis()->SetRangeUser(0, 20);
	  }
        }

	if (plot) {
	  ofstream out, out1, out2;
	  if (isratio==0) {
	    gSystem->mkdir((path + "/paper/" + version + "/xsecs_unfolding" + dirbSel + "/").c_str(), kTRUE);
	    c1->SaveAs((path + "/paper/" + version + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.pdf").c_str());
	    out.open((path + "/paper/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.dat").c_str());
	    out1.open((path + "/paper/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.txt").c_str());
	    out2.open((path + "/paper/" + version + "/" + "/xsecs_unfolding" + dirbSel + "/" + title_b + "_xsecs_unfolding.tex").c_str());
	  }
	  if (isratio==1) {
	    gSystem->mkdir((path + "/paper/" + version + "/ratios_unfolding" + dirbSel + "/").c_str(), kTRUE);
	    c1->SaveAs((path + "/paper/" + version + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.pdf").c_str());
	    out.open((path + "/paper/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.dat").c_str());
	    out1.open((path + "/paper/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.txt").c_str());
	    out2.open((path + "/paper/" + version + "/" + "/ratios_unfolding" + dirbSel + "/" + title_b + "_ratio_unfolding.tex").c_str());
	  }
	  if (isratio==0 && drawInclusive) {
	    out << h_data_tot->GetName();
	    out << endl;
	    out << std::setw(25) << "total";
	    out << std::setw(12) << "total";
	    out << std::setw(12) << "total";
	    out << endl;
	    out << std::setw(25) << "stat";
	    out << std::setw(12) << "syst";
	    out << std::setw(12) << "error";
	    out << std::setw(8) << "%";
	    out << endl;
	    for (int i=0;i<=h_data_tot->GetNbinsX()+1;i++) {
	      out << std::fixed;
	      out << std::setw(2) << i;
	      out << " ";
	      out << std::setprecision(6);
	      out << std::setw(10) << h_data_tot->GetBinContent(i);
	      out << " +- ";
	      out << std::setw(8) << h_data_stat->GetBinError(i);
	      out << " +- ";
	      out << std::setw(8) << h_data_syst->GetBinError(i);
	      out << " => ";
	      out << std::setw(8) << h_data_tot->GetBinError(i);
	      out << " => ";
	      out << std::setprecision(1);
	      out << std::setw(4) << TMath::Abs(100.*(h_data_stat->GetBinContent(i)==0 ? 0 : h_data_tot->GetBinError(i)/h_data_stat->GetBinContent(i)));
	      out << endl;
	    }
	  }
	  out << h_data_b_tot->GetName();
	  out << endl;
	  out << std::setw(25) << "total";
	  out << std::setw(12) << "total";
	  out << std::setw(12) << "total";
	  out << endl;
	  out << std::setw(25) << "stat";
	  out << std::setw(12) << "syst";
	  out << std::setw(12) << "error";
	  out << std::setw(8) << "%";
	  out << endl;
	  for (int i=0;i<=h_data_b_tot->GetNbinsX()+1;i++) {
	    out << std::fixed;
	    out << std::setw(2) << i;
	    out << " ";
	    out << std::setprecision(6);
	    out << std::setw(10) << h_data_b_tot->GetBinContent(i);
	    out << " +- ";
	    out << std::setw(8) << h_data_b_stat->GetBinError(i);
	    out << " +- ";
	    out << std::setw(8) << h_data_b_syst->GetBinError(i);
	    out << " => ";
	    out << std::setw(8) << h_data_b_tot->GetBinError(i);
	    out << " => ";
	    out << std::setprecision(1);
	    out << std::setw(4) << TMath::Abs(100.*(h_data_b_stat->GetBinContent(i)==0 ? 0 : h_data_b_tot->GetBinError(i)/h_data_b_stat->GetBinContent(i)));
	    out << endl;
	  }
	  out.close();
	  if (isratio==0 && drawInclusive) {
	    out1 << h_data_tot->GetName() << " - RELATIVE ERRORS";
	    out1 << endl;
	    out1 << std::setw(7) << "total";
	    out1 << std::setw(8) << "total";
	    out1 << std::setw(8) << "total";
	    out1 << endl;
	    out1 << std::setw(7) << "stat";
	    out1 << std::setw(8) << "syst";
	    out1 << std::setw(8) << "error";
	    out1 << endl;
	    for (int i=0;i<=h_data_tot->GetNbinsX()+1;i++) {
	      double val = TMath::Abs(100.*(h_data_tot->GetBinContent(i)==0 ? 0 : 1./h_data_tot->GetBinContent(i)));
	      out1 << std::fixed;
	      out1 << std::setw(2) << i;
	      out1 << " ";
	      out1 << std::setprecision(1);
	      out1 << std::setw(4) << h_data_stat->GetBinError(i)*val;
	      out1 << " +- ";
	      out1 << std::setw(4) << h_data_syst->GetBinError(i)*val;
	      out1 << " => ";
	      out1 << std::setw(4) << h_data_tot->GetBinError(i)*val;
	      out1 << endl;
	    }
	  }
	  out1 << h_data_b_tot->GetName() << " - RELATIVE ERRORS";
	  out1 << endl;
	  out1 << std::setw(7) << "total";
	  out1 << std::setw(8) << "total";
	  out1 << std::setw(8) << "total";
	  out1 << endl;
	  out1 << std::setw(7) << "stat";
	  out1 << std::setw(8) << "syst";
	  out1 << std::setw(8) << "error";
	  out1 << endl;
	  for (int i=0;i<=h_data_b_tot->GetNbinsX()+1;i++) {
	    double val = TMath::Abs(100.*(h_data_b_tot->GetBinContent(i)==0 ? 0 : 1./h_data_b_tot->GetBinContent(i)));
	    out1 << std::fixed;
	    out1 << std::setw(2) << i;
	    out1 << " ";
	    out1 << std::setprecision(1);
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
	    out2 << std::setw(8) << "\\textbf{total} &" ;
	    out2 << std::setw(8) << "\\textbf{total} &" ;
	    out2 << std::setw(8) << "\\textbf{total} &" ;
	    out2 << endl;
	    out2 << std::setw(8) << "\\textbf{stat} & ";
	    out2 << std::setw(8) << "\\textbf{syst} & ";
	    out2 << std::setw(8) << "\\textbf{error & ";
	    out2 << endl;
	    /*for (int i=0;i<=h_data->GetNbinsX()+1;i++) {
	      double val = TMath::Abs(100.*(h_data->GetBinContent(i)==0 ? 0 : 1./h_data->GetBinContent(i)));
	      out2 << std::fixed;
	      out2 << std::setw(2) << i;
	      out2 << " ";
	      out2 << std::setprecision(1);
	      out2 << " & ";
	      out2 << std::setw(4) << h_data_stat->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << h_data_syst->GetBinError(i)*val;
	      out2 << " & ";
	      out2 << std::setw(4) << h_data_tot->GetBinError(i)*val;
	      out2 << std::setw(4) << "\\tabularnewline" << "   " << "\\hline";
	      out2 << endl;
	    }*/
	  }/*
	  //out2 << h_data_b_tot->GetName() << " - RELATIVE ERRORS";
	  //out2 << endl;
	  out2 << std::setw(8) << "\\textbf{total} &";
	  out2 << std::setw(8) << "\\textbf{total} &";
	  out2 << std::setw(8) << "\\textbf{total} &";
	  out2 << endl;
	  out2 << std::setw(8) << "\\textbf{stat} &";
	  out2 << std::setw(8) << "\\textbf{syst} &";
	  out2 << std::setw(8) << "\\textbf{error} &";
	  out2 << endl;*/
	  for (int i=1;i<=h_data_b_tot->GetNbinsX()+1;i++) {
	    double val = TMath::Abs(100.*(h_data_b_tot->GetBinContent(i)==0 ? 0 : 1./h_data_b_tot->GetBinContent(i)));
	    out2 << std::fixed;
	    out2 << std::setw(2) << i;
	    out2 << " & ";
	    out2 << std::setprecision(1);
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
