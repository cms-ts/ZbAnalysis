#include "TLatex.h"

#include <TROOT.h>
#include "TStyle.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TLine.h"
#include "TObject.h"
#include <iostream>
#include <sstream>
#include "LumiLabel.C"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TTree.h"

#include <string.h>

void DataMCComp(string& s="", int plot=0){

bool isMu=true;
double Lumi2012A_muon = 792.472; //Jul13
double Lumi2012A_ele = 808.472;
double Lumi2012A;

int ilepton = 0; // if 1 --> electron; else muon

if(ilepton == 1)
	Lumi2012A = Lumi2012A_ele;
	else
	Lumi2012A = Lumi2012A_muon;

//////////////////////// DY

double Ngen_dy = 4577352;
double Xsec_dy = 3503.7; // NLO-corrected
double norm1 = (Lumi2012A / (Ngen_dy/ Xsec_dy)); // xsec corrected to NLO

///////////////////////

//////////////////////// tt

double Ngen_tt = 1360000;
double Xsec_tt = 225.197; //NLO
double norm2 = ( (Lumi2012A * Xsec_tt) / Ngen_tt);

///////////////////////

//////////////////////// ZZ

double Ngen_zz = 2686130;
double Xsec_zz = 8.059; // NLO with CTEQ
double norm3 = ( (Lumi2012A * Xsec_zz) / Ngen_zz);

///////////////////////

//////////////////////// WZ

double Ngen_wz = 2480000;
double Xsec_wz = 33.21;    // NLO with CTEQ
double norm4 = ( (Lumi2012A * Xsec_wz) / Ngen_wz);

///////////////////////

//////////////////////// QCD

double Ngen_qcd = 2000000;
double Xsec_qcd = 3.64E8; //search the NLO !
double norm5 = ( (Lumi2012A * Xsec_qcd) / Ngen_qcd);

///////////////////////

//////////////////////// WW

double Ngen_ww = 2480000;
double Xsec_ww = 54.838; // NLO with CTEQ
double norm6 = ( (Lumi2012A * Xsec_ww) / Ngen_ww); //search the NLO

///////////////////////

//////////////////////// WW

double Ngen_wj = 58117034;
double Xsec_wj = 31200; // NLO with CTEQ
double norm7 = ( (Lumi2012A * Xsec_wj) / Ngen_wj); //search the NLO

///////////////////////




	string title;
	if (s.empty()) string title = "w_jetmultiplicity";
	else title = s;

	//string title = "w_jet_pt";
	//string title = "w_ele_pt";
	//string title = "w_muon_pt";
	//string title = "w_mm_inv";
	//string title = "w_ee_inv";
	//string title = "w_secondvtx_N";
	//string title = "w_pu_weights";
	//string title = "recoVTX";
	//string title = "recoVTXw"

	if (title=="w_jetmultiplicity" || title=="w_jet_pt" || title=="recoVTX" || title=="recoVTXw" || title=="w_secondvtx_N" || title=="h_tracks" || title=="w_tracks" || title=="w_MET" || title == "w_bjetmultiplicity" || title == "w_bleading_pt") {
		norm1*=0.5;
		norm2*=0.5;
		norm3*=0.5;
		norm4*=0.5;
		norm5*=0.5;
		norm6*=0.5;
		norm7*=0.5;
	}

	if (ilepton == 1)
	  TFile *data = TFile::Open("/home/candelis/gpfs/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/test/v00/DoubleElectron_2012A_13Jul12.root"); //data file
	else
	  TFile *data = TFile::Open("/home/candelis/gpfs/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/test/v00/DoubleMu_2012A_13Jul12.root"); //data file

	TFile *mc1 = TFile::Open("/home/candelis/gpfs/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/test/v00/DYJetsToLL.root");
	TFile *mc2 = TFile::Open("/home/candelis/gpfs/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/test/v00/TTbar.root");
	TFile *mc3 = TFile::Open("/home/candelis/gpfs/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/test/v00/ZZ.root");
	TFile *mc4 = TFile::Open("/home/candelis/gpfs/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/test/v00/WZ.root");
	TFile *mc5 = TFile::Open("/home/candelis/gpfs/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/test/v00/QCD.root");
	TFile *mc6 = TFile::Open("/home/candelis/gpfs/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/test/v00/WW.root");
	TFile *mc7 = TFile::Open("/home/candelis/gpfs/CMSSW_5_3_7_patch4/src/ZbAnalysis/ZbSkim/test/v00/Wj.root");

	data->cd("demo");
	TH1F* h_data = (TH1F*)gDirectory->Get(title.c_str());

        mc1->cd("demo");
	TH1F* h_mc1 = (TH1F*)gDirectory->Get(title.c_str());
	TH1F* h_csv = (TH1F*)gDirectory->Get("bquarks");
         
	mc2->cd("demo");
	TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());

        mc3->cd("demo");
	TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());

        mc4->cd("demo");
	TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());

        mc5->cd("demo");
	TH1F* h_mc5 = (TH1F*)gDirectory->Get(title.c_str());
        
	mc6->cd("demo");
	TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());
	
	mc7->cd("demo");
	TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());

	THStack *hs = new THStack("hs","");

	h_mc1 -> SetLineColor(kBlack);
	h_mc1 -> SetFillColor(kYellow-4);
        //h_mc1 -> SetFillStyle(3004);
        
	h_mc4 -> SetLineColor(kGray+3);
	h_mc4 -> SetFillColor(kGray+3);
        //h_mc4 -> SetFillStyle(3004);

	
        h_mc2 -> SetLineColor(kBlack);
	h_mc2 -> SetFillColor(kBlue);
        //h_mc2 -> SetFillStyle(3004);

	h_mc3 -> SetLineColor(kBlack);
	h_mc3 -> SetFillColor(kGray+2);
        //h_mc3 -> SetFillStyle(3004);
	
        h_mc6 -> SetFillColor(kRed+2);
        h_mc6 -> SetLineColor(kGray+3);
        //h_mc6 -> SetFillStyle(3004);
        
	h_mc7 -> SetFillColor(kGray);
        h_mc7 -> SetLineColor(kGray+3);
        //h_mc6 -> SetFillStyle(3004);

	h_mc1->Scale(norm1); 	
	h_mc2->Scale(norm2);
	h_mc3->Scale(norm3);
	h_mc4->Scale(norm4);
	h_mc5->Scale(norm5);
	h_mc6->Scale(norm6);
        h_mc7->Scale(norm7);
	h_csv->Scale(norm1);	

	hs->Add(h_mc7);
	hs->Add(h_mc6);
	hs->Add(h_mc5);
	hs->Add(h_mc4);
	hs->Add(h_mc3);
	hs->Add(h_mc2);
	hs->Add(h_mc1);

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
	c1->cd();

	TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
	pad1->SetBottomMargin(0.001);
	pad1->Draw();
	pad1->cd();
        pad1->SetLogy();
	
	TH1F *h3=h_data->DrawCopy();
	h3->SetMinimum(-100);
	
	hs->Draw();
	hs->GetYaxis()->SetTitle("Events");
 	hs->GetXaxis()->SetLabelSize(0.08);
        hs->GetXaxis()->SetTitleOffset (0.7);
	hs->SetMinimum(8);
	
	h_data->Draw("EPSAMES");
	h_data->SetMarkerColor(kBlack);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize (1.0);
	//h_data->SetMaximum(1.2);

	if(title== "w_secondvtx_N"){		
		h_data->Draw("EP");
		h_mc1->Draw("same");
		h_mc1->SetFillColor(kYellow-4);
		h_csv->SetLineColor(kMagenta);
		h_csv->SetFillColor(kMagenta);
		h_csv->Draw("same");
		h_csv->SetFillStyle(3004);
	}

	leg = new TLegend(0.6, 0.60, 0.82, 0.78);
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	if(title!= "w_secondvtx_N"){
	leg->AddEntry(h_data,"DATA","p");
	leg->AddEntry(h_mc1,"Z+jets","f");
	leg->AddEntry(h_mc2,"t#bar{t}","f");
	leg->AddEntry(h_mc3,"ZZ","f");
	leg->AddEntry(h_mc4,"WZ","f");
	leg->AddEntry(h_mc5,"QCD","f");
	leg->AddEntry(h_mc6,"WW","f");
	leg->AddEntry(h_mc7,"W+jets", "f");
	}
	if(title== "w_secondvtx_N"){
	leg->AddEntry(h_csv,"b quarks in Z+jets","f");
	leg->AddEntry(h_mc1, "Z+jets", "f");
	leg->AddEntry(h_data, "DATA", "f");
	}
	leg->Draw();

	pad1->Update();
	c1->Update();

	c1->cd();


	h_ratio = (TH1D*) h_data->Clone("h_ratio");

	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);		  	  
	pad2->Draw();
        pad2->cd();
	h_ratio->Sumw2();
      	h_ratio->SetStats(0);
	h_ratio->SetTitle("");

	if(title=="w_jetmultiplicity"){
	h_ratio->GetXaxis ()->SetTitle("Number of jets");
	}else if(title=="w_mm_inv"){
	h_ratio->GetXaxis ()->SetTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
	}else if(title=="w_ee_inv"){
	h_ratio->GetXaxis ()->SetTitle("e^{+}e^{-} invariant mass [GeV/c^{2}]");
	}else if(title=="w_muon_pt"){
	h_ratio->GetXaxis ()->SetTitle("muon p_{T} [GeV/c]");
	}else if(title=="w_ele_pt"){
	h_ratio->GetXaxis ()->SetTitle("electron p_{T} [GeV/c]");
	}else if(title=="w_jet_pt"){
	h_ratio->GetXaxis ()->SetTitle("jet p_{T} [GeV/c]");
	}else if(title=="w_secondvtx_N"){
	h_ratio->GetXaxis ()->SetTitle("CSV discriminator");
	}else if(title=="w_recoVTX"){
	h_ratio->GetXaxis ()->SetTitle("Number of offline vertices");
	}else if(title=="w_muon_pt"){
	h_ratio->GetXaxis ()->SetTitle("muon p_{T} [GeV/c]");
	}else if(title=="w_MET"){
	h_ratio->GetXaxis ()->SetTitle("MET [GeV/c]");
	}

	h_ratio->GetXaxis ()->SetTitleOffset (0.7);
 	h_ratio->GetXaxis ()->SetTitleSize (0.1);
	h_ratio->GetXaxis ()->SetLabelFont (42);
 	h_ratio->GetXaxis()->SetLabelSize(0.08);
	h_ratio->GetXaxis ()->SetTitleFont (42);	  
	h_ratio->GetYaxis()->SetTitle("Data/MC");
	//h_ratio->GetYaxis()->SetNdivisions(5);
	h_ratio->GetYaxis()->SetTitleSize(0.09);
	h_ratio->GetYaxis()->SetLabelSize(0.05);
	h_ratio->GetYaxis()->SetRangeUser(-0.2, 1.9);
	h_ratio->GetYaxis()->SetTitleOffset(0.3);
	h_ratio->Divide(h_mc1);
	h_ratio->SetMarkerStyle(20);
	h_ratio->Draw("ep");

	  TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
	  OLine->SetLineColor(kRed);
	  OLine->SetLineStyle(2);
	  OLine->Draw();

	c1->cd();

 	TLatex *latexLabel=CMSPrel(0.8,"",0.12,0.85); // make fancy label
	latexLabel->Draw("same");

	if (plot) c1->SaveAs(("plots/" + title + ".pdf").c_str());

}

