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

void JECscaling(){

	TFile *data      = TFile::Open("/home/candelis/gpfs/CMSSW_5_3_9/src/ZbAnalysis/ZbSkim/test/data/v05/DoubleMu_2012_merge.root");
	data->cd("demo_mm");
	TH1F  *jet1 = (TH1F*)gDirectory->Get("w_first_jet_pt");
	data->cd("demo_mm_up");
	TH1F  *jet2 = (TH1F*)gDirectory->Get("w_first_jet_pt");
	data->cd("demo_mm_down");
        TH1F  *jet3 = (TH1F*)gDirectory->Get("w_first_jet_pt");
	
	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
	c1->cd();
	c1->SetLogy();


	jet1->SetLineColor(kBlack);
	jet1->SetMarkerStyle(20);
  	jet1->SetMarkerSize (1.0);
	jet2->SetLineColor(kBlue);
	jet3->SetLineColor(kRed);
	
	jet1->Draw("EPX0");
	jet2->Draw("LSAME");
	jet3->Draw("LSAME");

	jet1->GetXaxis()->SetTitle("jet p_{T} [GeV]");
	jet1->GetYaxis()->SetTitle("Events");


	leg = new TLegend(0.68, 0.58, 0.88, 0.88);
	leg->SetBorderSize(0);
	leg->SetEntrySeparation(0.01);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->AddEntry(jet1,"jet p_{T}","p");
	leg->AddEntry(jet2,"jet p_{T} JEC up","l");
	leg->AddEntry(jet3,"jet p_{T} JEC down","l");
	leg->Draw();

	c1->Update();



}

