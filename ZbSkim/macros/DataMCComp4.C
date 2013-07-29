#include "LumiLabel.C"
#include "LumiInfo_v09.h"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";

void DataMCComp4(string& title="", int plot=0, int ilepton=1) {

	gSystem->Load("libRooUnfold");

	if (ilepton<1 || ilepton>2) {
	  ilepton = 1 + ilepton % 2;
        }

	if (title.empty()) title = "w_jetmultiplicity";

	if (ilepton==1) {
	  if (title.find("muon")!=string::npos) return;
	  if (title.find("mm")!=string::npos) return;
	}
	if (ilepton==2) {
	  if (title.find("ele")!=string::npos) return;
	  if (title.find("ee")!=string::npos) return;
	}

	string file = title;

	if (file.find("_b_")!=string::npos) {
	  file.erase(file.find("_b_"), 2);
	}

        if (file.find("_bjet_")!=string::npos) {
	  file.erase(file.find("_bjet_")+1, 1);
	}

	TFile *data;
	if (ilepton==1) data = TFile::Open((path + "/electrons/" + version + "/xsecs/" + file + "_xsecs.root").c_str());
	if (ilepton==2) data = TFile::Open((path + "/muons/" + version + "/xsecs/" + file + "_xsecs.root").c_str());

	TFile *mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	TFile *mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	//TFile *mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_sherpa_gen.root").c_str());

	if (ilepton==1) {
	  mc1->cd("demoEleGen");
	  TH1F* h_truth = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd("demoEle");
	  TH1F* h_reco  = (TH1F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoEleDump");
	  TH2F* h_mtx  = (TH2F*)gDirectory->Get(title.c_str());
	}
	if (ilepton==2) {
	  mc1->cd("demoMuoGen");
	  TH1F* h_truth = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd("demoMuo");
	  TH1F* h_reco  = (TH1F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoMuoDump");
	  TH2F* h_mtx  = (TH2F*)gDirectory->Get(title.c_str());
	}

        RooUnfoldResponse response (h_reco, h_truth, h_mtx);
	response.UseOverflow();
// n=100: numero di mc toys usati per la propagazione dell'incertezza della 
// matrice di covarianza dopo l'unfolding.
// default n=1000, cambia di molto la velocita' dell'algo
	RooUnfoldSvd unfold (&response, h_reco, 7, 100);
	TH1F* h_unf= (TH1F*) unfold.Hreco();
       
	double N = (h_truth->Integral())/(h_unf->Integral());
        h_unf->Scale(N);

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
	c1->cd();
	h_unf->Draw();
	h_truth->Draw("SAME");
	h_reco->Draw("SAME");
	h_unf->SetLineColor(kBlack);
	h_truth->SetLineColor(kRed);
	h_reco->SetLineColor(kGreen);

	TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
	c2->cd();
	h_mtx->Draw();

	if (plot) {
	  if (ilepton==1) {
	    gSystem->mkdir((path + "/electrons/" + version + "/unfolding/").c_str(), kTRUE);
	    c2->SaveAs((path + "/electrons/" + version + "/unfolding/" + title + "_unfolding.pdf").c_str());
	    TFile f((path + "/electrons/" + version + "/unfolding/" + title + "_unfolding.root").c_str(),"RECREATE");
	    h_unf->Write(title.c_str());
	    f.Close();
	  }
	  if (ilepton==2) {
	    gSystem->mkdir((path + "/muons/" + version + "/unfolding/").c_str(), kTRUE);
	    c2->SaveAs((path + "/muons/" + version + "/unfolding/" + title + "_unfolding.pdf").c_str());
	    TFile f((path + "/muons/" + version + "/unfolding/" + title + "_unfolding.root").c_str(),"RECREATE");
	    h_unf->Write(title.c_str());
	    f.Close();
	  }
	}

}
