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

	if (file.find("_b")==string::npos) {
          if (file.find("_jet_")!=string::npos) {
            file.insert(file.find("_jet_")+1, "b");
          } else {
            file = file + "_b";
          }
        }

	TFile *data;
	if (ilepton==1) data = TFile::Open((path + "/electrons/" + version + "/xsecs/" + file + "_xsecs.root").c_str());
	if (ilepton==2) data = TFile::Open((path + "/muons/" + version + "/xsecs/" + file + "_xsecs.root").c_str());

	data->cd();
	TH1F* h_data_reco = (TH1F*)gDirectory->Get(title.c_str());

	TFile *mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	TFile *mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	//TFile *mc2 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());

	if (ilepton==1) {
	  mc1->cd("demoEleGen");
	  TH1F* h_mc1_truth = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd("demoEle");
	  TH1F* h_mc1_reco  = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd("demoEleDump");
	  TH2F* h_mc1_mtx  = (TH2F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoEleGen");
	  TH1F* h_mc2_truth = (TH1F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoEle");
	  TH1F* h_mc2_reco  = (TH1F*)gDirectory->Get(title.c_str());
	}
	if (ilepton==2) {
	  mc1->cd("demoMuoGen");
	  TH1F* h_mc1_truth = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd("demoMuo");
	  TH1F* h_mc1_reco  = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd("demoMuoDump");
	  TH2F* h_mc1_mtx  = (TH2F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoMuo");
	  TH1F* h_mc2_reco  = (TH1F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoMuoGen");
	  TH1F* h_mc2_truth = (TH1F*)gDirectory->Get(title.c_str());
	}

        RooUnfoldResponse response (h_mc1_reco, h_mc1_truth, h_mc1_mtx);

	if (title!="w_pt_Z_ee_b"&&title!="w_pt_Z_mm_b") {
	  response.UseOverflow();
	}

// n=100: numero di mc toys usati per la propagazione dell'incertezza della 
// matrice di covarianza dopo l'unfolding.
// default n=1000, cambia di molto la velocita' dell'algo
	RooUnfoldSvd unfold (&response, h_mc2_reco, 7, 100);

	//RooUnfoldBayes unfold (&response, h_mc2_reco, 7);

	TH1F* h_mc2_unf= (TH1F*) unfold.Hreco();
       
	double N = (h_mc2_truth->Integral())/(h_mc2_unf->Integral());
        h_mc2_unf->Scale(N);

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);
	c1->cd();
	c1->SetLogy();
	h_mc2_unf->Draw();
	h_mc2_truth->Draw("SAME");
	h_mc2_reco->Draw("SAME");
	h_mc1_truth->Draw("SAME");
	h_mc1_reco->Draw("SAME");

	h_mc2_unf->SetLineColor(kBlack);
	h_mc2_truth->SetLineColor(kRed);
	h_mc2_reco->SetLineColor(kGreen);
	h_mc1_truth->SetLineColor(kRed);
	h_mc1_reco->SetLineColor(kGreen);
	h_mc1_truth->SetLineStyle(2);
	h_mc1_reco->SetLineStyle(2);

	TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
	c2->cd();
	h_mc1_mtx->Draw();

	if (plot) {
	  if (ilepton==1) {
	    gSystem->mkdir((path + "/electrons/" + version + "/unfolding/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electrons/" + version + "/unfolding/" + title + "_unfolding.pdf").c_str());
	    TFile f((path + "/electrons/" + version + "/unfolding/" + title + "_unfolding.root").c_str(),"RECREATE");
	    h_mc2_unf->Write(title.c_str());
	    f.Close();
	  }
	  if (ilepton==2) {
	    gSystem->mkdir((path + "/muons/" + version + "/unfolding/").c_str(), kTRUE);
	    c1->SaveAs((path + "/muons/" + version + "/unfolding/" + title + "_unfolding.pdf").c_str());
	    TFile f((path + "/muons/" + version + "/unfolding/" + title + "_unfolding.root").c_str(),"RECREATE");
	    h_mc2_unf->Write(title.c_str());
	    f.Close();
	  }
	}

}
