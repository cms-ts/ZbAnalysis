#include "LumiLabel.C"
#include "LumiInfo_v09.h"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";

void DataMCComp4(string& title="", int plot=0, int ilepton=1, int imode=0) {

// imode = -1; // identity test using pattuples
// imode =  0; // identity test using MadGraph
// imode =  1; // closure test using MadGraph + Sherpa
// imode =  2; // unfolding data

	gSystem->Load("libRooUnfold");

	if (ilepton<1 || ilepton>2) {
	  ilepton = 1 + ilepton % 2;
        }

	double Lumi2012;

	if (ilepton==1) Lumi2012 = Lumi2012_ele;
	if (ilepton==2) Lumi2012 = Lumi2012_muon;

	double norm1 = ( (Lumi2012 * Xsec_dy) / Ngen_dy);
	double norm1_1 = ( (Lumi2012 * Xsec_dy_1) / Ngen_dy_1);

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

	TFile *mc1;
	if (imode==-1) mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_patgen.root").c_str());
	if (imode>= 0) mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());

	TFile *mc2;
	if (imode==-1) mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_patgen.root").c_str());
	if (imode== 0) mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	if (imode== 1) mc2 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());
	if (imode== 2) mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());

	TH1F* h_data_reco;
	data->cd();
	h_data_reco = (TH1F*)gDirectory->Get((title+"_raw").c_str());

	TH1F* h_mc1_truth;
	TH1F* h_mc1_reco;
	TH2F* h_mc1_matrix;
	TH1F* h_mc2_truth;
	TH1F* h_mc2_reco;

	if (ilepton==1) {
	  mc1->cd("demoEleGen");
	  h_mc1_truth  = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd("demoEle");
	  h_mc1_reco   = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd("demoEleDump");
	  h_mc1_matrix = (TH2F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoEleGen");
	  h_mc2_truth  = (TH1F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoEle");
	  h_mc2_reco   = (TH1F*)gDirectory->Get(title.c_str());
	}
	if (ilepton==2) {
	  mc1->cd("demoMuoGen");
	  h_mc1_truth   = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd("demoMuo");
	  h_mc1_reco    = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd("demoMuoDump");
	  h_mc1_matrix  = (TH2F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoMuoGen");
	  h_mc2_truth   = (TH1F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoMuo");
	  h_mc2_reco    = (TH1F*)gDirectory->Get(title.c_str());
	}

	h_mc1_truth->Scale(norm1);
	h_mc1_reco->Scale(norm1);
	h_mc1_matrix->Scale(norm1);

	if (imode<=0) {
	  h_mc2_truth->Scale(norm1);
	  h_mc2_reco->Scale(norm1);
	}
	if (imode==1) {
	  h_mc2_truth->Scale(norm1_1);
	  h_mc2_reco->Scale(norm1_1);
	}
	if (imode==2) {
	  h_mc2_truth->Scale(norm1);
	  h_mc2_reco->Scale(norm1);
	}

        RooUnfoldResponse response (h_mc1_reco, h_mc1_truth, h_mc1_matrix);

	response.UseOverflow();

	int kreg = 7; // default 0 -> nbins/2
	int ntoys = 100; // default 1000

	RooUnfoldSvd unfold_mc (&response, h_mc2_reco, kreg, ntoys);
	RooUnfoldSvd unfold_data (&response, h_data_reco, kreg, ntoys);

	if (imode<=0) unfold_mc.PrintTable(cout, h_mc1_truth);

	TH1F* h_mc2_unf;

	TH1F* h_data_unf;

	TCanvas* c1 = new TCanvas("c", "c", 800, 600);

	c1->cd();
        TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
        pad1->SetBottomMargin(0.001);
        pad1->Draw();
        pad1->cd();
        pad1->SetLogy();

	if (imode<=1) {
	  h_mc2_unf = (TH1F*) unfold_mc.Hreco();

	  float val = TMath::Max(h_mc2_unf->GetMaximum(), h_mc2_reco->GetMaximum());
	  val = TMath::Max(val, h_mc2_truth->GetMaximum());
	  val = TMath::Max(val, h_mc1_reco->GetMaximum());
	  val = TMath::Max(val, h_mc1_truth->GetMaximum());

	  h_mc2_unf->SetMaximum(1.5*val);

	  h_mc2_unf->Draw();
	  h_mc2_reco->Draw("SAME");
	  h_mc2_truth->Draw("SAME");

	  h_mc2_unf->SetLineColor(kBlack);
	  h_mc2_reco->SetLineColor(kGreen);
	  h_mc2_truth->SetLineColor(kRed);

	  h_mc1_reco->Draw("SAME");
	  h_mc1_truth->Draw("SAME");

	  h_mc1_reco->SetLineColor(kGreen);
	  h_mc1_truth->SetLineColor(kRed);
	  h_mc1_reco->SetLineStyle(2);
	  h_mc1_truth->SetLineStyle(2);
	}

	if (imode==2) {
	  c1->cd();
          TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
          pad1->SetBottomMargin(0.001);
          pad1->Draw();
          pad1->cd();
          pad1->SetLogy();

	  h_data_unf = (TH1F*) unfold_data.Hreco();

	  float val = TMath::Max(h_data_unf->GetMaximum(),h_data_reco->GetMaximum());
	  val = TMath::Max(val, h_mc1_reco->GetMaximum());
	  val = TMath::Max(val, h_mc1_truth->GetMaximum());

	  h_data_unf->SetMaximum(1.5*val);

	  h_data_unf->Draw();
	  h_data_reco->Draw("SAME");

	  h_mc1_reco->Draw("SAME");
	  h_mc1_truth->Draw("SAME");

	  h_data_unf->SetLineColor(kBlack);
	  h_data_reco->SetLineColor(kGreen);

	  h_mc1_reco->SetLineColor(kGreen);
	  h_mc1_truth->SetLineColor(kRed);
	  h_mc1_reco->SetLineStyle(2);
	  h_mc1_truth->SetLineStyle(2);
	}

        pad1->Update();
        c1->Update();

        c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        pad2->Draw();
        pad2->cd();

        TH1F *h_ratio;
        if (imode<=1) h_ratio = (TH1F*) h_mc2_unf->Clone();
        if (imode==2) h_ratio = (TH1F*) h_data_unf->Clone();

        h_ratio->SetTitle("");
        h_ratio->SetStats(0);

        h_ratio->GetXaxis()->SetTitleOffset(0.9);
        h_ratio->GetXaxis()->SetTitleSize(0.1);
        h_ratio->GetXaxis()->SetLabelFont(42);
        h_ratio->GetXaxis()->SetLabelSize(0.08);
        h_ratio->GetXaxis()->SetTitleFont(42);
        h_ratio->GetYaxis()->SetTitle("ratio");
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.09);
        h_ratio->GetYaxis()->SetLabelSize(0.08);
        h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
        h_ratio->GetYaxis()->SetTitleOffset(0.4);
        h_ratio->Divide(h_mc2_truth);
        //h_ratio->SetMarkerStyle(20);
        h_ratio->Draw();

        TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
        OLine->SetLineColor(kRed);
        OLine->SetLineWidth(1);
        OLine->Draw();

        c1->cd();

	if (plot) {
	  if (ilepton==1) {
	    gSystem->mkdir((path + "/electrons/" + version + "/unfolding/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electrons/" + version + "/unfolding/" + title + "_unfolding.pdf").c_str());
	    TFile f((path + "/electrons/" + version + "/unfolding/" + title + "_unfolding.root").c_str(),"RECREATE");
	    if (imode<=1) h_mc2_unf->Write(title.c_str());
	    if (imode==2) h_data_unf->Write(title.c_str());
	    f.Close();
	  }
	  if (ilepton==2) {
	    gSystem->mkdir((path + "/muons/" + version + "/unfolding/").c_str(), kTRUE);
	    c1->SaveAs((path + "/muons/" + version + "/unfolding/" + title + "_unfolding.pdf").c_str());
	    TFile f((path + "/muons/" + version + "/unfolding/" + title + "_unfolding.root").c_str(),"RECREATE");
	    if (imode<=1) h_mc2_unf->Write(title.c_str());
	    if (imode==2) h_data_unf->Write(title.c_str());
	    f.Close();
	  }
	}

}
