#include "LumiLabel.C"
#include "LumiInfo_v10.h"

#include "fixrange.C"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";

void DataMCComp4(string& title="", int plot=0, int ilepton=1, int imode=3, int method=0) {

// imode = -1; // identity test using MadGraph PAT
// imode =  0; // identity test using MadGraph GEN
// imode =  1; // closure test using MadGraph + Sherpa
// imode =  2; // closure test using MadGraph + Powheg
// imode =  3; // unfolding data with MadGraph
// imode =  4; // unfolding data with Sherpa
// imode =  5; // unfolding data with Powheg

// method = 0; // use SVD
// method = 1; // use Bayes
// method = 2; // use TUnfold
// method = 3; // use BinByBin

	gSystem->Load("libRooUnfold");

// use fit results for c_b

        double c_b=1.0;
	if (imode>=3) {
	  if (ilepton==1) c_b = 0.770;
	  if (ilepton==2) c_b = 0.747;
	}

	double Lumi2012;

	if (ilepton==1) Lumi2012 = Lumi2012_ele;
	if (ilepton==2) Lumi2012 = Lumi2012_muon;

	double norm1 = ( (Lumi2012 * Xsec_dy) / Ngen_dy);
	double norm1_1 = ( (Lumi2012 * Xsec_dy_1) / Ngen_dy_1);
	double norm1_2;
	if (ilepton==1) norm1_2 = ( (Lumi2012 * Xsec_dy_2) / Ngen_dy_2_ee);
	if (ilepton==2) norm1_2 = ( (Lumi2012 * Xsec_dy_2) / Ngen_dy_2_mm);

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
	if (imode== 0) mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	if (imode== 1) mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	if (imode== 2) mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	if (imode== 3) mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	if (imode== 4) mc1 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());
	if (imode== 5) {
	  if (ilepton==1) mc1 = TFile::Open((path + "/" + version + "/" + "DYToEE_powheg_gen.root").c_str());
	  if (ilepton==2) mc1 = TFile::Open((path + "/" + version + "/" + "DYToMuMu_powheg_gen.root").c_str());
	}

	TFile *mc2;
	if (imode==-1) mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_patgen.root").c_str());
	if (imode== 0) mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	if (imode== 1) mc2 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());
	if (imode== 2) {
	  if (ilepton==1) mc2 = TFile::Open((path + "/" + version + "/" + "DYToEE_powheg_gen.root").c_str());
	  if (ilepton==2) mc2 = TFile::Open((path + "/" + version + "/" + "DYToMuMu_powheg_gen.root").c_str());
	}
	if (imode== 3) mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	if (imode== 4) {
	  if (ilepton==1) mc2 = TFile::Open((path + "/" + version + "/" + "DYToEE_powheg_gen.root").c_str());
	  if (ilepton==2) mc2 = TFile::Open((path + "/" + version + "/" + "DYToMuMu_powheg_gen.root").c_str());
	}
	if (imode== 5) mc2 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());

	TH1F* h_data_reco;
	data->cd();
	h_data_reco = (TH1F*)gDirectory->Get((title+"_raw").c_str());

	string title_b = title;

	if (title.find("_b")!=string::npos) {
	  title_b = "b"+title.substr(1);
	}

	TH1F* h_mc1_truth;
	TH1F* h_mc1_reco;
	TH2F* h_mc1_matrix;
	TH1F* h_mc2_truth;
	TH1F* h_mc2_reco;

	if (ilepton==1) {
	  mc1->cd("demoEleGen");
	  h_mc1_truth  = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd("demoEle");
	  h_mc1_reco   = (TH1F*)gDirectory->Get(title_b.c_str());
	  mc1->cd("demoEleDump");
	  h_mc1_matrix = (TH2F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoEleGen");
	  h_mc2_truth  = (TH1F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoEle");
	  h_mc2_reco   = (TH1F*)gDirectory->Get(title_b.c_str());
	}
	if (ilepton==2) {
	  mc1->cd("demoMuoGen");
	  h_mc1_truth   = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd("demoMuo");
	  h_mc1_reco    = (TH1F*)gDirectory->Get(title_b.c_str());
	  mc1->cd("demoMuoDump");
	  h_mc1_matrix  = (TH2F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoMuoGen");
	  h_mc2_truth   = (TH1F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoMuo");
	  h_mc2_reco    = (TH1F*)gDirectory->Get(title_b.c_str());
	}

	h_data_reco->Sumw2();
	h_mc1_truth->Sumw2();
	h_mc1_reco->Sumw2();
	h_mc1_matrix->Sumw2();
	h_mc2_truth->Sumw2();
	h_mc2_reco->Sumw2();

	h_data_reco = fixrange(h_data_reco);
	h_mc1_truth = fixrange(h_mc1_truth);
	h_mc1_reco = fixrange(h_mc1_reco);
	h_mc1_matrix = fixrange(h_mc1_matrix);
	h_mc2_truth = fixrange(h_mc2_truth);
	h_mc2_reco = fixrange(h_mc2_reco);

	RooUnfoldResponse response (h_mc1_reco, h_mc1_truth, h_mc1_matrix);
	response.UseOverflow();

	h_mc1_truth->Scale(norm1);
	h_mc1_reco->Scale(norm1);
	if (imode==4) {
	  h_mc1_truth->Scale(norm1_1/norm1);
	  h_mc1_reco->Scale(norm1_1/norm1);
	}
	if (imode==5) {
	  h_mc1_truth->Scale(norm1_2/norm1);
	  h_mc1_reco->Scale(norm1_2/norm1);
	}
	h_mc2_truth->Scale(norm1);
	h_mc2_reco->Scale(norm1);
	if (imode==1) {
	  h_mc2_truth->Scale(norm1_1/norm1);
	  h_mc2_reco->Scale(norm1_1/norm1);
	}
	if (imode==2) {
	  h_mc2_truth->Scale(norm1_2/norm1);
	  h_mc2_reco->Scale(norm1_2/norm1);
	}

	if (title.find("_b")!=string::npos) {
	  h_mc1_truth->Scale(c_b);
	  h_mc1_reco->Scale(c_b);
	  h_mc2_truth->Scale(c_b);
	  h_mc2_reco->Scale(c_b);
	}

	if (method==0) {
	  int kreg = 0; // default 0 -> nbins/2
	  int ntoys = 100; // default 1000
	  RooUnfoldSvd unfold_mc (&response, h_mc2_reco, kreg, ntoys);
	  RooUnfoldSvd unfold_data (&response, h_data_reco, kreg, ntoys);
	}

	if (method==1) {
	  int niter = 4; // default 4 -> number of iterations
	  RooUnfoldBayes unfold_mc (&response, h_mc2_reco, niter);
	  RooUnfoldBayes unfold_data (&response, h_data_reco, niter);
	}

	if (method==2) {
	  RooUnfoldTUnfold unfold_mc (&response, h_mc2_reco);
	  RooUnfoldTUnfold unfold_data (&response, h_data_reco);
	}

	if (method==3) {
	  RooUnfoldBinByBin unfold_mc (&response, h_mc2_reco);
	  RooUnfoldBinByBin unfold_data (&response, h_data_reco);
	}

	if (imode<=0) unfold_mc.PrintTable(cout, h_mc1_truth);

	TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);

	c1->cd();
        TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
        pad1->SetBottomMargin(0.001);
        pad1->Draw();
        pad1->cd();
        pad1->SetLogy();

	TH1F* h_mc2_unf;
	if (imode<=2) {
	  h_mc2_unf = (TH1F*) unfold_mc.Hreco();

	  h_mc2_unf->SetMinimum(TMath::Max(1.0, 0.1*h_mc2_reco->GetMinimum()));

	  float val = TMath::Max(h_mc2_unf->GetMaximum(), h_mc2_reco->GetMaximum());
	  val = TMath::Max(val, h_mc2_truth->GetMaximum());
	  val = TMath::Max(val, h_mc1_reco->GetMaximum());
	  val = TMath::Max(val, h_mc1_truth->GetMaximum());

	  h_mc2_unf->SetMaximum(1.5*val);

	  h_mc2_unf->Draw("HIST");
	  h_mc2_reco->Draw("HISTSAME");
	  h_mc2_truth->Draw("HISTSAME");

	  h_mc2_unf->SetLineColor(kBlack);
	  h_mc2_reco->SetLineColor(kGreen);
	  h_mc2_truth->SetLineColor(kRed);

	  h_mc1_reco->Draw("HISTSAME");
	  h_mc1_truth->Draw("HISTSAME");

	  h_mc1_reco->SetLineColor(kGreen);
	  h_mc1_truth->SetLineColor(kRed);
	  h_mc1_reco->SetLineStyle(2);
	  h_mc1_truth->SetLineStyle(2);
	}

	TH1F* h_data_unf;
	if (imode>=3) {
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

	  h_data_reco->SetMaximum(1.5*val);

          h_data_reco->SetLineColor(kGreen);
          h_data_reco->SetMarkerColor(kGreen);
          h_data_reco->SetMarkerStyle(20);
          h_data_reco->SetMarkerSize(0.7);
          h_data_reco->Draw("EPX0");

	  h_data_unf->SetLineColor(kBlack);
	  h_data_unf->SetMarkerColor(kBlack);
	  h_data_unf->SetMarkerStyle(20);
	  h_data_unf->SetMarkerSize(0.7);
	  h_data_unf->Draw("EPX0SAME");

          h_mc1_truth->SetLineColor(kRed);
          h_mc1_truth->SetLineStyle(2);
          h_mc1_truth->Draw("HISTSAME");

	  h_mc1_reco->SetLineColor(kGreen);
	  h_mc1_reco->SetLineStyle(2);
	  h_mc1_reco->Draw("HISTSAME");
	}

	TH1F *tmp;
	if (imode<=2) tmp = h_mc2_reco;
	if (imode>=3) tmp = h_data_reco;
	tmp->SetTitle("");
	tmp->SetStats(0);
	tmp->GetYaxis()->SetTitle("Events");

	if (imode<=2) tmp = h_mc2_unf;
	if (imode>=3) tmp = h_data_unf;
	if (title=="w_first_jet_pt") {
	  tmp->GetXaxis()->SetTitle("leading jet p_{T} [GeV/c]");
	} else if (title=="w_first_jet_eta") {
	  tmp->GetXaxis()->SetTitle("leading jet #eta");
	} else if (title=="w_pt_Z_ee"||title=="w_pt_Z_mm") {
	  tmp->GetXaxis()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_delta_phi_ee"||title=="w_delta_phi_mm") {
	  tmp->GetXaxis()->SetTitle("#Delta #phi(jZ) [rad]");
	} else if (title=="w_Ht") {
          tmp->GetXaxis()->SetTitle("H_{T} [GeV/c]");
	}
	if (title=="w_first_bjet_pt") {
	  tmp->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV/c]");
	} else if (title=="w_first_bjet_eta") {
	  tmp->GetXaxis()->SetTitle("leading b-jet #eta");
	} else if (title=="w_pt_Z_ee_b"||title=="w_pt_Z_mm_b") {
	  tmp->GetXaxis()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_delta_phi_ee_b"||title=="w_delta_phi_mm_b") {
	  tmp->GetXaxis()->SetTitle("#Delta #phi(bZ) [rad]");
	} else if (title=="w_Ht_b") {
          tmp->GetXaxis()->SetTitle("H_{T} [GeV/c]");
	}

        TLegend *leg = new TLegend(0.42, 0.580, 0.68, 0.88);
        leg->SetBorderSize(0);
        leg->SetEntrySeparation(0.01);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);

        if (imode<=3) {
          leg->AddEntry(h_mc1_reco,"MADGRAPH reco","l");
          leg->AddEntry(h_mc1_truth,"MADGRAPH truth","l");
	}
        if (imode==4) {
          leg->AddEntry(h_mc1_reco,"SHERPA reco","l");
          leg->AddEntry(h_mc1_truth,"SHERPA truth","l");
	}
        if (imode==5) {
          leg->AddEntry(h_mc1_reco,"POWHEG reco","l");
          leg->AddEntry(h_mc1_truth,"POWHEG truth","l");
	}
        if (imode<=0) leg->AddEntry(h_mc2_unf,"MADGRAPH unfold","l");
        if (imode==1) {
          leg->AddEntry(h_mc2_reco,"SHERPA reco","l");
          leg->AddEntry(h_mc2_truth,"SHERPA truth","l");
          leg->AddEntry(h_mc2_unf,"SHERPA unfold","l");
        }
        if (imode==2) {
          leg->AddEntry(h_mc2_reco,"POWHEG reco","l");
          leg->AddEntry(h_mc2_truth,"POWHEG truth","l");
          leg->AddEntry(h_mc2_unf,"POWHEG unfold","l");
        }
        if (imode>=3) {
          leg->AddEntry(h_data_reco,"DATA reco","p");
          leg->AddEntry(h_data_unf,"DATA unfold","p");
        }

        leg->Draw();

        pad1->Update();
        c1->Update();

        c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        pad2->Draw();
        pad2->cd();

        TH1F *h_ratio;
        if (imode<=2) h_ratio = (TH1F*) h_mc2_unf->Clone();
        if (imode>=3) h_ratio = (TH1F*) h_data_unf->Clone();

	h_ratio->SetTitle("");
        h_ratio->SetStats(0);

        h_ratio->GetXaxis()->SetTitleOffset(0.9);
        h_ratio->GetXaxis()->SetTitleSize(0.1);
        h_ratio->GetXaxis()->SetLabelFont(42);
        h_ratio->GetXaxis()->SetLabelSize(0.08);
        h_ratio->GetXaxis()->SetTitleFont(42);
        h_ratio->GetYaxis()->SetTitle("reco / truth");
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.09);
        h_ratio->GetYaxis()->SetLabelSize(0.08);
        h_ratio->GetYaxis()->SetRangeUser(-0.2, 2.2);
        if (imode<=0) h_ratio->GetYaxis()->SetRangeUser(0.95, 1.05);
        h_ratio->GetYaxis()->SetTitleOffset(0.4);
        if (imode<=2) h_ratio->Divide(h_mc2_truth);
        if (imode>=3) {
          h_ratio->Divide(h_mc1_truth);
          h_ratio->SetMarkerStyle(20);
        }
        h_ratio->SetMarkerSize(0.7);
        h_ratio->Draw("EPX0");
        if (imode<=2) h_ratio->Draw("EP0SAME");

        TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
        OLine->SetLineColor(kRed);
        OLine->SetLineWidth(1);
        OLine->Draw();

	TCanvas* c2;
	if (method==0) {
          c2 = new TCanvas("c2", "c2", 800, 600);
	  c2->cd();
	  c2->SetLogy();
	  TH1D *d;
	  if (imode<=2) d = unfold_mc.Impl()->GetD();
	  if (imode>=3) d = unfold_data.Impl()->GetD();
	  d->DrawCopy();
	}

	if (plot) {
	  if (ilepton==1) {
	    gSystem->mkdir((path + "/electrons/" + version + "/unfolding/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electrons/" + version + "/unfolding/" + title + "_unfolding.pdf").c_str());
	    if (c2) c2->SaveAs((path + "/electrons/" + version + "/unfolding/" + title + "_unfolding_check.pdf").c_str());
	    if (imode>=3) {
	      TFile f((path + "/electrons/" + version + "/unfolding/" + title + "_unfolding.root").c_str(),"RECREATE");
	      h_data_unf->Write(title.c_str());
	      f.Close();
	    }
	  }
	  if (ilepton==2) {
	    gSystem->mkdir((path + "/muons/" + version + "/unfolding/").c_str(), kTRUE);
	    c1->SaveAs((path + "/muons/" + version + "/unfolding/" + title + "_unfolding.pdf").c_str());
	    if (c2) c2->SaveAs((path + "/muons/" + version + "/unfolding/" + title + "_unfolding_check.pdf").c_str());
	    if (imode>=3) {
	      TFile f((path + "/muons/" + version + "/unfolding/" + title + "_unfolding.root").c_str(),"RECREATE");
	      h_data_unf->Write(title.c_str());
	      f.Close();
	    }
	  }
	}

}
