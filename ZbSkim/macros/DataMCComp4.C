#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v11.h"

#include "RooUnfold.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldParms.h"
#include "TSVDUnfold_local.h"

#include "fixrange.C"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";

void DataMCComp4(int irun=0, string title="", int plot=0, int ilepton=1, int imode=3, int method=0) {

bool verbose = false;
// bool verbose = true;

// imode = -1; // identity test using MadGraph PAT
// imode =  0; // identity test using MadGraph GEN
// imode =  1; // closure test using MadGraph + Sherpa
// imode =  2; // closure test using MadGraph + Powheg
// imode =  3; // unfolding data with MadGraph
// imode =  4; // unfolding data with Sherpa
// imode =  5; // unfolding data with Powheg
// imode =  6; // unfolding data with MadGraph 4FS

// method = 0; // use SVD
// method = 1; // use Bayes
// method = 2; // use BinByBin

string subdir="0";
string postfix="";
if (irun==1) {             // irun==1 => JEC Up
  subdir="1";
  postfix="Up";
}
if (irun==2) {             // irun==2 => JEC Down
  subdir="2";
  postfix="Down";
}
if (irun==3) {             // irun==3 => PU Up
  subdir="3";
  postfix="Pup";
}
if (irun==4) {             // irun==4 => PU Down
  subdir="4";
  postfix="Pum";
}
if (irun==5) {             // irun==5 => top bkg
  subdir="5";
  postfix="";  
}
if (irun==6) {             // irun==6 => b purity
  subdir="6";
  postfix="";   
}
if (irun==7) {             // irun==7 => unfolding
  subdir="7";
  postfix="";   
}
if (irun==8) {             // irun==8 => unfolding with Sherpa
  subdir="8";
  postfix="";   
}
if (irun==9) {             // irun==9 => unfolding with Powheg
  subdir="9";
  postfix="";   
}
if (irun==10) {            // irun==10 => bkg systematics
  subdir="10";
  postfix="";
}
if (irun==11) {            // irun==11 => JER Up
  subdir="11";
  postfix="JerUp";
}
if (irun==12) {            // irun==12 => JER Down
  subdir="12";
  postfix="JerDown";
}
if (irun==13) {            // irun==13 => bkg statistics
  subdir="13";
  postfix="";
}
if (irun==77) {            // irun==77 => unfolding with MadGraph 4FS
  subdir="77";
  postfix="";
}
if (irun==88) {            // irun==88 => deltaR
  subdir="88";
  postfix="DR";
}
if (irun==99) {            // irun==99 => pur
  subdir="99";
  postfix="Pur";
}

        if (irun==8) imode = 4;
        if (irun==9) imode = 5;
        if (irun==77) imode = 6;

	if (imode<=2 && subdir!="0") return;

	//if (!verbose) gErrorIgnoreLevel = kError;

	/* purity */

        double c_b=1.0;
        double ec_b=0.0;
	double c_c=1.0;
	double ec_c=0.0;
	double c_uds=1.0;
	double ec_uds=0.0;
       
	ifstream in3;
	if (imode>=3) {
	  if (ilepton==1) {
	    in3.open((path + "/electrons/" + version + "/" + subdir + "/distributions/" + "w_BJP_doFit" + ".dat").c_str());
	  }
	  if (ilepton==2) {
	    in3.open((path + "/muons/" + version + "/" + subdir + "/distributions/" + "w_BJP_doFit" + ".dat").c_str());
	  }
	  in3 >> c_uds >> ec_uds;
	  in3 >> c_b >> ec_b;
	  in3 >> c_c >> ec_c;
	  in3.close();
	}

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

	string file = title;

	if (file.find("_b")==string::npos) {
	  if (file.find("_jet_")!=string::npos) {
	    file.insert(file.find("_jet_")+1, "b");
	  } else {
	    file = file + "_b";
	  }
	}

	TFile* data=0;
	if (ilepton==1) data = TFile::Open((path + "/electrons/" + version + "/" + subdir + "/xsecs/" + file + "_xsecs.root").c_str());
	if (ilepton==2) data = TFile::Open((path + "/muons/" + version + "/" + subdir + "/xsecs/" + file + "_xsecs.root").c_str());

	TFile* mc1=0;
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
	if (imode== 6) mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL2_gen.root").c_str());

	TFile* mc2=0;
	if (imode==-1) mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_patgen.root").c_str());
	if (imode== 0) mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	if (imode== 1) mc2 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());
	if (imode== 2) {
	  if (ilepton==1) mc2 = TFile::Open((path + "/" + version + "/" + "DYToEE_powheg_gen.root").c_str());
	  if (ilepton==2) mc2 = TFile::Open((path + "/" + version + "/" + "DYToMuMu_powheg_gen.root").c_str());
	}
	if (imode== 3) mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	if (imode== 4) mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	if (imode== 5) mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	if (imode== 6) mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());

	TH1F* h_data_reco;
	data->cd();
	h_data_reco = (TH1F*)gDirectory->Get((title+"_raw").c_str());

	string title_b = title;

	if (title.find("_b")!=string::npos) {
	  title_b = "b"+title.substr(1);
	}

	TH1F* h_mc1_truth=0;
	TH1F* h_mc1_reco=0;
	TH2F* h_mc1_matrix=0;
	TH1F* h_mc2_truth=0;
	TH1F* h_mc2_reco=0;

	if (ilepton==1) {
	  mc1->cd("demoEleGen");
	  h_mc1_truth  = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd(("demoEle"+postfix).c_str());
	  h_mc1_reco   = (TH1F*)gDirectory->Get(title_b.c_str());
	  mc1->cd(("demoEleDump"+postfix).c_str());
	  h_mc1_matrix = (TH2F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoEleGen");
	  h_mc2_truth  = (TH1F*)gDirectory->Get(title.c_str());
	  mc2->cd(("demoEle"+postfix).c_str());
	  h_mc2_reco   = (TH1F*)gDirectory->Get(title_b.c_str());
	}
	if (ilepton==2) {
	  mc1->cd("demoMuoGen");
	  h_mc1_truth   = (TH1F*)gDirectory->Get(title.c_str());
	  mc1->cd(("demoMuo"+postfix).c_str());
	  h_mc1_reco    = (TH1F*)gDirectory->Get(title_b.c_str());
	  mc1->cd(("demoMuoDump"+postfix).c_str());
	  h_mc1_matrix  = (TH2F*)gDirectory->Get(title.c_str());
	  mc2->cd("demoMuoGen");
	  h_mc2_truth   = (TH1F*)gDirectory->Get(title.c_str());
	  mc2->cd(("demoMuo"+postfix).c_str());
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

// FIX: SHERPA
if (ilepton==2) {
  if (imode==1) h_mc2_reco->Scale(0.9488);
  if (imode==4) {
    h_mc1_reco->Scale(0.9488);
    h_mc1_matrix->Scale(0.9488);
  }
}
// FIX: SHERPA

	RooUnfoldResponse response(h_mc1_reco, h_mc1_truth, h_mc1_matrix);
	response.UseOverflow(kTRUE);
	if (verbose) response.Print();

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
	if (imode==6) {
	  h_mc1_truth->Scale(norm1_3/norm1);
	  h_mc1_reco->Scale(norm1_3/norm1);
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

	RooUnfold::ErrorTreatment err;
	RooUnfold* unfold_mc=0;
	RooUnfold* unfold_data=0;

	int kreg = 0; // default 0 -> nbins/2

	kreg = response.GetNbinsMeasured()/2;
	if (title.find("jet_eta")!=string::npos) kreg = 4;

	if (method==0) {
	  unfold_mc = new RooUnfoldSvd(&response, h_mc2_reco, kreg);
	  unfold_data = new RooUnfoldSvd(&response, h_data_reco, kreg);
	}

	int niter = 4; // default 4 -> number of iterations

	if (method==1) {
	  unfold_mc = new RooUnfoldBayes(&response, h_mc2_reco, niter);
	  unfold_data = new RooUnfoldBayes(&response, h_data_reco, niter);
	}

	if (method==2) {
	  unfold_mc = new RooUnfoldBinByBin(&response, h_mc2_reco);
	  unfold_data = new RooUnfoldBinByBin(&response, h_data_reco);
	}

	if (!verbose) {
	  unfold_mc->SetVerbose(-1);
	  unfold_data->SetVerbose(-1);
	}

	int ntoys = 50; // default 50
	if (irun==7) ntoys = 100;
	if (irun==8) ntoys = 100;
	if (irun==9) ntoys = 100;
	if (irun==77) ntoys = 100;
	unfold_mc->SetNToys(ntoys);
	unfold_data->SetNToys(ntoys);

	int dosys = 0; // default 0 -> 0=stat, 1=stat+sys, 2=sys only
	if (irun==7) dosys = 1;
	if (irun==8) dosys = 1;
	if (irun==9) dosys = 1;
	if (irun==77) dosys = 1;
	unfold_mc->IncludeSystematics(dosys);
	unfold_data->IncludeSystematics(dosys);

	if (imode<=2) {
	  unfold_mc->Print();
	  err = RooUnfold::kErrors;
	  unfold_mc->PrintTable(cout, h_mc2_truth, err);
	}
	if (imode>=3) {
	  unfold_data->Print();
	  err = RooUnfold::kErrors;
	  unfold_data->PrintTable(cout, 0, err);
	}

	TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
	c1->cd();
        TPad* pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
        pad1->SetBottomMargin(0.001);
        pad1->Draw();
        pad1->cd();
        pad1->SetLogy();

	TH1F* h_mc2_unfold=0;
	if (imode<=2) {
	  err = RooUnfold::kErrors;
	  h_mc2_unfold = (TH1F*) unfold_mc->Hreco(err);

	  double vmin = TMath::Max(1.0, 0.1*h_mc2_reco->GetMinimum());
	  h_mc2_unfold->SetMinimum(vmin);

	  double vmax = TMath::Max(0.0, h_mc2_unfold->GetMaximum());
	  vmax = TMath::Max(vmax, h_mc2_reco->GetMaximum());
	  vmax = TMath::Max(vmax, h_mc2_truth->GetMaximum());
	  vmax = TMath::Max(vmax, h_mc1_reco->GetMaximum());
	  vmax = TMath::Max(vmax, h_mc1_truth->GetMaximum());
	  h_mc2_unfold->SetMaximum(1.5*vmax);

	  h_mc2_unfold->SetStats(0);

	  h_mc2_unfold->Draw("HIST");
	  h_mc2_reco->Draw("HISTSAME");
	  h_mc2_truth->Draw("HISTSAME");

	  h_mc2_unfold->SetLineColor(kBlack);
	  h_mc2_reco->SetLineColor(kGreen);
	  h_mc2_truth->SetLineColor(kRed);

	  h_mc1_reco->Draw("HISTSAME");
	  h_mc1_truth->Draw("HISTSAME");

	  h_mc1_reco->SetLineColor(kGreen);
	  h_mc1_truth->SetLineColor(kRed);
	  h_mc1_reco->SetLineStyle(2);
	  h_mc1_truth->SetLineStyle(2);
	}

	TH1F* h_data_unfold=0;
	if (imode>=3) {
	  err = RooUnfold::kErrors;
	  h_data_unfold = (TH1F*) unfold_data->Hreco(err);

	  double vmax = TMath::Max(0.0, h_data_unfold->GetMaximum());
	  vmax = TMath::Max(vmax, h_data_reco->GetMaximum());
	  vmax = TMath::Max(vmax, h_mc1_reco->GetMaximum());
	  vmax = TMath::Max(vmax, h_mc1_truth->GetMaximum());
	  h_data_reco->SetMaximum(1.5*vmax);

	  h_data_reco->SetStats(0);

          h_data_reco->SetLineColor(kGreen);
          h_data_reco->SetMarkerColor(kGreen);
          h_data_reco->SetMarkerStyle(20);
          h_data_reco->SetMarkerSize(0.7);
          h_data_reco->Draw("EPX0");

	  h_data_unfold->SetLineColor(kBlack);
	  h_data_unfold->SetMarkerColor(kBlack);
	  h_data_unfold->SetMarkerStyle(20);
	  h_data_unfold->SetMarkerSize(0.7);
	  h_data_unfold->Draw("EPX0SAME");

          h_mc1_truth->SetLineColor(kRed);
          h_mc1_truth->SetLineStyle(2);
          h_mc1_truth->Draw("HISTSAME");

	  h_mc1_reco->SetLineColor(kGreen);
	  h_mc1_reco->SetLineStyle(2);
	  h_mc1_reco->Draw("HISTSAME");
	}

	TH1F* tmp;
	if (imode<=2) tmp = h_mc2_reco;
	if (imode>=3) tmp = h_data_reco;
	tmp->SetTitle("");
	tmp->SetStats(0);
	tmp->GetYaxis()->SetTitle("Events");

	if (imode<=2) tmp = h_mc2_unfold;
	if (imode>=3) tmp = h_data_unfold;
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

        TLegend* leg1 = new TLegend(0.42, 0.580, 0.68, 0.88);
        leg1->SetBorderSize(0);
        leg1->SetEntrySeparation(0.01);
        leg1->SetFillColor(0);
        leg1->SetFillStyle(0);

        if (imode<=3) {
          leg1->AddEntry(h_mc1_reco,"MADGRAPH reco","l");
          leg1->AddEntry(h_mc1_truth,"MADGRAPH truth","l");
	}
        if (imode==4) {
          leg1->AddEntry(h_mc1_reco,"SHERPA reco","l");
          leg1->AddEntry(h_mc1_truth,"SHERPA truth","l");
	}
        if (imode==5) {
          leg1->AddEntry(h_mc1_reco,"POWHEG reco","l");
          leg1->AddEntry(h_mc1_truth,"POWHEG truth","l");
	}
        if (imode<=0) leg1->AddEntry(h_mc2_unfold,"MADGRAPH unfold","l");
        if (imode==1) {
          leg1->AddEntry(h_mc2_reco,"SHERPA reco","l");
          leg1->AddEntry(h_mc2_truth,"SHERPA truth","l");
          leg1->AddEntry(h_mc2_unfold,"SHERPA unfold","l");
        }
        if (imode==2) {
          leg1->AddEntry(h_mc2_reco,"POWHEG reco","l");
          leg1->AddEntry(h_mc2_truth,"POWHEG truth","l");
          leg1->AddEntry(h_mc2_unfold,"POWHEG unfold","l");
        }
        if (imode>=3) {
          leg1->AddEntry(h_data_reco,"DATA reco","p");
          leg1->AddEntry(h_data_unfold,"DATA unfold","p");
        }

        leg1->Draw();

        TLatex* t = new TLatex();
        t->SetTextSize(0.05);
        t->SetTextFont(42);
        t->SetLineWidth(2);
        t->SetNDC();
        if (method==0) t->DrawLatex(0.13,0.85,"SVD");
        if (method==1) t->DrawLatex(0.13,0.85,"Bayes");
        if (method==2) t->DrawLatex(0.13,0.85,"BinByBin");

        pad1->Update();
        c1->Update();

        c1->cd();

        TPad* pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        pad2->Draw();
        pad2->cd();

        TH1F* h_ratio;
        if (imode<=2) h_ratio = (TH1F*)h_mc2_unfold->Clone();
        if (imode>=3) h_ratio = (TH1F*)h_data_unfold->Clone();

	h_ratio->SetTitle("");
        h_ratio->SetStats(0);

        h_ratio->GetXaxis()->SetTitleOffset(0.9);
        h_ratio->GetXaxis()->SetTitleSize(0.1);
        h_ratio->GetXaxis()->SetLabelFont(42);
        h_ratio->GetXaxis()->SetLabelSize(0.08);
        h_ratio->GetXaxis()->SetTitleFont(42);
        h_ratio->GetYaxis()->SetTitle("unfold / truth");
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.09);
        h_ratio->GetYaxis()->SetLabelSize(0.08);
        h_ratio->GetYaxis()->SetRangeUser(0.3, 1.7);
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

        TLine* OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
        OLine->SetLineColor(kRed);
        OLine->SetLineWidth(1);
        OLine->Draw();

	TCanvas* c2=0;
	if (method==0) {
          c2 = new TCanvas("c2", "c2", 800, 600);
	  c2->cd();
	  c2->SetLogy();
	  TH1D* d;
	  if (imode<=2) d = ((RooUnfoldSvd*)unfold_mc)->Impl()->GetD();
	  if (imode>=3) d = ((RooUnfoldSvd*)unfold_data)->Impl()->GetD();
	  d->DrawCopy();
	  if (method==0) {
	    TMarker *marker = new TMarker(kreg,d->GetYaxis()->GetXmin(),20);
	    marker->SetMarkerColor(kRed);
	    marker->Draw();
	  }
	}

	TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);
	c3->cd();
	c3->SetLogz();
	TH2F* h_response = (TH2F*)response.Hresponse()->Clone();
	h_response->SetStats(0);
	h_response->SetTitle("Response matrix: (x,y)=(measured,truth)");
	h_response->GetXaxis()->SetTitle(tmp->GetXaxis()->GetTitle());
	h_response->GetYaxis()->SetTitle(tmp->GetXaxis()->GetTitle());
	h_response->Draw("COLZ");

        if (method==0) t->DrawLatex(0.13,0.85,"SVD");
        if (method==1) t->DrawLatex(0.13,0.85,"Bayes");
        if (method==2) t->DrawLatex(0.13,0.85,"BinByBin");

	int nv = response.GetNbinsMeasured();
	if (response.UseOverflowStatus()) nv = nv + 2;
	TVectorD err_err(nv);
	TVectorD err_res(nv);
	if (imode<=2) {
	  err_err = unfold_mc->ErecoV(RooUnfold::kErrors);
	  err_res = unfold_mc->ErecoV(RooUnfold::kCovToy);
	}
	if (imode>=3) {
	  err_err = unfold_data->ErecoV(RooUnfold::kErrors);
	  err_res = unfold_data->ErecoV(RooUnfold::kCovToy);
	}

	int ntx = h_mc2_truth->GetNbinsX();
	float xlo = h_mc2_truth->GetXaxis()->GetXmin();
	float xhi = h_mc2_truth->GetXaxis()->GetXmax();
	TH1F* h_err_err = new TH1F("unferr", "Unfolding errors", ntx, xlo, xhi);
	TH1F* h_err_res = new TH1F("toyerr", "Toy MC RMS",       ntx, xlo, xhi);
	for (int i= 0; i<ntx; i++) {
	  h_err_err->SetBinContent(i+1,err_err[i+response.UseOverflowStatus()]);
	  h_err_res->SetBinContent(i+1,err_res[i+response.UseOverflowStatus()]);
	}
	h_err_err->SetMarkerColor(kBlue);
	h_err_err->SetLineColor(kBlue);
	h_err_err->SetMarkerStyle(24);
	h_err_err->SetMinimum(0);
	h_err_res->SetMarkerColor(kRed);
	h_err_res->SetMarkerStyle(4);
	h_err_res->SetMinimum(0);

	TCanvas* c4 = new TCanvas("c4", "c4", 800, 600);
	c4->cd();
	if (imode<=2) {
	  h_err_err->Divide(h_mc2_unfold);
	  h_err_res->Divide(h_mc2_unfold);
	}
	if (imode>=3) {
	  h_err_err->Divide(h_data_unfold);
	  h_err_res->Divide(h_data_unfold);
	}
	h_err_err->Scale(100);
	h_err_res->Scale(100);
	double vmax = TMath::Max(0.0, h_err_err->GetMaximum());
	vmax = TMath::Max(vmax, h_err_res->GetMaximum());
	h_err_err->SetMaximum(1.5*vmax);
	h_err_err->SetStats(0);
	h_err_err->GetXaxis()->SetTitle(tmp->GetXaxis()->GetTitle());
	h_err_err->GetYaxis()->SetTitle("Error [%]");
	h_err_err->Draw("HIST");
	h_err_res->Draw("PSAME");

	TLegend* leg2 = new TLegend (0.70, 0.75, 0.89, 0.89);
	leg2->SetBorderSize(0);
	leg2->SetEntrySeparation(0.01);
	leg2->SetFillColor(0);
	leg2->SetFillStyle(0);

	leg2->AddEntry(h_err_err, "Unfolding errors", "L");
	leg2->AddEntry(h_err_res, "Toy MC RMS",       "P");

	leg2->Draw();

        if (method==0) t->DrawLatex(0.13,0.85,"SVD");
        if (method==1) t->DrawLatex(0.13,0.85,"Bayes");
        if (method==2) t->DrawLatex(0.13,0.85,"BinByBin");

	TH2F* h_err_cov;
	err = RooUnfold::kCovariance;
	if (imode<=2) h_err_cov = new TH2F(TMatrix(unfold_mc->Ereco(err)));
	if (imode>=3) h_err_cov = new TH2F(TMatrix(unfold_data->Ereco(err)));

	TCanvas* c5 = new TCanvas("c5", "c5", 800, 600);
	c5->cd();
	h_err_cov->SetStats(0);
	h_err_cov->SetTitle("Covariance matrix (including under/overflows)");
	h_err_cov->GetXaxis()->SetTitle(tmp->GetXaxis()->GetTitle());
	h_err_cov->GetYaxis()->SetTitle(tmp->GetXaxis()->GetTitle());
	gStyle->SetHistMinimumZero();
	h_err_cov->SetMarkerSize(0.85);
	h_err_cov->GetXaxis()->SetNdivisions(020);
	h_err_cov->GetYaxis()->SetNdivisions(020);
	gStyle->SetGridStyle(1);
	gStyle->SetPaintTextFormat("+11.1e");
	h_err_cov->Draw("TEXT30");

        if (method==0) t->DrawLatex(0.01,0.95,"SVD");
        if (method==1) t->DrawLatex(0.01,0.95,"Bayes");
        if (method==2) t->DrawLatex(0.01,0.95,"BinByBin");

	RooUnfoldParms* parms;
	err = RooUnfold::kErrors;
	if (imode<=2) parms = new RooUnfoldParms(unfold_mc, err, h_mc2_truth);
	if (imode>=3) parms = new RooUnfoldParms(unfold_data, err);

	float maxparm=0.0;
	if (method==0) maxparm = response.GetNbinsMeasured();
	if (method==1) maxparm = 50.0;
	if (method==2) maxparm = 1.0;
	parms->SetMinParm(1.0);
	parms->SetMaxParm(maxparm);
	parms->SetStepSizeParm(1.0);

	TProfile* hParmChi2 = parms->GetChi2();
	TProfile* hParmErr = parms->GetRMSError();
	TProfile* hParmRes = parms->GetMeanResiduals();
	TH1* hParmRms = parms->GetRMSResiduals();

	hParmChi2->SetStats(0);
	hParmErr->SetStats(0);
	hParmRes->SetStats(0);
	hParmRms->SetStats(0);

	TCanvas* c6 = new TCanvas("c6", "c6", 800, 600);
	c6->cd();
	c6->Divide(2, 2);
	c6->cd(1);
	hParmChi2->Draw("P");
	hParmChi2->GetXaxis()->SetTitle("regparm");
        if (method==0) t->DrawLatex(0.13,0.85,"SVD");
        if (method==1) t->DrawLatex(0.13,0.85,"Bayes");
        if (method==2) t->DrawLatex(0.13,0.85,"BinByBin");
	c6->cd(2);
	hParmErr->Draw("P");
	hParmErr->GetXaxis()->SetTitle("regparm");
        if (method==0) t->DrawLatex(0.13,0.85,"SVD");
        if (method==1) t->DrawLatex(0.13,0.85,"Bayes");
        if (method==2) t->DrawLatex(0.13,0.85,"BinByBin");
	c6->cd(3);
	hParmRes->Draw("P");
	hParmRes->GetXaxis()->SetTitle("regparm");
        if (method==0) t->DrawLatex(0.13,0.85,"SVD");
        if (method==1) t->DrawLatex(0.13,0.85,"Bayes");
        if (method==2) t->DrawLatex(0.13,0.85,"BinByBin");
	c6->cd(4);
	hParmRms->Draw("P");
	hParmRms->GetXaxis()->SetTitle("regparm");
        if (method==0) t->DrawLatex(0.13,0.85,"SVD");
        if (method==1) t->DrawLatex(0.13,0.85,"Bayes");
        if (method==2) t->DrawLatex(0.13,0.85,"BinByBin");

	if (imode==0) title = title + "_identity_madgraph";
	if (imode==1) title = title + "_closure_sherpa";
	if (imode==2) title = title + "_closure_powheg";

	if (plot) {
	  ofstream out;
	  if (ilepton==1) {
	    gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/unfolding/").c_str(), kTRUE);
	    c1->SaveAs((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding.pdf").c_str());
	    if (c2) c2->SaveAs((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding_check.pdf").c_str());
	    if (c3) c3->SaveAs((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding_response.pdf").c_str());
	    if (c4) c4->SaveAs((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding_errors.pdf").c_str());
	    if (c5) c5->SaveAs((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding_covariance.pdf").c_str());
	    if (c6) c6->SaveAs((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding_scan.pdf").c_str());
	    if (imode>=3) {
	      TFile f((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding.root").c_str(),"RECREATE");
	      h_data_unfold->Write(title.c_str());
	      f.Close();
	      out.open((path + "/electrons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding.dat").c_str());
	    }
	  }
	  if (ilepton==2) {
	    gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/unfolding/").c_str(), kTRUE);
	    c1->SaveAs((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding.pdf").c_str());
	    if (c2) c2->SaveAs((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding_check.pdf").c_str());
	    if (c3) c3->SaveAs((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding_response.pdf").c_str());
	    if (c4) c4->SaveAs((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding_errors.pdf").c_str());
	    if (c5) c5->SaveAs((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding_covariance.pdf").c_str());
	    if (c6) c6->SaveAs((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding_scan.pdf").c_str());
	    if (imode>=3) {
	      TFile f((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding.root").c_str(),"RECREATE");
	      h_data_unfold->Write(title.c_str());
	      f.Close();
	      out.open((path + "/muons/" + version + "/" + subdir + "/unfolding/" + title + "_unfolding.dat").c_str());
	    }
	  }
	  if (imode>=3) {
	    out << std::fixed << std::setw( 11 ) << std::setprecision( 2 );
	    out << h_data_reco->Integral(0,h_data_reco->GetNbinsX()+1) << endl;
	    out << std::fixed << std::setw( 11 ) << std::setprecision( 2 );
	    out << h_data_unfold->Integral(0,h_data_unfold->GetNbinsX()+1) << endl;
	    out << std::fixed << std::setw( 11 ) << std::setprecision( 2 );
	    out << h_mc1_truth->Integral(0,h_mc1_truth->GetNbinsX()+1) << endl;
	    out.close();
	  }
	}

}
