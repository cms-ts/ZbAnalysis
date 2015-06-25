#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v14.h"

#include "fixrange.C"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";
//string path = "/gpfs/cms/users/lalicata/work/test/data/";

void DataMCComp3(int irun=0, string title="", int plot=0, int ilepton=1, int numB=0) {

int useDY = 0; // use MadGraph
//int useDY = 1; // use Sherpa
//int useDY = 2; // use Powheg
//int useDY = 3; // use MG_aMC@NLO+P8

//int useWeights = 0; // do not use weights for numB=0
int useWeights = 1; // use weights for numB=0

string subdir="0";
string postfix="";
string dirbSel="";
string bSel="";
string genPostfix="";

bool bbBkg = false;
bool bbSig = false;

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
if (irun==14) {            // irun==14 => unfolding with MadGraph 4FS
  subdir="14";
  postfix="";
}
if (irun==15) {            // irun==15 => unfolding with data weight
  subdir="15";
  postfix="";
}
if (irun==16) {            // irun==16 => unfolding with MadGraph aMC@NLO
  subdir="16";
  postfix="";
}
if (irun==17) {            // irun==17 => templates from data
  subdir="17";
  postfix="";
}
if (irun==18) {            // irun==18 => templates from Sherpa
  subdir="18";
  postfix="";
}
if (irun==19) {            // irun==19 => templates from MadGraph aMC@NLO
  subdir="19";
  postfix="";
}
if (numB==1) {
  postfix = postfix + "1b";
  dirbSel = "_1b";
  bSel = "Z + (= 1) b-jet";
  genPostfix = "1b";
}
if (numB==2) {
  postfix = postfix + "2b";
  dirbSel = "_2b";
  bSel = "Z + (#geq 2) b-jet";
  genPostfix = "2b";
}

if (numB==1) bbBkg = true;
if (numB==2) bbSig = true;

	if (irun==16) useDY = 3;

	if (title.empty()) title = "w_jetmultiplicity";

        if (ilepton==1) {
          if (title.find("muon")!=string::npos) return;
          if (title.find("mm")!=string::npos) return;
        }
        if (ilepton==2) {
          if (title.find("ele")!=string::npos) return;
          if (title.find("ee")!=string::npos) return;
        }

	TFile* mc1 = 0;
	TFile* mc2 = 0;
	if (useDY==0) {
	  mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	  mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen.root").c_str());
	  if (useWeights && numB==0) {
	    mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen_weights.root").c_str());
	    mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_gen_weights.root").c_str());
	  }
	}
	if (useDY==1) {
	  mc1 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());
	  mc2 = TFile::Open((path + "/" + version + "/" + "DYJets_sherpa_gen.root").c_str());
	}
	if (useDY==2) {
	  if (ilepton==1) {
	    mc1 = TFile::Open((path + "/" + version + "/" + "DYToEE_powheg_gen.root").c_str());
	    mc2 = TFile::Open((path + "/" + version + "/" + "DYToEE_powheg_gen.root").c_str());
	  }
	  if (ilepton==2) {
	    mc1 = TFile::Open((path + "/" + version + "/" + "DYToMuMu_powheg_gen.root").c_str());
	    mc2 = TFile::Open((path + "/" + version + "/" + "DYToMuMu_powheg_gen.root").c_str());
	  }
	}
	if (useDY==3) {
	  mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_aMC_gen.root").c_str());
	  mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_aMC_gen.root").c_str());
	  if (useWeights && numB==0) {
	    mc1 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_aMC_gen_weights.root").c_str());
	    mc2 = TFile::Open((path + "/" + version + "/" + "DYJetsToLL_aMC_gen_weights.root").c_str());
	  }
	}

/* efficiency:  e_Z / e_Zb = e_Z / e_Z_1 * e_Z_b */

int itype = 0; // e_Z and e_Zb = e_Z_1 * e_Z_b
//int itype = 1; // e_Z_1
//int itype = 2; // e_Z_b

	string title_b = title;

        //if (title.find("_b")!=string::npos) {
	if (title.find("_b")!=string::npos || numB==1 || numB==2) {
	  if (itype==0) title_b = "b"+title.substr(1);
	  if (itype==2) title_b = "b"+title.substr(1);
	}

        TH1F* h_reco_bbBkg = 0;
        TH1F* h_reco_bbSig = 0;

	if (ilepton==1&&itype==0) mc1->cd(("demoEle"+postfix).c_str());
	if (ilepton==2&&itype==0) mc1->cd(("demoMuo"+postfix).c_str());
	if (ilepton==1&&itype==1) mc1->cd(("demoEleBtag"+genPostfix).c_str());
	if (ilepton==2&&itype==1) mc1->cd(("demoMuoBtag"+genPostfix).c_str());
	if (ilepton==1&&itype==2) mc1->cd(("demoEle"+postfix).c_str());
	if (ilepton==2&&itype==2) mc1->cd(("demoMuo"+postfix).c_str());
	TH1F* h_reco = (TH1F*)gDirectory->Get(title_b.c_str());
        if (bbBkg) h_reco_bbBkg = (TH1F*)gDirectory->Get(("bbBkg"+title.substr(1)).c_str());
	if (bbSig) h_reco_bbSig = (TH1F*)gDirectory->Get(("bbSig"+title.substr(1)).c_str());
	if (ilepton==1&&itype==0) mc2->cd(("demoEleGen"+genPostfix).c_str());
	if (ilepton==2&&itype==0) mc2->cd(("demoMuoGen"+genPostfix).c_str());
	if (ilepton==1&&itype==1) mc2->cd(("demoEleGen"+genPostfix).c_str());
	if (ilepton==2&&itype==1) mc2->cd(("demoMuoGen"+genPostfix).c_str());
	if (ilepton==1&&itype==2) mc2->cd(("demoEleBtag"+genPostfix).c_str());
	if (ilepton==2&&itype==2) mc2->cd(("demoMuoBtag"+genPostfix).c_str());
	TH1F* h_gen = (TH1F*)gDirectory->Get(title.c_str());

	if (bbBkg) {
	  for (int i=0; i<=h_reco->GetNbinsX()+1; i++) {
	    double c = h_reco->GetBinContent(i);
	    c = c - h_reco_bbBkg->GetBinContent(i);
	    if (c<0) c = 0.0;
	    h_reco->SetBinContent(i,c);
	    double e = TMath::Power(h_reco->GetBinError(i),2);
	    e = e - TMath::Power(h_reco_bbBkg->GetBinError(i),2);
if (e<0) {
cout << "e31 " << i << " " << h_reco->GetName() << " " << c << " " << e << endl;
}
	    if (e<0) e = 0.0;
	    h_reco->SetBinError(i,TMath::Sqrt(e));
	  }
	}
	if (bbSig) h_reco = h_reco_bbSig;

	h_reco = fixrange(h_reco, numB);
	h_gen = fixrange(h_gen, numB);

	double N = 1.0;
	double errN = 0.0;
	N = h_reco->IntegralAndError(0,h_reco->GetNbinsX()+1,errN) / h_gen->Integral(0,h_gen->GetNbinsX()+1);
	errN = errN / h_gen->Integral(0,h_gen->GetNbinsX()+1);

//	TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
//	c1->cd();

//	h_gen->Draw("EPX0");
//	h_gen->SetMarkerStyle(20);
//	h_reco->Draw("histsame");

//    	c1->Update();

        TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
	c2->cd();

	h_reco->SetTitle("");
	h_reco->SetStats(0);

	if (title=="w_first_jet_pt") {
	  h_reco->GetXaxis()->SetTitle("leading jet p_{T} [GeV/c]");
	} else if (title=="w_first_jet_eta") {
	  h_reco->GetXaxis()->SetTitle("leading jet #eta");
	} else if (title=="w_pt_Z_ee"||title=="w_pt_Z_mm") {
	  h_reco->GetXaxis()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_delta_phi_ee"||title=="w_delta_phi_mm") {
	  h_reco->GetXaxis()->SetTitle("#Delta #phi(Zj) [rad]");
	} else if (title=="w_Ht") {
          h_reco->GetXaxis()->SetTitle("H_{T} [GeV/c]");
	}
	if (title=="w_first_bjet_pt") {
	  h_reco->GetXaxis()->SetTitle("leading b-jet p_{T} [GeV/c]");
	  h_reco->GetXaxis()->SetRangeUser(30, 300);
	  if (numB==2) h_reco->GetXaxis()->SetRangeUser(30, 200);
	  h_reco->GetYaxis()->SetRangeUser(0, 1);
	} else if (title=="w_second_bjet_pt") {
	  h_reco->GetXaxis()->SetTitle("sub-leading b-jet p_{T} [GeV/c]");
	  h_reco->GetXaxis()->SetRangeUser(30, 300);
	  if (numB==2) h_reco->GetXaxis()->SetRangeUser(30, 120);
	  h_reco->GetYaxis()->SetRangeUser(0, 1);
	} else if (title=="w_first_bjet_eta") {
	  h_reco->GetXaxis()->SetTitle("leading b-jet #eta");
	  h_reco->GetYaxis()->SetRangeUser(0, 1);
	} else if (title=="w_pt_Z_ee_b"||title=="w_pt_Z_mm_b") {
	  h_reco->GetXaxis()->SetTitle("Z boson p_{T} [GeV/c]");
	  h_reco->GetXaxis()->SetRangeUser(0, 300);
	  if (numB==2) h_reco->GetXaxis()->SetRangeUser(0, 230);
	  h_reco->GetYaxis()->SetRangeUser(0, 1);
	} else if (title=="w_delta_phi_ee_b"||title=="w_delta_phi_mm_b") {
	  h_reco->GetXaxis()->SetTitle("#Delta #phi(Zb) [rad]");
	  h_reco->GetYaxis()->SetRangeUser(0, 1);
	} else if (title=="w_Ht_b") {
          h_reco->GetXaxis()->SetTitle("H_{T} [GeV/c]");
	  h_reco->GetXaxis()->SetRangeUser(30, 500);
	  if (numB==2) h_reco->GetXaxis()->SetRangeUser(30, 400);
	  h_reco->GetYaxis()->SetRangeUser(0, 1);
	}

	h_reco->GetXaxis()->SetTitleOffset(0.95);
	h_reco->GetXaxis()->SetTitleSize(0.04);
	h_reco->GetXaxis()->SetLabelFont(42);
	h_reco->GetXaxis()->SetLabelSize(0.04);
	h_reco->GetXaxis()->SetTitleFont(42);
	if (title_b == title) {
	  h_reco->GetYaxis()->SetTitle("#epsilon = N_{Z}^{RECO} / N_{Z}^{GEN}");
	} else {
	  h_reco->GetYaxis()->SetTitle("#epsilon = N_{Z+b}^{RECO} / N_{Z+b}^{GEN}");
	}
	h_reco->GetYaxis()->SetNdivisions(505);
	h_reco->GetYaxis()->SetTitleSize(0.04);
	h_reco->GetYaxis()->SetLabelSize(0.04);
	h_reco->GetYaxis()->SetRangeUser(0.0, 1.5);
	h_reco->GetYaxis()->SetTitleOffset(1.04);
	h_reco->SetMarkerStyle(20);

	h_reco->Divide(h_gen);
	h_reco->Draw("EPX0");

        TLatex *Label = new TLatex();
	Label->SetTextSize(0.0275);
	Label->SetTextFont(42);
	Label->SetLineWidth(2);
	Label->SetNDC();
	char buff[100];
	sprintf(buff, "< #epsilon > = #frac{#int RECO}{#int GEN} = %5.3f #pm %5.3f", N, errN);
	Label->DrawLatex(0.50, 0.68, buff);

	c2->Update();

	if (useDY==1) subdir = subdir + "_sherpa";
	if (useDY==2) subdir = subdir + "_powheg";

	if (plot) {
	  ofstream out;
	  if (ilepton==1) {
	    gSystem->mkdir((path + "/electrons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/").c_str(), kTRUE);
	    c2->SaveAs((path + "/electrons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + title + "_efficiency.pdf").c_str());
	    TFile f((path + "/electrons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + title + "_efficiency.root").c_str(),"RECREATE");
	    h_reco->Write(title.c_str());
	    f.Close();
	    out.open((path + "/electrons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + title + "_efficiency.dat").c_str());
	  }
	  if (ilepton==2) {
	    gSystem->mkdir((path + "/muons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/").c_str(), kTRUE);
	    c2->SaveAs((path + "/muons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + title + "_efficiency.pdf").c_str());
	    TFile f((path + "/muons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + title + "_efficiency.root").c_str(),"RECREATE");
	    h_reco->Write(title.c_str());
	    f.Close();
	    out.open((path + "/muons/" + version + "/" + subdir + "/efficiency" + dirbSel + "/" + title + "_efficiency.dat").c_str());
	  }
	  out << N << " " << errN << endl;
	  out.close();
	}
}

