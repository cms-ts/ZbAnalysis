#include "DataMCComp.h"
#include "LumiLabel.C"
#include "LumiInfo_v14.h"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";
//string path = "/gpfs/cms/users/lalicata/work/test/data/";

void DataMCComp6(int irun=0, string title="", int plot=0, int isratio=1, int numB=0) {

string bSel="";
string dirbSel="";
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
if (irun==66) {            // irun==66 => unfolding with data weight
  subdir="66";
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
  postfix="";
}
if (numB==1) {
  postfix = postfix + "1b";
  dirbSel="_1b";
  bSel="Z + (= 1) b-jet";
}
if (numB==2) {
  postfix = postfix + "2b";
  dirbSel="_2b";
  bSel="Z + (#geq 2) b-jet";
}

        string title_ee = title;
        string title_ee_b = title;
        string title_mm = title;
        string title_mm_b = title;


        if (numB==0) {
          if (title=="w_pt_Z"||title=="w_mass_Zj"||title=="w_delta_phi") {
            title_ee = title_ee + "_ee";
            title_ee_b = title_ee_b + "_ee";
            title_mm = title_mm + "_mm";
            title_mm_b = title_mm_b + "_mm";
          }

          if (title.find("_bjet_")!=string::npos) {
            title_ee.erase(title_ee.find("_bjet_")+1, 1);
            title_mm.erase(title_mm.find("_bjet_")+1, 1);
          } else {
            title_ee_b = title_ee + "_b";
            title_mm_b = title_mm + "_b";
          }
        }

        if (numB==1 || numB==2) {
          if (title=="w_pt_Z" || title=="w_delta_phi") {
            title_ee = title_ee + "_ee_b";
            title_ee_b = title_ee_b + "_ee_b";
            title_mm = title_mm + "_mm_b";
            title_mm_b = title_mm_b + "_mm_b";
          }
          if (title=="w_Ht") {
            title_ee = title_ee + "_b";
            title_ee_b = title_ee_b + "_b";
            title_mm = title_mm + "_b";
            title_mm_b = title_mm_b + "_b";
          }
        }


        TFile f_ee((path + "/electrons/" + version + "/" + subdir + "/unfolding" + dirbSel + "/" + title_ee + "_unfolding.root").c_str());
        TFile f_ee_b((path + "/electrons/" + version + "/" + subdir + "/unfolding" + dirbSel + "/" + title_ee_b + "_unfolding.root").c_str());
        TH1F* h_data_ee = (TH1F*)f_ee.Get(title_ee.c_str())->Clone();
        TH1F* h_data_ee_b = (TH1F*)f_ee_b.Get(title_ee_b.c_str())->Clone();
        h_data_ee->SetDirectory(0);
        h_data_ee_b->SetDirectory(0);
        f_ee.Close();
        f_ee_b.Close();

        TFile f_mm((path + "/muons/" + version + "/" + subdir + "/unfolding" + dirbSel + "/" + title_mm + "_unfolding.root").c_str());
        TFile f_mm_b((path + "/muons/" + version + "/" + subdir + "/unfolding" + dirbSel + "/" + title_mm_b + "_unfolding.root").c_str());
        TH1F* h_data_mm = (TH1F*)f_mm.Get(title_mm.c_str())->Clone();
        TH1F* h_data_mm_b = (TH1F*)f_mm_b.Get(title_mm_b.c_str())->Clone();
        h_data_mm->SetDirectory(0);
        h_data_mm_b->SetDirectory(0);
        f_mm.Close();
        f_mm_b.Close();

        if (numB!=0) {
          h_data_ee = (TH1F*)h_data_ee->Clone();
          h_data_mm = (TH1F*)h_data_mm->Clone();
        }

        h_data_ee->SetStats(0);
        h_data_ee_b->SetStats(0);
        h_data_mm->SetStats(0);
        h_data_mm_b->SetStats(0);

	h_data_ee->Scale(1./Lumi2012_ele, "width");
	h_data_ee_b->Scale(1./Lumi2012_ele, "width");
	h_data_mm->Scale(1./Lumi2012_muon, "width");
	h_data_mm_b->Scale(1./Lumi2012_muon, "width");

        TCanvas* c1 = new TCanvas("c", "c", 800, 600);
        c1->cd();

        TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
        pad1->SetBottomMargin(0.001);
        pad1->Draw();
        pad1->cd();
        if (isratio==0) pad1->SetLogy();

        h_data_ee->SetTitle("");
        h_data_ee_b->SetTitle("");

        if (isratio==0) {
          h_data_ee->GetYaxis()->SetTitle("#sigma [pb]");
        }
        if (isratio==1) {
          h_data_ee_b->GetYaxis()->SetTitle("#sigma(Z+b) / #sigma(Z+j) [%]");
        }

        h_data_ee->SetMarkerStyle(20);
        h_data_ee->SetMarkerSize(0.7);
        h_data_ee->SetMarkerColor(kRed);
        h_data_ee_b->SetMarkerStyle(24);
        h_data_ee_b->SetMarkerSize(0.7);
        h_data_ee_b->SetMarkerColor(kRed);
        h_data_mm->SetMarkerStyle(20);
        h_data_mm->SetMarkerSize(0.7);
        h_data_mm->SetMarkerColor(kBlue);
        h_data_mm_b->SetMarkerStyle(24);
        h_data_mm_b->SetMarkerSize(0.7);
        h_data_mm_b->SetMarkerColor(kBlue);

        if (isratio==0) {
          h_data_ee->SetMinimum(TMath::Max(0.000002,0.25*h_data_ee_b->GetBinContent(h_data_ee_b->GetMinimumBin())));
          h_data_ee->Draw("EPX");
          h_data_ee_b->Draw("EPXSAME");
          h_data_mm->Draw("EPXSAME");
          h_data_mm_b->Draw("EPXSAME");
        }
        if (isratio==1) {
          h_data_ee_b->Divide(h_data_ee);
          h_data_mm_b->Divide(h_data_mm);
          h_data_ee_b->Draw("EPX");
          h_data_mm_b->Draw("EPXSAME");
        }

        TLegend *leg = new TLegend(0.62, 0.580, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetEntrySeparation(0.01);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);

        if (isratio==0) {
          leg->AddEntry(h_data_ee,"Z(#rightarrow ee) DATA","p");
        }
        leg->AddEntry(h_data_ee_b,"Z(#rightarrow ee)+b DATA","p");
        if (isratio==0) {
          leg->AddEntry(h_data_mm,"Z(#rightarrow #mu#mu) DATA","p");
        }
        leg->AddEntry(h_data_mm_b,"Z(#rightarrow #mu#mu)+b DATA","p");

        leg->Draw();

        pad1->Update();
        c1->Update();

        c1->cd();

        TLatex *latexLabel = CMSPrel((Lumi2012_ele+Lumi2012_muon)/2./1000.,"",0.15,0.94);
        latexLabel->Draw("same");

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        pad2->Draw();
        pad2->cd();

        TH1F* h_ratio = (TH1F*)h_data_ee->Clone("h_ratio");
        TH1F* h_ratio_b = (TH1F*)h_data_ee_b->Clone("h_ratio_b");
        if (isratio==0) {
          h_ratio->Divide(h_data_mm);
          h_ratio_b->Divide(h_data_mm_b);
        }
        if (isratio==1) {
          h_ratio_b->Divide(h_data_mm_b);
        }

        h_ratio_b->SetTitle("");

        h_ratio_b->GetXaxis()->SetTitleOffset(0.9);
        h_ratio_b->GetXaxis()->SetTitleSize(0.1);
        h_ratio_b->GetXaxis()->SetLabelFont(42);
        h_ratio_b->GetXaxis()->SetLabelSize(0.08);
        h_ratio_b->GetXaxis()->SetTitleFont(42);
        h_ratio_b->GetYaxis()->SetTitle("Electrons/Muons");
        h_ratio_b->GetYaxis()->SetNdivisions(505);
        h_ratio_b->GetYaxis()->SetTitleSize(0.09);
        h_ratio_b->GetYaxis()->SetLabelSize(0.08);
        h_ratio_b->GetYaxis()->SetRangeUser(0.5, 1.5);
        h_ratio_b->GetYaxis()->SetTitleOffset(0.4);

        h_ratio->SetMarkerStyle(20);
        h_ratio->SetMarkerColor(kBlack);
        h_ratio_b->SetMarkerStyle(24);
        h_ratio_b->SetMarkerColor(kBlack);

        h_ratio_b->Draw("EPX");
        if (isratio==0) {
          h_ratio->Draw("EPXSAME");
        }

        TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
        OLine->SetLineColor(kGreen);
        OLine->SetLineWidth(1);
        OLine->Draw();

        c1->cd();

	if (plot) {
          if (isratio==0) {
            gSystem->mkdir((path + "/combined/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/").c_str(), kTRUE);
            c1->SaveAs((path + "/combined/" + version + "/" + subdir + "/xsecs_unfolding" + dirbSel + "/" + title + "_xsecs_unfolding.pdf").c_str());
          }
          if (isratio==1) {
            gSystem->mkdir((path + "/combined/" + version + "/" + subdir + "/ratios_unfolding" + dirbSel + "/").c_str(), kTRUE);
            c1->SaveAs((path + "/combined/" + version + "/" + subdir + "/ratios_unfolding" + dirbSel + "/" + title + "_ratios_unfolding.pdf").c_str());
          }
	}

}
