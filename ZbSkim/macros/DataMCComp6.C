#include "LumiLabel.C"
#include "LumiInfo_v09.h"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/";

void DataMCComp6(string& title="", int plot=0, int isratio=1) {

        string title_ee = title;
        string title_ee_b = title;
        string title_mm = title;
        string title_mm_b = title;

        if (title=="w_pt_Z"||title=="w_delta_phi") {
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

        TH1F* h_data;
        TH1F* h_data_b;

        TFile f_ee((path + "/electrons/" + version + "/unfolding/" + title_ee + "_unfolding.root").c_str());
        TFile f_ee_b((path + "/electrons/" + version + "/unfolding/" + title_ee_b + "_unfolding.root").c_str());
        h_data_ee = (TH1F*)f_ee.Get(title_ee.c_str())->Clone();
        h_data_ee_b = (TH1F*)f_ee_b.Get(title_ee_b.c_str())->Clone();
        h_data_ee->SetDirectory(0);
        h_data_ee_b->SetDirectory(0);
        f_ee.Close();
        f_ee_b.Close();

        TFile f_mm((path + "/muons/" + version + "/unfolding/" + title_mm + "_unfolding.root").c_str());
        TFile f_mm_b((path + "/muons/" + version + "/unfolding/" + title_mm_b + "_unfolding.root").c_str());
        h_data_mm = (TH1F*)f_mm.Get(title_mm.c_str())->Clone();
        h_data_mm_b = (TH1F*)f_mm_b.Get(title_mm_b.c_str())->Clone();
        h_data_mm->SetDirectory(0);
        h_data_mm_b->SetDirectory(0);
        f_mm.Close();
        f_mm_b.Close();

        h_data_ee->SetStats(0);
        h_data_ee_b->SetStats(0);
        h_data_mm->SetStats(0);
        h_data_mm_b->SetStats(0);

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

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.3);
        pad2->Draw();
        pad2->cd();

        TH1F* h_ratio = h_data_ee->Clone("h_ratio");
        TH1F* h_ratio_b = h_data_ee_b->Clone("h_ratio_b");
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

        if (isratio==0) {
          gSystem->mkdir((path + "/combined/" + version + "/xsecs/").c_str(), kTRUE);
          c1->SaveAs((path + "/combined/" + version + "/xsecs/" + title + "_xsecs.pdf").c_str());
        }
        if (isratio==1) {
          gSystem->mkdir((path + "/combined/" + version + "/ratios/").c_str(), kTRUE);
          c1->SaveAs((path + "/combined/" + version + "/ratios/" + title + "_ratios.pdf").c_str());
        }

}