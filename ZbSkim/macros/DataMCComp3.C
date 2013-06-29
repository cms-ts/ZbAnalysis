#include "LumiLabel.C"
#include "LumiInfo_v09.h"

string path = "/gpfs/cms/users/candelis/work/ZbSkim/test/data/" + version + "/";

void DataMCComp3(string& title="", int plot=0, int ilepton=1) {

        if (ilepton<1 || ilepton>2) {
        ilepton = 1 + ilepton % 2;
	}

	if (title.empty()) title = "w_jetmultiplicity";
	TFile *mc1 = TFile::Open((path+"DYJetsToLL_gen.root").c_str());
	TFile *mc2 = TFile::Open((path+"DYJetsToLL_gen.root").c_str());

        if (ilepton==1) {
          if (title.find("muon")!=string::npos) return;
          if (title.find("mm")!=string::npos) return;
        }
        if (ilepton==2) {
          if (title.find("ele")!=string::npos) return;
          if (title.find("ee")!=string::npos) return;
        }
        
	if (ilepton==1) mc1->cd("demo_ee");
	if (ilepton==2) mc1->cd("demo_mm");  
	TH1F* h_reco = (TH1F*)gDirectory->Get(title.c_str());
	if (ilepton==1) mc2->cd("demo_ee_gen");
	if (ilepton==2) mc2->cd("demo_mm_gen");
	TH1F* h_gen = (TH1F*)gDirectory->Get(title.c_str());

//	TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
//	c1->cd();
          
	h_reco->Sumw2();
	h_gen->Sumw2();
//	double N = h_reco->Integral() / h_gen->Integral();
//	double errN = sqrt(h_reco->Integral()) / h_gen->Integral();
	double N = h_reco->GetEffectiveEntries() / h_gen->GetEffectiveEntries();
	double errN = sqrt(h_reco->GetEffectiveEntries()) / h_gen->GetEffectiveEntries();
//	h_gen->Draw("EPX0");
//	h_gen->SetMarkerStyle(20);
//	h_reco->Draw("histsame");

//    	c1->Update();

        TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
	c2->cd();

	h_reco->SetTitle("");
	h_reco->SetStats(0);
	
	if (title=="w_first_jet_pt") {
	h_reco->GetXaxis ()->SetTitle("leading jet p_{T} [GeV/c]");
	} else if (title=="w_first_jet_eta") {
	h_reco->GetXaxis ()->SetTitle("leading jet #eta");
	} else if (title=="w_pt_Z_ee"||title=="w_pt_Z_mm") {
	h_reco->GetXaxis ()->SetTitle("Z boson p_{T} [GeV/c]");
	} else if (title=="w_delta_phi_ee"||title=="w_delta_phi_mm") {
	h_reco->GetXaxis ()->SetTitle("#Delta #phi(Zj) [rad]");
	} else if (title=="w_Ht") {
        h_reco->GetXaxis ()->SetTitle("H_{T} [GeV/c]");
	}
	if (title=="w_first_bjet_pt") {
	h_reco->GetXaxis ()->SetTitle("leading b-jet p_{T} [GeV/c]");
	h_reco->GetXaxis()->SetRangeUser(0, 200);
	h_reco->GetYaxis()->SetRangeUser(0, 1);
	} else if (title=="w_first_bjet_eta") {
	h_reco->GetYaxis()->SetRangeUser(0, 1);
	h_reco->GetXaxis ()->SetTitle("leading b-jet #eta");
	} else if (title=="w_pt_Z_ee_b"||title=="w_pt_Z_mm_b") {
	h_reco->GetYaxis()->SetRangeUser(0, 1);
	h_reco->GetXaxis()->SetRangeUser(0, 200);
	h_reco->GetXaxis ()->SetTitle("Z boson p_{T} + b-jets [GeV/c]");
	} else if (title=="w_delta_phi_ee_b"||title=="w_delta_phi_mm_b") {
	h_reco->GetYaxis()->SetRangeUser(0, 1);
	h_reco->GetXaxis ()->SetTitle("#Delta #phi(Zb) [rad]");
	} else if (title=="w_Ht_b") {
	h_reco->GetYaxis()->SetRangeUser(0, 1);
	h_reco->GetXaxis()->SetRangeUser(0, 200);
        h_reco->GetXaxis ()->SetTitle("H_{T} [GeV/c]");
	}
	
	//h_reco->GetXaxis()->SetTitle("variable X");
	h_reco->GetXaxis()->SetTitleOffset(0.95);
	h_reco->GetXaxis()->SetTitleSize(0.04);
	h_reco->GetXaxis()->SetLabelFont(42);
	h_reco->GetXaxis()->SetLabelSize(0.04);
	h_reco->GetXaxis()->SetTitleFont(42);
	h_reco->GetYaxis()->SetTitle("#epsilon = N_{Z+b}^{RECO} / N_{Z+b}^{GEN}");
	h_reco->GetYaxis()->SetNdivisions(505);
	h_reco->GetYaxis()->SetTitleSize(0.04);
	h_reco->GetYaxis()->SetLabelSize(0.04);
	//h_reco->GetYaxis()->SetRangeUser(0.5, 1.5);
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
        //for(int i=0; i<h_reco->GetXaxis()->GetNbins(); i++){
	//	cout<<"bin"<<i<<"="<<h_reco->GetBinContent(i)<<endl;
	//}  

	/*TGraphAsymmErrors eff (h_gen, h_reco);
	eff.SetTitle("");
	eff.GetXaxis()->SetTitle("variable X");
	eff.GetXaxis()->SetTitleOffset(0.9);
	eff.GetXaxis()->SetTitleSize(0.04);
	eff.GetXaxis()->SetLabelFont(42);
	eff.GetXaxis()->SetLabelSize(0.04);
	eff.GetXaxis()->SetTitleFont(42);
	eff.GetYaxis()->SetTitle("Efficiency");
	eff.GetYaxis()->SetNdivisions(505);
	eff.GetYaxis()->SetTitleSize(0.04);
	eff.GetYaxis()->SetLabelSize(0.04);
	eff.GetYaxis()->SetRangeUser(0.5, 1.5);
	eff.GetYaxis()->SetTitleOffset(0.9);
	eff.SetMarkerStyle(20);
	eff.Draw();*/
	 
	 if(plot){
	 if (ilepton==1) {
	 gSystem->mkdir(("electrons/" + version + "/efficiency/").c_str());
	 c2->SaveAs(("electrons/" + version + "/efficiency" + "/" + title + "_efficiency" + ".pdf").c_str());
	 }
	 if (ilepton==2) {
	 gSystem->mkdir(("muons/" + version + "/efficiency/").c_str());
	 c2->SaveAs(("muons/" + version + "/efficiency" + "/" + title + "_efficiency" + ".pdf").c_str());
	 }
	 }

}