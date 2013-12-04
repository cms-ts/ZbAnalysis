#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

TH1F* fixrange(TH1F* old) {

  float x1, x2;
  string name = old->GetName();

  if (name.find("Ht")!=string::npos) {
    x1 = 30.;
    x2 = 250.;
  } else if (name.find("jet_pt")!=string::npos) {
    x1 = 30.;
    x2 = 200.;
  } else if (name.find("pt_Z")!=string::npos) {
    x1 = 0.;
    x2 = 200.;
  } else {
    x1 = old->GetXaxis()->GetBinCenter(1);
    x2 = old->GetXaxis()->GetBinCenter(old->GetNbinsX());
  }

  int nx = old->GetXaxis()->FindBin(x2)-old->GetXaxis()->FindBin(x1)+1;

  x1 = old->GetXaxis()->GetBinLowEdge(old->GetXaxis()->FindBin(x1));
  x2 = old->GetXaxis()->GetBinUpEdge(old->GetXaxis()->FindBin(x2));

  TH1F* tmp = new TH1F("tmp",old->GetTitle(),nx,x1,x2);
  tmp->Sumw2();

  tmp->GetXaxis()->SetTitle(old->GetXaxis()->GetTitle());
  tmp->GetYaxis()->SetTitle(old->GetYaxis()->GetTitle());

  for (int i=0;i<=old->GetNbinsX()+1;i++) {
    int ii = tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i));
    float c1 = tmp->GetBinContent(ii);
    float e1 = tmp->GetBinError(ii);
    float c2 = old->GetBinContent(i);
    float e2 = old->GetBinError(i);

    tmp->SetBinContent(ii,c1+c2);
    tmp->SetBinError(ii,TMath::Sqrt(e1*e1+e2*e2));
  }

  tmp->SetEntries(old->GetEntries());

  old->Delete();
  tmp->SetName(name.c_str());

  return tmp;
}

TH2F* fixrange(TH2F* old) {

  float x1, x2;
  string name = old->GetName();

  if (name.find("Ht")!=string::npos) {
    x1 = 30.;
    x2 = 250.;
  } else if (name.find("jet_pt")!=string::npos) {
    x1 = 30.;
    x2 = 200.;
  } else if (name.find("pt_Z")!=string::npos) {
    x1 = 0.;
    x2 = 200.;
  } else {
    x1 = old->GetXaxis()->GetBinCenter(1);
    x2 = old->GetXaxis()->GetBinCenter(old->GetNbinsX());
  }

  float y1=x1;
  float y2=x2;

  int nx = old->GetXaxis()->FindBin(x2)-old->GetXaxis()->FindBin(x1)+1;
  int ny = old->GetYaxis()->FindBin(y2)-old->GetYaxis()->FindBin(y1)+1;

  x1 = old->GetXaxis()->GetBinLowEdge(old->GetXaxis()->FindBin(x1));
  x2 = old->GetXaxis()->GetBinUpEdge(old->GetXaxis()->FindBin(x2));
  y1 = old->GetYaxis()->GetBinLowEdge(old->GetYaxis()->FindBin(y1));
  y2 = old->GetYaxis()->GetBinUpEdge(old->GetYaxis()->FindBin(y2));

  TH2F* tmp = new TH2F("tmp",old->GetTitle(),nx,x1,x2,ny,y1,y2);
  tmp->Sumw2();

  tmp->GetXaxis()->SetTitle(old->GetXaxis()->GetTitle());
  tmp->GetYaxis()->SetTitle(old->GetYaxis()->GetTitle());

  for (int i=0;i<=old->GetNbinsX()+1;i++) {
    for (int j=0;j<=old->GetNbinsY()+1;j++) {
      int ii = tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i));
      int jj = tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j));
      float c1 = tmp->GetBinContent(ii,jj);
      float e1 = tmp->GetBinError(ii,jj);
      float c2 = old->GetBinContent(i,j);
      float e2 = old->GetBinError(i,j);
//
// underflows and overflows *must* be discarded before the unfolding:
// inefficiencies and fakes/background are taken from the additional entries
// in the 'reco' and 'truth' histograms, respectively
//
      if (ii==0) continue;
      if (ii==tmp->GetNbinsX()+1) continue;
      if (jj==0) continue;
      if (jj==tmp->GetNbinsY()+1) continue;

      tmp->SetBinContent(ii,jj,c1+c2);
      tmp->SetBinError(ii,jj,TMath::Sqrt(e1*e1+e2*e2));
    }
  }

  tmp->SetEntries(old->GetEntries());

  old->Delete();
  tmp->SetName(name.c_str());

  return tmp;
}
