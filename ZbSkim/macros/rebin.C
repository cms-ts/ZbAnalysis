#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

TH1F* rebin(TH1F* old, int numB) {

  string name = old->GetName();

  int nb=0;
  float nbins[100];

  if (numB!=2 && (name.find("first_jet_pt")!=string::npos||name.find("first_bjet_pt")!=string::npos)) {
    nb = 14;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
    nbins[10] = 1;
    nbins[11] = 2;
    nbins[12] = 3;
    nbins[13] = 5;
  } else if (numB!=2 && name.find("pt_Z")!=string::npos) {
    nb = 19;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
    nbins[10] = 1;
    nbins[11] = 1;
    nbins[12] = 1;
    nbins[13] = 1;
    nbins[14] = 2;
    nbins[15] = 2;
    nbins[16] = 3;
    nbins[17] = 4;
    nbins[18] = 6;
  } else if (numB!=2 && name.find("Ht")!=string::npos) {
    nb = 15;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
    nbins[10] = 1;
    nbins[11] = 2;
    nbins[12] = 3;
    nbins[13] = 4;
    nbins[14] = 5;
  } else if (numB==2 && name.find("first_bjet_pt")!=string::npos) {
    nb = 12;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
    nbins[10] = 1;
    nbins[11] = 2;
  } else if (numB==2 && name.find("second_bjet_pt")!=string::npos) {
    nb = 10;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
  } else if (numB==2 && name.find("pt_Z")!=string::npos) {
    nb = 18;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
    nbins[10] = 1;
    nbins[11] = 1;
    nbins[12] = 1;
    nbins[13] = 1;
    nbins[14] = 1;
    nbins[15] = 1;
    nbins[16] = 3;
    nbins[17] = 5;
  } else if (numB==2 && (name.find("Zbb_mass")!=string::npos||name.find("eebb_mass")!=string::npos||name.find("mmbb_mass")!=string::npos||name.find("embb_mass")!=string::npos)) {
    nb = 9;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 2;
    nbins[8] = 3;
  } else if (numB==2 && name.find("bb_mass")!=string::npos) {
    nb = 11;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
    nbins[10] = 3;
  } else if (numB==2 && name.find("A_")!=string::npos) {
    nb = 8;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 3;
  } else {
    return old;
  }

  int n=0;
  for (int i=0; i<nb; i++) {
    n=n+nbins[i];
  }
  if (n!=old->GetNbinsX()) {
    cout << "error: " << name << " - " << n << " " << old->GetNbinsX() << endl;
    return old;
  }

  int s=0;
  const int nbs=nb+1;
  Double_t xbins[nbs];
  for (int i=0; i<=nb; i++) {
    if (i>0) s=s+nbins[i-1];
    xbins[i] = old->GetXaxis()->GetBinUpEdge(s);
  }

  TH1F* tmp = new TH1F("tmp", old->GetTitle(),nb,xbins);

  for (int i=1; i<=old->GetXaxis()->GetNbins()+1; i++) {
    float c1 = old->GetBinContent(i);
    float e1 = old->GetBinError(i);
    float c1b = tmp->GetBinContent(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)));
    float e1b = tmp->GetBinError(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)));
    tmp->SetBinContent(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), c1+c1b);
    tmp->SetBinError(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), TMath::Sqrt(e1*e1+e1b*e1b));
  }

  old->Delete();
  tmp->SetName(name.c_str());

  return tmp;
}

TH2F* rebin(TH2F* old, int numB) {

  string name = old->GetName();

  int nb=0;
  float nbins[100];

  if (numB!=2 && (name.find("first_jet_pt")!=string::npos||name.find("first_bjet_pt")!=string::npos)) {
    nb = 14;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
    nbins[10] = 1;
    nbins[11] = 2;
    nbins[12] = 3;
    nbins[13] = 5;
  } else if (numB!=2 && name.find("pt_Z")!=string::npos) {
    nb = 19;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
    nbins[10] = 1;
    nbins[11] = 1;
    nbins[12] = 1;
    nbins[13] = 1;
    nbins[14] = 2;
    nbins[15] = 2;
    nbins[16] = 3;
    nbins[17] = 4;
    nbins[18] = 6;
  } else if (numB!=2 && name.find("Ht")!=string::npos) {
    nb = 15;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
    nbins[10] = 1;
    nbins[11] = 2;
    nbins[12] = 3;
    nbins[13] = 4;
    nbins[14] = 5;
  } else if (numB==2 && name.find("first_bjet_pt")!=string::npos) {
    nb = 12;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
    nbins[10] = 1;
    nbins[11] = 2;
  } else if (numB==2 && name.find("second_bjet_pt")!=string::npos) {
    nb = 10;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
  } else if (numB==2 && name.find("pt_Z")!=string::npos) {
    nb = 18;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
    nbins[10] = 1;
    nbins[11] = 1;
    nbins[12] = 1;
    nbins[13] = 1;
    nbins[14] = 1;
    nbins[15] = 1;
    nbins[16] = 3;
    nbins[17] = 5;
  } else if (numB==2 && (name.find("Zbb_mass")!=string::npos||name.find("eebb_mass")!=string::npos||name.find("mmbb_mass")!=string::npos||name.find("embb_mass")!=string::npos)) {
    nb = 9;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 2;
    nbins[8] = 3;
  } else if (numB==2 && name.find("bb_mass")!=string::npos) {
    nb = 11;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 1;
    nbins[8] = 1;
    nbins[9] = 1;
    nbins[10] = 3;
  } else if (numB==2 && name.find("A_")!=string::npos) {
    nb = 8;
    nbins[0] = 1;
    nbins[1] = 1;
    nbins[2] = 1;
    nbins[3] = 1;
    nbins[4] = 1;
    nbins[5] = 1;
    nbins[6] = 1;
    nbins[7] = 3;
  } else {
    return old;
  }

  int n=0;
  for (int i=0; i<nb; i++) {
    n=n+nbins[i];
  }
  if (n!=old->GetNbinsX()) {
    cout << "error: " << n << " " << old->GetNbinsX() << endl;
    return old;
  }

  int s=0;
  const int nbs=nb+1;
  Double_t xbins[nbs];
  for (int i=0; i<=nb; i++) {
    if (i>0) s=s+nbins[i-1];
    xbins[i] = old->GetXaxis()->GetBinUpEdge(s);
  }

  TH2F* tmp = new TH2F("tmp", old->GetTitle(),nb,xbins,nb,xbins);

  for (int i=1; i<=old->GetXaxis()->GetNbins(); i++) {
    for (int j=1; j<=old->GetYaxis()->GetNbins(); j++) {
      float c2 = old->GetBinContent(i,j);
      float e2 = old->GetBinError(i,j);
      float c2b = tmp->GetBinContent(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)),tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j)));
      float e2b = tmp->GetBinError(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)),tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j)));
      tmp->SetBinContent(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j)), c2+c2b);
      tmp->SetBinError(tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i)), tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j)), TMath::Sqrt(e2*e2+e2b*e2b));
    }
  }

  old->Delete();
  tmp->SetName(name.c_str());

  return tmp;
}

