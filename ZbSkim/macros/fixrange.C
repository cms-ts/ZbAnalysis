TH1F* fixrange(TH1F *old) {

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
    return old;
  }

  int nx = old->GetXaxis()->FindBin(x2)-old->GetXaxis()->FindBin(x1)+1;

  x1 = old->GetXaxis()->GetBinLowEdge(old->GetXaxis()->FindBin(x1));
  x2 = old->GetXaxis()->GetBinUpEdge(old->GetXaxis()->FindBin(x2));

  TH1F *tmp = new TH1F("tmp",old->GetTitle(),nx,x1,x2);
  tmp->Sumw2();

  for (int i=0;i<=old->GetNbinsX()+1;i++) {
    int ii = tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i));
    float c = tmp->GetBinContent(ii)+old->GetBinContent(i);
    float e = TMath::Sqrt(tmp->GetBinError(ii)**2+old->GetBinError(i)**2);
    tmp->SetBinContent(ii,c);
    tmp->SetBinError(ii,e);
  }

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
    return old;
  }

  float y1=x1;
  float y2=x2;

  int nx = old->GetXaxis()->FindBin(x2)-old->GetXaxis()->FindBin(x1)+1;
  int ny = old->GetYaxis()->FindBin(y2)-old->GetYaxis()->FindBin(y1)+1;

  x1 = old->GetXaxis()->GetBinLowEdge(old->GetXaxis()->FindBin(x1));
  x2 = old->GetXaxis()->GetBinUpEdge(old->GetXaxis()->FindBin(x2));
  y1 = old->GetYaxis()->GetBinLowEdge(old->GetYaxis()->FindBin(y1));
  y2 = old->GetYaxis()->GetBinUpEdge(old->GetYaxis()->FindBin(y2));

  TH2F *tmp = new TH2F("tmp",old->GetTitle(),nx,x1,x2,ny,y1,y2);
  tmp->Sumw2();

  for (int i=0;i<=old->GetNbinsX()+1;i++) {
    for (int j=0;j<=old->GetNbinsY()+1;j++) {
      int ii = tmp->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i));
      int jj = tmp->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j));
      float c = tmp->GetBinContent(ii,jj)+old->GetBinContent(i,j);
      float e = TMath::Sqrt(tmp->GetBinError(ii,jj)**2+old->GetBinError(i,j)**2);
      tmp->SetBinContent(ii,jj,c);
      tmp->SetBinError(ii,jj,e);
    }
  }

  old->Delete();
  tmp->SetName(name.c_str());

  return tmp;
}
