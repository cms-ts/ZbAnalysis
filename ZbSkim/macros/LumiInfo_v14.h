string version = "v14";

//////////////////// Integrated luminosity in 1/picobarn

double Lumi2012_ele  = 19789.0; // full 2012 electrons luminosity
double Lumi2012_muon = 19751.0; // full 2012 muons luminosity
double Lumi2012_ele_muon = 19780.0; // full 2012 electrons+muons luminosity

// single period lumi, same for electrons and muons dataset in 1/picobarn

double Lumi_2012A_22Jan =  889.362;
double Lumi_2012B_22Jan = 4429.000;
double Lumi_2012C_22Jan = 7152.000;
double Lumi_2012D_22Jan = 7318.000;

//////////////////////// DY MadGraph

double Ngen_dy = 30459503;
double Xsec_dy = 3531.8; // NNLO
double eXsec_dy = 0;

//////////////////////// DY MadGraph_aMC@NLO + Pythia 8

//double Ngen_dy_amc = 28944033;
// weights extracted from h_gen_weights->GetBinContent(2)
double Ngen_dy_amc = 2.59285811200000000e+11;
double Xsec_dy_amc = 3424.0; // NNLO
double eXsec_dy_amc = 0;

//////////////////////// DY Sherpa

//double Ngen_dy_1 = 127014144;
// weights extracted from h_gen_weights->GetBinContent(2)
double Ngen_dy_1 = 4.37323040000000000e+07; 
double Xsec_dy_1 = 3531.8; // NNLO
double eXsec_dy_1 = 0;

//////////////////////// DY Powheg

double Ngen_dy_2_ee = 3297045; // EE
double Ngen_dy_2_mm = 3283850; // MM
double Xsec_dy_2 = 5745.25/3.0; // NNLO
double eXsec_dy_2 = 0;

//////////////////////// DY MadGraph 2

double Ngen_dy_3 = 7551580;
double Xsec_dy_3 = 3.0*20.3; // NLO
double eXsec_dy_3 = 0.15*Xsec_dy_3;

//////////////////////// DY Powheg 2

double Ngen_dy_4_ee = 100*10000; // EE
double Ngen_dy_4_mm = 100*10000; // MM
double Xsec_dy_4 = 333.866; // NNLO
double eXsec_dy_4 = 0;

//////////////////////// DY Powheg MiNLO

//double Ngen_dy_5_ee = 40*25000; // EE
//double Ngen_dy_5_mm = 40*25000; // MM
// weights extracted from h_gen_weights->GetBinContent(2)
double Ngen_dy_5_ee = 1.17808870400000000e+09; // EE
double Ngen_dy_5_mm = 1.17592755200000000e+09; // MM
double Xsec_dy_5 = 1156.41; // NNLO
double eXsec_dy_5 = 0;

//////////////////////// TTbar

double Ngen_tt = 6923750;
double Xsec_tt = 225.197; // NLO
double eXsec_tt = 0.07*Xsec_tt;

//////////////////////// ZZ

double Ngen_zz = 9799908;
double Xsec_zz = 8.059; // NLO with CTEQ
double eXsec_zz = 0.15*Xsec_zz; 

//////////////////////// WZ

double Ngen_wz = 10000283;
double Xsec_wz = 33.21;    // NLO with CTEQ
double eXsec_wz = 0.15*Xsec_wz; 

//////////////////////// QCD

double Ngen_qcd = 7529312;
double Xsec_qcd = 3.64E8; // search the NLO !
double eXsec_qcd = 0.15*Xsec_qcd; 

//////////////////////// WW

double Ngen_ww = 10000431;
double Xsec_ww = 54.838; // NLO with CTEQ
double eXsec_ww = 0.15*Xsec_ww; 

//////////////////////// Wj

double Ngen_wj = 57709905;
double Xsec_wj = 36703.2; // NNLO
double eXsec_wj = 0.15*Xsec_wj; 

//////////////////////// single top s-channel

double Ngen_tS = 139974;
double Xsec_tS = 1.76; // NNO
double eXsec_tS = 0.15*Xsec_tS; 

//////////////////////// single top t-channel

double Ngen_tT = 1935072;
double Xsec_tT = 30.7; // NNLO
double eXsec_tT = 0.15*Xsec_tT; 

//////////////////////// single top tW-channel

double Ngen_tW = 493460;
double Xsec_tW = 11.1; // NNLO
double eXsec_tW = 0.15*Xsec_tW; 

//////////////////////// single antitop s-channel

double Ngen_tSb = 259961;
double Xsec_tSb = 3.79; // NNO
double eXsec_tSb = 0.15*Xsec_tSb; 

//////////////////////// single antitop t-channel

double Ngen_tTb = 3758227;
double Xsec_tTb = 56.4; // NNLO
double eXsec_tTb = 0.15*Xsec_tTb; 

//////////////////////// single antitop tW-channel

double Ngen_tWb = 497658;
double Xsec_tWb = 11.1; // NNLO
double eXsec_tWb = 0.15*Xsec_tWb; 

///////////////////////

