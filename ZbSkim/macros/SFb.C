#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>


// Calculation of the SFb by means of the POG prescription for 2012 data https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_payload_Moriond13.txt
// for the CSVT btagger in ranges of pT and |eta| < 2.4

void SFb(){

       	double SFb;
	int i,j;
	
	float ptmin[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
	float ptmax[] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
	
	float pT[]     = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};

	float SFb_error[] = {
		0.0515703,
		0.0264008,
		0.0272757,
		0.0275565,
		0.0248745,
		0.0218456,
		0.0253845,
		0.0239588,
		0.0271791,
		0.0273912,
		0.0379822,
		0.0411624,
		0.0786307,
		0.0866832,
		0.0942053,
		0.102403 };

	for(i=0; i<17; i++){
		float x = (ptmin[i]+ptmax[i])/2.;
		SFb = (0.927563+(1.55479e-05*x))+(-1.90666e-07*(x*x));
		if(i<16){
//			std::cout<< ptmin[i]<<"<"<<"pT"<<"<"<< ptmax[i] <<"         "<< "SFb = " << SFb <<" +- "<< SFb_error[i] << std::endl;
			std::cout<< ptmin[i]<<"    "<<ptmax[i]<<"      -2.4       2.4       "<<SFb<<"  "<<SFb_error[i]<<" "<<SFb_error[i]<<std::endl;
		}
	}
}

