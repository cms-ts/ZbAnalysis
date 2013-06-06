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
		0.0567059,
     		0.0266907,
    		0.0263491,
		0.0342831,
		0.0303327,
 		0.024608,
		0.0333786,
		0.0317642,
		0.031102,
		0.0295603,
		0.0474663,
		0.0503182,
		0.0580424,
		0.0575776,
		0.0769779,
		0.0898199 };

	for(i=0; i<17; i++){
		
		SFb =0.869965*((1.+(0.0335062*ptmax[i]))/(1.+(0.0304598*ptmax[i])));
		if(i<16){
			std::cout<< ptmin[i]<<"<"<<"pT"<<"<"<< ptmax[i] <<"         "<< "SFb = " << SFb <<" +- "<< SFb_error[i] << std::endl;
		}
	}
}
       	
