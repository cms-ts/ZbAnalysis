#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>


class run_lumi{
public:
	class record{
	public:
		double run;
		double lumi;
	record():run(0),lumi(0){
      	std::cout<<"record: Setting everything to zero !!!\n";
    	}
	record(double run_i, double lumi_i):
		run(run_i),lumi(lumi_i){}

	bool trovato(double run_number){
		return (run_number==run);
	}
	
	};

		

run_lumi(std::string filename){
	std::ifstream file(filename.c_str());
	double run_i, lumi_i;
	while(file>>run_i>>lumi_i)
	recd.push_back(record(run_i,lumi_i));
}

double luminosity(double run_number){
	for(unsigned int i=0; i!=recd.size(); i++){
		if((recd[i]).trovato(run_number)){
		return (recd[i].lumi);
		}
	}
	return 1;
}


private:
  run_lumi(){};
  std::vector<record> recd;

};

