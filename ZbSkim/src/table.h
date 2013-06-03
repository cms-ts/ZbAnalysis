/*
   A utility to read information from table indexed by four columns

   ptLow,ptHi,etaLow,etaHi

   Information called Val() is taken from fifth colum:

   effi

   The information called Err() is the average of the sixth and seventh columns

   errLow, errHi

*/

#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

class table
{
public:
  class record
  {
  public:
    double ptLow;
    double ptHi;
    double etaLow;
    double etaHi;
    double effi;
    double erlw;
    double erhi;
    record() :
      ptLow (0), ptHi (0), etaLow (0), etaHi (0), effi (0), erlw (0), erhi (0)
    {
      std::cout << "record: Setting everything to zero !!!\n";
    }
    record (double pt1, double pt2,
	    double eta1, double eta2,
	    double eff, double erl, double erh) :
      ptLow (pt1), ptHi (pt2), etaLow (eta1), etaHi (eta2), effi (eff), erlw (erl), erhi (erh)
    {
    }
    bool belongTo (double pt, double eta)
    {
      //std::cout<<ptHi<<"\t"<<pt<<"\t"<<ptLow<<"\t";
      //std::cout<<etaHi<<"\t"<<eta<<"\t"<<etaLow<<"\n";
      //std::cout<<((pt<ptHi&&pt>=ptLow)&&(eta<etaHi&&eta>=etaLow))<<std::endl;
      return (pt < ptHi && pt >= ptLow) && (eta < etaHi && eta >= etaLow);
    }
  };

public:
  table (std::string filename)
  {
    std::ifstream file (filename.c_str ());
    double pt1, pt2, eta1, eta2, effi, errLw, errHi;
    while (file >> pt1 >> pt2 >> eta1 >> eta2 >> effi >> errLw >> errHi)
      {
	recd.push_back (record (pt1, pt2, eta1, eta2, effi, errLw, errHi));
      }
  }

  double Err (double pt, double eta)
  {
    double hiPtBin = 0;
    for (unsigned int i = 0; i != recd.size (); i++)
      {
	if ((recd[i]).belongTo (pt, eta))
	  {
	    return 0.5 * (recd[i].erlw + recd[i].erhi);
	  }
	if ((recd[i]).belongTo (100, eta))
	  hiPtBin = 0.5 * (recd[i].erlw + recd[i].erhi);
      }
    return hiPtBin;
  }

  double Val (double pt, double eta)
  {
    double hiPtBin = 0;
    for (unsigned int i = 0; i != recd.size (); i++)
      {
	if ((recd[i]).belongTo (pt, eta))
	  {
	    return recd[i].effi;
	  }
	if ((recd[i]).belongTo (100, eta))
	  hiPtBin = recd[i].effi;
      }
    return hiPtBin;
  }

private:
  table ()
  {
  };
  std::vector < record > recd;
};
