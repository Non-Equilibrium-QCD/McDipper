/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */

#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_dht.h>



// ---------------------- INCLUDE PARAMETERS ---------------------------- //

// ---------------------- INCLUDE OBJECTS ---------------------------- //
#include "include/config.h"
#include "include/nucleus.h"
#include "include/event.h"
#include "include/pdfs_lhapdf.h"

// -------------------------- HERE STARTS THE MAIN -------------------------------- //


int main (int argc, char **argv) {

  if (argc!=3) {
    std::cerr << "FAILURE: No config file found!" << std::endl;
    std::cerr << argv[0] << " -i path/to/config.yaml" << std::endl;
    exit(EXIT_FAILURE) ;
  }
  
  std::string path_to_config(argv[2]);
  Config config(path_to_config);

  PDFs quark_dist(&config);

  int Nx=100; 
  double xmax= 0.99;
  double xmin= 1e-6;
  double dy=log(xmax/xmin)/(Nx-1.);

  double Qmin_p = 0.01;
  double Qmax_p = 100.0;
  int NQ_p = 1000;
  double dQ = std::log(Qmax_p/Qmin_p)/(NQ_p-1.);
  
  std::ofstream xuV_f;
  std::ofstream xdV_f;
  std::ostringstream xuVname;
  std::ostringstream xdVname;
  xuVname << "u_V_dist.txt"  ;
  xdVname << "d_V_dist.txt"  ;
  std::cout<< "Qmax=" <<quark_dist.get_QMax()<<std::endl;
  xuV_f.open(xuVname.str());
  for(int ix=0;ix<Nx;ix++){
    double x= xmin*exp(dy*ix);
    xuV_f<< x; 
    for (int iQ = 0; iQ < NQ_p; iQ++)
    {
      double Qi=Qmin_p*exp(iQ*dQ);
      xuV_f<<  "\t"  << quark_dist.xfpart(QuarkID::u,x, Qi)-quark_dist.xfpart(QuarkID::ubar,x, Qi);
    }
    xuV_f<< std::endl;
  }

  xuV_f.close();
  xdV_f.open(xdVname.str());
  for(int ix=0;ix<Nx;ix++){
    double x= xmin*exp(dy*ix);
    xdV_f<< x; 
    for (int iQ = 0; iQ < NQ_p; iQ++)
    {
       double Qi=Qmin_p*exp(iQ*dQ);
      xdV_f<<  "\t"  << quark_dist.xfpart(QuarkID::d,x, Qi)-quark_dist.xfpart(QuarkID::dbar,x, Qi);
    }
    xdV_f<< std::endl;
  }
  xdV_f.close();
  
  // Event EventGen(config);
  // EventGen.MakeEventByEvent();
  
  // LHAPDF::setVerbosity(1);
	return 0;
}
