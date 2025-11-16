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
#include "include/config_mv.h"
#include "include/nucleus.h"
#include "include/event.h"
#include "include/pdfs_lhapdf.h"
#include "include/dipole_mv.h"
// -------------------------- HERE STARTS THE MAIN -------------------------------- //

int main (int argc, char **argv) {

  if (argc!=3) {
    std::cerr << "FAILURE: No config file found!" << std::endl;
    std::cerr << argv[0] << " -i path/to/config.yaml" << std::endl;
    exit(EXIT_FAILURE) ;
  }
  
  std::string path_to_config(argv[2]);
  Config config(path_to_config);
  MVDipole mvdip(&config);
  
  int nr = 101;
  double rmin = 0.0001;
  double rmax = 10;
  double dq = log(rmax/rmin)/(nr-1.);

  double Yt = 2.8;
  double x=0.01*exp(-Yt);
  std::cout<< x << std::endl;
  for (size_t i = 0; i < nr; i++)
  {
    double r = rmin * exp(i*dq);
    // std::cout<< r << '\t' << mvdip.FundamentalDipole(x,r,0.5)<< '\t' << mvdip.FundamentalDipole(x,r,1.0)<< '\t' << mvdip.FundamentalDipole(x,r,2.0)<< '\t' << mvdip.FundamentalDipole(x,r,8) << std::endl;
    std::cout<< r << '\t' << mvdip.AdjointDipole(x,r,0.5)<< '\t' << mvdip.AdjointDipole(x,r,1.0)<< '\t' << mvdip.AdjointDipole(x,r,2.0)<< '\t' << mvdip.AdjointDipole(x,r,8) << std::endl;

  }
  
  
  // Event EventGen(config);
  // EventGen.MakeEventByEvent();

  // LHAPDF::setVerbosity(1);
	return 0;
}
