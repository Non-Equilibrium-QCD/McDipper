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
  
  Event EventGen(config);
  EventGen.MakeEventByEvent();
  
  LHAPDF::setVerbosity(1);
	return 0;
}
