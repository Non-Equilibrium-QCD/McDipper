/* Copyright (c) 2023 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef GENPARS_H
#define GENPARS_H

#include <math.h>

namespace gen_pars{
  /* General parameters namespace. */
  const int NC = 3;
  const int dA = NC*NC-1;
  const double CA = double(NC);
  const double CF = double(dA)/(2.*double(NC));
  const int NF =3;

  const double qu=2./3.; const double qu2= qu*qu;
  const double qd=-1./3.; const double qd2= qd*qd;
  const double qs=-1./3.; const double qs2= qs*qs;

  const double qB = 1./3.;
  const double alpha_S = 0.3;

  const double hbarc= 0.1973;
  const double fm_to_GeVm1 =1/hbarc;
  const double fm2_to_GeVm2 =fm_to_GeVm1*fm_to_GeVm1;

  const double GeV_to_fmm1 =1/hbarc;
  const double GeV2_to_fmm2 =GeV_to_fmm1*GeV_to_fmm1;

  const double GeVm1_to_fm =hbarc;
  const double GeVm2_to_fm2 =GeVm1_to_fm*GeVm1_to_fm;

  const double fmm1_to_GeV =hbarc;
  const double fmm2_to_GeV2 =fmm1_to_GeV*fmm1_to_GeV;

  const double mb_to_fm2 = 0.1;
  const double fm2_to_mb = 1/mb_to_fm2;

  //Computation Parameters

  const double TMax=10;
  const double TMin=0.0;
  const int NT=101;
  const double dT = (TMax-TMin)/(NT-1);
  const double T_tolerance = dT/2.;
  const double Q_tolerance = 0.;

  const double epsabs= 0.0;
  const double epsrel= 0.001;
  const double limit=50000;
  const double PMIN=0.01;
  const double PMAX=50;
  const double NP = 1001;
  const double DP = (PMAX-PMIN)/(NP-1.);
  const int skip=2;

  const size_t routine = 6;

  //DipoleComputation
  const double pref_glue = dA*pow(2*M_PI,-4)/alpha_S/NC;
  const double XdipMax = 1.0;
  const double XdipMin = 0.0;

}
#endif /* gen_pars */
 