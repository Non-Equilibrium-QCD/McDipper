/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef PARSIPSAT_H
#define PARSIPSAT_H

namespace IPsat_pars{
  const int HankelTransMode=1; /// Chooses the way to perform the Hankel transform, the naive integration is preferred to the fast Hankel transform

  const int NY = 101;
  const int NR =101;
  const double x0=0.67;

  const double rmin=2.0e-4;
  const double rmin_tol=1.0e-4;
  const double rmax=3.0;
  const double du =log(rmax/rmin)/(NR-1.);
  const double umax=log(rmax/rmin);
  const int N_UMAX_INT=200;

  const double rsf=2.9;


  const double Ymin = 0.81;
  const double Ymax = 9.99;

  const double dYextra= 0.05;
  const double drextra= 2e-6;

  const double dR= 0.025;

  const double RMaxTol = 0.2;

  const double KTMAX=30;
  const double KTMIN=0.1;
  const double KTDIST=3;


  const double k_dip_min = 0.005; //In GeV!
  const double k_dip_max = 51; // In GeV!
  const int k_dip_points = 121;
  const double k_dip_dk=log(k_dip_max/k_dip_min)/(k_dip_points-1.);

  const double y_dip_min = 1.0;
  const double y_dip_max = 30.0;
  const int y_dip_points = 121;
  const double y_dip_dy=(y_dip_max-y_dip_min)/(y_dip_points-1.);

  const double T_dip_min = 0.;
  const double T_dip_max = 10.1;
  const int T_dip_points= 101;
  const double T_dip_dT = (T_dip_max-T_dip_min)/(T_dip_points-1.);

  const double x0_scaling=0.01;  // Sets in the geometrical scaling

}
#endif /* params_ipsat */
 