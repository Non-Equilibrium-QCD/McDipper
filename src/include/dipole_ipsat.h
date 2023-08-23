/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef DIP_IPSAT_H
#define DIP_IPSAT_H
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "params_ipsat.h"
#include "Hankel.h"

enum class Rep : int { Fundamental = 0, Adjoint = 1};

class Dipole{
	/* This is the (IP-Sat) Dipole class. It reads in input from the */
	public:
	Dipole();
	Dipole(int set, bool Verbose);
	virtual ~Dipole();
	// --------------------------------------- EVALUATION OF PDFS --------------------------------------- //
	// G Functions -> rx-space
	double G(double x, double r);
	double GPrime(double x, double r);
	double GDoublePrime(double x, double r);

	// Dipole Functions rx-space
	// double D(double x, double r, double T){return exp(-T*G(x, r));}
	double DDoublePrime(double x, double r, double T);
	double DDoublePrimeAdjoint(double x, double r, double T);
	double AdjointDipole(double x, double r, double T);
	double FundamentalDipole(double x, double r, double T);
	//Moment functions
	double FundamentalDipole_k(double x, double k, double T);
	double AdjointDipole_k(double x, double k, double T);
	double MomentDipole(int n, double x, double T, Rep rep);
	double EffectiveRadius(int n, double x, double T, Rep rep);

	//// Transform Dipole functions. Set of function taking D(r)-> D(k);
	void Make_Dipole_Grids();
	void Make_Dipole_Interpolators();
	void Get_Momentum_Dipoles();
	void read_in_dipoles();
	void Make_Momentum_Dipoles();
	void Transform_Dipole(Rep dipole_rep);
	void Transform_Dipole_Naive(Rep dipole_rep);

	void set_new_xscaling(double xscaling_new){xscaling=xscaling_new;}
	double get_xscaling(){return xscaling;}
	void get_Q2();


	// Import Functions
	bool import_dipole(Rep dipole_rep);

	// Auxiliary Functions: Variables
	double get_Y(double x){return log(1/x);}
	double get_U(double r){return log(r/IPsat_pars::rmin);}
	double get_R(double u){return IPsat_pars::rmin*exp(u);}
	double get_X(double Y){return exp(-Y);}

	double get_XMIN(){return  X_MIN;}
	double get_XMAX(){return  X_MAX;}

	//output
	void write_config();
	void printProgress(double percentage1,double percentage2);
	void print_Dk_config();
	void dump_momentum_Dipole(double T);
	void dump_transformed_norm(double T);
	void Transform_Dipole_Naive_Test(Rep rep,double T_t);
	void make_test_output();



	private:

	bool DipVerbose;
	int DipSet;
	int HankelMode;

	int DipSet_conf;
	int HankelMode_conf;

	double aa;
	double xscaling;
	// gsl_interp_accel *accPDF[pdf_pars::NQQ];
	// gsl_spline *PDF[pdf_pars::NQQ];

	gsl_spline2d * G_spl ;
	gsl_spline2d * Gp_spl ;
	gsl_spline2d * Gpp_spl ;

	gsl_interp_accel * xacc ;
	gsl_interp_accel * yacc ;

	gsl_interp_accel * xaccP ;
	gsl_interp_accel * yaccP ;

	gsl_interp_accel * xaccPP ;
	gsl_interp_accel * yaccPP ;

	int N1,N2,N3;
	double DK,DY,DT;
	double K_MIN,Y_MIN,T_MIN;
	double K_MAX,Y_MAX,T_MAX;
	double X_MIN,X_MAX;


	double * k_grid_homo;
	double * q_grid_homo;
	double * Y_grid;
	double * T_grid;

	gsl_spline2d ** DF_spl ; //spline in 2d for (k,T)
	gsl_spline2d ** DA_spl ; //spline in 2d for (k,T), the Y direction  is made by

	gsl_interp_accel ** xaccA ;
	gsl_interp_accel ** yaccA ;

	gsl_interp_accel ** xaccF ;
	gsl_interp_accel ** yaccF ;

	gsl_spline2d * Q2F_spl ; //spline in 2d for (k,T)
	gsl_spline2d * Q2A_spl ; //spline in 2d for (k,T), the Y direction  is made by

	gsl_interp_accel * xaccQ2A ;
	gsl_interp_accel * yaccQ2A ;

	gsl_interp_accel * xaccQ2F ;
	gsl_interp_accel * yaccQ2F ;

	// Data Managing
	std::string path_to_tabs;
	std::string path_to_set;

	double MinFactor_conf;
	double MaxFactor_conf;

	double MinFactor;
	double MaxFactor;

	double kmin_conf,kmax_conf;
	double Ymin_conf,Ymax_conf;
	double Tmin_conf,Tmax_conf;
	int Nk_conf,NY_conf,NT_conf;
	double dk_conf,dY_conf,dT_conf;

	bool computed_in_run=false;
	bool is_initialized=false;


	// Tools//
	double smooth(double r, double r0){return (1+tanh((r-r0)/aa))/2.;}
	bool IsPathExist(const std::string &s);
	void check_for_tabs_folders();
	std::string get_set_name(int n);
	void get_name_and_value(std::string testline,std::string &name, std::string &value);
	void process_dip_parameters(std::string testline);
	void make_momentum_parameters(bool in_run);
	void quadratic_extrapolation(double x1,double f1,double x2,double f2,double x3,double f3,double &a,double &b,double &c);

	//checkers
	bool check_dipole_sets();
	bool check_config();
	void read_in_config();
	bool is_in_range(double x, double x1,double x2);



 };

 struct DipFT{
   double p;
   double x;
   double T;
   Dipole * dip;
 };

#endif /* dipole_ipsat */
 