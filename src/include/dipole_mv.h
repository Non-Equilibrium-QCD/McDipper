/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef DIP_MV_H
#define DIP_MV_H
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
#include "config_mv.h"
#include "config.h"

enum class MVLargeX : int { Frozen=0,Simple = 1, TMD = 2, PDFMatch = 3 };
enum class MVRep : int { Fundamental = 0, Adjoint = 1};

namespace MV_pars{
  const int HankelTransMode=1; /// Chooses the way to perform the Hankel transform, the naive integration is preferred to the fast Hankel transform
  const double xMax =0.8;
  const int large_X_modes= 20; 
}

class MVDipole{
	/* This is the (IP-Sat) Dipole class. It reads in input from the */
	public:
	MVDipole();
	// MVDipole(std::string path,  int LargeXExt, int NLX,double XMax, bool Verbose);
	MVDipole(Config * conf);
	void read_config_in();
	virtual ~MVDipole();
	// --------------------------------------- EVALUATION OF PDFS --------------------------------------- //
	// G Functions -> rx-space
	void initialize();
	void make_pos_dipole_interpolators();
	void initialize_position_space_dipole_interpolators();
	void initialize_Qs2_interpolators();
	
	///////////////////// POS-SPACE FUNCTIONS //////////////////////
	
	double AdjointDipole(double x, double r, double T);
	double FundamentalDipole(double x, double r, double T);
	double FundamentalQ2(double x, double T);
	double AdjointQ2(double x, double T);
	double FundamentaldQ2dY(double x, double T);
	double AdjointdQ2dY(double x, double T);

	///////////////////// MOM-SPACE FUNCTIONS //////////////////////
	double FundamentalDipole_k(double x, double k, double T);
	double AdjointDipole_k(double x, double k, double T);
	double AdjointDipole_k_thread_local(double x, double k, double T);
	double MomentDipole(int n, double x, double T, MVRep rep);
	double EffectiveRadius(int n, double x, double T, MVRep rep);

	////////////////////// TRANS. FUNCTIONS ///////////////////////
	void Make_Dipole_Interpolators();
	void Get_Momentum_Dipoles();
	void read_in_dipoles();
	void Make_Momentum_Dipoles();
	void Transform_Dipole();
	void Transform_Dipole_Naive();

	//////////////////// IMPORT FUNCTIONS ////////////////////////

		// Import Functions
	bool import_dipole();
	//////////////////////// CONFIG TOOLS ////////////////////////

	// Tools
	// void write_config();
	void printProgress(double percentage1,double percentage2);

	double getKmax(){return mvconf.getKmax();}
	double getKmin(){return mvconf.getKmin();}
	

	/////////////////////// RETRIEVING TOOLS ////////////////////////
	double getTmax(){return mvconf.getTmax();}
	double getTmin(){return mvconf.getTmin();}
	int getNT(){return mvconf.getNT();}
	double getdT(){return mvconf.getdT();}
	double T_at(int i){return mvconf.T_at(i);}

	double getYminLX(){return YminLX;}
	double getYmax(){return mvconf.getYmax();}
	double getYmin(){return mvconf.getYmin();}
	double getXmin(){return mvconf.getXmin();}
	double getXmax(){return xmax;}
	double getMaxEta(){return 0.5*log(xmax/mvconf.getXmin());}

	int indexY(double Y);


	/////////////////////// OUTPUT TOOLS ////////////////////////

	void test_output_fixed_T(double T);
	void dump_test_norm();

	private:
		std::string path;
		VerboseLevel DipVerbose;
		int HankelMode;
		MVLargeX LargeXType;

		MVConfig mvconf;
		
		int DipSet_conf;
		int HankelMode_conf;
		double * LargeXParameters;
		int NLargeXParameters;
		int NFixedParameters=4;

		int NLargeX;
		int NYTot;
		double xmax;
		double YminLX;
		double dYLX;

		// double aa;
		// double xscaling;
		// gsl_interp_accel *accPDF[pdf_pars::NQQ];
		// gsl_spline *PDF[pdf_pars::NQQ];

		gsl_spline2d ** SF_spl;
		gsl_spline2d ** SA_spl;

		gsl_interp_accel ** xaccSA ;
		gsl_interp_accel ** yaccSA ;

		gsl_interp_accel ** xaccSF ;
		gsl_interp_accel ** yaccSF ;


		// double * k_grid_homo;
		// double * q_grid_homo;
		// double * Y_grid;
		// double * T_grid;

		gsl_spline2d ** DF_spl ; //spline in 2d for (k,T)
		gsl_spline2d ** DA_spl ; //spline in 2d for (k,T), the Y direction  is made by

		gsl_interp_accel ** xaccDA ;
		gsl_interp_accel ** yaccDA ;

		gsl_interp_accel ** xaccDF ;
		gsl_interp_accel ** yaccDF ;

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


		bool computed_in_run=false;
		bool is_initialized=false;


		// Tools//
		// exptrapolation tools 
		double get_LargeX_Dipole(double x, double r, double T, MVRep rep);
		void get_LargeX_Type(int large_x_extrapolator);
		double get_Y_at(int i);

		bool IsPathExist(const std::string &s);
		void check_for_tabs_folders();
		std::string get_set_name(int n);
		// void get_name_and_value(std::string testline,std::string &name, std::string &value);
		// void process_config_parameters(std::string testline);
		// void make_momentum_parameters(bool in_run);
		void quadratic_extrapolation(double x1,double f1,double x2,double f2,double x3,double f3,double &a,double &b,double &c);

		//checkers
		bool check_dipole_sets();
		bool check_config();
		void read_in_config();
		bool is_in_range(double x, double x1,double x2);

		void set_largeX_values();

	 	double MVExponent(double r, double T,double Q2x, MVRep rep );

		

 };

 struct DipMVFT{
   double p;
   double x;
   double T;
   MVDipole * dip;
 };

#endif /* dipole_ipsat */
 


