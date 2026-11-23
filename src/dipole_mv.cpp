

/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */
#include <algorithm>
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <math.h>
#include <sstream>
#include <sys/stat.h>
#include <filesystem>
#include <fstream>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_multifit.h>

namespace fs = std::filesystem;

// BESSEL FUNCTIONS //
#ifndef BesselJ0
#define BesselJ0(x) jn(0,x)
#endif

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include "include/dipole_mv.h"
#include "include/params_gen.h"

// To do: 
// [ ] Finish evaluation function
// [ ] Test Dipole and Q2 interpolators.
// [ ] Large X matching? 
// [ ] Adapt Fourier transform routine to MVDipole
// [ ] Interpolator MVDipole (momentum space)


// HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL HANKEL 

namespace HankelDipoleMV{

	double BesselTol=1e-6;
	double EpsRel=1e-5;
	// Absolute tolerance of each interval, as a fraction of the largest |running sum|. With epsabs = 0, qag tries to reach
	// EpsRel on far-tail intervals that contribute ~1e-15 of the result (spline-interpolated S ~ 1e-20) and hits round-off
	// or the subdivision limit, costing seconds per interval. 1e-11: changes results by < 3e-10 (Claude/mv_io_check/).
	double EpsAbsFrac=1e-11;
	int BessQuad =6;
	int GenQuad =6;

	static const int NumberOfIntegrationPoints=64000;

	// n-th positive zero (n >= 1) of J_order, order = 0 or 1
	double BesselZero(int order, int n){
		return (order==0) ? gsl_sf_bessel_zero_J0(n) : gsl_sf_bessel_zero_J1(n);
	}

	// Int_0^inf dr F(r), with F oscillating like J_order(k r): integrated interval by interval between
	// consecutive zeros of J_order(k r); stops when the last 100 intervals add less than BesselTol of the sum.
	double BesselIntegrate(double k, gsl_function F, int order){
			// SETUP GSL INTEGRATOR //
			gsl_integration_workspace *Workspace=gsl_integration_workspace_alloc(NumberOfIntegrationPoints);
			// CALCULATE INTEGRAL BY PARTITION //
			int max_steps=96000; double a=0.0, b=0.0; double Value=0.0,Error; double Sum=0.0; double PartialSum=0.0; double MaxSum=0.0;
			for(int n=0; n< max_steps; n++){
					if(n%100==0 && n>0){
							if(std::abs(PartialSum)<BesselTol*std::abs(Sum) || (PartialSum==0.0 && Sum==0.0)){
									gsl_integration_workspace_free (Workspace);
									return Sum;
							}
							PartialSum=0.0;
					}

					a = b;
					b = BesselZero(order,n+1);

					gsl_integration_qag(&F,a/k,b/k,EpsAbsFrac*MaxSum,EpsRel,NumberOfIntegrationPoints,BessQuad,Workspace,&Value,&Error);

					if(gsl_isnan(Value)){
							std::cerr << "#ERROR IN NUMERICAL INTEGRATION " << a << " " << b <<  " " << k << " " << Value <<  " " << gsl_isnan(Value) << std::endl;
							gsl_integration_workspace_free (Workspace);
							return NAN;
					}
					Sum+=Value;  PartialSum+=Value;
					MaxSum=std::max(MaxSum,std::abs(Sum));
			}

			std::cerr << "#WARNING SLOW CONVERGENCE FOR k=" << k << " VALUE IS " << Sum << " AND INCREMENTS IS " << Value << " AFTER " << max_steps << " STEPS"  << std::endl;

			gsl_integration_workspace_free (Workspace);
			return Sum;

	}

	double DipoleAdjFT_Integrand(double r,void * pars){
		DipMVFT * params= ((struct DipMVFT *) pars);
		return ( 2.0*M_PI*r*BesselJ0( (params->p)*r ) ) * (params->dip)->AdjointDipole(params->x,r,params->T);

	}

	double DipoleFunFT_Integrand(double r,void * pars){
		DipMVFT * params= ((struct DipMVFT *) pars);
		return ( 2.0*M_PI*r*BesselJ0( (params->p)*r ) ) * (params->dip)->FundamentalDipole(params->x,r,params->T);
	}

	double DipoleFunFT(DipMVFT * Variables){
		// CALCULATE qg TMDS //
		gsl_function GSL_FundDipole_Integrand;
		GSL_FundDipole_Integrand.function=&DipoleFunFT_Integrand;
		GSL_FundDipole_Integrand.params=Variables;
		double DipMom=BesselIntegrate(Variables->p,GSL_FundDipole_Integrand,0);
		return DipMom;
	}


	double DipoleAdjFT(DipMVFT * Variables){
		// CALCULATE qg TMDS //
		gsl_function GSL_AdjDipole_Integrand;
		GSL_AdjDipole_Integrand.function=&DipoleAdjFT_Integrand;
		GSL_AdjDipole_Integrand.params=Variables;
		double DipMom=BesselIntegrate(Variables->p,GSL_AdjDipole_Integrand,0);
		return DipMom;
	}

	// J1 form (integration by parts of the J0 form, using d/dx[x J1(x)] = x J0(x)):
	//   D(k) = (2 pi)^2 S_inf delta^2(k) - (2 pi / k) Int_0^inf dr r J1(k r) dS/dr
	// Only the k > 0 part is computed. The boundary term vanishes because S - S_inf -> 0 at large r,
	// so the constant S_inf (S(rmax) for a BK table) never enters the integral.
	double DipoleFunFT_J1_Integrand(double r,void * pars){
		DipMVFT * params= ((struct DipMVFT *) pars);
		return r*gsl_sf_bessel_J1( (params->p)*r ) * (params->dip)->FundamentalDipole_dr(params->x,r,params->T);
	}

	double DipoleAdjFT_J1_Integrand(double r,void * pars){
		DipMVFT * params= ((struct DipMVFT *) pars);
		return r*gsl_sf_bessel_J1( (params->p)*r ) * (params->dip)->AdjointDipole_dr(params->x,r,params->T);
	}

	double DipoleFunFT_J1(DipMVFT * Variables){
		gsl_function F;
		F.function=&DipoleFunFT_J1_Integrand;
		F.params=Variables;
		return -(2.0*M_PI/Variables->p) * BesselIntegrate(Variables->p,F,1);
	}

	double DipoleAdjFT_J1(DipMVFT * Variables){
		gsl_function F;
		F.function=&DipoleAdjFT_J1_Integrand;
		F.params=Variables;
		return -(2.0*M_PI/Variables->p) * BesselIntegrate(Variables->p,F,1);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





MVDipole::MVDipole(){}

// MVDipole::MVDipole(std::string PATH, int LargeXExt,int NLX,double XMax, bool Verbose){
// 	DipVerbose=Verbose;
// 	path= PATH;
// 	mvconf=MVConfig(path);
// 	if(DipVerbose){mvconf.terminal_setup_output();}
// 	initialize();
// 	get_LargeX_Type(LargeXExt);
// 	NLargeX=NLX;
// 	xmax=XMax;
// 	set_largeX_values();
// 	NTot= mvconf.getNY() + NLargeX;

// 	// HERE I HAVE TO INCLUDE ALL THE OTHER PARAMETERS

// 	Get_Momentum_Dipoles();
// }

MVDipole::MVDipole(Config * conf){
	DipVerbose=conf->get_Verbose();
	path= conf->get_modelPath();
	if(DipVerbose>VerboseLevel::None){std::cout<<"[MVDipole]: BK table: "<<path<<std::endl;}

	mvconf=MVConfig(path);
	if(mvconf.getAdjoint()!="Casimir" && mvconf.getAdjoint()!="Squared"){
		std::cerr<<"[MVDipole]: config.yaml in "<<path<<" must give BK-Parameters: Adjoint = Casimir or Squared (got '"
		         <<mvconf.getAdjoint()<<"')"<<std::endl;
		exit(EXIT_FAILURE);
	}
	{
		auto in = [](const std::string& s, std::initializer_list<const char*> L){
			for(auto v: L){ if(s==v) return true; } return false; };
		if(!in(mvconf.getIC(),{"MV-Logarithm","MV-Regulated"})
		 || !in(mvconf.getKtype(),{"LO","Daughter","Parent","Smallest","Balitsky"})
		 || !in(mvconf.getCoupling(),{"Smooth","Frozen"})){
			std::cerr<<"[MVDipole]: config.yaml in "<<path<<" has missing/unknown IC, Ktype or Coupling (got IC='"
			         <<mvconf.getIC()<<"', Ktype='"<<mvconf.getKtype()<<"', Coupling='"<<mvconf.getCoupling()<<"')"<<std::endl;
			exit(EXIT_FAILURE);
		}
	}
	if(DipVerbose>VerboseLevel::None){mvconf.terminal_setup_output();}
	initialize();
	get_LargeX_Type(conf->get_ModelParams(1));
	NLargeX=int(conf->get_ModelParams(2));
	xmax=conf->get_ModelParams(3);
	if(NLargeX>0 && std::fabs(mvconf.getYmin())>1e-9){
		// the large-x extrapolation is matched to the table at x = x0, i.e. Y = 0
		std::cerr<<"[MVDipole]: the table in "<<path<<" starts at Ymin = "<<mvconf.getYmin()
		         <<", but the large-x extrapolation (Npoints = "<<NLargeX<<") needs a table starting at x = x0 (Ymin = 0)."
		         <<" Use a table with Ymin = 0 or set Npoints: 0."<<std::endl;
		exit(EXIT_FAILURE);
	}
	set_largeX_values();
 
	NYTot= mvconf.getNY() + NLargeX;

	// We will get the mode dependent parameters.  (NFixedParameters=4)
	NLargeXParameters=conf->get_mpars_number()-NFixedParameters;
	LargeXParameters=new double[NLargeXParameters];
	for (size_t i = 0; i < NLargeXParameters; i++) {LargeXParameters[i]=conf->get_ModelParams(NFixedParameters+i);}
	if(mvconf.getIC()=="MV-Regulated" && LargeXType==MVLargeX::Simple && NLargeX>0){
		if(NLargeXParameters<1){
			std::cerr<<"[MVDipole]: MV-Regulated large x (LargeX: 1) needs beta as the first large-x parameter"<<std::endl;
			exit(EXIT_FAILURE);
		}
		initialize_MVReg_LargeX();
	}
	// std::cerr << LargeXParameters[0] <<std::endl;
	// Settings that enter the momentum table: stored in and compared with tabs/MV-Dipole/SetN/config.yaml
	mvconf.set_McDipper_settings(int(conf->get_ModelParams(1)), NLargeX, xmax,
		std::vector<double>(LargeXParameters, LargeXParameters+NLargeXParameters), MV_pars::HankelTransMode);
	Get_Momentum_Dipoles();
	
}


MVDipole::~MVDipole(){
  	for (int i = 0; i < mvconf.getNY(); i++) {
		gsl_spline2d_free(SA_spl[i]);
		gsl_spline2d_free(SF_spl[i]);
		gsl_interp_accel_free (xaccSA[i]);
		gsl_interp_accel_free (yaccSA[i]);
		gsl_interp_accel_free (xaccSF[i]);
		gsl_interp_accel_free (yaccSF[i]);
	}
	if(LamF_spl!=nullptr){gsl_spline2d_free(LamF_spl);}
  }

  void MVDipole::make_pos_dipole_interpolators(){
	SF_spl = new gsl_spline2d*[mvconf.getNY()] ;
	SA_spl = new gsl_spline2d*[mvconf.getNY()] ;

	xaccSA = new gsl_interp_accel*[mvconf.getNY()] ;
	yaccSA = new gsl_interp_accel*[mvconf.getNY()] ;

	xaccSF = new gsl_interp_accel*[mvconf.getNY()] ;
	yaccSF = new gsl_interp_accel*[mvconf.getNY()] ;
	if(DipVerbose>VerboseLevel::None){std::cout<<"[MVDipole]: Interpolator Grid created for Dipoles" <<std::endl;}

}


  void MVDipole::initialize(){
	make_pos_dipole_interpolators();
	//initialize SF(x,T,r) and SA(x,T,r)
	initialize_position_space_dipole_interpolators();
	// initialize QSA2(x,T) and QSF2(x,T,r)
	initialize_Qs2_interpolators();
	
  }

  // Checks that a grid value read from a BK-EV-NIP table matches the grid implied by the MV config.
  // The evaluators assume uniform T, Y (and log-r) grids built from config.yaml, so any mismatch scrambles the interpolation.
  static void check_grid_value(const char * table, const char * var, int idx, double read, double expected){
	const double relTol = 1e-5;
	const double absTol = 1e-9;
	if (std::fabs(read-expected) > absTol + relTol*std::max(std::fabs(read),std::fabs(expected))){
		std::cerr << "[MVDipole]: Grid mismatch in " << table << ": " << var << "[" << idx << "] = " << read
		          << " in file, but config.yaml gives " << expected << ". Table and config are inconsistent." << std::endl;
		exit(1);
	}
  }

  void MVDipole::initialize_position_space_dipole_interpolators(){
	// File written by BK-EV-NIP (BK::RunModel) with columns  T  Y  r[fm]  S_F  S_A
	// and loop order T (outer) -> Y -> r (inner).
	std::ostringstream Gname;
	Gname << path<< "/mv_dipoles.dat";
	FILE *table1 = fopen(Gname.str().c_str(),"r");
	if(table1==NULL){std::cerr<<"[MVDipole]: Could not open " << Gname.str() << std::endl; exit(1);}

	const int NY = mvconf.getNY();
	const int NT = mvconf.getNT();
	const int Nr = mvconf.getNr();

	// Buffers in file order: index (iT*NY + iY)*Nr + ir
	std::vector<double> SFbuf(static_cast<size_t>(NT)*NY*Nr);
	std::vector<double> SAbuf(static_cast<size_t>(NT)*NY*Nr);

	std::vector<double> T(NT), q(Nr);

	double y_t,r_t,SA_t,SF_t,T_t;

	for (int iT=0; iT < NT; iT++){
		for (int iY=0; iY < NY; iY++){
			for (int ir=0; ir < Nr; ir++){
				if (fscanf(table1,"%lf %lf %lf %lf %lf", &T_t, &y_t, &r_t, &SF_t, &SA_t) != 5){
					std::cerr<<"[MVDipole]: Error reading MV-Dipole table at (iT,iY,ir)=("<<iT<<","<<iY<<","<<ir<<")"<<std::endl;
					exit(1);
				}
				check_grid_value("mv_dipoles.dat","T",iT,T_t,mvconf.T_at(iT));
				check_grid_value("mv_dipoles.dat","Y",iY,y_t,mvconf.Y_at(iY));
				check_grid_value("mv_dipoles.dat","r",ir,r_t,mvconf.r_at(ir));

				size_t idx = (static_cast<size_t>(iT)*NY + iY)*Nr + ir;
				SFbuf[idx]=SF_t;
				SAbuf[idx]=SA_t;

				if(iY==0 && ir==0){T[iT]=T_t;}
				if(iT==0 && iY==0){q[ir]=mvconf.getq(r_t);}
			}
		}
	}
	// The table must end here: trailing data means the config grid is smaller than the table
	if (fscanf(table1,"%lf", &T_t) == 1){
		std::cerr<<"[MVDipole]: mv_dipoles.dat has more entries than NT*NY*Nr = "<<NT*NY*Nr<<" from config.yaml"<<std::endl;
		exit(1);
	}
	fclose(table1);

	// One (T, q) spline per Y slice. gsl_spline2d layout: z[j*xsize + i] with x=T (i), y=q (j)
	std::vector<double> GAarr(static_cast<size_t>(NT)*Nr);
	std::vector<double> GFarr(static_cast<size_t>(NT)*Nr);

	for (int iY=0; iY < NY; iY++){
		for (int iT=0; iT < NT; iT++){
			for (int ir=0; ir < Nr; ir++){
				size_t idx = (static_cast<size_t>(iT)*NY + iY)*Nr + ir;
				GAarr[NT*ir + iT]=SAbuf[idx];
				GFarr[NT*ir + iT]=SFbuf[idx];
			}
		}
		const gsl_interp2d_type * TA= gsl_interp2d_bicubic;
		SA_spl[iY] = gsl_spline2d_alloc(TA,NT, Nr);
		xaccSA[iY] = gsl_interp_accel_alloc();
		yaccSA[iY] = gsl_interp_accel_alloc();
		gsl_spline2d_init(SA_spl[iY],  T.data(), q.data(), GAarr.data(), NT, Nr);

		const gsl_interp2d_type * TF = gsl_interp2d_bicubic;
		SF_spl[iY] = gsl_spline2d_alloc(TF,NT, Nr);
		xaccSF[iY] = gsl_interp_accel_alloc();
		yaccSF[iY] = gsl_interp_accel_alloc();
		gsl_spline2d_init(SF_spl[iY], T.data(), q.data(), GFarr.data(), NT, Nr);
	}

	if(DipVerbose>VerboseLevel::None){std::cerr<<"[MVDipole]: Dipole interpolator succesfully created! "<<std::endl;}

  }

  void MVDipole::initialize_Qs2_interpolators(){
	// File written by BK-EV-NIP with columns  T  Y  Qs2_F[GeV^2]  Qs2_A[GeV^2]
	// and loop order T (outer) -> Y (inner).
	std::ostringstream Qname;
	Qname << path<< "/mv_Q2.dat";
	FILE *table1 = fopen(Qname.str().c_str(),"r");
	if(table1==NULL){std::cerr<<"[MVDipole]: Could not open " << Qname.str() << std::endl; exit(1);}

	const int NY = mvconf.getNY();
	const int NT = mvconf.getNT();

	// gsl_spline2d layout: z[j*xsize + i] with x=Y (i), y=T (j)
	std::vector<double> QA2arr(static_cast<size_t>(NY)*NT);
	std::vector<double> QF2arr(static_cast<size_t>(NY)*NT);

	std::vector<double> Y(NY), T(NT);

	double y_t,QA2_t,QF2_t,T_t;

	for (int iT=0; iT < NT; iT++){
		for (int iY=0; iY < NY; iY++){
			if (fscanf(table1,"%lf %lf %lf %lf", &T_t, &y_t, &QF2_t, &QA2_t) != 4){
				std::cerr<<"[MVDipole]: Error reading MV-Q2 table at (iT,iY)=("<<iT<<","<<iY<<")"<<std::endl;
				exit(1);
			}
			check_grid_value("mv_Q2.dat","T",iT,T_t,mvconf.T_at(iT));
			check_grid_value("mv_Q2.dat","Y",iY,y_t,mvconf.Y_at(iY));

			QA2arr[NY*iT + iY]=QA2_t;
			QF2arr[NY*iT + iY]=QF2_t;

			if(iY==0){T[iT]=T_t;}
			if(iT==0){Y[iY]=y_t;}
		}
	}
	if (fscanf(table1,"%lf", &T_t) == 1){
		std::cerr<<"[MVDipole]: mv_Q2.dat has more entries than NT*NY = "<<NT*NY<<" from config.yaml"<<std::endl;
		exit(1);
	}
	fclose(table1);

	const gsl_interp2d_type *TQA= gsl_interp2d_bicubic;
	Q2A_spl = gsl_spline2d_alloc(TQA, NY, NT);
	xaccQ2A = gsl_interp_accel_alloc();
	yaccQ2A = gsl_interp_accel_alloc();
	gsl_spline2d_init(Q2A_spl, Y.data(), T.data(), QA2arr.data(), NY, NT);

	const gsl_interp2d_type *TQF= gsl_interp2d_bicubic;
	Q2F_spl = gsl_spline2d_alloc(TQF, NY, NT);
	xaccQ2F = gsl_interp_accel_alloc();
	yaccQ2F = gsl_interp_accel_alloc();
	gsl_spline2d_init(Q2F_spl, Y.data(), T.data(), QF2arr.data(), NY, NT);

	if(DipVerbose>VerboseLevel::None){std::cerr<<"[MVDipole]: Qs2 interpolator succesfully created. "<<std::endl;}
}
//////////////////////////////////////////////////////////
//////       Position-Space Calling Methods       ////////
//////////////////////////////////////////////////////////

double MVDipole::FundamentalDipole(double x, double r, double T){
	double dipole=0.0;
	double Y = mvconf.getY(x);
	double q = mvconf.getq(r);
	if(Y<mvconf.getYmin()){
		//large x extrapolation
		dipole=get_LargeX_Dipole(x,r,T,MVRep::Fundamental);
	}
	else if(Y>mvconf.getYmax()){
		//small-x scaling
		dipole = 0.0;
	}
	else{
		if(r<mvconf.getRmin()){
		//small r extrapolation// For now, let's just freeze it 
			dipole=FundamentalDipole(x, mvconf.getRmin(),T);
		}
		else if(r>mvconf.getRmax()){
		//large-r extrapolation// For now, let's just freeze it 
			// dipole=FundamentalDipole(x, mvconf.getRmax(),T);
			// no freeze, set to 0 
			dipole =0.0;
		}
		else{
			// Bulk interpolation
			int jY;
			double Y1,Y2;
			double tmp1,tmp2,tmp;


			jY = static_cast<int>(std::floor( (Y-mvconf.getYmin())/mvconf.getdY() ));
			// std::cerr<<jY<<"\t"<< Y<<std::endl;	
			if(jY<mvconf.getNY()-1){
				Y1= mvconf.getYmin() + jY*mvconf.getdY();
				Y2=mvconf.getYmin() +  (jY+1)*mvconf.getdY();

				tmp1 = gsl_spline2d_eval(SF_spl[jY],T,q,NULL, NULL);
				tmp2 = gsl_spline2d_eval(SF_spl[jY+1],T,q,NULL, NULL);
				tmp= tmp1 + (tmp2-tmp1) * (Y-Y1)/(Y2-Y1) ;  // lin interpolation in Y
				dipole=tmp;
			}
			else if(jY==mvconf.getNY()-1){
				Y1= mvconf.getYmin() + jY*mvconf.getdY();
				dipole= gsl_spline2d_eval(SF_spl[jY],T,q,NULL, NULL);
			}

			// else()
			
			
		}

	}
	return dipole;
}


double MVDipole::AdjointDipole(double x, double r, double T){
	double dipole=0.0;
	double Y = mvconf.getY(x);
	double q = mvconf.getq(r);
	if(Y<mvconf.getYmin()){
		//large x extrapolation
		dipole=get_LargeX_Dipole(x,r,T,MVRep::Adjoint);
	}
	else if(Y>mvconf.getYmax()){
		//small-x scaling
		dipole = 0.0;
	}
	else{
		if(r<mvconf.getRmin()){
		//small r extrapolation// For now, let's just freeze it 
			dipole=AdjointDipole(x, mvconf.getRmin(),T);
		}
		else if(r>mvconf.getRmax()){
		//large-r extrapolation// For now, let's just freeze it 
			// dipole=AdjointDipole(x, mvconf.getRmax(),T);
		// no freeze, set to 0 
			dipole =0.0;
		}
		else{
			// Bulk interpolation
			int jY;
			double Y1,Y2;
			double tmp1,tmp2,tmp;

			jY = static_cast<int>(std::floor( (Y-mvconf.getYmin())/mvconf.getdY() ));

			if(jY<mvconf.getNY()-1){
				Y1 = mvconf.getYmin() + jY*mvconf.getdY();
				Y2 = mvconf.getYmin() +  (jY+1)*mvconf.getdY();

				tmp1 = gsl_spline2d_eval(SA_spl[jY],T,q,NULL, NULL);
				tmp2 = gsl_spline2d_eval(SA_spl[jY+1],T,q,NULL, NULL);
				tmp= tmp1 + (tmp2-tmp1) * (Y-Y1)/(Y2-Y1) ;  // lin interpolation in Y
				dipole=tmp;
			}
			else if(jY==mvconf.getNY()-1){
				Y1= mvconf.getYmin() + jY*mvconf.getdY();
				dipole= gsl_spline2d_eval(SA_spl[jY],T,q,NULL, NULL);
			}

			
			
		}

	}
	
	return dipole;
}


// dS/dr [fm^-1], used by the J1 Hankel transform.
// Same branches and Y-interpolation as FundamentalDipole/AdjointDipole. Outside [rmin, rmax] S is taken
// constant (frozen below rmin, S(rmax) above), so dS/dr = 0 there.
double MVDipole::FundamentalDipole_dr(double x, double r, double T){return Dipole_dr(x,r,T,MVRep::Fundamental);}
double MVDipole::AdjointDipole_dr(double x, double r, double T){return Dipole_dr(x,r,T,MVRep::Adjoint);}

double MVDipole::Dipole_dr(double x, double r, double T, MVRep rep){
	const double Y = mvconf.getY(x);
	const bool fund = (rep==MVRep::Fundamental);
	const bool largeX = (Y<mvconf.getYmin());   // as in FundamentalDipole/AdjointDipole: x = x0 is taken from the BK table
	if(largeX){return get_LargeX_Dipole_dr(x,r,T,rep);}
	if(Y>mvconf.getYmax()){return 0.0;}
	if(r<mvconf.getRmin() || r>mvconf.getRmax()){return 0.0;}

	gsl_spline2d     ** spl  = fund ? SF_spl : SA_spl;
	// no GSL accelerators (shared mutable state): the transform runs in parallel

	const double q  = mvconf.getq(r);
	const int    jY = static_cast<int>(std::floor( (Y-mvconf.getYmin())/mvconf.getdY() ));
	double dSdq;
	if(jY < mvconf.getNY()-1){
		const double Y1 = mvconf.Y_at(jY);
		const double Y2 = mvconf.Y_at(jY+1);
		const double d1 = gsl_spline2d_eval_deriv_y(spl[jY],  T,q,NULL,NULL);
		const double d2 = gsl_spline2d_eval_deriv_y(spl[jY+1],T,q,NULL,NULL);
		dSdq = d1 + (d2-d1)*(Y-Y1)/(Y2-Y1);  // lin interpolation in Y
	}
	else{
		dSdq = gsl_spline2d_eval_deriv_y(spl[jY],T,q,NULL,NULL);
	}
	return dSdq/r;  // q = ln(r/rmin) -> dq/dr = 1/r
}


// Evaluating Q2

double MVDipole::FundamentalQ2(double x, double T){
	double dipole=0.0;
	double Y = mvconf.getY(x);

	if(Y<mvconf.getYmin()){dipole = 0.0;}
	else if(Y>mvconf.getYmax()){dipole = 0.0;}
	else if(T<mvconf.getTmin()){dipole = 0.0;}
	else if(T>mvconf.getTmax()){dipole = 0.0;}
	else{dipole=gsl_spline2d_eval(Q2F_spl,Y,T,NULL, NULL);}
	return dipole;
}


double MVDipole::AdjointQ2(double x, double T){
	double dipole=0.0;
	double Y = mvconf.getY(x);

	if(Y<mvconf.getYmin()){dipole = 0.0;}
	else if(Y>mvconf.getYmax()){dipole = 0.0;}
	else if(T<mvconf.getTmin()){dipole = 0.0;}
	else if(T>mvconf.getTmax()){dipole = 0.0;}
	else{dipole=gsl_spline2d_eval(Q2A_spl,Y,T,NULL, NULL);}
	return dipole;
}
// rapidity derivatives of Q2
double MVDipole::FundamentaldQ2dY(double x, double T){
	double dipole=0.0;
	double Y = mvconf.getY(x);

	if(Y<mvconf.getYmin()){dipole = 0.0;}
	else if(Y>mvconf.getYmax()){dipole = 0.0;}
	else if(T<mvconf.getTmin()){dipole = 0.0;}
	else if(T>mvconf.getTmax()){dipole = 0.0;}
	else{dipole=gsl_spline2d_eval_deriv_x(Q2F_spl,Y,T,NULL, NULL);}
	return dipole;
}

double MVDipole::AdjointdQ2dY(double x, double T){
	double dipole=0.0;
	double Y = mvconf.getY(x);

	if(Y<mvconf.getYmin()){dipole = 0.0;}
	else if(Y>mvconf.getYmax()){dipole = 0.0;}
	else if(T<mvconf.getTmin()){dipole = 0.0;}
	else if(T>mvconf.getTmax()){dipole = 0.0;}
	else{dipole=gsl_spline2d_eval_deriv_x(Q2A_spl,Y,T,NULL, NULL);}
	return dipole;
}


double MVDipole::MVExponent(double r, double T, double Q2x,MVRep rep ){

	//NOTICE THAT Q2 has to be given in fm!!
	double Ci =0.0;
	if(rep==MVRep::Adjoint){Ci=mvconf.getAdjointPower();}   // S_A = S_F^Ci, the same relation as in the BK table
	else if(rep==MVRep::Fundamental){Ci=1.;}
	double Lr=log(mvconf.getEc()*M_E+ 1.0/(r*mvconf.getLambdaMV_in_fm()) );
	// MV-Logarithm form with Q^2(x) as the scale: (sigma0 T/4) (r^2 Q^2(x))^gamma L(r); equals the IC at x = x0
	double expo = 0.25*(mvconf.getSigma0()*T) * pow(Q2x*r*r, mvconf.getGamma()) * Lr;
	return Ci*expo;
}

double MVDipole::MVExponent_dr(double r, double T, double Q2x, MVRep rep ){
	// d/dr of MVExponent. With E = Ci (sigma0 T/4) (Q2 r^2)^gamma L and L = ln(ec e + 1/(r Lambda)):
	// dL/dr = -1/[r (1 + ec e r Lambda)]  =>  dE/dr = (E/r) [2 gamma - 1/(L (1 + ec e r Lambda))]
	const double ecE = mvconf.getEc()*M_E;
	const double rL  = r*mvconf.getLambdaMV_in_fm();
	const double L   = log(ecE + 1.0/rL);
	return (MVExponent(r,T,Q2x,rep)/r) * (2.0*mvconf.getGamma() - 1.0/(L*(1.0 + ecE*rL)));
}


double MVDipole::get_LargeX_Q2(double x, double T, MVRep rep){
	// Q^2(x) in fm^-2 used in the large-x extrapolation; returns -1 for modes not implemented yet.
	double Q2X=-1.0;
	switch (LargeXType)
	{
		case MVLargeX::Frozen:
			Q2X=mvconf.getQs02_in_fm();
			break;
		case MVLargeX::Simple:
			{				// lets start with the large-x part 
				// LargeXParameters[0] = beta
				// We can use everything from here onwards in GeV as only dimensionless radios contribute
				double Q2x0=0.0;
				double dQ2x0dY=0.0;
				if(rep==MVRep::Adjoint){
					Q2x0= MVDipole::AdjointQ2(mvconf.getX0(), T)*gen_pars::GeV2_to_fmm2;
					dQ2x0dY= MVDipole::AdjointdQ2dY(mvconf.getX0(), T)*gen_pars::GeV2_to_fmm2;

				}
				else if(rep==MVRep::Fundamental){
					Q2x0= MVDipole::FundamentalQ2(mvconf.getX0(), T)*gen_pars::GeV2_to_fmm2;
					dQ2x0dY= MVDipole::FundamentaldQ2dY(mvconf.getX0(), T)*gen_pars::GeV2_to_fmm2;
				}
				// The matching divides by Qs^2(x0, T); it is 0 where the BK table has no
				// saturation crossing at Y = 0 (Qs finder wrote 0), and lambda would be NaN.
				if(!(Q2x0 > 0.0)){
					std::cerr<<"[MVDipole]: large-x Simple mode needs Qs^2(x0, T) > 0, but the BK table gives "<<Q2x0
					         <<" fm^-2 at T = "<<T<<" fm^-2 (no saturation crossing at Y = 0).\n"
					         <<"             Use a T range where the dipole saturates at x0, or set Npoints: 0."<<std::endl;
					exit(EXIT_FAILURE);
				}
				double qbarx0 =sqrt(Q2x0) /(sqrt(2.)*mvconf.getLambdaMV_in_fm());
				double Ls = log(qbarx0 + M_E*mvconf.getEc());
				// now let's get the matching at x0
				double lambda = (dQ2x0dY/Q2x0) *( 1- (1./(2*Ls*mvconf.getGamma()))*(qbarx0/(qbarx0 + M_E*mvconf.getEc())) );
				lambda -= LargeXParameters[0] * mvconf.getX0()/(1-mvconf.getX0());

				Q2X=mvconf.getQs02_in_fm() *pow( (1-x)/(1-mvconf.getX0()),LargeXParameters[0] ) * pow(mvconf.getX0()/x, lambda);
				// all these formulas are defined in the extended documentation.
			}
			break;
		case MVLargeX::TMD:
			break;
		case MVLargeX::PDFMatch:
			break;
		default:
			break;
	}
	return Q2X;
}

double MVDipole::get_LargeX_Dipole(double x, double r, double T, MVRep rep){
	// MV-Regulated: Simple mode in MVReg_LargeX_Exponent. Other modes are not implemented: no scattering for x > x0, S = 1
	if(mvconf.getIC()=="MV-Regulated"){
		if(LamF_spl==nullptr || r<=0.0){return 1.0;}
		return std::exp(-MVReg_LargeX_Exponent(x,r,T,rep,NULL));
	}
	const double Q2X = get_LargeX_Q2(x,T,rep);
	if(Q2X<0){return 0.0;}  // TMD, PDFMatch: not implemented
	return exp(- MVExponent(r, T, Q2X, rep ));
}

double MVDipole::get_LargeX_Dipole_dr(double x, double r, double T, MVRep rep){
	if(mvconf.getIC()=="MV-Regulated"){
		if(LamF_spl==nullptr || r<=0.0){return 0.0;}   // S = 1  =>  D(k>0) = 0
		double dEdr;
		const double E = MVReg_LargeX_Exponent(x,r,T,rep,&dEdr);
		return -std::exp(-E)*dEdr;
	}
	// S = exp(-E)  =>  dS/dr = -S dE/dr
	const double Q2X = get_LargeX_Q2(x,T,rep);
	if(Q2X<0){return 0.0;}
	return -exp(- MVExponent(r, T, Q2X, rep )) * MVExponent_dr(r, T, Q2X, rep);
}


// MV-Regulated large x (x > x0), Simple mode:
//   E(x,r,T) = -ln S = Ci R(x,r,T) E_IC(r,T),   E_IC = (2 pi kappa T/m^2) [1 - m r K1(m r)],
//   R = ((1-x)/(1-x0))^beta (x0/x)^lambda(r,T),  lambda = d ln(-ln S_F)/dY at Y = 0 (BK table) - beta x0/(1-x0).
// E and dE/dY are then continuous at x0 for every r; D(k) -> 0 smoothly as x -> 1. beta = LargeXParameters[0] (as MV-Log).
// S_A = S_F^Ci in the table, so one lambda serves both representations. lambda is frozen outside [rmin, rmax].
void MVDipole::initialize_MVReg_LargeX(){
	const int NT = mvconf.getNT();
	const int Nr = mvconf.getNr();
	if(mvconf.getNY()<3){
		std::cerr<<"[MVDipole]: MV-Regulated large x needs at least 3 Y rows in the BK table (NY = "<<mvconf.getNY()<<")"<<std::endl;
		exit(EXIT_FAILURE);
	}
	const double dY = mvconf.getdY();
	std::vector<double> T(NT), q(Nr), lam(static_cast<size_t>(NT)*Nr);
	for (int iT=0; iT<NT; iT++){T[iT]=mvconf.T_at(iT);}
	for (int ir=0; ir<Nr; ir++){q[ir]=mvconf.getq(mvconf.r_at(ir));}
	int nBad=0;
	for (int iT=0; iT<NT; iT++){
		for (int ir=0; ir<Nr; ir++){
			// 3-point forward difference of ln E in Y on the rows Y = 0, dY, 2 dY
			double lnE[3]; bool ok=true;
			for (int j=0; j<3; j++){
				const double S = gsl_spline2d_eval(SF_spl[j],T[iT],q[ir],NULL,NULL);
				if(!(S>0.0 && S<1.0)){ok=false; break;}
				lnE[j] = std::log(-std::log(S));
			}
			double l = 0.0;   // lambda = 0 (no x-dependence beyond the (1-x)^beta factor) where ln(-ln S) is undefined
			if(ok){l = (-3.0*lnE[0] + 4.0*lnE[1] - lnE[2])/(2.0*dY);}
			else{nBad++;}
			lam[static_cast<size_t>(NT)*ir + iT] = l;   // gsl_spline2d layout z[j*xsize + i], x = T (i), y = q (j)
		}
	}
	if(nBad>0){std::cerr<<"[MVDipole]: MV-Regulated large x: S_F not in (0,1) at "<<nBad<<" (T,r) points, lambda set to 0 there"<<std::endl;}
	LamF_spl = gsl_spline2d_alloc(gsl_interp2d_bicubic, NT, Nr);
	gsl_spline2d_init(LamF_spl, T.data(), q.data(), lam.data(), NT, Nr);
}

double MVDipole::MVReg_LargeX_Exponent(double x, double r, double T, MVRep rep, double * dEdr){
	// E = -ln S for x > x0 and, if dEdr != NULL, dE/dr [fm^-1]. Needs r > 0.
	const double Ci   = (rep==MVRep::Adjoint) ? mvconf.getAdjointPower() : 1.0;
	const double m    = mvconf.getMIC()*gen_pars::GeV_to_fmm1;
	const double A    = 2.0*M_PI*mvconf.getKappa()*T/(m*m);
	const double u    = m*r;
	const double beta = LargeXParameters[0];
	const double x0   = mvconf.getX0();
	const double L    = std::log(x0/x);   // < 0 for x > x0
	const double rc   = std::min(std::max(r,mvconf.getRmin()),mvconf.getRmax());
	const double qc   = mvconf.getq(rc);
	const double lam  = gsl_spline2d_eval(LamF_spl,T,qc,NULL,NULL) - beta*x0/(1.0-x0);
	const double R    = std::pow((1.0-x)/(1.0-x0),beta) * std::exp(lam*L);
	const double f    = (u<700.0) ? 1.0 - u*gsl_sf_bessel_K1(u) : 1.0;   // K1 underflows beyond u ~ 700
	const double E    = Ci*R*A*f;
	if(dEdr!=NULL){
		const double dfdr   = (u<700.0) ? m*u*gsl_sf_bessel_K0(u) : 0.0;   // d/du [u K1(u)] = -u K0(u)
		const double dlamdr = (r>mvconf.getRmin() && r<mvconf.getRmax()) ? gsl_spline2d_eval_deriv_y(LamF_spl,T,qc,NULL,NULL)/r : 0.0;
		*dEdr = Ci*R*A*dfdr + E*L*dlamdr;   // dR/dr = R ln(x0/x) dlambda/dr
	}
	return E;
}

void MVDipole::get_LargeX_Type(int large_x_extrapolator){
	switch (large_x_extrapolator)
	{
	case 0:
		LargeXType = MVLargeX::Frozen;
		break;
	case 1:
		LargeXType = MVLargeX::Simple;
		break;
	case 2:
		LargeXType = MVLargeX::TMD;
		break;
	case 3:
		LargeXType = MVLargeX::PDFMatch;
		break;
	default:
		break;
	}
}

//////////////////////////////////////////////////////////
//////            Construction Methods            ////////
//////////////////////////////////////////////////////////

/// JUST COPY THIS FROM IPSAT!!!

//////////////////////////////////////////////////////////
//////            Transformation Methods          ////////
//////////////////////////////////////////////////////////


void MVDipole::Make_Dipole_Interpolators(){
	DF_spl = new gsl_spline2d*[ NYTot ] ;
	DA_spl = new gsl_spline2d*[ NYTot ] ;


	if(DipVerbose>VerboseLevel::None){std::cout<<"\n[MVDipole]: Interpolator Grid created for Dipoles" <<std::endl;}

}

void MVDipole::Get_Momentum_Dipoles(){
	//Check wether it is pre-computed
	check_for_tabs_folders();
	computed_in_run = !(check_dipole_sets());
	// If not computed already -> compute new momentum space dipoles
	if(computed_in_run){
		Make_Momentum_Dipoles();
	}
	// Output config
	// print_Dk_config();  This NEEDS TO BE DONE TOO. , but can be directly performed in the mvconfig files
	// Read in dipole tabulations and interpolate
	read_in_dipoles();
	// Make Output for selected dipoles
	// make_test_output(); laterzzzz
}

void MVDipole::read_in_dipoles(){
	Make_Dipole_Interpolators();
	bool is_imported = import_dipole() ;
	if(!is_imported){std::cerr<<"[MVDipole]: Error while importing dipole!"<<std::endl;exit(EXIT_FAILURE);}
}

void MVDipole::Make_Momentum_Dipoles(){
  	fs::create_directories(path_to_set);
	MinFactor=HankelParameters::MinFactor;
	MaxFactor=HankelParameters::MinFactor;
	// Transform
	// Debug output of S(r) at fixed T into the set folder (off: not needed for the tables,
	// and it reads the splines outside the T range for tables that do not cover these T).
	// test_output_fixed_T(0.5);
	// test_output_fixed_T(1);
	// test_output_fixed_T(4);
	// test_output_fixed_T(8);
	if(MV_pars::HankelTransMode==0){
		if(DipVerbose>VerboseLevel::None){std::cout<<"[MVDipole]: Transforming IP-Sat Dipoles using the Fast Hankel Transform Method" <<std::endl;}
		Transform_Dipole();
	}
	else if(MV_pars::HankelTransMode==1){
		if(DipVerbose>VerboseLevel::None){std::cerr<<"[MVDipole]: Transforming MV Dipoles: J0 integration of S between Bessel zeros" <<std::endl;}
		Transform_Dipole_Naive();
	}
	else if(MV_pars::HankelTransMode==2){
		if(DipVerbose>VerboseLevel::None){std::cerr<<"[MVDipole]: Transforming MV Dipoles: J1 integration of dS/dr between Bessel zeros" <<std::endl;}
		Transform_Dipole_Naive();
	}
	// Written only after the transform has finished: a set without config.yaml is incomplete and skipped
	mvconf.write_config(path_to_set);

}

void MVDipole::Transform_Dipole(){
	// double DipF[N1];

	double T_t,x_t;
	double Y_t;
	double DipA_k=0;
	double DipF_k=0;
	double RDomA,RDomF;


	std::ofstream dip_f;
	std::ostringstream path_to_dip;
	path_to_dip << path_to_set<< "/MomentumDipoles.dat" ;

	if(DipVerbose>VerboseLevel::None){
		std::cout<< "[MVDipole]: Transforming Dipoles to momentum space"<< std::endl;
	}
	if(DipVerbose==VerboseLevel::Dynamic){
		std::cout<< "Total (Y,T,k)                                                          Local (T,k)"<< std::endl; 
	}
	dip_f.open(path_to_dip.str());

	for (int iy = - NLargeX; iy < mvconf.getNY(); iy++){
		Y_t = get_Y_at(iy);
		x_t = mvconf.getx(Y_t);
			
		for (int iT = 0; iT < mvconf.getNT(); iT++) {

			T_t=mvconf.T_at(iT);

			Hankel HankA(0,false);
			Hankel HankF(0,false);
			// Set Hankel Environment around the characteristic scale of the the dipole
			RDomA= sqrt(2./AdjointQ2(x_t,T_t));
			RDomF= sqrt(2./FundamentalQ2(x_t,T_t));
			HankA.Initialize_Domain_w_CharScale(RDomA);
			HankF.Initialize_Domain_w_CharScale(RDomF);

			if(iy==0 && iT==0){
				if(DipVerbose>VerboseLevel::None){ 
					std::cout<< "[MVDipole]: Transforming for fundamental function using "<< HankF.get_Npoints() << " points "<<std::endl; 
					std::cout<< "[MVDipole]: Transforming for adjoint function using "<< HankA.get_Npoints() << " points "<<std::endl; 
				}
			}
			//Set points
			for (int ix = 0; ix < HankF.get_Npoints(); ix++) {HankF.set_FX(ix,FundamentalDipole(x_t,HankF.get_X(ix),T_t)) ;}
			for (int ix = 0; ix < HankA.get_Npoints(); ix++) {HankA.set_FX(ix,AdjointDipole(x_t,HankA.get_X(ix),T_t)) ;}
			// Transform
			HankA.Transform();
			HankF.Transform();
			//retrieval

			double kA_int_grid[HankA.get_Npoints()];
			double kF_int_grid[HankF.get_Npoints()];
			double DAk_int_grid[HankA.get_Npoints()];
			double DFk_int_grid[HankF.get_Npoints()];

			// Important to get the right UNITS!!
			// k needs conversion fm^{-1}-> GeV
			// Dk needs conversion fm^{2}-> GeV^{-2}
			for (int ik = 0; ik < HankA.get_Npoints(); ik++) {
				kA_int_grid[ik]=HankA.get_K(ik)*gen_pars::fmm1_to_GeV;
				DAk_int_grid[ik]=2*M_PI*HankA.get_FK(ik)*gen_pars::fm2_to_GeVm2;
			}
			for (int ik = 0; ik < HankF.get_Npoints(); ik++) {
				kF_int_grid[ik]=HankF.get_K(ik)*gen_pars::fmm1_to_GeV;
				DFk_int_grid[ik]=2*M_PI*HankF.get_FK(ik)*gen_pars::fm2_to_GeVm2;
			}

			// interpolate 1D
			gsl_interp_accel *accA_t= gsl_interp_accel_alloc ();
			gsl_spline *splineA_t = gsl_spline_alloc (gsl_interp_cspline, HankA.get_Npoints());

			gsl_interp_accel *accF_t= gsl_interp_accel_alloc ();
			gsl_spline *splineF_t = gsl_spline_alloc (gsl_interp_cspline, HankF.get_Npoints());

			gsl_spline_init (splineA_t, kA_int_grid, DAk_int_grid, HankA.get_Npoints());
			gsl_spline_init (splineF_t, kF_int_grid, DFk_int_grid, HankF.get_Npoints());


			// Here extrapolation goes on
			for (int iu = 0; iu < mvconf.getNK(); iu++) {
				double k_t = mvconf.k_at(iu);
				double u_t = mvconf.u_at(iu);

				//Adjoint 
				if(k_t>mvconf.k_at(HankA.get_Npoints()-1)){ 
					DipA_k =0;
				}
				else if(k_t<kA_int_grid[0]){
					double k1,k2,k3;
					double a1,a2,a3;//adjoint
					double a,b,c;

					k1 = kA_int_grid[0]; k2 = kA_int_grid[1]; k3 = kA_int_grid[2];
					a1 = DAk_int_grid[0];a2 = DAk_int_grid[1];a3 = DAk_int_grid[2];
					quadratic_extrapolation(k1,a1,k2,a2,k3,a3,a,b,c);
					DipA_k =a*k_t*k_t + b*k_t+ c;
				}
				else{DipA_k =gsl_spline_eval (splineA_t,k_t, accA_t);}
				
				//Fundamental
				if(k_t>mvconf.k_at(HankA.get_Npoints()-1)){ 
					DipF_k =0;
				}
				else if(k_t<kA_int_grid[0]){
					double k1,k2,k3;
					double f1,f2,f3;//fundamenta
					double a,b,c;

					k1 = kF_int_grid[0]; k2 = kF_int_grid[1]; k3 = kF_int_grid[2];
					f1 = DFk_int_grid[0];f2 = DFk_int_grid[1];f3 = DFk_int_grid[2];
					quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a,b,c);
					DipF_k =a*k_t*k_t + b*k_t+ c;
				}
				else{DipF_k =gsl_spline_eval (splineF_t,k_t, accF_t);}



				dip_f<< T_t<<"\t"<< Y_t<<"\t"<< u_t<<"\t" << DipF_k <<"\t" << DipA_k << std::endl;

				if(DipVerbose==VerboseLevel::Dynamic){
					if(remainder(iy,gen_pars::skip)==0){
						double percentage_done2 = double(iy)/double(IPsat_pars::y_dip_points);
						double percentage_done1 = double(iT)/double(IPsat_pars::T_dip_points);
						printProgress(percentage_done2,percentage_done1);}
				}
			}
			gsl_spline_free (splineF_t);
			gsl_interp_accel_free (accF_t);			
			gsl_spline_free (splineA_t);
			gsl_interp_accel_free (accA_t);
		}
	}
	dip_f.close();
	if(DipVerbose==VerboseLevel::Dynamic){
		printProgress(1,1);
		std::cout<< std::endl;
	}

}

void MVDipole::Transform_Dipole_Naive(){
	std::ofstream dip_f;
	std::ostringstream path_to_dip;
	path_to_dip << path_to_set<< "/MomentumDipoles.dat" ;

	if(DipVerbose>VerboseLevel::None){
		std::cout<< "[MVDipole]: Transforming Dipoles to momentum space"<< std::endl;
	}
	dip_f.open(path_to_dip.str());
	dip_f << std::setprecision(12);   // was 6 digits; same as the BK-EV-NIP tables

	const int NT = mvconf.getNT();
	const int NK = mvconf.getNK();
	std::vector<double> DF(long(NT)*NK), DA(long(NT)*NK);

	for (int iy = - NLargeX; iy < mvconf.getNY(); iy++){
		const double Y_t = get_Y_at(iy);
		const double x_t = mvconf.getx(Y_t);

		// Every (T,k) point of one Y row is an independent integral over a table
		// that is only read, so the row parallelises directly. Results go to arrays
		// and are written afterwards in the original order, so the file does not
		// depend on the number of threads.
		#pragma omp parallel for collapse(2) schedule(dynamic)
		for (int iT = 0; iT < NT; iT++) {
			for (int ik = 0; ik < NK; ik++) {
				const double T_t = mvconf.T_at(iT);
				const double k_t = mvconf.k_at(ik);
				double DipF_k = 0, DipA_k = 0;
				if(T_t!=0){
					DipMVFT HTPars={k_t*gen_pars::GeV_to_fmm1,x_t,T_t,this};
					if(MV_pars::HankelTransMode==2){
						DipF_k = HankelDipoleMV::DipoleFunFT_J1(&HTPars)*gen_pars::fm2_to_GeVm2;
						DipA_k = HankelDipoleMV::DipoleAdjFT_J1(&HTPars)*gen_pars::fm2_to_GeVm2;
					}
					else{
						DipF_k = HankelDipoleMV::DipoleFunFT(&HTPars)*gen_pars::fm2_to_GeVm2;
						DipA_k = HankelDipoleMV::DipoleAdjFT(&HTPars)*gen_pars::fm2_to_GeVm2;
					}
				}
				DF[long(iT)*NK+ik] = DipF_k;
				DA[long(iT)*NK+ik] = DipA_k;
			}
		}

		for (int iT = 0; iT < NT; iT++) {
			for (int ik = 0; ik < NK; ik++) {
				dip_f<< Y_t<<"\t"<<mvconf.T_at(iT)<<"\t"<< mvconf.u_at(ik)<<"\t" << DF[long(iT)*NK+ik] <<"\t" << DA[long(iT)*NK+ik] << std::endl;
			}
		}
		if(DipVerbose==VerboseLevel::Dynamic){
			printProgress(double(iy+NLargeX+1)/double(NYTot), 1.0);
		}
	}
	dip_f.close();
	if(DipVerbose==VerboseLevel::Dynamic){
		std::cout<< std::endl;
	}
}


//////////////////////////////////////////////////////////
//////         Transformation Aide-Methods        ////////
//////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////
//////              Importing Methods             ////////
//////////////////////////////////////////////////////////

bool MVDipole::import_dipole(){
	bool is_read=false;
	std::ostringstream tablename;

	tablename << path_to_set <<"/MomentumDipoles.dat";

	if(DipVerbose>VerboseLevel::None){std::cout<<"[MVDipole]: Importing momentum-space dipoles. SRCFILE=" << tablename.str()<< std::endl ;}

	double DipAArray[mvconf.getNK()*mvconf.getNT()];
	double DipFArray[mvconf.getNK()*mvconf.getNT()];

	double Tarray[mvconf.getNT()];
	double Uarray[mvconf.getNK()];
	// for (int iT = 0; iT < mvconf.getNT(); iT++) {Tarray[iT]= mvconf.T_at(iT);}
	// for (int iU = 0; iU < mvconf.getNK(); iU++) {Uarray[iU]= mvconf.u_at(iU);}

	// THIS NEEDS 
	FILE *table = fopen(tablename.str().c_str(),"r");
	if(table==NULL){std::cerr<<"[MVDipole]: Could not open " << tablename.str() << std::endl; exit(1);}
	double T_t,y_t,u_t,DipA_t,DipF_t;
	// grid values in the file have 6 significant digits
	auto close = [](double a, double b){ return std::fabs(a-b) <= 1e-6 || std::fabs(a-b) <= 1e-5*std::max(std::fabs(a),std::fabs(b)); };
	auto stop = [&](const std::string& msg){
		std::cerr<<"[MVDipole]: "<<tablename.str()<<": "<<msg<<"\n             The set does not match its config; delete "<<path_to_set<<" and rerun."<<std::endl;
		exit(EXIT_FAILURE); };
	long line=0;
	const bool useBicubic = (mvconf.getNK()>=4 && mvconf.getNT()>=4);
	if(!useBicubic){std::cerr<<"[MVDipole]: WARNING: NT or NK < 4, momentum dipoles interpolated bilinearly in (u, T)"<<std::endl;}
	for (int iY = 0; iY < NYTot; iY++) {
		for (int iT = 0; iT < mvconf.getNT(); iT++) {		
			for (int iU = 0; iU < mvconf.getNK(); iU++) {
				line++;
				if (fscanf(table,"%lf %lf %lf %lf %lf",  &y_t, &T_t, &u_t, &DipF_t, &DipA_t) != 5){stop("cannot read line "+std::to_string(line)+" (file shorter than the grid?)");}
				if(!close(y_t,get_Y_at(iY-NLargeX)) || !close(T_t,mvconf.T_at(iT)) || !close(u_t,mvconf.u_at(iU))){
					std::ostringstream m; m<<"grid mismatch at line "<<line<<": (Y,T,u) = ("<<y_t<<", "<<T_t<<", "<<u_t<<"), expected ("
					  <<get_Y_at(iY-NLargeX)<<", "<<mvconf.T_at(iT)<<", "<<mvconf.u_at(iU)<<")";
					stop(m.str());
				}
				if(!std::isfinite(DipF_t) || !std::isfinite(DipA_t)){stop("non-finite dipole at line "+std::to_string(line));}
				DipFArray[mvconf.getNK()*iT + iU]=DipF_t;
				DipAArray[mvconf.getNK()*iT + iU]=DipA_t;
				if(iU==0 && iY ==0 ){Tarray[iT]=T_t;}
				if(iY==0 && iT ==0 ){Uarray[iU]=u_t;}
			}
		}
		
		// bicubic in (u, T); GSL needs >= 4 points per axis, otherwise bilinear
		const gsl_interp2d_type *TF= useBicubic ? gsl_interp2d_bicubic : gsl_interp2d_bilinear;
	    DF_spl[iY] = gsl_spline2d_alloc(TF,mvconf.getNK(),mvconf.getNT());
	    // no GSL accelerators: they are mutable state shared between threads (gluon table is OpenMP-parallel)
	    gsl_spline2d_init(DF_spl[iY], Uarray, Tarray, DipFArray, mvconf.getNK(),mvconf.getNT());
		
		const gsl_interp2d_type *TA= useBicubic ? gsl_interp2d_bicubic : gsl_interp2d_bilinear;
	    DA_spl[iY] = gsl_spline2d_alloc(TA,mvconf.getNK(),mvconf.getNT());
	    gsl_spline2d_init(DA_spl[iY], Uarray, Tarray, DipAArray,mvconf.getNK(),mvconf.getNT());
	                                                                                                                                                                                                                                      }
	{ double extra; if(fscanf(table,"%lf",&extra)==1){stop("more lines than the grid ("+std::to_string(line)+" expected)");} }
	fclose(table);
 	is_read=true;
	if(DipVerbose>VerboseLevel::None){std::cout<<"[MVDipole]: Momentum-space dipole ---> Imported!" << std::endl;}
  	return is_read;
}

// Not sure if this will be needed this time though

//////////////////////////////////////////////////////////
//////       Momentum-Space Calling Methods       ////////
//////////////////////////////////////////////////////////

	// int NLargeX;
	// 	int NYTot;
	// 	double xmax;
	// 	double YminLX;
	// 	double dYLX;
int MVDipole::indexY(double Y){	
    int j;
    if (Y < mvconf.getYmin()) {
        // Left side: YminLX <= Y < Ymin
        j = static_cast<int>( std::floor((Y - YminLX) / dYLX) );
    } else {
        // Right side: Ymin <= Y <= Ymax (rows counted from Ymin, not from 0)
        j = NLargeX + static_cast<int>( std::floor((Y - mvconf.getYmin()) / mvconf.getdY() ) );
    }
    // rounding (e.g. dY written with 6 digits) can give the last row; the caller also reads row j+1
    return std::max(0, std::min(j, NYTot-2));
}

double MVDipole::FundamentalDipole_k(double x, double k, double T){
	double tmp=0;
	double Y_t = mvconf.getY(x);
	if( x>=0.999 ){return 0.0;}
	else if ( k<0 || k >= mvconf.getKmax() ){return 0.0;}
	else if ( T > mvconf.getTmax() ){return 0.0;}
	else if( Y_t < YminLX ){return 0;}
	else if( NLargeX==0 && Y_t < mvconf.getYmin() ){return 0;}   // no large-x rows: D = 0 for x > x0 (placeholder, see CLAUDE.md)
	else if( Y_t >= mvconf.getYmax() ){return 0;}
	else{
		int jY;
		double Y1,Y2;
		double tmp1,tmp2;

		double u_t;

		jY = indexY(Y_t);
		Y1= get_Y_at(jY-NLargeX);
		Y2= get_Y_at(jY+1-NLargeX);

		if(k <= mvconf.getKmin()){
			double U1,U2,U3,k1,k2,k3;
			double f1,f2,f3;
			double a1,b1,c1,a2,b2,c2;
			U1=mvconf.u_at(0);U2=mvconf.u_at(1),U3=mvconf.u_at(2);
			k1 = mvconf.k_at(0);k2 = mvconf.k_at(1);k3 = mvconf.k_at(2);
			if (T > mvconf.getTmin()){
				f1 = gsl_spline2d_eval(DF_spl[jY],U1,T,NULL, NULL);
				f2 = gsl_spline2d_eval(DF_spl[jY],U2,T,NULL, NULL);
				f3 = gsl_spline2d_eval(DF_spl[jY],U3,T,NULL, NULL);
			}
			else{
				f1 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DF_spl[jY],U1,mvconf.getTmin(),NULL, NULL);
				f2 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DF_spl[jY],U2,mvconf.getTmin(),NULL, NULL);
				f3 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DF_spl[jY],U3,mvconf.getTmin(),NULL, NULL);
			}
			quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a1,b1,c1);

			if (T > mvconf.getTmin()){
				f1 = gsl_spline2d_eval(DF_spl[jY+1],U1,T,NULL, NULL);
				f2 = gsl_spline2d_eval(DF_spl[jY+1],U2,T,NULL, NULL);
				f3 = gsl_spline2d_eval(DF_spl[jY+1],U3,T,NULL, NULL);
			}
			else{
				f1 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DF_spl[jY+1],U1,mvconf.getTmin(),NULL, NULL);
				f2 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DF_spl[jY+1],U2,mvconf.getTmin(),NULL, NULL);
				f3 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DF_spl[jY+1],U3,mvconf.getTmin(),NULL, NULL);
			}

			quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a2,b2,c2);
			tmp1= a1*k*k + b1*k+ c1;
			tmp2= a2*k*k + b2*k+ c2;
		}

		else{
			u_t = mvconf.getu(k);
			
			if (T > mvconf.getTmin()){
				tmp1 = gsl_spline2d_eval(DF_spl[jY],u_t,T,NULL, NULL);
				tmp2 = gsl_spline2d_eval(DF_spl[jY+1],u_t,T,NULL, NULL);
			}
			else{
				tmp1 = (T/mvconf.getTmin()) *gsl_spline2d_eval(DF_spl[jY],u_t,mvconf.getTmin(),NULL, NULL);
				tmp2 = (T/mvconf.getTmin()) *gsl_spline2d_eval(DF_spl[jY+1],u_t,mvconf.getTmin(),NULL, NULL);
			}
		}
		tmp = tmp1 + (tmp2-tmp1) * (Y_t-Y1)/(Y2-Y1) ;  // lin interpolation in Q

	}
	return tmp;
}
double MVDipole::AdjointDipole_k_thread_local(double x, double k, double T){
	double tmp=0;
	double Y_t = mvconf.getY(x);
	if( x>=0.999 ){return 0.0;}
	else if ( k<0 || k >= mvconf.getKmax() ){return 0.0;}
	else if ( T > mvconf.getTmax() ){return 0.0;}
	else if( Y_t < YminLX ){return 0;}
	else if( NLargeX==0 && Y_t < mvconf.getYmin() ){return 0;}   // no large-x rows: D = 0 for x > x0 (placeholder, see CLAUDE.md)
	else if( Y_t > mvconf.getYmax() ){return 0;}
	else{
		int jY;
		double Y1,Y2;
		double tmp1,tmp2;

		double u_t;

		jY = indexY(Y_t);
		Y1= mvconf.Y_at(jY);
		Y2= mvconf.Y_at(jY+1);
		thread_local gsl_interp_accel* xacc_local = gsl_interp_accel_alloc();
        thread_local gsl_interp_accel* yacc_local = gsl_interp_accel_alloc();
		thread_local gsl_interp_accel* xacc2_local = gsl_interp_accel_alloc();
        thread_local gsl_interp_accel* yacc2_local = gsl_interp_accel_alloc();

		if(k <= mvconf.getKmin()){
			double U1,U2,U3,k1,k2,k3;
			double f1,f2,f3;
			double a1,b1,c1,a2,b2,c2;
			U1=mvconf.u_at(0);U2=mvconf.u_at(1),U3=mvconf.u_at(2);
			k1 = mvconf.k_at(0);k2 = mvconf.k_at(1);k3 = mvconf.k_at(2);
			if (T > mvconf.getTmin()){
				f1 = gsl_spline2d_eval(DA_spl[jY],U1,T,xacc_local, yacc_local);
				f2 = gsl_spline2d_eval(DA_spl[jY],U2,T,xacc_local, yacc_local);
				f3 = gsl_spline2d_eval(DA_spl[jY],U3,T,xacc_local, yacc_local);
			}
			else{
				f1 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DA_spl[jY],U1,mvconf.getTmin(),xacc_local, yacc_local);
				f2 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DA_spl[jY],U2,mvconf.getTmin(),xacc_local, yacc_local);
				f3 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DA_spl[jY],U3,mvconf.getTmin(),xacc_local, yacc_local);
			}
			quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a1,b1,c1);

			if (T > mvconf.getTmin()){
				f1 = gsl_spline2d_eval(DA_spl[jY+1],U1,T,xacc2_local, yacc2_local);
				f2 = gsl_spline2d_eval(DA_spl[jY+1],U2,T,xacc2_local, yacc2_local);
				f3 = gsl_spline2d_eval(DA_spl[jY+1],U3,T,xacc2_local, yacc2_local);
			}
			else{
				f1 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DA_spl[jY+1],U1,mvconf.getTmin(),xacc2_local, yacc2_local);
				f2 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DA_spl[jY+1],U2,mvconf.getTmin(),xacc2_local, yacc2_local);
				f3 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DA_spl[jY+1],U3,mvconf.getTmin(),xacc2_local, yacc2_local);
			}

			quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a2,b2,c2);
			tmp1= a1*k*k + b1*k+ c1;
			tmp2= a2*k*k + b2*k+ c2;
		}

		else{
			u_t = mvconf.getu(k);
			
			if (T > mvconf.getTmin()){
				tmp1 = gsl_spline2d_eval(DA_spl[jY],u_t,T,NULL, NULL);
				tmp2 = gsl_spline2d_eval(DA_spl[jY+1],u_t,T,NULL, NULL);
			}
			else{
				tmp1 = (T/mvconf.getTmin()) *gsl_spline2d_eval(DA_spl[jY],u_t,mvconf.getTmin(),xacc_local, yacc_local);
				tmp2 = (T/mvconf.getTmin()) *gsl_spline2d_eval(DA_spl[jY+1],u_t,mvconf.getTmin(),xacc2_local, yacc2_local);
			}
		}
		tmp = tmp1 + (tmp2-tmp1) * (Y_t-Y1)/(Y2-Y1) ;  // lin interpolation in Q

	}
	return tmp;
}

double MVDipole::AdjointDipole_k(double x, double k, double T){
	double tmp=0;
	double Y_t = mvconf.getY(x);
	if( x>=0.999 ){return 0.0;}
	else if ( k<0 || k >= mvconf.getKmax() ){return 0.0;}
	else if ( T > mvconf.getTmax() ){return 0.0;}
	else if( Y_t < YminLX ){return 0;}
	else if( NLargeX==0 && Y_t < mvconf.getYmin() ){return 0;}   // no large-x rows: D = 0 for x > x0 (placeholder, see CLAUDE.md)
	else if( Y_t >= mvconf.getYmax() ){return 0;}
	else{
		int jY;
		double Y1,Y2;
		double tmp1,tmp2;

		double u_t;

		jY = indexY(Y_t);
		// Y1= mvconf.Y_at(jY);
		// Y2= mvconf.Y_at(jY+1);
		Y1= get_Y_at(jY-NLargeX);
		Y2= get_Y_at(jY+1-NLargeX);
		
		if(k <= mvconf.getKmin()){
			double U1,U2,U3,k1,k2,k3;
			double f1,f2,f3;
			double a1,b1,c1,a2,b2,c2;
			U1=mvconf.u_at(0);U2=mvconf.u_at(1),U3=mvconf.u_at(2);
			k1 = mvconf.k_at(0);k2 = mvconf.k_at(1);k3 = mvconf.k_at(2);
			if (T > mvconf.getTmin()){
				f1 = gsl_spline2d_eval(DA_spl[jY],U1,T,NULL, NULL);
				f2 = gsl_spline2d_eval(DA_spl[jY],U2,T,NULL, NULL);
				f3 = gsl_spline2d_eval(DA_spl[jY],U3,T,NULL, NULL);
			}
			else{
				f1 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DA_spl[jY],U1,mvconf.getTmin(),NULL, NULL);
				f2 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DA_spl[jY],U2,mvconf.getTmin(),NULL, NULL);
				f3 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DA_spl[jY],U3,mvconf.getTmin(),NULL, NULL);
			}
			quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a1,b1,c1);

			if (T > mvconf.getTmin()){
				f1 = gsl_spline2d_eval(DA_spl[jY+1],U1,T,NULL, NULL);
				f2 = gsl_spline2d_eval(DA_spl[jY+1],U2,T,NULL, NULL);
				f3 = gsl_spline2d_eval(DA_spl[jY+1],U3,T,NULL, NULL);
			}
			else{
				f1 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DA_spl[jY+1],U1,mvconf.getTmin(),NULL, NULL);
				f2 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DA_spl[jY+1],U2,mvconf.getTmin(),NULL, NULL);
				f3 = (T/mvconf.getTmin()) * gsl_spline2d_eval(DA_spl[jY+1],U3,mvconf.getTmin(),NULL, NULL);
			}

			quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a2,b2,c2);
			tmp1= a1*k*k + b1*k+ c1;
			tmp2= a2*k*k + b2*k+ c2;
		}

		else{
			u_t = mvconf.getu(k);
			
			if (T > mvconf.getTmin()){
				tmp1 = gsl_spline2d_eval(DA_spl[jY],u_t,T,NULL, NULL);
				tmp2 = gsl_spline2d_eval(DA_spl[jY+1],u_t,T,NULL, NULL);
			}
			else{
				tmp1 = (T/mvconf.getTmin()) *gsl_spline2d_eval(DA_spl[jY],u_t,mvconf.getTmin(),NULL, NULL);
				tmp2 = (T/mvconf.getTmin()) *gsl_spline2d_eval(DA_spl[jY+1],u_t,mvconf.getTmin(),NULL, NULL);
			}
		}
		tmp = tmp1 + (tmp2-tmp1) * (Y_t-Y1)/(Y2-Y1) ;  // lin interpolation in Q

	}
	return tmp;
}


//////////////////////////////////////////////////////////
//////              Output Methods             ////////
//////////////////////////////////////////////////////////

void MVDipole::printProgress(double percentage1, double percentage2) {
	int val1 = (int) (percentage1 * 100);
	int lpad1 = (int) (percentage1 * PBWIDTH);
	int rpad1 = PBWIDTH - lpad1;

	int val2 = (int) (percentage2 * 100);
	int lpad2 = (int) (percentage2 * PBWIDTH);
	int rpad2 = PBWIDTH - lpad2;

	printf("\r%3d%% [%.*s%*s]", val1, lpad1, PBSTR, rpad1, "");
	printf("   %3d%% [%.*s%*s]", val2, lpad2, PBSTR, rpad2, "");
	// printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
	fflush(stdout);
}

//////////////////////////////////////////////////////////
//////                 Other Tools                ////////
//////////////////////////////////////////////////////////



bool MVDipole::IsPathExist(const std::string &s){
	struct stat buffer;
	return (stat (s.c_str(), &buffer) == 0);
}


void MVDipole::check_for_tabs_folders(){
  //Check the existance of the tabs-folder
  std::ostringstream tabspath;
  fs::path PWD = fs::current_path();
  tabspath << PWD.string() << "/tabs" ;
  if(!IsPathExist(tabspath.str())){ fs::create_directories(tabspath.str()); }
	//Create the tabs folder if not existing
	//Check the existance of the IPSat-Dipole folder
	std::ostringstream mvpath;
	mvpath << PWD.string()<< "/tabs/MV-Dipole" ;
	path_to_tabs = mvpath.str();
	if(!IsPathExist(path_to_tabs)){ fs::create_directories( path_to_tabs);}
	//Create the IPSat-Dipole folder if not existing
}

std::string MVDipole::get_set_name(int n){
  std::ostringstream temppath;
  temppath << path_to_tabs<< "/Set"<< n;
  return temppath.str();
}

bool MVDipole::check_dipole_sets(){
  bool is_computed =false;
  int set=0;// counter for the sets (which are mumbered as Setn, n={0,...,N})
  path_to_set=get_set_name(set);
  bool is_path=IsPathExist(path_to_set);
  while (is_path) {
	if(!fs::exists(path_to_set+"/config.yaml")){
		// incomplete set (transform interrupted): never reused
		std::cerr<<"[MVDipole]: WARNING: "<<path_to_set<<" has no config.yaml (interrupted transform?); skipping it"<<std::endl;
		set++;
		path_to_set=get_set_name(set);
		is_path=IsPathExist(path_to_set);
		continue;
	}
	MVConfig tab_config = MVConfig(path_to_set);
	is_computed = mvconf.compare(tab_config);// Checks configfiles match
	if(is_computed){break;}
	else{
		set++;
		path_to_set=get_set_name(set);// goes to next Set and
		is_path=IsPathExist(path_to_set); // checks if exists
    }
  }

  return is_computed;
}


//////////////////////////////////////////////////////////
//////           Extrapolation Tools              ////////
//////////////////////////////////////////////////////////



void MVDipole::quadratic_extrapolation(double x1,double f1,double x2,double f2,double x3,double f3,double& a,double& b,double& c){
	double DET = (x1-x2)*(x1-x3)*(x2-x3);
	a = ( f3*(x1 - x2) + f1*(x2 - x3) + f2*(-x1 + x3) ) /DET;
	b = ( f3*(-pow(x1,2) + pow(x2,2)) + f2*(pow(x1,2) - pow(x3,2)) + f1*(-pow(x2,2) + pow(x3,2)) ) /DET;
	c = ( f3*x1*(x1 - x2)*x2 + f1*x2*(x2 - x3)*x3 + f2*x1*x3*(-x1 + x3) )/DET;
}

void MVDipole::set_largeX_values(){
	YminLX = log(mvconf.getX0()/xmax);
	dYLX = (NLargeX>0) ? (mvconf.getYmin()-YminLX)/double(NLargeX) : 0.0;   // unused if NLargeX == 0
}

double MVDipole::get_Y_at(int i){
	if(i<0){return mvconf.getYmin() + i*dYLX;}
	else{return mvconf.Y_at(i);}
}



//////////////////////////////////////////////////////////
//////               Output Tools                 ////////
//////////////////////////////////////////////////////////

void MVDipole::test_output_fixed_T(double T){
	double x_t,Y_t;
	double r_t;
	std::ostringstream path_to_test;
	path_to_test << path_to_set<< "/position_dipoles_test_fixed_T_"<< T << ".dat" ;
	// std::cerr << path_to_dip.str() << std::endl;
	std::ofstream test_f(path_to_test.str());
	if (!test_f.is_open()) {
   		 std::cerr << "ERROR: Cannot open file: " << path_to_test.str() << std::endl;
	}

	for (int iy = - NLargeX; iy < mvconf.getNY(); iy++){
		Y_t = get_Y_at(iy);
		x_t = mvconf.getx(Y_t);
		for (int ir=0; ir < mvconf.getNr(); ir++)
		{
			r_t = mvconf.r_at(ir);
			test_f<< Y_t << "\t"<< r_t << "\t" << FundamentalDipole(x_t, r_t, T) << "\t"<< AdjointDipole(x_t, r_t, T) << std::endl;
		}
	}
	test_f.close();
}	



void MVDipole::dump_test_norm(){
	std::ofstream dip_f;
	std::ostringstream path_to_file;
	path_to_file << path_to_set<< "/Dipole_test.dat" ;
	dip_f.open(path_to_file.str());
	int Nx_t = 121;
	int NT_t = 121;
	int Nk_t = 121;
	double Y_MIN_OUT = YminLX;
	double Y_MAX_OUT = mvconf.getYmax() ;
	double Tmin=0.1;
	double Tmax=10;
	double dT = (Tmax-Tmin)/(NT_t-1.);
	double kmin=0.01;
	double kmax=50;
	double dk = log(kmax/kmin)/(Nk_t-1.);
	std::cout << Y_MIN_OUT << "\t" << Y_MAX_OUT << std::endl;
	double T_t,kt_t;
	double Dip_A_k,Dip_F_k;
	for (int iT = 0; iT < NT_t; iT++) {
		T_t=Tmin + iT * dT;
		
		for (int iy = 0; iy < Nx_t; iy++) {
			double y_t = Y_MIN_OUT + iy*(Y_MAX_OUT-Y_MIN_OUT)/(Nx_t-1.);
			double x_t = mvconf.getx(y_t);
			for (int ik = 0; ik < Nk_t; ik++) {
				kt_t = kmin*exp(ik*dk);
				Dip_F_k = FundamentalDipole_k(x_t,kt_t, T_t);
				Dip_A_k = AdjointDipole_k(x_t,kt_t, T_t);

				dip_f<< T_t<<"\t"<< x_t<<"\t"<< kt_t <<"\t" << Dip_F_k<<"\t" << Dip_A_k << std::endl;
			}
		}
	}
	dip_f.close();
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



