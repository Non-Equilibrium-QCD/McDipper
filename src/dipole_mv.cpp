

/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */
#include <algorithm>
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
	int BessQuad =6;
	int GenQuad =6;

	static const int NumberOfIntegrationPoints=64000;
	double BesselIntegrate(double k,gsl_function F){
			// SETUP GSL INTEGRATOR //
			gsl_integration_workspace *Workspace=gsl_integration_workspace_alloc(NumberOfIntegrationPoints);
			// CALCULATE INTEGRAL BY PARTITION //
			int max_steps=96000; double a,b; double Value,Error; double Sum=0.0; double PartialSum=0.0;
			for(int n=0; n< max_steps; n++){
					if(n%100==0){
							if(std::abs(PartialSum)<BesselTol*std::abs(Sum)){
									gsl_integration_workspace_free (Workspace);
									return Sum;
							}
							else if(n>0 && PartialSum==0.0 && Sum==0.0){
								return Sum ;
							}
							else{
									PartialSum=0.0;
							}
					}

					if(n==0){
							a = 0.0;
							b = gsl_sf_bessel_zero_J0(1);
					}
					else{
							a = b;
							b = gsl_sf_bessel_zero_J0(n+1);
					}


					gsl_integration_qag(&F,a/k,b/k,0,EpsRel,NumberOfIntegrationPoints,BessQuad,Workspace,&Value,&Error);

					if(gsl_isnan(Value)){
							std::cerr << "#ERROR IN NUMERICAL INTEGRATION " << a << " " << b <<  " " << k << " " << Value <<  " " << gsl_isnan(Value) << std::endl;
							return NAN;
					}
					Sum+=Value;  PartialSum+=Value;
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
		double DipMom=BesselIntegrate(Variables->p,GSL_FundDipole_Integrand);
		return DipMom;
	}


	double DipoleAdjFT(DipMVFT * Variables){
		// CALCULATE qg TMDS //
		gsl_function GSL_AdjDipole_Integrand;
		GSL_AdjDipole_Integrand.function=&DipoleAdjFT_Integrand;
		GSL_AdjDipole_Integrand.params=Variables;
		double DipMom=BesselIntegrate(Variables->p,GSL_AdjDipole_Integrand);
		return DipMom;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





MVDipole::MVDipole(){}

MVDipole::MVDipole(std::string PATH, int LargeXExt,int NLX,double XMax, bool Verbose){
	DipVerbose=Verbose;
	path= PATH;
	mvconf=MVConfig(path);
	if(DipVerbose){mvconf.terminal_setup_output();}
	initialize();
	get_LargeX_Type(LargeXExt);
	NLargeX=NLX;
	xmax=XMax;
	set_largeX_values();
	NTot= mvconf.getNY() + NLargeX;

	Get_Momentum_Dipoles();
}

MVDipole::MVDipole(Config * conf){
	DipVerbose=conf->get_Verbose();
	path= conf->get_modelPath();
	mvconf=MVConfig(path);
	if(DipVerbose){mvconf.terminal_setup_output();}
	initialize();
	get_LargeX_Type(conf->get_ModelParams(1));
	NLargeX=int(conf->get_ModelParams(2));
	xmax=conf->get_ModelParams(3);
	set_largeX_values();

	NTot= mvconf.getNY() + NLargeX;
	
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
  }

  void MVDipole::make_pos_dipole_interpolators(){
	SF_spl = new gsl_spline2d*[mvconf.getNY()] ;
	SA_spl = new gsl_spline2d*[mvconf.getNY()] ;

	xaccSA = new gsl_interp_accel*[mvconf.getNY()] ;
	yaccSA = new gsl_interp_accel*[mvconf.getNY()] ;

	xaccSF = new gsl_interp_accel*[mvconf.getNY()] ;
	yaccSF = new gsl_interp_accel*[mvconf.getNY()] ;
	if(DipVerbose){std::cout<<"[MVDipole]: Interpolator Grid created for Dipoles" <<std::endl;}

}


  void MVDipole::initialize(){
	make_pos_dipole_interpolators();
	//initialize SF(x,T,r) and SA(x,T,r)
	initialize_position_space_dipole_interpolators();
	// initialize QSA2(x,T) and QSF2(x,T,r)
	initialize_Qs2_interpolators();
	
  }

  void MVDipole::initialize_position_space_dipole_interpolators(){
	std::ostringstream Gname;
	Gname << path<< "/mv_dipoles.dat";
	FILE *table1 = fopen(Gname.str().c_str(),"r");

	// double Gf[IPsat_pars::NY][IPsat_pars::NR];
  	double GAarr[mvconf.getNr()*mvconf.getNT()];
	double GFarr[mvconf.getNr()*mvconf.getNT()];


	double Y[mvconf.getNY()];
	double r[mvconf.getNr()];
	double q[mvconf.getNr()];
	double T[mvconf.getNT()];


	double y_t,r_t,SA_t,SF_t,T_t;

	for (int iY=0; iY < mvconf.getNY(); iY++)
	{
		for (int iT=0; iT < mvconf.getNT(); iT++)
		{
			for (int ir=0; ir < mvconf.getNr(); ir++)
			{
				
				if (fscanf(table1,"%lf %lf %lf %lf %lf", &y_t, &T_t, &r_t, &SF_t, &SA_t) != 5){printf("[MVDipole]: Error reading MV-Dipole table!\n");exit(1);}
				GAarr[mvconf.getNr()*iT + ir]=-log(SA_t);
				GFarr[mvconf.getNr()*iT + ir]=-log(SA_t);
				
				if(ir==0 && iY ==0 ){T[iT]=T_t;}
				if(ir==0 && iT ==0 ){Y[iY]=y_t;}
				if(iY==0 && iT ==0 ){
					r[ir]=r_t;
					q[ir]=mvconf.getq(r_t);
				}
			}
		}
		const gsl_interp2d_type * TA= gsl_interp2d_bilinear;
		SA_spl[iY] = gsl_spline2d_alloc(TA,mvconf.getNr(), mvconf.getNT());
		xaccSA[iY] = gsl_interp_accel_alloc();
		yaccSA[iY] = gsl_interp_accel_alloc();
		gsl_spline2d_init(SA_spl[iY], q, T, GAarr, mvconf.getNr(), mvconf.getNT());

		const gsl_interp2d_type * TF = gsl_interp2d_bilinear;
		SF_spl[iY] = gsl_spline2d_alloc(TF,mvconf.getNr(), mvconf.getNT());
		xaccSF[iY] = gsl_interp_accel_alloc();
		yaccSF[iY] = gsl_interp_accel_alloc();
		gsl_spline2d_init(SF_spl[iY], q,T, GFarr, mvconf.getNr(), mvconf.getNT());
	}

	fclose(table1);
	if(DipVerbose){std::cerr<<"[MVDipole]: Dipole interpolator succesfully created! "<<std::endl;}

  }

  void MVDipole::initialize_Qs2_interpolators(){
	std::ostringstream Qname;
	Qname << path<< "/mv_Q2.dat";
	FILE *table1 = fopen(Qname.str().c_str(),"r");

	// double Gf[IPsat_pars::NY][IPsat_pars::NR];
  	double QA2arr[mvconf.getNY()*mvconf.getNT()];
	double QF2arr[mvconf.getNY()*mvconf.getNT()];


	double Y[mvconf.getNY()];
	double T[mvconf.getNT()];


	double y_t,QA2_t,QF2_t,T_t;

	for (int iY=0; iY < mvconf.getNY(); iY++)
	{
		for (int iT=0; iT < mvconf.getNT(); iT++)
		{
			if (fscanf(table1,"%lf %lf %lf %lf", &y_t, &T_t, &QF2_t, &QA2_t) != 4){printf("[MVDipole]: Error reading MV-Q2 table!\n");exit(1);}
			QA2arr[mvconf.getNT()*iY + iT]=QA2_t;
			QF2arr[mvconf.getNT()*iY + iT]=QF2_t;
			
			if(iY ==0 ){T[iT]=T_t;}
			if(iT ==0 ){Y[iY]=y_t;}
		}
	}
	fclose(table1);

	// Now we get interpolate the array
	const gsl_interp2d_type *TQA= gsl_interp2d_bicubic;
	Q2A_spl = gsl_spline2d_alloc(TQA, mvconf.getNT(), mvconf.getNY());
	xaccQ2A = gsl_interp_accel_alloc();
	yaccQ2A = gsl_interp_accel_alloc();
	gsl_spline2d_init(Q2A_spl, T, Y, QA2arr,mvconf.getNT(), mvconf.getNY());

	// std::cout<< "Adjoint ready" << std::endl;

	const gsl_interp2d_type *TQF= gsl_interp2d_bicubic;
	Q2F_spl = gsl_spline2d_alloc(TQF, mvconf.getNT(), mvconf.getNY());
	xaccQ2F = gsl_interp_accel_alloc();
	yaccQ2F = gsl_interp_accel_alloc();
	gsl_spline2d_init(Q2F_spl, T, Y, QF2arr, mvconf.getNT(), mvconf.getNY());
	// std::cout<< "Fundamental ready" << std: :endl;

	if(DipVerbose){std::cerr<<"[MVDipole]: Qs2 interpolator succesfully created. "<<std::endl;}
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


			jY = (int) ( (T-mvconf.getYmin())/mvconf.getdY() );
			Y1= mvconf.getYmin() + jY*mvconf.getdY();
			Y2=mvconf.getYmin() +  (jY+1)*mvconf.getdY();

			tmp1 = gsl_spline2d_eval(SF_spl[jY],q,T,xaccSF[jY], yaccSF[jY]);
			tmp2 = gsl_spline2d_eval(SF_spl[jY+1],q,T,xaccSF[jY+1], yaccSF[jY+1]);
			tmp= tmp1 + (tmp2-tmp1) * (Y-Y1)/(Y2-Y1) ;  // lin interpolation in Y
			dipole=exp(-tmp);
			
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

			jY = (int) ( (T-mvconf.getYmin())/mvconf.getdY() );
			Y1= mvconf.getYmin() + jY*mvconf.getdY();
			Y2=mvconf.getYmin() +  (jY+1)*mvconf.getdY();

			tmp1 = gsl_spline2d_eval(SA_spl[jY],q,T,xaccSA[jY], yaccSA[jY]);
			tmp2 = gsl_spline2d_eval(SA_spl[jY+1],q,T,xaccSA[jY+1], yaccSA[jY+1]);
			tmp= tmp1 + (tmp2-tmp1) * (Y-Y1)/(Y2-Y1) ;  // lin interpolation in Y
			dipole=exp(-tmp);
			
		}

	}
	
	return dipole;
}

double MVDipole::FundamentalQ2(double x, double T){
	double dipole=0.0;
	double Y = mvconf.getY(x);

	if(Y<mvconf.getYmin()){dipole = 0.0;}
	else if(Y>mvconf.getYmax()){dipole = 0.0;}
	else if(T<mvconf.getTmin()){dipole = 0.0;}
	else if(T>mvconf.getTmax()){dipole = 0.0;}
	else{dipole=gsl_spline2d_eval(Q2F_spl,T,Y,xaccQ2F, yaccQ2F);}
	return dipole;
}

double MVDipole::AdjointQ2(double x, double T){
	double dipole=0.0;
	double Y = mvconf.getY(x);

	if(Y<mvconf.getYmin()){dipole = 0.0;}
	else if(Y>mvconf.getYmax()){dipole = 0.0;}
	else if(T<mvconf.getTmin()){dipole = 0.0;}
	else if(T>mvconf.getTmax()){dipole = 0.0;}
	else{dipole=gsl_spline2d_eval(Q2A_spl,T,Y,xaccQ2A, yaccQ2A);}
	return dipole;
}



double MVDipole::get_LargeX_Dipole(double x, double r, double T, MVRep rep){
	double largexdip=0.0;
	double Ci ;
	if(rep==MVRep::Adjoint){Ci=2.;}
	else if(rep==MVRep::Fundamental){Ci=1.;}
	switch (XMatch)
	{
		case MVLargeX::Frozen:
			{	
				double expo,lambda,delta;
				expo = 0.25*(mvconf.getSigma0()*T);
				expo *= mvconf.getQs02_in_fm()*pow(r,2.);  // Here large-x extrapolation is missing 
				expo *= log(mvconf.getEc()*M_E+ 1.0/(r*mvconf.getLambdaMV_in_fm()) );
				largexdip=exp(-Ci*expo);
			}
			break;
		case MVLargeX::Simple:
			{
				double expo,lambda,delta;
				expo = 0.25*(mvconf.getSigma0()*T)*pow(r,2.);
				expo *= mvconf.getQs02_in_fm();  // Here large-x extrapolation is missing 
				expo *= log(mvconf.getEc()*M_E+ 1.0/(r*mvconf.getLambdaMV_in_fm()) );
				largexdip=exp(-Ci*expo);
			}
			break;
		case MVLargeX::TMD:
			break;
		case MVLargeX::PDFMatch:
			break;
		default:
			break;
	}
	return largexdip;
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
	DF_spl = new gsl_spline2d*[mvconf.getNY() ] ;
	DA_spl = new gsl_spline2d*[mvconf.getNY() ] ;

	xaccDA = new gsl_interp_accel*[mvconf.getNY() ] ;
	yaccDA = new gsl_interp_accel*[mvconf.getNY() ] ;

	xaccDF = new gsl_interp_accel*[mvconf.getNY() ] ;
	yaccDF = new gsl_interp_accel*[mvconf.getNY() ] ;
	if(DipVerbose){std::cout<<"\nInterpolator Grid created for Dipoles" <<std::endl;}

}

void MVDipole::Get_Momentum_Dipoles(){
	//Check wether it is pre-computed
	check_for_tabs_folders();
	computed_in_run = !(check_dipole_sets());
	// If not computed already -> compute new momentum space dipoles
	if(computed_in_run){Make_Momentum_Dipoles();}
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
	mvconf.write_config(path_to_set);
	// Transform
	if(MV_pars::HankelTransMode==0){
		if(DipVerbose){std::cout<<"[MVDipole]: Transforming IP-Sat Dipoles using the Fast Hankel Transform Method" <<std::endl;}
		Transform_Dipole();
	}
	else if(MV_pars::HankelTransMode==1){
		if(DipVerbose){std::cerr<<"[MVDipole]: Transforming IP-Sat Dipoles using the Bessel-Zero method integration method" <<std::endl;}
		Transform_Dipole_Naive();
	}

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

	if(DipVerbose){
		std::cout<< "[MVDipole]: Transforming Dipoles to momentum space"<< std::endl;
		std::cout<< "Total (T,Y,k)                                                          Local (Y,k)"<< std::endl; 
	}
	dip_f.open(path_to_dip.str());

	for (int iy = 0; iy < mvconf.getNY(); iy++) {
		Y_t = mvconf.Y_at(iy);
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
				if(DipVerbose){ 
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

				if(DipVerbose){
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
	if(DipVerbose){
		printProgress(1,1);
		std::cout<< std::endl;
	}

}

void MVDipole::Transform_Dipole_Naive(){
	// double DipF[N1];
	double DipA_k = 0;
	double DipF_k = 0;
	double T_t,x_t,Y_t;
	double k_t,u_t;
	std::ofstream dip_f;
	std::ostringstream path_to_dip;
	path_to_dip << path_to_set<< "/MomentumDipoles.dat" ;

	if(DipVerbose){
		std::cout<< "[MVDipole]: Transforming Dipoles to momentum space"<< std::endl;
		std::cout<< "Total (T,Y,k)                                                          Local (Y,k)"<< std::endl; 
	}
	dip_f.open(path_to_dip.str());

	for (int iy = 0; iy < mvconf.getNY(); iy++) {
		Y_t = mvconf.Y_at(iy);
		x_t = mvconf.getx(Y_t);
			
		for (int iT = 0; iT < mvconf.getNT(); iT++) {

			T_t=mvconf.T_at(iT);

			for (int ik = 0; ik < mvconf.getNK(); ik++) {
				k_t = mvconf.k_at(ik);
				u_t = mvconf.u_at(ik);
				if(T_t==0){
					dip_f<< Y_t<<"\t"<< T_t<<"\t"<< u_t <<"\t" << 0<<"\t" << 0 << std::endl;
				}
				else{
					DipMVFT HTPars={k_t*gen_pars::GeV_to_fmm1,x_t,T_t,this};
					DipF_k = HankelDipoleMV::DipoleFunFT(&HTPars)*gen_pars::fm2_to_GeVm2;
					DipA_k = HankelDipoleMV::DipoleAdjFT(&HTPars)*gen_pars::fm2_to_GeVm2;

					dip_f<< Y_t<<"\t"<<T_t<<"\t"<< u_t<<"\t" << DipF_k <<"\t" << DipA_k << std::endl;
				}
			}
			if(DipVerbose){
				if(remainder(iy,gen_pars::skip)==0){
					double percentage_done2 = double(iy)/double(IPsat_pars::y_dip_points);
					double percentage_done1 = double(iT)/double(IPsat_pars::T_dip_points);
					printProgress(percentage_done2,percentage_done1);}
			}
		}
	}
	dip_f.close();
	if(DipVerbose){
		printProgress(1,1);
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

	if(DipVerbose){std::cout<<"[MVDipole]: Importing momentum-space dipoles. SRCFILE=" << tablename.str()<< std::endl ;}

	double DipAArray[mvconf.getNK()*mvconf.getNT()];
	double DipFArray[mvconf.getNK()*mvconf.getNT()];

	double Tarray[mvconf.getNT()];
	double Uarray[mvconf.getNK()];
	for (int iT = 0; iT < mvconf.getNT(); iT++) {Tarray[iT]= mvconf.T_at(iT);}
	for (int iU = 0; iU < mvconf.getNK(); iU++) {Uarray[iU]= mvconf.u_at(iU);}

	//
	FILE *table = fopen(tablename.str().c_str(),"r");
	double T_t,y_t,q_t,DipA_t,DipF_t;
	for (int iY = 0; iY <NTot; iY++) {
		for (int iT = 0; iT < mvconf.getNT(); iT++) {		
			for (int iU = 0; iU < mvconf.getNK(); iU++) {

				if (fscanf(table,"%lf %lf %lf %lf %lf", &T_t, &y_t, &q_t, &DipF_t, &DipA_t) != 5){printf("Error reading distribution table!\n");exit(1);}
				DipFArray[mvconf.getNK()*iT + iU]=-log(DipF_t);
				DipAArray[mvconf.getNK()*iT + iU]=-log(DipA_t);
			}
		}
		
		const gsl_interp2d_type *TF= gsl_interp2d_bilinear;
	    DF_spl[iY] = gsl_spline2d_alloc(TF,mvconf.getNK(),mvconf.getNT());
	    xaccDF[iY] = gsl_interp_accel_alloc();
	    yaccDF[iY] = gsl_interp_accel_alloc();
	    gsl_spline2d_init(DF_spl[iY], Uarray, Tarray, DipFArray, mvconf.getNK(),mvconf.getNT());
		
		const gsl_interp2d_type *TA= gsl_interp2d_bilinear;
	    DA_spl[iY] = gsl_spline2d_alloc(TA,mvconf.getNK(),mvconf.getNT());
	    xaccDA[iY] = gsl_interp_accel_alloc();
	    yaccDA[iY] = gsl_interp_accel_alloc();
	    gsl_spline2d_init(DA_spl[iY], Uarray, Tarray, DipAArray,mvconf.getNK(),mvconf.getNT());
	                                                                                                                                                                                                                                      }
	fclose(table);
 	is_read=true;
	if(DipVerbose){std::cout<<"[MVDipole]: Momentum-space dipole ---> Imported!" << std::endl;}
  	return is_read;
}

// Not sure if this will be needed this time though

//////////////////////////////////////////////////////////
//////       Momentum-Space Calling Methods       ////////
//////////////////////////////////////////////////////////

double MVDipole::FundamentalDipole_k(double x, double k, double T){
	double tmp=0;
	double Y_t = mvconf.getY(x);
	if( x>=0.999 ){return 0.0;}
	else if ( k<0 || k >= mvconf.getKmax() ){return 0.0;}
	else if ( T > mvconf.getTmax() ){return 0.0;}
	else if( Y_t < YminLX ){return 0;}
	else if( Y_t > mvconf.getYmax() ){return 0;}
	else{
		// CHECK T-SCALING!
		if (T < mvconf.getTmin()){
			double T1= mvconf.T_at(0);
			double T2= mvconf.T_at(1);
			double Q12 = gsl_spline2d_eval(Q2F_spl,T1,Y_t,xaccQ2F, yaccQ2F);
			double Q22 = gsl_spline2d_eval(Q2F_spl,T2,Y_t,xaccQ2F, yaccQ2F);
			double a=(-(Q22*T1) + Q12*T2)/(T1*(T1 - T2)*T2);
			double b= (Q22*pow(T1,2) - Q12*pow(T2,2))/(pow(T1,2)*T2 - T1*pow(T2,2));
			// double a = Q12*pow(T1,-2.);
			double Q2T = a * pow(T,2.) + b * T;
			tmp= (Q12/Q2T)*FundamentalDipole_k(x, k * sqrt(Q12/Q2T), T1);
		}
		else{
			int jY;
			double Y1,Y2;
			double tmp1,tmp2;

			double Q_t;

			jY = (int) ( (Y_t-mvconf.getYmin())/mvconf.getdY() );
			Y1= mvconf.Y_at(jY);
			Y2= mvconf.Y_at(jY+1);

			if(k <= mvconf.getKmin()){
				double U1,U2,U3,k1,k2,k3;
				double f1,f2,f3;
				double a1,b1,c1,a2,b2,c2;
				U1=mvconf.u_at(0);U2=mvconf.u_at(1),U3=mvconf.u_at(2);
				k1 = mvconf.k_at(0);k2 = mvconf.k_at(1);k3 = mvconf.k_at(2);
				f1 = gsl_spline2d_eval(DF_spl[jY],U1,T,xaccDF[jY], yaccDF[jY]);
				f2 = gsl_spline2d_eval(DF_spl[jY],U2,T,xaccDF[jY], yaccDF[jY]);
				f3 = gsl_spline2d_eval(DF_spl[jY],U3,T,xaccDF[jY], yaccDF[jY]);
				quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a1,b1,c1);

				f1 = gsl_spline2d_eval(DF_spl[jY+1],U1,T,xaccDF[jY], yaccDF[jY]);
				f2 = gsl_spline2d_eval(DF_spl[jY+1],U2,T,xaccDF[jY], yaccDF[jY]);
				f3 = gsl_spline2d_eval(DF_spl[jY+1],U3,T,xaccDF[jY], yaccDF[jY]);
				quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a2,b2,c2);
				tmp1= a1*k*k + b1*k+ c1;
				tmp2= a2*k*k + b2*k+ c2;
			}
			else{
				Q_t = mvconf.getu(k);
				tmp1 = gsl_spline2d_eval(DF_spl[jY],Q_t,Y_t,xaccDF[jY], yaccDF[jY]);
				tmp2 = gsl_spline2d_eval(DF_spl[jY+1],Q_t,Y_t,xaccDF[jY+1], yaccDF[jY+1]);
			}
			tmp = tmp1 + (tmp2-tmp1) * (Y_t-Y1)/(Y2-Y1) ;  // lin interpolation in Q
		}
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
	else if( Y_t > mvconf.getYmax() ){return 0;}
	else{
		// CHECK T-SCALING!
		if (T < mvconf.getTmin()){
			double T1= mvconf.T_at(0);
			double T2= mvconf.T_at(1);
			double Q12 = gsl_spline2d_eval(Q2A_spl,T1,Y_t,xaccQ2A, yaccQ2A);
			double Q22 = gsl_spline2d_eval(Q2A_spl,T2,Y_t,xaccQ2A, yaccQ2A);
			double a=(-(Q22*T1) + Q12*T2)/(T1*(T1 - T2)*T2);
			double b= (Q22*pow(T1,2) - Q12*pow(T2,2))/(pow(T1,2)*T2 - T1*pow(T2,2));
			// double a = Q12*pow(T1,-2.);
			double Q2T = a * pow(T,2.) + b * T;
			tmp= (Q12/Q2T)*AdjointDipole_k(x, k * sqrt(Q12/Q2T), T1);
		}
		else{
			int jY;
			double Y1,Y2;
			double tmp1,tmp2;

			double Q_t;

			jY = (int) ( (Y_t-mvconf.getYmin())/mvconf.getdY() );
			Y1= mvconf.Y_at(jY);
			Y2= mvconf.Y_at(jY+1);

			if(k <= mvconf.getKmin()){
				double U1,U2,U3,k1,k2,k3;
				double f1,f2,f3;
				double a1,b1,c1,a2,b2,c2;
				U1=mvconf.u_at(0);U2=mvconf.u_at(1),U3=mvconf.u_at(2);
				k1 = mvconf.k_at(0);k2 = mvconf.k_at(1);k3 = mvconf.k_at(2);
				f1 = gsl_spline2d_eval(DA_spl[jY],U1,T,xaccDA[jY], yaccDA[jY]);
				f2 = gsl_spline2d_eval(DA_spl[jY],U2,T,xaccDA[jY], yaccDA[jY]);
				f3 = gsl_spline2d_eval(DA_spl[jY],U3,T,xaccDA[jY], yaccDA[jY]);
				quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a1,b1,c1);

				f1 = gsl_spline2d_eval(DA_spl[jY+1],U1,T,xaccDA[jY], yaccDA[jY]);
				f2 = gsl_spline2d_eval(DA_spl[jY+1],U2,T,xaccDA[jY], yaccDA[jY]);
				f3 = gsl_spline2d_eval(DA_spl[jY+1],U3,T,xaccDA[jY], yaccDA[jY]);
				quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a2,b2,c2);
				tmp1= a1*k*k + b1*k+ c1;
				tmp2= a2*k*k + b2*k+ c2;
			}
			else{
				Q_t = mvconf.getu(k);
				tmp1 = gsl_spline2d_eval(DA_spl[jY],Q_t,Y_t,xaccDA[jY], yaccDA[jY]);
				tmp2 = gsl_spline2d_eval(DA_spl[jY+1],Q_t,Y_t,xaccDA[jY+1], yaccDA[jY+1]);
			}
			tmp = tmp1 + (tmp2-tmp1) * (Y_t-Y1)/(Y2-Y1) ;  // lin interpolation in Q
		}
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
	dYLX = (mvconf.getYmin()-YminLX)/double(NLargeX);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



