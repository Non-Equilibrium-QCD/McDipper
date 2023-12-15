

/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */

#include <iostream>
#include <math.h>
#include <sstream>
#include <sys/stat.h>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_multifit.h>

namespace fs = std::filesystem;

// BESSEL FUNCTIONS //
#ifndef BesselJ0
#define BesselJ0(x) jn(0,x)
#endif

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include "include/dipole_ipsat.h"
#include "include/params_gen.h"



namespace HankelDipole{

	double BesselTol=1e-5;
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
		DipFT * params= ((struct DipFT *) pars);
		return ( 2.0*M_PI*r*BesselJ0( (params->p)*r ) ) * (params->dip)->AdjointDipole(params->x,r,params->T);

	}

	double DipoleFunFT_Integrand(double r,void * pars){
		DipFT * params= ((struct DipFT *) pars);
		return ( 2.0*M_PI*r*BesselJ0( (params->p)*r ) ) * (params->dip)->FundamentalDipole(params->x,r,params->T);
	}

	double DipoleFunFT(DipFT * Variables){
		// CALCULATE qg TMDS //
		gsl_function GSL_FundDipole_Integrand;
		GSL_FundDipole_Integrand.function=&DipoleFunFT_Integrand;
		GSL_FundDipole_Integrand.params=Variables;
		double DipMom=BesselIntegrate(Variables->p,GSL_FundDipole_Integrand);
		return DipMom;
	}


	double DipoleAdjFT(DipFT * Variables){
		// CALCULATE qg TMDS //
		gsl_function GSL_AdjDipole_Integrand;
		GSL_AdjDipole_Integrand.function=&DipoleAdjFT_Integrand;
		GSL_AdjDipole_Integrand.params=Variables;
		double DipMom=BesselIntegrate(Variables->p,GSL_AdjDipole_Integrand);
		return DipMom;
	}
}


Dipole::Dipole(){}

Dipole::Dipole(int set, bool Verbose){
	DipVerbose=Verbose;
	DipSet= set;
	HankelMode=IPsat_pars::HankelTransMode;

	xscaling = IPsat_pars::x0_scaling;

	char Gname[256];
	switch (set) {
			case 1:
				sprintf(Gname,"ipsat/alpha_S_xg_set1.dat");break;
			case 2:
				sprintf(Gname,"ipsat/alpha_S_xg_set2.dat");break;
			default:
        std::cerr<<"Error: Parameter set not implemented. Exiting.";exit(EXIT_FAILURE);
		}

	FILE *table1 = fopen(Gname,"r");
	// double Gf[IPsat_pars::NY][IPsat_pars::NR];
  	double Gf[IPsat_pars::NY*IPsat_pars::NR];
	double dGfdu[IPsat_pars::NY*IPsat_pars::NR];
	double d2Gfdu2[IPsat_pars::NY*IPsat_pars::NR];

	double Y[IPsat_pars::NY];
	[[maybe_unused]] double rY[IPsat_pars::NY][IPsat_pars::NR];
	[[maybe_unused]] double r[IPsat_pars::NR];

	[[maybe_unused]] double uY[IPsat_pars::NY][IPsat_pars::NR];
	double u[IPsat_pars::NR];

	double y_t,r_t,G_t,Gp_t,Gpp_t,u_t;

	for (int iY=0; iY < IPsat_pars::NY; iY++)
	{
		for (int ir=0; ir < IPsat_pars::NR; ir++)
		{
			if (fscanf(table1,"%lf %lf %lf %lf %lf", &y_t, &r_t, &G_t, &Gp_t, &Gpp_t) != 5){printf("Error reading table 1!\n");exit(1);}
			Y[iY]=y_t;
			if(ir==0){
				if ((r_t-IPsat_pars::rmin)/r_t > IPsat_pars::rmin_tol){std::cerr<<"Minimum relative radius, r_0 ="<< r_t<< ", not matching parameters, rmin=" <<IPsat_pars::rmin << ". Check IP-Sat input!" << set <<std::endl;exit(1);}
			}

			u_t=get_U(r_t);

			rY[iY][ir]=r_t;
			uY[iY][ir]=u_t;

			Gf[IPsat_pars::NY*ir+iY] = G_t;
			dGfdu[IPsat_pars::NY*ir+iY] = Gp_t;
			d2Gfdu2[IPsat_pars::NY*ir+iY] = Gpp_t;
			if(iY==0){
				r[ir]=r_t;
				u[ir]=u_t;
			};
		}
	}

	int ny = sizeof(Y) / sizeof(Y[0]);
	int nu = sizeof(u) / sizeof(u[0]);

	const gsl_interp2d_type *T= gsl_interp2d_bicubic;
	G_spl = gsl_spline2d_alloc(T, ny, nu);
	xacc = gsl_interp_accel_alloc();
	yacc = gsl_interp_accel_alloc();
	gsl_spline2d_init(G_spl, Y, u, Gf, ny, nu);

	const gsl_interp2d_type *TP= gsl_interp2d_bicubic;
	Gp_spl = gsl_spline2d_alloc(TP, ny, nu);
	xaccP = gsl_interp_accel_alloc();
	yaccP = gsl_interp_accel_alloc();
	gsl_spline2d_init(Gp_spl, Y, u, dGfdu, ny, nu);

	const gsl_interp2d_type *TPP= gsl_interp2d_bicubic;
	Gpp_spl = gsl_spline2d_alloc(TPP, ny, nu);
	xaccPP = gsl_interp_accel_alloc();
	yaccPP = gsl_interp_accel_alloc();
	gsl_spline2d_init(Gpp_spl, Y, u, d2Gfdu2, ny, nu);

	fclose(table1);

	if(DipVerbose){std::cerr<<"Dipole succesfully created for Parameter Set  " << set <<std::endl;}

	Get_Momentum_Dipoles();


	//Initialisation of internal variables.
}

Dipole::~Dipole(){
  	gsl_spline2d_free (G_spl);
  	gsl_interp_accel_free (xacc);
		gsl_interp_accel_free (yacc);

		gsl_spline2d_free (Gp_spl);
  	gsl_interp_accel_free (xaccP);
		gsl_interp_accel_free (yaccP);

		gsl_spline2d_free (Gpp_spl);
  	gsl_interp_accel_free (xaccPP);
		gsl_interp_accel_free (yaccPP);

		for (int i = 0; i < IPsat_pars::NY; i++) {
	    gsl_spline2d_free(DA_spl[i]);
	    gsl_spline2d_free(DF_spl[i]);
	    gsl_interp_accel_free (xaccA[i]);
	    gsl_interp_accel_free (yaccA[i]);
	    gsl_interp_accel_free (xaccF[i]);
	    gsl_interp_accel_free (yaccF[i]);
		}
  }

//////////////////////////////////////////////////////////
//////       Position-Space Calling Methods       ////////
//////////////////////////////////////////////////////////

double Dipole::G(double x, double r){
double Y_t= get_Y(x);
double u_t= get_U(r);

double G_t=0;
if(x<0){std::cerr<<"Dipole.cpp: Error! Accesed x is x>1 or x<0!"<<std::endl;exit(1);}
if(x>1){return 0;}
if(r<0){std::cerr<<"Dipole.cpp: Error! Accesed r<0!" <<std::endl;exit(1);}
else{
// There are 9 cases for the regimes
	if(r < IPsat_pars::rmin){
		/* Here we deal with small-r cases. The second derivative is kept continuous to keep the Hankel integral well behaved.
		In this sector, logarithmic extrapolation is used to approximate the small-r dynamics. G"(r,Y)=G"(Y)*(A + B*log(r0/r)) */
		double u1= 0;
		double u2= IPsat_pars::du;

		double Ay,By;
		double GPP1,GPP2;

		if(Y_t < IPsat_pars::Ymin){
			/* Here we deal with small-Y cases (x~1). = In this sector, pow-law extrapolation is used, assuming G(r,0)= 0;
			This means G"(r,Y)=G"(r, Ymin)*(Y/Ymin)**B */
			double Y1= IPsat_pars::Ymin;
			double Y2= IPsat_pars::Ymin+IPsat_pars::dYextra;
			double GPP1Y1= gsl_spline2d_eval(Gpp_spl,Y1,u1,xaccPP, yaccPP);
			double GPP1Y2= gsl_spline2d_eval(Gpp_spl,Y2,u1,xaccPP, yaccPP);
			double GPP2Y1= gsl_spline2d_eval(Gpp_spl,Y1,u2,xaccPP, yaccPP);
			double GPP2Y2= gsl_spline2d_eval(Gpp_spl,Y2,u2,xaccPP, yaccPP);

			double B1= log(GPP1Y2/GPP1Y1)/log(Y2/Y1);
			double B2= log(GPP2Y2/GPP2Y1)/log(Y2/Y1);
			GPP1 =  GPP1Y1* pow(Y_t/IPsat_pars::Ymin,B1) ;
			GPP2 =  GPP2Y1* pow(Y_t/IPsat_pars::Ymin,B2) ;


		}
		else if(Y_t > IPsat_pars::Ymax){
			/* Here we deal with large-Y cases (x~0). = Based on the dynamics small-x BK dynamics at fixed Q^2 we
			assume a behaviour of log(1/x). This means G"(r,Y)=A(r) + B(r)*log(x0/x) =  A(r) + B(r)*(Y-Ymax)  */
			double Y1= IPsat_pars::Ymax-IPsat_pars::dYextra;
			double Y2= IPsat_pars::Ymax;

			double GPP1Y1 = gsl_spline2d_eval(Gpp_spl,Y1,u1,xaccPP, yaccPP);
			double GPP1Y2 = gsl_spline2d_eval(Gpp_spl,Y2,u1,xaccPP, yaccPP);

			double GPP2Y1 = gsl_spline2d_eval(Gpp_spl,Y1,u2,xaccPP, yaccPP);
			double GPP2Y2 = gsl_spline2d_eval(Gpp_spl,Y2,u2,xaccPP, yaccPP);

			GPP1 = (GPP1Y2-GPP1Y1)*(Y_t - Y2)/(Y2-Y1) + GPP1Y2;
			GPP2 = (GPP2Y2-GPP2Y1)*(Y_t - Y2)/(Y2-Y1) + GPP2Y2;

		}
		else {/* Standard case, included in the data. Interpolation is naively used here. */
			GPP1 = gsl_spline2d_eval(Gpp_spl,Y_t,u1,xaccPP, yaccPP);
			GPP2 = gsl_spline2d_eval(Gpp_spl,Y_t,u2,xaccPP, yaccPP);
		}

		Ay= GPP1;
		By= (GPP1-GPP2)/u2;
		G_t = r*r*(2*Ay-2*By*u_t + 3*By)/4.;

	}
	else if(r >= IPsat_pars::rmax){
		/* Here we deal with large-r cases. The second derivative is kept continuous to keep the Hankel integral well behaved.
		In this sector, constant extrapolation is used to approximate the small-r dynamics. G"(r,Y)=K*G"(Y) */
		if(Y_t < IPsat_pars::Ymin){
			/* Here we deal with small-Y cases (x~1). = In this sector, pow-law extrapolation is used, assuming G(r,0)= 0;
			This means G"(r,Y)=G"(r, Ymin)*(Y/Ymin)**B */
			double Y1= IPsat_pars::Ymin;
			double Y2= IPsat_pars::Ymin+IPsat_pars::dYextra;
			double G1= gsl_spline2d_eval(Gpp_spl,Y1,IPsat_pars::umax,xaccPP, yaccPP);
			double G2= gsl_spline2d_eval(Gpp_spl,Y2,IPsat_pars::umax,xaccPP, yaccPP);
			double B= log(G2/G1)/log(Y2/Y1);
			G_t =  G1*(r*r/2)* pow(Y_t/IPsat_pars::Ymin,B) ;

		}
		else if(Y_t > IPsat_pars::Ymax){
			/* Here we deal with large-Y cases (x~0). = Based on the dynamics small-x BK dynamics at fixed Q^2 we
			assume a behaviour of log(1/x). This means G"(r,Y)=A(r) + B(r)*log(x0/x) =  A(r) + B(r)*(Y-Ymax)  */
			double Y1= IPsat_pars::Ymax-IPsat_pars::dYextra;
			double Y2= IPsat_pars::Ymax;
			double G2 = gsl_spline2d_eval(Gpp_spl,Y2,IPsat_pars::umax,xaccPP,yaccPP);
			double G1 = gsl_spline2d_eval(Gpp_spl,Y1,IPsat_pars::umax,xaccPP,yaccPP) ;
			G_t = (r*r/2)*( (G2-G1)*(Y_t - Y2)/(Y2-Y1) + G2 ) ;
		}
		else {/* Standard case, included in the data. Interpolation is naively used here. */
			G_t = (r*r/2)*gsl_spline2d_eval(Gpp_spl,Y_t,IPsat_pars::umax,xacc, yacc);
		}
	}
	else {
	/* Here we deal with r cases included in the data */
	if(Y_t < IPsat_pars::Ymin){
		/* Here we deal with small-Y cases (x~1). = In this sector, pow-law extrapolation is used, assuming G(r,0)= 0;
		This means G"(r,Y)=G"(r, Ymin)*(Y/Ymin)**B */
		double Y1= IPsat_pars::Ymin;
		double Y2= IPsat_pars::Ymin+IPsat_pars::dYextra;
		double G1=gsl_spline2d_eval(G_spl,Y1,u_t,xacc, yacc);
		double G2=gsl_spline2d_eval(G_spl,Y2,u_t,xacc, yacc);
		double B= log(G2/G1)/log(Y2/Y1);
		G_t =  G1* pow(Y_t/IPsat_pars::Ymin,B) ;
	}
	else if(Y_t > IPsat_pars::Ymax){
		/* Here we deal with large-Y cases (x~0). = Based on the dynamics small-x BK dynamics at fixed Q^2 we
		assume a behaviour of log(1/x). This means G"(r,Y)=A(r) + B(r)*log(x0/x) =  A(r) + B(r)*(Y-Ymax)  */
		double Y1= IPsat_pars::Ymax-IPsat_pars::dYextra;
		double Y2= IPsat_pars::Ymax;
		double G2 = gsl_spline2d_eval(G_spl,Y2,u_t,xacc, yacc);
		double G1 = gsl_spline2d_eval(G_spl,Y1,u_t,xacc,yacc) ;
		G_t = (G2-G1)*(Y_t - Y2)/(Y2-Y1) + G2;
	}
	else {/* Standard case, included in the data. Interpolation is naively used here. */
		G_t = gsl_spline2d_eval(G_spl,Y_t,u_t,xacc, yacc);
	}
}
}
return G_t;
}
double Dipole::GPrime(double x, double r){
double Y_t= get_Y(x);
double u_t= get_U(r);

double Gp_t=0;
if(x<0){std::cerr<<"Dipole.cpp: Error! Accesed x is x>1 or x<0!"<<std::endl;exit(1);}
if(x>1){return 0;}
if(r<0){std::cerr<<"Dipole.cpp: Error! Accesed r<0!" <<std::endl;exit(1);}
else{
// There are 9 cases for the regimes
	if(r < IPsat_pars::rmin){
		/* Here we deal with small-r cases. The second derivative is kept continuous to keep the Hankel integral well behaved.
		In this sector, logarithmic extrapolation is used to approximate the small-r dynamics. G"(r,Y)=G"(Y)*(A + B*log(r0/r)) */

		double u1= 0;
		double u2= IPsat_pars::du;

		double Ay,By;
		double GPP1,GPP2;

		if(Y_t < IPsat_pars::Ymin){
			/* Here we deal with small-Y cases (x~1). = In this sector, pow-law extrapolation is used, assuming G(r,0)= 0;
			This means G"(r,Y)=G"(r, Ymin)*(Y/Ymin)**B */
			double Y1= IPsat_pars::Ymin;
			double Y2= IPsat_pars::Ymin+IPsat_pars::dYextra;
			double GPP1Y1= gsl_spline2d_eval(Gpp_spl,Y1,u1,xaccPP, yaccPP);
			double GPP1Y2= gsl_spline2d_eval(Gpp_spl,Y2,u1,xaccPP, yaccPP);
			double GPP2Y1= gsl_spline2d_eval(Gpp_spl,Y1,u2,xaccPP, yaccPP);
			double GPP2Y2= gsl_spline2d_eval(Gpp_spl,Y2,u2,xaccPP, yaccPP);

			double B1= log(GPP1Y2/GPP1Y1)/log(Y2/Y1);
			double B2= log(GPP2Y2/GPP2Y1)/log(Y2/Y1);
			GPP1 =  GPP1Y1* pow(Y_t/IPsat_pars::Ymin,B1) ;
			GPP2 =  GPP2Y1* pow(Y_t/IPsat_pars::Ymin,B2) ;

		}
		else if(Y_t > IPsat_pars::Ymax){
			/* Here we deal with large-Y cases (x~0). = Based on the dynamics small-x BK dynamics at fixed Q^2 we
			assume a behaviour of log(1/x). This means G"(r,Y)=A(r) + B(r)*log(x0/x) =  A(r) + B(r)*(Y-Ymax)  */
			double Y1= IPsat_pars::Ymax-IPsat_pars::dYextra;
			double Y2= IPsat_pars::Ymax;

			double GPP1Y1 = gsl_spline2d_eval(Gpp_spl,Y1,u1,xaccPP, yaccPP);
			double GPP1Y2 = gsl_spline2d_eval(Gpp_spl,Y2,u1,xaccPP, yaccPP);

			double GPP2Y1 = gsl_spline2d_eval(Gpp_spl,Y1,u2,xaccPP, yaccPP);
			double GPP2Y2 = gsl_spline2d_eval(Gpp_spl,Y2,u2,xaccPP, yaccPP);

			GPP1 = (GPP1Y2-GPP1Y1)*(Y_t - Y2)/(Y2-Y1) + GPP1Y2;
			GPP2 = (GPP2Y2-GPP2Y1)*(Y_t - Y2)/(Y2-Y1) + GPP2Y2;

		}
		else {/* Standard case, included in the data. Interpolation is naively used here. */
			GPP1 = gsl_spline2d_eval(Gpp_spl,Y_t,u1,xaccPP, yaccPP);
			GPP2 = gsl_spline2d_eval(Gpp_spl,Y_t,u2,xaccPP, yaccPP);
		}

		Ay= GPP1;
		By= (GPP1-GPP2)/u2;
		Gp_t = r*(Ay-By*u_t+By);//Ay - u_t *By;
	}
	else if(r >= IPsat_pars::rmax){
		/* Here we deal with large-r cases. The second derivative is kept continuous to keep the Hankel integral well behaved.
		In this sector, constant extrapolation is used to approximate the small-r dynamics. G"(r,Y)=K*G"(Y) */
		if(Y_t < IPsat_pars::Ymin){
			/* Here we deal with small-Y cases (x~1). = In this sector, pow-law extrapolation is used, assuming G(r,0)= 0;
			This means G"(r,Y)=G"(r, Ymin)*(Y/Ymin)**B */
			double Y1= IPsat_pars::Ymin;
			double Y2= IPsat_pars::Ymin+IPsat_pars::dYextra;
			double G1= gsl_spline2d_eval(Gpp_spl,Y1,IPsat_pars::umax,xaccPP, yaccPP);
			double G2= gsl_spline2d_eval(Gpp_spl,Y2,IPsat_pars::umax,xaccPP, yaccPP);
			double B= log(G2/G1)/log(Y2/Y1);
			Gp_t =  r*G1* pow(Y_t/IPsat_pars::Ymin,B) ;

		}
		else if(Y_t > IPsat_pars::Ymax){
			/* Here we deal with large-Y cases (x~0). = Based on the dynamics small-x BK dynamics at fixed Q^2 we
			assume a behaviour of log(1/x). This means G"(r,Y)=A(r) + B(r)*log(x0/x) =  A(r) + B(r)*(Y-Ymax)  */
			double Y1= IPsat_pars::Ymax-IPsat_pars::dYextra;
			double Y2= IPsat_pars::Ymax;
			double G2 = gsl_spline2d_eval(Gpp_spl,Y2,IPsat_pars::umax,xaccPP,yaccPP);
			double G1 = gsl_spline2d_eval(Gpp_spl,Y1,IPsat_pars::umax,xaccPP,yaccPP) ;
			Gp_t = r*( (G2-G1)*(Y_t - Y2)/(Y2-Y1) + G2 ) ;
		}
		else {/* Standard case, included in the data. Interpolation is naively used here. */
			Gp_t = r*gsl_spline2d_eval(Gpp_spl,Y_t,IPsat_pars::umax,xacc, yacc);
		}
	}
	else {
	/* Here we deal with r cases included in the data */
	if(Y_t < IPsat_pars::Ymin){
		/* Here we deal with small-Y cases (x~1). = In this sector, pow-law extrapolation is used, assuming G(r,0)= 0;
		This means G"(r,Y)=G"(r, Ymin)*(Y/Ymin)**B */
		double Y1= IPsat_pars::Ymin;
		double Y2= IPsat_pars::Ymin+IPsat_pars::dYextra;
		double G1=gsl_spline2d_eval(Gp_spl,Y1,u_t,xaccP, yaccP);
		double G2=gsl_spline2d_eval(Gp_spl,Y2,u_t,xaccP, yaccP);
		double B= log(G2/G1)/log(Y2/Y1);
		Gp_t =  G1* pow(Y_t/IPsat_pars::Ymin,B) ;
	}
	else if(Y_t > IPsat_pars::Ymax){
		/* Here we deal with large-Y cases (x~0). = Based on the dynamics small-x BK dynamics at fixed Q^2 we
		assume a behaviour of log(1/x). This means G"(r,Y)=A(r) + B(r)*log(x0/x) =  A(r) + B(r)*(Y-Ymax)  */
		double Y1= IPsat_pars::Ymax-IPsat_pars::dYextra;
		double Y2= IPsat_pars::Ymax;
		double G2 = gsl_spline2d_eval(Gp_spl,Y2,u_t,xaccP, yaccP);
		double G1 = gsl_spline2d_eval(Gp_spl,Y1,u_t,xaccP,yaccP) ;
		Gp_t = (G2-G1)*(Y_t - Y2)/(Y2-Y1) + G2;
	}
	else {/* Standard case, included in the data. Interpolation is naively used here. */
		Gp_t = gsl_spline2d_eval(Gp_spl,Y_t,u_t,xaccP, yaccP);
	}
}
}
return Gp_t;

}
double Dipole::GDoublePrime(double x, double r){
	double Y_t= get_Y(x);
	double u_t= get_U(r);

	double Gpp_t=0;
	if(x<0){std::cerr<<"Dipole.cpp: Error! Accesed x is x>1 or x<0!"<<std::endl;exit(1);}
	if(x>1){return 0;}
	if(r<0){std::cerr<<"Dipole.cpp: Error! Accesed r<0!" <<std::endl;exit(1);}
	else{
	// There are 9 cases for the regimes
		if(r < IPsat_pars::rmin){
			/* Here we deal with small-r cases. The second derivative is kept continuous to keep the Hankel integral well behaved.
			In this sector, logarithmic extrapolation is used to approximate the small-r dynamics. G"(r,Y)=G"(Y)*(A + B*log(r0/r)) */

			double u1= 0;
			double u2= IPsat_pars::du;

			double Ay,By;
			double GPP1,GPP2;

			if(Y_t < IPsat_pars::Ymin){
				/* Here we deal with small-Y cases (x~1). = In this sector, pow-law extrapolation is used, assuming G(r,0)= 0;
				This means G"(r,Y)=G"(r, Ymin)*(Y/Ymin)**B */
				double Y1= IPsat_pars::Ymin;
				double Y2= IPsat_pars::Ymin+IPsat_pars::dYextra;
				double GPP1Y1= gsl_spline2d_eval(Gpp_spl,Y1,u1,xaccPP, yaccPP);
				double GPP1Y2= gsl_spline2d_eval(Gpp_spl,Y2,u1,xaccPP, yaccPP);
				double GPP2Y1= gsl_spline2d_eval(Gpp_spl,Y1,u2,xaccPP, yaccPP);
				double GPP2Y2= gsl_spline2d_eval(Gpp_spl,Y2,u2,xaccPP, yaccPP);

				double B1= log(GPP1Y2/GPP1Y1)/log(Y2/Y1);
				double B2= log(GPP2Y2/GPP2Y1)/log(Y2/Y1);
				GPP1 =  GPP1Y1* pow(Y_t/IPsat_pars::Ymin,B1) ;
				GPP2 =  GPP2Y1* pow(Y_t/IPsat_pars::Ymin,B2) ;

			}
			else if(Y_t > IPsat_pars::Ymax){
				/* Here we deal with large-Y cases (x~0). = Based on the dynamics small-x BK dynamics at fixed Q^2 we
				assume a behaviour of log(1/x). This means G"(r,Y)=A(r) + B(r)*log(x0/x) =  A(r) + B(r)*(Y-Ymax)  */
				double Y1= IPsat_pars::Ymax-IPsat_pars::dYextra;
				double Y2= IPsat_pars::Ymax;

				double GPP1Y1 = gsl_spline2d_eval(Gpp_spl,Y1,u1,xaccPP, yaccPP);
				double GPP1Y2 = gsl_spline2d_eval(Gpp_spl,Y2,u1,xaccPP, yaccPP);

				double GPP2Y1 = gsl_spline2d_eval(Gpp_spl,Y1,u2,xaccPP, yaccPP);
				double GPP2Y2 = gsl_spline2d_eval(Gpp_spl,Y2,u2,xaccPP, yaccPP);

				GPP1 = (GPP1Y2-GPP1Y1)*(Y_t - Y2)/(Y2-Y1) + GPP1Y2;
				GPP2 = (GPP2Y2-GPP2Y1)*(Y_t - Y2)/(Y2-Y1) + GPP2Y2;

			}
			else {/* Standard case, included in the data. Interpolation is naively used here. */
				GPP1 = gsl_spline2d_eval(Gpp_spl,Y_t,u1,xaccPP, yaccPP);
				GPP2 = gsl_spline2d_eval(Gpp_spl,Y_t,u2,xaccPP, yaccPP);
			}

			Ay= GPP1;
			By= (GPP1-GPP2)/u2;
			Gpp_t = Ay - u_t *By;//r*(Ay-By*u_t+By);

		}
		else if(r >= IPsat_pars::rmax){
			/* Here we deal with large-r cases. The second derivative is kept continuous to keep the Hankel integral well behaved.
			In this sector, constant extrapolation is used to approximate the small-r dynamics. G"(r,Y)=K*G"(Y) */
			if(Y_t < IPsat_pars::Ymin){
				/* Here we deal with small-Y cases (x~1). = In this sector, pow-law extrapolation is used, assuming G(r,0)= 0;
				This means G"(r,Y)=G"(r, Ymin)*(Y/Ymin)**B */
				double Y1= IPsat_pars::Ymin;
				double Y2= IPsat_pars::Ymin+IPsat_pars::dYextra;
				double G1= gsl_spline2d_eval(Gpp_spl,Y1,IPsat_pars::umax,xaccPP, yaccPP);
				double G2= gsl_spline2d_eval(Gpp_spl,Y2,IPsat_pars::umax,xaccPP, yaccPP);
				double B= log(G2/G1)/log(Y2/Y1);
				Gpp_t =  G1* pow(Y_t/IPsat_pars::Ymin,B) ;

			}
			else if(Y_t > IPsat_pars::Ymax){
				/* Here we deal with large-Y cases (x~0). = Based on the dynamics small-x BK dynamics at fixed Q^2 we
				assume a behaviour of log(1/x). This means G"(r,Y)=A(r) + B(r)*log(x0/x) =  A(r) + B(r)*(Y-Ymax)  */
				double Y1= IPsat_pars::Ymax-IPsat_pars::dYextra;
				double Y2= IPsat_pars::Ymax;
				double G2 = gsl_spline2d_eval(Gpp_spl,Y2,IPsat_pars::umax,xaccPP,yaccPP);
				double G1 = gsl_spline2d_eval(Gpp_spl,Y1,IPsat_pars::umax,xaccPP,yaccPP) ;
				Gpp_t = ( (G2-G1)*(Y_t - Y2)/(Y2-Y1) + G2 ) ;
			}
			else {/* Standard case, included in the data. Interpolation is naively used here. */
				Gpp_t = gsl_spline2d_eval(Gpp_spl,Y_t,IPsat_pars::umax,xacc, yacc);
			}
		}
		else {
		/* Here we deal with r cases included in the data */
		if(Y_t < IPsat_pars::Ymin){
			/* Here we deal with small-Y cases (x~1). = In this sector, pow-law extrapolation is used, assuming G(r,0)= 0;
			This means G"(r,Y)=G"(r, Ymin)*(Y/Ymin)**B */
			double Y1= IPsat_pars::Ymin;
			double Y2= IPsat_pars::Ymin+IPsat_pars::dYextra;
			double G1=gsl_spline2d_eval(Gpp_spl,Y1,u_t,xaccPP, yaccPP);
			double G2=gsl_spline2d_eval(Gpp_spl,Y2,u_t,xaccPP, yaccPP);
			double B= log(G2/G1)/log(Y2/Y1);
			Gpp_t =  G1* pow(Y_t/IPsat_pars::Ymin,B) ;
		}
		else if(Y_t > IPsat_pars::Ymax){
			/* Here we deal with large-Y cases (x~0). = Based on the dynamics small-x BK dynamics at fixed Q^2 we
			assume a behaviour of log(1/x). This means G"(r,Y)=A(r) + B(r)*log(x0/x) =  A(r) + B(r)*(Y-Ymax)  */
			double Y1= IPsat_pars::Ymax-IPsat_pars::dYextra;
			double Y2= IPsat_pars::Ymax;
			double G2 = gsl_spline2d_eval(Gpp_spl,Y2,u_t,xaccPP, yaccPP);
			double G1 = gsl_spline2d_eval(Gpp_spl,Y1,u_t,xaccPP,yaccPP) ;
			Gpp_t = (G2-G1)*(Y_t - Y2)/(Y2-Y1) + G2;
		}
		else {/* Standard case, included in the data. Interpolation is naively used here. */
			Gpp_t = gsl_spline2d_eval(Gpp_spl,Y_t,u_t,xaccPP, yaccPP);
		}
	}
	}
	return Gpp_t;
}
double Dipole::DDoublePrime(double x, double r, double T){
		double Gp = GPrime(x,r);
		double Gpp = GDoublePrime(x,r);
		return (1/r)*T * exp(-G(x, r)*T) * ( r*T*Gp*Gp - Gp - r*Gpp ) ;
}


double Dipole::DDoublePrimeAdjoint(double x, double r, double T){
		double GA= (gen_pars::CA/gen_pars::CF) * G(x, r);
		double GAp = (gen_pars::CA/gen_pars::CF) * GPrime(x,r);
		double GApp = (gen_pars::CA/gen_pars::CF) * GDoublePrime(x,r);
		return (1/r)*T * exp(-GA*T) * ( r*T*GAp*GAp - GAp - r*GApp ) ;
}

double Dipole::FundamentalDipole(double x, double r, double T){
	return exp( -G(x, r)*T );
}
double Dipole::AdjointDipole(double x, double r, double T){
	return exp( -(gen_pars::CA/gen_pars::CF) * G(x, r)*T );
}

//////////////////////////////////////////////////////////
//////            Construction Methods            ////////
//////////////////////////////////////////////////////////

void Dipole::Get_Momentum_Dipoles(){
	//Check wether it is pre-computed
	check_for_tabs_folders();
	computed_in_run = !(check_dipole_sets());
	//Make momentum parameters and grid
	make_momentum_parameters(computed_in_run);
	Make_Dipole_Grids();
	// If not computed already -> compute new momentum space dipoles
	if(computed_in_run){Make_Momentum_Dipoles();}
	// Output config
	print_Dk_config();
	// Read in dipole tabulations and interpolate
	read_in_dipoles();
	// Get Q2 for x-scaling at large-x
	get_Q2();
	// Make Output for selected dipoles
	make_test_output();
}



void Dipole::Make_Dipole_Grids(){

	// First Make the K-grid
	k_grid_homo= new double[N1];
	q_grid_homo= new double[N1];
	for (int ik = 0; ik < N1; ik++) {
			q_grid_homo[ik]=ik*DK;
			k_grid_homo[ik]=K_MIN*exp(ik*DK);
	}

	Y_grid= new double[N2];
	for (int iy = 0; iy < N2; iy++) {Y_grid[iy]=Y_MIN + iy*DY;}

	T_grid= new double[N3];
	for (int iT = 0; iT < N3; iT++) {T_grid[iT]=T_MIN + iT*DT;}

	if(DipVerbose){std::cout<<"---> Grid created for Dipoles with Nk = " << N1<< " and  NY = " << N2<<std::endl;}
	is_initialized=true;
}

void Dipole::Make_Dipole_Interpolators(){
	DF_spl = new gsl_spline2d*[N3] ;
	DA_spl = new gsl_spline2d*[N3] ;

	xaccA = new gsl_interp_accel*[N3] ;
	yaccA = new gsl_interp_accel*[N3] ;

	xaccF = new gsl_interp_accel*[N3] ;
	yaccF = new gsl_interp_accel*[N3] ;
	if(DipVerbose){std::cout<<"\nInterpolator Grid created for Dipoles" <<std::endl;}

}

void Dipole::read_in_dipoles(){
	Make_Dipole_Interpolators();
	bool is_imported = import_dipole(Rep::Fundamental) && import_dipole(Rep::Adjoint) ;
	if(!is_imported){std::cerr<<"Error while importing dipole!"<<std::endl;exit(EXIT_FAILURE);}
}

void Dipole::Make_Momentum_Dipoles(){
  fs::create_directories(path_to_set);
	MinFactor=HankelParameters::MinFactor;
	MaxFactor=HankelParameters::MinFactor;
	write_config();
	//Transform
	if(IPsat_pars::HankelTransMode==0){
		if(DipVerbose){std::cout<<"Transforming IP-Sat Dipoles using the Fast Hankel Transform Method" <<std::endl;}
		Transform_Dipole(Rep::Fundamental);
		Transform_Dipole(Rep::Adjoint);
	}
	else if(IPsat_pars::HankelTransMode==1){
		if(DipVerbose){std::cerr<<"Transforming IP-Sat Dipoles using the Bessel-Zero integration method" <<std::endl;}
		Transform_Dipole_Naive(Rep::Fundamental);
		Transform_Dipole_Naive(Rep::Adjoint);
	}

}

//////////////////////////////////////////////////////////
//////            Transformation Methods          ////////
//////////////////////////////////////////////////////////

void Dipole::Transform_Dipole(Rep rep){
	// double DipF[N1];

	double T_t,RDom_t,x_t;
	double Dip_k=0;

	Hankel Hank(0,false);

	std::ofstream dip_f;
  std::ostringstream path_to_dip;
	if(rep==Rep::Fundamental){path_to_dip << path_to_set<< "/FundamentalDipole.dat" ;}
	else if(rep==Rep::Adjoint){path_to_dip << path_to_set<< "/AdjointDipole.dat" ;}

	if(DipVerbose){
		if(rep==Rep::Fundamental){std::cout<< "Transforming Fundamental Dipole"<< std::endl;}
		else if(rep==Rep::Adjoint){std::cout<< "Transforming Adjoint Dipole "<< std::endl;}
		std::cout<< "Total (T,Y,k)                                                          Local (Y,k)"<< std::endl;  ;
	}  
	dip_f.open(path_to_dip.str());

	for (int iT = 0; iT < IPsat_pars::T_dip_points; iT++) {
		T_t=IPsat_pars::T_dip_min + iT * IPsat_pars::T_dip_dT;

		for (int iy = 0; iy < N2; iy++) {
			x_t=get_X(Y_grid[iy]);
			// Set Hankel Environment around the characteristic scale of the the dipole
			RDom_t= EffectiveRadius(2,x_t,T_t,rep);
			Hank.Initialize_Domain_w_CharScale(RDom_t);
			if(iy==0 && iT==0){
				if(DipVerbose){ std::cout<< "Transforming for discretized function using "<< Hank.get_Npoints() << " points "<<std::endl; }
			}
			//Set points
			for (int ix = 0; ix < Hank.get_Npoints(); ix++) {
				if(rep==Rep::Fundamental){
					Hank.set_FX(ix,FundamentalDipole(x_t,Hank.get_X(ix),T_t)) ;
				}else if(rep==Rep::Adjoint){
					Hank.set_FX(ix,AdjointDipole(x_t,Hank.get_X(ix),T_t)) ;
				}
			}
			// Transform
			Hank.Transform();
			//retrieval

			double k_int_grid[Hank.get_Npoints()];
			double Dk_int_grid[Hank.get_Npoints()];

			for (int ik = 0; ik < Hank.get_Npoints(); ik++) {
				// Important to get the right UNITS!!
				// k needs conversion fm^{-1}-> GeV
				// Dk needs conversion fm^{2}-> GeV^{-2}
				k_int_grid[ik]=Hank.get_K(ik)*gen_pars::fmm1_to_GeV;
				Dk_int_grid[ik]=2*M_PI*Hank.get_FK(ik)*gen_pars::fm2_to_GeVm2;
			}

			// interpolate 1D
			gsl_interp_accel *acc_t= gsl_interp_accel_alloc ();
			gsl_spline *spline_t = gsl_spline_alloc (gsl_interp_cspline, Hank.get_Npoints());

			gsl_spline_init (spline_t, k_int_grid, Dk_int_grid, Hank.get_Npoints());


			// Here extrapolation goes on
			for (int iq = 0; iq < N1; iq++) {

				if(k_grid_homo[iq]>k_int_grid[Hank.get_Npoints()-1]){ Dip_k =0;}
				else if(k_grid_homo[iq]<k_int_grid[0]){
					double k1,k2,k3;
					double f1,f2,f3;
					double a,b,c;
					double k_t= k_grid_homo[iq];
					k1 = k_int_grid[0]; k2 = k_int_grid[1]; k3 = k_int_grid[2];
					f1 = Dk_int_grid[0];f2 = Dk_int_grid[1];f3 = Dk_int_grid[2];
					quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a,b,c);
					Dip_k =a*k_t*k_t + b*k_t+ c;
				}
				else{Dip_k =gsl_spline_eval (spline_t,k_grid_homo[iq], acc_t);}

				dip_f<< T_t<<"\t"<< Y_grid[iy]<<"\t"<< q_grid_homo[iq]<<"\t" << Dip_k << std::endl;

				if(DipVerbose){
					if(remainder(iy,gen_pars::skip)==0){
						double percentage_done1 = double(iy)/double(IPsat_pars::y_dip_points);
						double percentage_done2 = double(iT)/double(IPsat_pars::T_dip_points);
						printProgress(percentage_done2,percentage_done1);}
				}
			}
			gsl_spline_free (spline_t);
    	gsl_interp_accel_free (acc_t);
		}
	}
	dip_f.close();
	if(DipVerbose){
		printProgress(1,1);
		std::cout<< std::endl;
	}
}

void Dipole::Transform_Dipole_Naive(Rep rep){
	// double DipF[N1];
	double Dip_k = 0;
	double T_t,x_t;

	std::ofstream dip_f;
  std::ostringstream path_to_dip;
	if(rep==Rep::Fundamental){path_to_dip << path_to_set<< "/FundamentalDipole.dat" ;}
	else if(rep==Rep::Adjoint){path_to_dip << path_to_set<< "/AdjointDipole.dat" ;}

	if(DipVerbose){
		if(rep==Rep::Fundamental){std::cout<< "Transforming Fundamental Dipole"<< std::endl;}
		else if(rep==Rep::Adjoint){std::cout<< "Transforming Adjoint Dipole "<< std::endl;}
		std::cout<< "Total (T,Y,k)                                                          Local (Y,k)"<< std::endl;  ;
	}
	dip_f.open(path_to_dip.str());

	for (int iT = 0; iT < IPsat_pars::T_dip_points; iT++) {
		T_t=IPsat_pars::T_dip_min + iT * IPsat_pars::T_dip_dT;

		for (int iy = 0; iy < N2; iy++) {
			x_t=get_X(Y_grid[iy]);
			for (int ik = 0; ik < N1; ik++) {
				if(T_t==0){
					dip_f<< T_t<<"\t"<< Y_grid[iy]<<"\t"<< q_grid_homo[ik]<<"\t" << 0 << std::endl;
				}
				else{
					DipFT HTPars={k_grid_homo[ik]*gen_pars::GeV_to_fmm1,x_t,T_t,this};
					if(rep==Rep::Fundamental){Dip_k = HankelDipole::DipoleFunFT(&HTPars)*gen_pars::fm2_to_GeVm2;}
					else if(rep==Rep::Adjoint){Dip_k = HankelDipole::DipoleAdjFT(&HTPars)*gen_pars::fm2_to_GeVm2;}

					dip_f<< T_t<<"\t"<< Y_grid[iy]<<"\t"<< q_grid_homo[ik]<<"\t" << Dip_k << std::endl;
				}
			}
			if(DipVerbose){
				if(remainder(iy,gen_pars::skip)==0){
					double percentage_done1 = double(iy)/double(IPsat_pars::y_dip_points);
					double percentage_done2 = double(iT)/double(IPsat_pars::T_dip_points);
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



void Dipole::Transform_Dipole_Naive_Test(Rep rep,double Tt){
	// double DipF[N1];
	double Dip_k = 0;
	double x_t;

	int Nx_t = 10;
	double x[ Nx_t ];
	double y[ Nx_t ];
	double Y_MIN_OUT = Y_MIN *1.01;
	double Y_MAX_OUT = Y_MAX *0.99;
	for (int iy = 0; iy < Nx_t; iy++) {
 		y[iy] = Y_MIN_OUT + iy*(Y_MAX_OUT-Y_MIN_OUT)/(Nx_t-1.);
		x[iy] = get_X(y[iy]);
	}


	// double A_ext, B_ext, C_ext;

	std::ofstream dip_f;
  std::ostringstream path_to_dip;
	if(rep==Rep::Fundamental){path_to_dip << path_to_set<< "/FundamentalDipole_test_T_"<< Tt<< ".dat" ;}
	else if(rep==Rep::Adjoint){path_to_dip << path_to_set<< "/AdjointDipole_test_T_"<< Tt<< ".dat" ;}

	if(DipVerbose){
		if(rep==Rep::Fundamental){std::cout<< "Transforming Fundamental Dipole"<< std::endl;}
		else if(rep==Rep::Adjoint){std::cout<< "Transforming Adjoint Dipole "<< std::endl;}
	}
	dip_f.open(path_to_dip.str());
	dip_f<< "# k";
	for (int ix = 0; ix < Nx_t; ix++) {dip_f<< "\t D(x="<<x[ix]<<")";}
	dip_f<< std::endl;

	for (int ik = 0; ik < N1; ik++) {
		dip_f<< K_MIN*exp(q_grid_homo[ik]);
		for (int iy = 0; iy < Nx_t; iy++) {
			x_t=get_X(y[iy]);
			DipFT HTPars={k_grid_homo[ik]*gen_pars::GeV_to_fmm1,x_t,Tt,this};
			if(rep==Rep::Fundamental){Dip_k = HankelDipole::DipoleFunFT(&HTPars)*gen_pars::fm2_to_GeVm2;}
			else if(rep==Rep::Adjoint){Dip_k = HankelDipole::DipoleAdjFT(&HTPars)*gen_pars::fm2_to_GeVm2;}
			dip_f<<"\t" << Dip_k ;
		}
		dip_f << std::endl;
	}

	dip_f.close();
}

//////////////////////////////////////////////////////////
//////         Transformation Aide-Methods        ////////
//////////////////////////////////////////////////////////

double Dipole::MomentDipole(int n, double x, double T, Rep rep){
	double r_t,u_t,c_n;
	double MomentSum=0;

	for (int iu = 0; iu < IPsat_pars::N_UMAX_INT; iu++) {
		 u_t= IPsat_pars::du*iu;
		 r_t= IPsat_pars::rmin*exp(u_t);
		 if(iu==0 || iu==IPsat_pars::N_UMAX_INT-1){c_n=1/2.;}
		 else{c_n=1.;}
		 if(rep==Rep::Fundamental){ MomentSum += c_n * pow(r_t,n+2) * FundamentalDipole(x,r_t,T);}
		 else if(rep==Rep::Adjoint){ MomentSum += c_n * pow(r_t,n+2) * AdjointDipole(x,r_t,T);}

	}
	return 2*M_PI*IPsat_pars::du*MomentSum;
}

double Dipole::EffectiveRadius(int n, double x, double T, Rep rep){
	return pow(MomentDipole(n,x,T,rep)/MomentDipole(0,x,T,rep),1./n);
}

//////////////////////////////////////////////////////////
//////       Momentum-Space Calling Methods       ////////
//////////////////////////////////////////////////////////

double Dipole::FundamentalDipole_k(double x, double k, double T){
	double tmp=0;
	if( x>=0.999 ){return 0.0;}
	else if ( k<0 || k >= K_MAX ){return 0.0;}
	else if ( T <= T_MIN || T >= T_MAX ){return 0.0;}
	else if(x <= X_MIN ){return 0;}
	else{
		double Y_t = get_Y(x);
		if(x > xscaling){
			double Y0 = get_Y(xscaling);
			double Q20= gsl_spline2d_eval(Q2F_spl,T,Y0,xaccQ2F, yaccQ2F);
			double Q2X; 
			if ( x > X_MAX ){
				double x1= get_X(Y_grid[0]);
				double x2= get_X(Y_grid[1]);
				double Q21= gsl_spline2d_eval(Q2F_spl,T,Y_grid[0],xaccQ2F, yaccQ2F);
				double Q22= gsl_spline2d_eval(Q2F_spl,T,Y_grid[1],xaccQ2F, yaccQ2F);
				double b= log(Q21/Q22)/ log(  (1-x1)/(1-x2) );
				double a= Q21 *pow(1-x1,-b);
				Q2X = a*pow(1-x,b);
			}
			else{
				Q2X = gsl_spline2d_eval(Q2A_spl,T,Y_t,xaccQ2A, yaccQ2A);
			}
			tmp= (Q20/Q2X)*FundamentalDipole_k(xscaling, k * sqrt(Q20/Q2X), T);
		}
		else{
			int jT;
			double T1,T2;
			double tmp1,tmp2;

			double Q_t;

			jT = (int) ( (T-T_MIN)/DT );
			T1= T_MIN + jT*DT;
			T2= T_MIN + (jT+1)*DT;

			if(k <= K_MIN){
				double Q1,Q2,Q3,k1,k2,k3;
				double f1,f2,f3;
				double a1,b1,c1,a2,b2,c2;
				Q1=0;Q2=DK,Q3=2*DK;
				k1 = K_MIN*exp(Q1);k2 = K_MIN*exp(Q2);k3 = K_MIN*exp(Q3);
				f1 = gsl_spline2d_eval(DF_spl[jT],Q1,Y_t,xaccF[jT], yaccF[jT]);
				f2 = gsl_spline2d_eval(DF_spl[jT],Q2,Y_t,xaccF[jT], yaccF[jT]);
				f3 = gsl_spline2d_eval(DF_spl[jT],Q3,Y_t,xaccF[jT], yaccF[jT]);
				quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a1,b1,c1);

				f1 = gsl_spline2d_eval(DF_spl[jT+1],Q1,Y_t,xaccF[jT], yaccF[jT]);
				f2 = gsl_spline2d_eval(DF_spl[jT+1],Q2,Y_t,xaccF[jT], yaccF[jT]);
				f3 = gsl_spline2d_eval(DF_spl[jT+1],Q3,Y_t,xaccF[jT], yaccF[jT]);
				quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a2,b2,c2);
				tmp1= a1*k*k + b1*k+ c1;
				tmp2= a2*k*k + b2*k+ c2;
			}
			else{
				Q_t = log(k/K_MIN);
				tmp1 = gsl_spline2d_eval(DF_spl[jT],Q_t,Y_t,xaccF[jT], yaccF[jT]);
				tmp2 = gsl_spline2d_eval(DF_spl[jT+1],Q_t,Y_t,xaccF[jT+1], yaccF[jT+1]);
			}
			tmp = tmp1 + (tmp2-tmp1) * (T-T1)/(T2-T1) ;  // lin interpolation in Q
		}
	}
	return tmp;
}

double Dipole::AdjointDipole_k(double x, double k, double T){
	double tmp=0;
	if( x>=0.999 ){tmp=0.0;}
	else if (  k<0 || k >= K_MAX ){tmp=0.0;}
	else if ( T <= T_MIN || T >= T_MAX ){tmp=0.0;}
	else if(x <= X_MIN ){tmp=0.0;}
	else{
		double Y_t = get_Y(x);
		if(x > xscaling){
			double Y0 = get_Y(xscaling);
			double Q20= gsl_spline2d_eval(Q2A_spl,T,Y0,xaccQ2A, yaccQ2A);
			double Q2X; 
			if ( x > X_MAX ){
				double x1= get_X(Y_grid[0]);
				double x2= get_X(Y_grid[1]);
				double Q21= gsl_spline2d_eval(Q2A_spl,T,Y_grid[0],xaccQ2A, yaccQ2A);
				double Q22= gsl_spline2d_eval(Q2A_spl,T,Y_grid[1],xaccQ2A, yaccQ2A);
				double b= log(Q21/Q22)/ log(  (1-x1)/(1-x2) );
				double a= Q21 *pow(1-x1,-b);
				Q2X = a*pow(1-x,b);
			}
			else{
				Q2X = gsl_spline2d_eval(Q2A_spl,T,Y_t,xaccQ2A, yaccQ2A);
			}
			tmp= (Q20/Q2X)*AdjointDipole_k(xscaling, k * sqrt(Q20/Q2X), T);
			// if(tmp!=tmp){std::cerr<<"TEST"<<x<<"\t"<<k<<"\t"<<T<< std::endl;}

		}
		else{
			int jT;
			double T1,T2;
			double tmp1,tmp2;

			double Q_t;

			jT = (int) ( (T-T_MIN)/DT );
			T1= T_MIN + jT*DT;
			T2= T_MIN + (jT+1)*DT;

			if(k <= K_MIN){
				double Q1,Q2,Q3,k1,k2,k3;
				double f1,f2,f3;
				double a1,b1,c1,a2,b2,c2;
				Q1=0;Q2=DK,Q3=2*DK;
				k1 = K_MIN*exp(Q1);k2 = K_MIN*exp(Q2);k3 = K_MIN*exp(Q3);
				f1 = gsl_spline2d_eval(DA_spl[jT],Q1,Y_t,xaccA[jT], yaccA[jT]);
				f2 = gsl_spline2d_eval(DA_spl[jT],Q2,Y_t,xaccA[jT], yaccA[jT]);
				f3 = gsl_spline2d_eval(DA_spl[jT],Q3,Y_t,xaccA[jT], yaccA[jT]);
				quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a1,b1,c1);

				f1 = gsl_spline2d_eval(DA_spl[jT+1],Q1,Y_t,xaccA[jT], yaccA[jT]);
				f2 = gsl_spline2d_eval(DA_spl[jT+1],Q2,Y_t,xaccA[jT], yaccA[jT]);
				f3 = gsl_spline2d_eval(DA_spl[jT+1],Q3,Y_t,xaccA[jT], yaccA[jT]);
				quadratic_extrapolation(k1,f1,k2,f2,k3,f3,a2,b2,c2);
				tmp1= a1*k*k + b1*k+ c1;
				tmp2= a2*k*k + b2*k+ c2;
			}
			else{
				Q_t = log(k/K_MIN);
				tmp1 = gsl_spline2d_eval(DA_spl[jT],Q_t,Y_t,xaccA[jT], yaccA[jT]);
				tmp2 = gsl_spline2d_eval(DA_spl[jT+1],Q_t,Y_t,xaccA[jT+1], yaccA[jT+1]);
			}
			tmp = tmp1 + (tmp2-tmp1) * (T-T1)/(T2-T1) ;  // lin interpolation in Q
		}
	}
	return tmp;
}

//////////////////////////////////////////////////////////
//////              Importing Methods             ////////
//////////////////////////////////////////////////////////

bool Dipole::import_dipole(Rep dipole_rep){
	bool is_read=false;

  std::ostringstream tablename;

	if(dipole_rep==Rep::Fundamental){tablename << path_to_set <<"/FundamentalDipole.dat" ;}
	else if(dipole_rep==Rep::Adjoint){ tablename << path_to_set <<"/AdjointDipole.dat" ;}

	if(DipVerbose){
		if(dipole_rep==Rep::Fundamental){std::cout<<"Importing Fundamental Dipole. SRCFILE=" << tablename.str()<< std::endl ;}
		else if(dipole_rep==Rep::Adjoint){ std::cout<<"Importing Adjoint Dipole. SRCFILE=" << tablename.str()<< std::endl ;}
	}

	double DipArray[N1*N2];

  double Yarray[N2];
	double Qarray[N1];
	for (int iY = 0; iY < N2; iY++) {Yarray[iY]= Y_MIN + iY*DY;}
	for (int iQ = 0; iQ < N1; iQ++) {Qarray[iQ]= iQ*DK;}

	//
	FILE *table = fopen(tablename.str().c_str(),"r");
  double T_t,y_t,q_t,Dip_t;
	for (int iT = 0; iT < N3; iT++) {
		for (int iY = 0; iY < N2; iY++) {
			for (int iq = 0; iq < N1; iq++) {
				if (fscanf(table,"%lf %lf %lf %lf", &T_t, &y_t, &q_t, &Dip_t) != 4){printf("Error reading distribution table!\n");exit(1);}
				DipArray[N1*iY + iq]=Dip_t;
			}
		}

		if(dipole_rep==Rep::Fundamental){
			const gsl_interp2d_type *T= gsl_interp2d_bilinear;
	    DF_spl[iT] = gsl_spline2d_alloc(T, N1, N2);
	    xaccF[iT] = gsl_interp_accel_alloc();
	    yaccF[iT] = gsl_interp_accel_alloc();
	    gsl_spline2d_init(DF_spl[iT], Qarray, Yarray, DipArray, N1, N2);
		}
		else if(dipole_rep==Rep::Adjoint){
			const gsl_interp2d_type *T= gsl_interp2d_bilinear;
	    DA_spl[iT] = gsl_spline2d_alloc(T, N1, N2);
	    xaccA[iT] = gsl_interp_accel_alloc();
	    yaccA[iT] = gsl_interp_accel_alloc();
	    gsl_spline2d_init(DA_spl[iT], Qarray, Yarray, DipArray, N1, N2);
		}
	}
	fclose(table);
 	is_read=true;
	if(DipVerbose){std::cout<<"---> Done "<< std::endl ;}
  return is_read;
}

//////////////////////////////
////////// TOOLS /////////////
//////////////////////////////

bool Dipole::IsPathExist(const std::string &s){
	struct stat buffer;
	return (stat (s.c_str(), &buffer) == 0);
}

void Dipole::check_for_tabs_folders(){
	//Check the existance of the tabs-folder
  std::ostringstream tabspath;
  fs::path PWD = fs::current_path();
  tabspath << PWD.string() << "/tabs" ;
  if(!IsPathExist(tabspath.str())){ fs::create_directories(tabspath.str()); }
	//Create the tabs folder if not existing
	//Check the existance of the IPSat-Dipole folder
	std::ostringstream ipsatpath;
	ipsatpath << PWD.string()<< "/tabs/IPSat-Dipole" ;
	path_to_tabs = ipsatpath.str();
	if(!IsPathExist(path_to_tabs)){ fs::create_directories( path_to_tabs);}
	//Create the IPSat-Dipole folder if not existing
}

std::string Dipole::get_set_name(int n){
  std::ostringstream temppath;
  temppath << path_to_tabs<< "/Set"<< n;
  return temppath.str();
}

void Dipole::get_name_and_value(std::string testline,std::string &name, std::string &value){
  std::string line_temp= testline;
  int pos = line_temp.find(":");
  name = line_temp.substr(0,pos);
  value = line_temp.substr(pos+1, line_temp.length()-pos );
  name.erase(std::remove_if(name.begin(),name.end(), ::isspace),name.end());
  value.erase(std::remove_if(value.begin(),value.end(), ::isspace),value.end());
}

//////////////////////////////
////////// CONFIG ////////////
//////////////////////////////

bool Dipole::check_dipole_sets(){
  bool is_computed =false;
  int set=0;// counter for the sets (which are mumbered as Setn, n={0,...,N})
  path_to_set=get_set_name(set);
  bool is_path=IsPathExist(path_to_set);
  while (is_path) {
		read_in_config();
    is_computed = check_config();// Checks configfiles match
    if(is_computed){break;}
    else{
      set++;
      path_to_set=get_set_name(set);// goes to next Set and
      is_path=IsPathExist(path_to_set); // checks if exists

    }
  }

  return is_computed;
}

void Dipole::process_dip_parameters(std::string testline){
  std::string name_t,value_t;
  get_name_and_value(testline,name_t,value_t);

 if(name_t=="ParameterSet"){DipSet_conf=std::stoi(value_t);}
 if(name_t=="TransMethod"){HankelMode_conf=std::stoi(value_t);}

  if(name_t=="kmin"){kmin_conf=std::stod(value_t);}
  if(name_t=="kmax"){kmax_conf=std::stod(value_t);}
	if(name_t=="NK"){Nk_conf=std::stoi(value_t);}

	if(name_t=="Ymin"){Ymin_conf=std::stod(value_t);}
  if(name_t=="Ymax"){Ymax_conf=std::stod(value_t);}
	if(name_t=="NY"){NY_conf=std::stoi(value_t);}

	if(name_t=="Tmin"){Tmin_conf=std::stod(value_t);}
  if(name_t=="Tmax"){Tmax_conf=std::stod(value_t);}
	if(name_t=="NT"){NT_conf=std::stoi(value_t);}

	if(name_t=="MinFactor"){MinFactor_conf=std::stod(value_t);}
	if(name_t=="MaxFactor"){MaxFactor_conf=std::stod(value_t);}

}

bool Dipole::check_config(){
	// This checks that the formerly computed dipole contains
	// the newly asked range for the dipole
  bool is_config = ( DipSet_conf==DipSet );
	is_config = is_config && ( HankelMode_conf==HankelMode);

	is_config = is_config && (kmin_conf<=IPsat_pars::k_dip_min);
	is_config = is_config && (kmax_conf>=IPsat_pars::k_dip_max);

	is_config = is_config && (Ymin_conf<=IPsat_pars::y_dip_min);
	is_config = is_config && (Ymax_conf>=IPsat_pars::y_dip_max);

	is_config = is_config && (Tmin_conf<=IPsat_pars::T_dip_min);
	is_config = is_config && (Tmax_conf>=IPsat_pars::T_dip_max);

	if(HankelMode_conf==0){
		is_config = is_config && (MinFactor_conf==HankelParameters::MinFactor);
		is_config = is_config && (MaxFactor_conf==HankelParameters::MaxFactor);
	}

  return is_config;
}

void Dipole::read_in_config(){
  std::string line;
  std::string header;

	std::ostringstream path_to_config;
  path_to_config << path_to_set<< "/config.yaml" ;

  std::ifstream dip_configfile(path_to_config.str());
  if (dip_configfile.is_open())
  {
    while ( getline (dip_configfile,line) )
    {
      if(line.length()>0){process_dip_parameters(line);}
      else{continue;}
    }
  }
  else{std::cerr<<"Dipole Error! Failed to open Dipole config-file!    Config path:   "<< path_to_config.str() <<std::endl; exit(EXIT_FAILURE); }
  dip_configfile.close();
  dk_conf=log(kmax_conf/kmin_conf)/(Nk_conf-1.);
	dY_conf=(Ymax_conf-Ymin_conf)/(NY_conf-1.);
	dT_conf=(Tmax_conf-Tmin_conf)/(NT_conf-1.);

}

void Dipole::make_momentum_parameters(bool in_run){
	if(in_run){
		N1=IPsat_pars::k_dip_points;//Number of r points before Hankel Transformation.
		N2=IPsat_pars::y_dip_points;
		N3=IPsat_pars::T_dip_points;//Number of T points

		DK= IPsat_pars::k_dip_dk;
		DY= IPsat_pars::y_dip_dy;
		DT= IPsat_pars::T_dip_dT;

		K_MIN=IPsat_pars::k_dip_min;
		Y_MIN=IPsat_pars::y_dip_min;
		T_MIN=IPsat_pars::T_dip_min;

		K_MAX=IPsat_pars::k_dip_max;
		Y_MAX=IPsat_pars::y_dip_max;
		T_MAX=IPsat_pars::T_dip_max;
	}else{
		read_in_config();
		N1=Nk_conf;//Number of r points before Hankel Transformation.
		N2=NY_conf;
		N3=NT_conf;//Number of T points

		DK= dk_conf;
		DY= dY_conf;
		DT= dT_conf;

		K_MIN=kmin_conf;
		Y_MIN=Ymin_conf;
		T_MIN=Tmin_conf;

		K_MAX=kmax_conf;
		Y_MAX=Ymax_conf;
		T_MAX=Tmax_conf;

		MinFactor=MinFactor_conf;
		MaxFactor=MaxFactor_conf;
	}
	X_MAX=exp(-Y_MIN);
	X_MIN=exp(-Y_MAX);

}


//////////////////////////////
////////// OUTPUT ////////////
//////////////////////////////

void Dipole::write_config(){
	std::ofstream dip_config_f;
  std::ostringstream path_to_config;
  path_to_config << path_to_set<< "/config.yaml" ;
	dip_config_f.open(path_to_config.str());
	dip_config_f<< "ParameterSet:   "<< DipSet<< std::endl;
	dip_config_f<< "TransMethod:   "<< HankelMode << std::endl;
	dip_config_f<< std::endl;
	dip_config_f<< "kmin: "<< IPsat_pars::k_dip_min<< std::endl;
	dip_config_f<< "kmax: "<< IPsat_pars::k_dip_max<< std::endl;
	dip_config_f<< "NK:   "<< IPsat_pars::k_dip_points<< std::endl;
	dip_config_f<< std::endl;
	dip_config_f<< "Ymin: "<< IPsat_pars::y_dip_min<< std::endl;
	dip_config_f<< "Ymax: "<< IPsat_pars::y_dip_max<< std::endl;
	dip_config_f<< "NY:   "<< IPsat_pars::y_dip_points<< std::endl;
	dip_config_f<< std::endl;
	dip_config_f<< "Tmin: "<< IPsat_pars::T_dip_min<< std::endl;
	dip_config_f<< "Tmax: "<< IPsat_pars::T_dip_max<< std::endl;
	dip_config_f<< "NT:   "<< IPsat_pars::T_dip_points<< std::endl;
	dip_config_f<< std::endl;
	if(HankelMode==0){
		dip_config_f<< "MinFactor: "<< HankelParameters::MinFactor<< std::endl;
		dip_config_f<< "MaxFactor: "<< HankelParameters::MaxFactor<< std::endl;
	}


	dip_config_f.close();
}

void Dipole::printProgress(double percentage1, double percentage2) {
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

void Dipole::print_Dk_config(){
	std::cout<< "ParameterSet:   "<< DipSet<< std::endl;
	std::cout<< "TransMethod:   "<< HankelMode<< std::endl;
	std::cout<< std::endl;
	std::cout<< "kmin: "<< K_MIN << std::endl;
	std::cout<< "kmax: "<< K_MAX << std::endl;
	std::cout<< "NK:   "<< N1<< std::endl;
	std::cout<< std::endl;
	std::cout<< "Ymin: "<< Y_MIN << std::endl;
	std::cout<< "Ymax: "<< Y_MAX << std::endl;
	std::cout<< "Xmin: "<< X_MIN << std::endl;
	std::cout<< "Xmax: "<< X_MAX << std::endl;
	std::cout<< "NY:   "<< N2<< std::endl;
	std::cout<< std::endl;
	std::cout<< "Tmin: "<< T_MIN<< std::endl;
	std::cout<< "Tmax: "<< T_MAX<< std::endl;
	std::cout<< "NT:   "<< N3<< std::endl;
	std::cout<< std::endl;
	if(HankelMode==0){
		std::cout<< "MinFactor: "<< MinFactor<< std::endl;
		std::cout<< "MaxFactor: "<< MaxFactor<< std::endl;
	}


};

void Dipole::make_test_output(){
	dump_momentum_Dipole(1);
	dump_momentum_Dipole(2);
	dump_momentum_Dipole(4);
	dump_transformed_norm(1);
	dump_transformed_norm(2);
	dump_transformed_norm(4);
}

void Dipole::dump_transformed_norm(double T){
	std::ofstream dip_A_f;
	std::ofstream dip_F_f;
	std::ostringstream path_to_file_A;
	std::ostringstream path_to_file_F;
	path_to_file_A << path_to_set<< "/Dipole_0_adj_T_"<< T<< "_trans_norm.dat" ;
	path_to_file_F << path_to_set<< "/Dipole_0_fund_T_"<< T<< "_trans_norm.dat" ;
	dip_A_f.open(path_to_file_A.str());
	dip_F_f.open(path_to_file_F.str());
	int Nx_t = 100;
	double Y_MIN_OUT = Y_MIN *1.01;
	double Y_MAX_OUT = Y_MAX *0.99;
	for (int iy = 0; iy <Nx_t; iy++) {
 		double y_t = Y_MIN_OUT + iy*(Y_MAX_OUT-Y_MIN_OUT)/(Nx_t-1.);
		double x_t = get_X(y_t);
		dip_A_f<< x_t << "\t" <<AdjointDipole_k(x_t,0, T) << std::endl  ;
		dip_F_f<< x_t << "\t" <<FundamentalDipole_k(x_t,0, T) << std::endl  ;
	}
	dip_A_f.close();
	dip_F_f.close();
}

void Dipole::dump_momentum_Dipole(double T){
	int Nx_t = 10;
	double x[ Nx_t ];
	double y[ Nx_t ];
	double Y_MIN_OUT = Y_MIN *1.01;
	double Y_MAX_OUT = Y_MAX *0.99;
	for (int iy = 0; iy < Nx_t; iy++) {
 		y[iy] = Y_MIN_OUT + iy*(Y_MAX_OUT-Y_MIN_OUT)/(Nx_t-1.);
		x[iy] = get_X(y[iy]);
	}

	std::ofstream dip_A_f;
	std::ofstream dip_F_f;
	std::ostringstream path_to_file_A;
	std::ostringstream path_to_file_F;
	path_to_file_A << path_to_set<< "/Dipole_k_adj_T_"<< T<< ".dat" ;
	path_to_file_F << path_to_set<< "/Dipole_k_fund_T_"<< T<< ".dat" ;
	dip_A_f.open(path_to_file_A.str());
	dip_F_f.open(path_to_file_F.str());

	dip_A_f<< "# k";
	dip_F_f<< "# k";

	for (int ix = 0; ix < Nx_t; ix++) {
		dip_A_f<< "\tQ2(x="<<x[ix]<<")\t D(x="<<x[ix]<<")";
		dip_F_f<< "\tQ2(x="<<x[ix]<<")\t D(x="<<x[ix]<<")";
	}
	dip_A_f<< std::endl;
	dip_F_f<< std::endl;

	for (int ik = 0; ik < N1; ik++) {
		double q_t = ik*DK;
		double k_t = K_MIN*exp(q_t );

		dip_A_f<< k_t ;
		dip_F_f<< k_t ;

		for (int j = 0; j < Nx_t; j++) {
			double Q2A_t= gsl_spline2d_eval(Q2A_spl,T,y[j],xaccQ2A, yaccQ2A);
			double Q2F_t= gsl_spline2d_eval(Q2F_spl,T,y[j],xaccQ2F, yaccQ2F);

			dip_A_f<< "\t"<< Q2A_t<< "\t"<< AdjointDipole_k(x[j], k_t, T);
			dip_F_f<< "\t"<< Q2F_t<< "\t"<< FundamentalDipole_k(x[j], k_t, T) ;
		}
		dip_A_f<< std::endl ;
		dip_F_f<< std::endl ;
	}
	dip_A_f.close();
	dip_F_f.close();
}


void Dipole::quadratic_extrapolation(double x1,double f1,double x2,double f2,double x3,double f3,double& a,double& b,double& c){
	double DET = (x1-x2)*(x1-x3)*(x2-x3);
	a = ( f3*(x1 - x2) + f1*(x2 - x3) + f2*(-x1 + x3) ) /DET;
	b = ( f3*(-pow(x1,2) + pow(x2,2)) + f2*(pow(x1,2) - pow(x3,2)) + f1*(-pow(x2,2) + pow(x3,2)) ) /DET;
	c = ( f3*x1*(x1 - x2)*x2 + f1*x2*(x2 - x3)*x3 + f2*x1*x3*(-x1 + x3) )/DET;
}

bool Dipole::is_in_range(double x, double x1,double x2){
	if(x1<x2){
		return ( (x>x1) && (x<x2) );
	}
	else if (x1>x2){
		return ( (x>x2) && (x<x1) );
	}
	return false;
}

void Dipole::get_Q2(){
	/* This function produces a 2D array for the Qbar2(x,T) quantity.
	This is obtained by assuiming a GBW form, and solving for r when
	D_{rep}(rbar,x)=exp[-0.25]. This list is used then to enforce
	geometrical scaling for x>xscaling */

	double Q2F[N3*N2];
	double Q2A[N3*N2];
	double peak_point = exp(-1./4.);

	std::ofstream Q2_scaling_f;
  std::ostringstream Q2_scaling_name;
  Q2_scaling_name << path_to_set <<"/Q2_scaling.dat";
	Q2_scaling_f.open(Q2_scaling_name.str());
	if(DipVerbose){
 	 std::cout<< Q2_scaling_name.str() << std::endl;}

	//Output to plot out the Scaling Q
	double T_t, x_t;
	int ix;
	for (int iT = 0; iT < N3; iT++) {
		T_t=T_MIN + iT * DT;
		for (int iy = 0; iy < N2; iy++) {
			x_t=get_X(Y_grid[iy]);

			// Find  Q2 for both representations
			bool is_adj_found=false;
			bool is_fund_found=false;

			for (int iu = 0; iu < IPsat_pars::NR - 1; iu++) {
				double r1 = IPsat_pars::rmin *exp(iu*IPsat_pars::du);
				double r2 = IPsat_pars::rmin *exp( (iu+1)*IPsat_pars::du);

				double dip_fund_1 = FundamentalDipole(x_t,r1, T_t);
				double dip_adj_1 = AdjointDipole(x_t,r1, T_t);

				double dip_fund_2 = FundamentalDipole(x_t,r2, T_t);
				double dip_adj_2 = AdjointDipole(x_t,r2, T_t);
				// fundamental
				if(T_t ==0){
					Q2F[N3*iy+iT] =0;
					Q2A[N3*iy+iT] =0;
					break;
				}
				else{
					if(is_in_range(peak_point,dip_fund_1,dip_fund_2)){
						is_fund_found=true;
						double rbar_fund = r1 + (peak_point -dip_fund_1)*(r2-r1)/(dip_fund_2-dip_fund_1);
						Q2F[N3*iy+iT] =pow(rbar_fund * gen_pars::fm_to_GeVm1 ,-2.);
					}
					// adjoint
					if(is_in_range(peak_point,dip_adj_1,dip_adj_2)){
						is_adj_found=true;
						double rbar_adj = r1 + (peak_point -dip_adj_1)*(r2-r1)/(dip_adj_2-dip_adj_1);
						Q2A[N3*iy+iT] =pow(rbar_adj * gen_pars::fm_to_GeVm1 ,-2.);
					}
					if(is_adj_found && is_fund_found){break;}
				}
			}
		}
	}

	//Output to plot out the Scaling Q
	for (int iT = 0; iT < N3; iT++) {
		T_t=T_MIN + iT * DT;
		for (int iy = 0; iy < N2; iy++) {
			ix= N2-1-iy;
			x_t=get_X(Y_grid[ix]);
			Q2_scaling_f<< T_t << "\t"<< x_t<< "\t"<< Q2F[N3*ix+iT]<< "\t"<< Q2A[N3*ix+iT] << std::endl;
		}
	}

	Q2_scaling_f.close();

	// Interpolate Q's
	const gsl_interp2d_type *TA= gsl_interp2d_bicubic;
	Q2A_spl = gsl_spline2d_alloc(TA, N3, N2);
	xaccQ2A = gsl_interp_accel_alloc();
	yaccQ2A = gsl_interp_accel_alloc();
	gsl_spline2d_init(Q2A_spl, T_grid, Y_grid, Q2A, N3, N2);

	const gsl_interp2d_type *TF= gsl_interp2d_bicubic;
	Q2F_spl = gsl_spline2d_alloc(TF, N3, N2);
	xaccQ2F = gsl_interp_accel_alloc();
	yaccQ2F = gsl_interp_accel_alloc();
	gsl_spline2d_init(Q2F_spl, T_grid, Y_grid, Q2F, N3, N2);


}