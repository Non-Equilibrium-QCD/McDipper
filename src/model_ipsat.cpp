/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */

#include <ctime>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <complex>
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <set>
#include <vector>
#include <functional>
#include <math.h>
#include <sstream>
#include <sys/stat.h>
#include <filesystem>

// GSL MATH ROUTINES //
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>

// GSL INTEGRATION ROUTINES //
#include <gsl/gsl_integration.h>

// GSL INTERPOLATION ROUTINES //
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "include/model_ipsat.h"
#include "include/params_gen.h"

#include "include/routines.h"

// BESSEL FUNCTIONS //
#ifndef BesselJ0
#define BesselJ0(x) jn(0,x)
#endif

namespace fs = std::filesystem;

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


namespace IPSatYields{

    int NPTS=3;
    double PTS[3]={0.0,M_PI,2.*M_PI};


    double Q02 =2.0;

    double BesselTol=1e-4;
    double EpsRel=1e-2;
    int BessQuad =4;
    int GenQuad =4;
    double EpsRelPhi=1e-4;

    int NPhi_m = 10;
    double hPhi= 2*M_PI/(3*NPhi_m) ;
    double epsilon =0.001;

    ///Cuhre implementation of gluon energy
    double phi_max = M_PI*0.99;
    double phi_min = -M_PI*0.99;
    double pcut = 0.1;
    int DIM_GLUE= 3;

    static const int NumberOfIntegrationPoints=64000;
    static const int VegasInitCalls=1e5;

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

    double DipoleFunFT(double p,DipFT * Variables){
      Variables->p=p;
      // CALCULATE qg TMDS //
      gsl_function GSL_FundDipole_Integrand;
      GSL_FundDipole_Integrand.function=&DipoleFunFT_Integrand;
      GSL_FundDipole_Integrand.params=Variables;
      double DipMom=BesselIntegrate(p,GSL_FundDipole_Integrand);
      return DipMom;
    }

    double DipoleAdjFT(double p,DipFT * Variables){
      Variables->p=p;
      // CALCULATE qg TMDS //
      gsl_function GSL_AdjDipole_Integrand;
      GSL_AdjDipole_Integrand.function=&DipoleAdjFT_Integrand;
      GSL_AdjDipole_Integrand.params=Variables;
      double DipMom=BesselIntegrate(p,GSL_AdjDipole_Integrand);
      return DipMom;
    }

    bool is_validX(double x, double xmin, double xmax){return (x>xmin) && (x<xmax); }


    inline double P(double q, double k, double phi ){ return sqrt(q*q + k*k + 2*q*k*cos(phi)); }

    bool check_if_zero(double y,double sqrts,double kmin,double kmax, double xmin, double xmax){

    	double x1min=kmin*exp( y)/sqrts;
    	double x2min=kmin*exp(-y)/sqrts;

      double x1max=kmax*exp( y)/sqrts;
    	double x2max=kmax*exp(-y)/sqrts;

      //check x1 is outside range if (x1min > xmax) then the whole region is above the dipole's area of validity
      // if instead (x1max < xmin), everything is extreme low-x, and once again we miss the region of valdity
    	bool is_null = (x1min > xmax) || (x1max < xmin);
      //check x2 is outside range
      is_null =  is_null || ((x2min > xmax) || (x2max < xmin)) ;
      // IF ANY of both are completely out of range, the integration gives 0, we are preventing this.
    	return !is_null;
      // We return the conditional validity is_valid=not(is_null) to comply with the structure above.
    }

    static int GluonEnergyDensIntegrand(const int *ndim, const double xx[],const int *ncomp, double ff[], void *params) {

          struct GluonParsIPSat * pars = (struct GluonParsIPSat *) params;

          double phi_t=(2.0 * M_PI)*xx[0];
      
          double q_t=(IPsat_pars::KTMAX-IPsat_pars::KTMIN)*xx[1]+ IPsat_pars::KTMIN;
          double k_t=(IPsat_pars::KTMAX-IPsat_pars::KTMIN)*xx[2]+ IPsat_pars::KTMIN;
          
          double p_t= pow( pow(pcut,2.) + pow(P(q_t,k_t,phi_t),2.) ,1/2.)  ;

          if(p_t==0){ff[0]=0;}
          else{
            double x1=p_t*exp( pars->y)/pars->sqrts;
            double x2=p_t*exp(-pars->y)/pars->sqrts;

              double D1 =(pars->dip)->AdjointDipole_k(x1,q_t,pars->T1);
              double D2 =(pars->dip)->AdjointDipole_k(x2,k_t,pars->T2);
              double result= 2*M_PI * pow(IPsat_pars::KTMAX-IPsat_pars::KTMIN,2.) * ( gen_pars::pref_glue/ ( 2.0*M_PI) ) * ( pow( q_t,3.)* pow(k_t,3.)/p_t) * D1*D2 *gen_pars::GeV2_to_fmm2 ;
              if(result>=0){
                ff[0] = result;
              }
              else{
                ff[0] = 0;
              }

          }
          return 0 ;
        }

        
    double GluonEnergyDens(GluonParsIPSat *  Variables){
      double sqrts_t= Variables->sqrts;
      double y_t= Variables->y;
      //
      double xmin_dip = (Variables->dip)->get_XMIN();
      double xmax_dip = (Variables->dip)->get_XMAX();

      double kmin= IPsat_pars::KTMIN;
      double kmax= IPsat_pars::KTMAX;

      double error, result ;
      //
      bool is_valid = (kmin<kmax);
      is_valid = is_valid && ( Variables->T1>0 && Variables->T2>0 );
      is_valid = is_valid && check_if_zero(y_t,sqrts_t,kmin,kmax,xmin_dip, xmax_dip);
      //
      if(is_valid){
        Routines::make_cuhre_1C(3, GluonEnergyDensIntegrand,Variables,result,error);
      }
      else{return 0;}
      return result;
    }


    double quark_n12(double p, void * params){
      QuarkParsIPSat * pars = ( (struct QuarkParsIPSat *) params);
      double x1=p*exp( pars->y)/pars->sqrts;
      double x2=p*exp(-pars->y)/pars->sqrts;


      double Dip_t = (pars->dip)->FundamentalDipole_k(x2,p,pars->T);
      double result;
      //Quark contribution
      if(x2>(pars->dip)->get_XMAX() || x2<(pars->dip)->get_XMIN() ){return 0;}
      else{
        result = pow(p,1+pars->k)* ( (pars->partons)->xfpart(pars->quark_id,x1,p) ) * Dip_t* pow(2*M_PI,-1);
        if(result>0){return result;}
        else{return 0;}
      }
    }

}

IPSat::IPSat(){}
IPSat::IPSat(Config ConfInput){
  gsl_set_error_handler_off();
  config=Config(ConfInput);
  p_set= int( config.get_ModelParams(0) );
  xscaling = config.get_ModelParams(1);

  if(config.get_Verbose()){
    std::cout<< std::endl;
    std::cout<<"                     ----->  Initializing IPSat model  <----- "<< std::endl;
    std::cout<< std::endl;
  }
  quark_dist = new PDFs(&ConfInput);

  if(config.get_Verbose()){std::cout<< std::endl;}
  Dip= new Dipole(p_set,config.get_Verbose());
  Dip->set_new_xscaling(xscaling);
  
	config.set_TMax(gen_pars::TMax);
	config.set_TMin(gen_pars::TMin);
	config.set_NT(gen_pars::NT);
	config.set_dT();

  std:: cout<< "The x0_scaling = "<<Dip->get_xscaling()<< std::endl;
} 

IPSat::~IPSat(){}

void IPSat::MakeTable(std::string path_to_set){			
	// MAKES ith SET DIRECTORY. Uncomment when config dump is done
	SETPATH = path_to_set;
  fs::create_directories(SETPATH);
	// Write new config to setpath
	config.set_dump(SETPATH);
	if(config.get_Verbose()){std::cout<<"New config written to "<<SETPATH  << std::endl;}

  if(config.get_Verbose()){std::cout<<"--> Tabulating conserved charges in the IP-Sat model framework"<<SETPATH  << std::endl;}
  make_gluon_energy();
	if(config.get_Verbose()){std::cout<<"\nGluon Energy written to"<<SETPATH  << std::endl;}

  make_baryon_stopping_all();
	if(config.get_Verbose()){std::cout<<"\n(All) Quark elements written to"<<SETPATH  << std::endl;}
}


void IPSat::make_gluon_energy(){

	std::ofstream density_f;
  std::ostringstream densityname;
  densityname << SETPATH <<"/"<< gluon_energy_table_name ;
  if(config.get_Verbose()){ std::cout<< "Total (y,T1,T2)                                                          Local (T1,T2)"<< std::endl;  ;}

	GluonParsIPSat parameters;
	parameters.sqrts= config.get_collEnergy();
  parameters.dip= Dip;
  
	double res=0;
  int counter = 0;
  density_f.open(densityname.str());
	for (size_t iy = 0; iy < config.get_NETA(); iy++) {
		double y_t = iy*config.get_dETA() + config.get_ETAMIN();
		parameters.y =y_t ;
		for (size_t i1 = 0; i1 < config.get_NT(); i1++) {
			double T1_t = i1*config.get_dT() + config.get_TMin();
			parameters.T1 = T1_t;
			for (size_t i2 = 0; i2 < config.get_NT(); i2++) {
				double T2_t = i2*config.get_dT() + config.get_TMin();
				parameters.T2 = T2_t;

        res = IPSatYields::GluonEnergyDens(&parameters);

				density_f<<parameters.y << "\t"<<parameters.T1 << "\t"<<parameters.T2 << "\t"<<res<< "\n";
        if(config.get_Verbose()){
					if(remainder(i1,gen_pars::skip)==0){
						double percentage_done1 = double(iy)/double(config.get_NETA());
						double percentage_done2 = double(i1)/double(config.get_NT());
						printProgress2(percentage_done1,percentage_done2);
          }
				}
        counter++;
			}
		}
	}

	if(config.get_Verbose()){printProgress2(1,1);}
}

void IPSat::make_baryon_stopping(int k, QuarkID qid, QuarkID aqid){
  // gsl_set_error_handler_off();
	std::ofstream density_f;
  std::ostringstream densityname;
  densityname << SETPATH <<"/"<< get_quark_construction_file(k, qid) ;
	QuarkParsIPSat parameters;
	parameters.sqrts= config.get_collEnergy();
	parameters.k=k;
	parameters.partons = quark_dist;
  parameters.dip= Dip;

	gsl_integration_workspace * w= gsl_integration_workspace_alloc ( gen_pars::limit);
	gsl_function F;
	F.function = &IPSatYields::quark_n12;//(double p, void * params);

	double res12q,err12q;
	double res21q,err21q;
	double res12aq,err12aq;
	double res21aq,err21aq;

  int status12q,status12aq,status21q,status21aq;

	density_f.open(densityname.str());
	for (size_t iy = 0; iy < config.get_NETA(); iy++) {
		double y_t = iy*config.get_dETA() + config.get_ETAMIN();

		for (size_t i1 = 0; i1 < config.get_NT(); i1++) {
			double T_t = i1*config.get_dT() + config.get_TMin();

      bool is_null12 = check_if_zero_F(y_t,T_t);
      bool is_null21 = check_if_zero_F(-y_t,T_t);

			if(is_null12){
				res12q=0;err12q=0;res12aq=0;err12aq=0;
			}
			else{
				parameters.y = y_t ;
				parameters.T = T_t;
				parameters.quark_id = qid;
				F.params = &parameters;
				status12q=gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12q,&err12q);
  

				parameters.quark_id = aqid;
        F.params = &parameters;
				status12aq=gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12aq,&err12aq);

			}

      if(is_null21){
				res21q=0;err21q=0;res21aq=0;err21aq=0;
			}
			else{
				parameters.y = -y_t ;
				parameters.T = T_t;
				parameters.quark_id = qid;
				F.params = &parameters;
				status21q=gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21q,&err21q);

				parameters.quark_id = aqid;
        F.params = &parameters;
				status21aq=gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21aq,&err21aq);
			}

      density_f<<y_t<< "\t"<<T_t<< "\t"<<res12q<< "\t"<<res12aq<< "\t"<<res21q<< "\t"<<res21aq<< "\n";
		}
		density_f<< std::endl;
		if(config.get_Verbose()){
			if(remainder(iy,gen_pars::skip)==0){double percentage_done = double(iy)/double(config.get_NETA());printProgress(percentage_done);}
		}

	}
	gsl_integration_workspace_free (w);
 	density_f.close();
	if(config.get_Verbose()){printProgress(1);std::cout<<std::endl;}
}

void IPSat::make_baryon_stopping_all(){
	make_baryon_stopping(0, QuarkID::u, QuarkID::ubar);
	make_baryon_stopping(1, QuarkID::u, QuarkID::ubar);
	if(config.get_Verbose()){std::cout<<"U-Quark elements written \n " << std::endl;}

	make_baryon_stopping(0, QuarkID::d, QuarkID::dbar);
	make_baryon_stopping(1, QuarkID::d, QuarkID::dbar);
	if(config.get_Verbose()){std::cout<<"D-Quark elements written \n " << std::endl;}

	make_baryon_stopping(0, QuarkID::s, QuarkID::sbar);
	make_baryon_stopping(1, QuarkID::s, QuarkID::sbar);
	if(config.get_Verbose()){std::cout<<"S-Quark elements written \n " << std::endl;}


}

bool IPSat::check_if_zero_F(double y,double T){

  double x1=gen_pars::PMIN*exp(y)/config.get_collEnergy();

	bool is_null = (T< gen_pars::T_tolerance) ;
  is_null = is_null || (x1 > gen_pars::XdipMax);


	return is_null;
}

void  IPSat::printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
void IPSat::printProgress2(double percentage1, double percentage2) {
    int val1 = (int) (percentage1 * 100);
    int lpad1 = (int) (percentage1 * PBWIDTH);
    int rpad1 = PBWIDTH - lpad1;

		int val2 = (int) (percentage2 * 100);
    int lpad2 = (int) (percentage2 * PBWIDTH);
    int rpad2 = PBWIDTH - lpad2;

		printf("\r%3d%% [%.*s%*s]", val1, lpad1, PBSTR, rpad1, "");
		printf("   %3d%% [%.*s%*s]", val2, lpad2, PBSTR, rpad2, "");
    fflush(stdout);
}

std::string IPSat::get_quark_construction_file(int k, QuarkID qid){
  std::ofstream quark_name_f;
  std::ostringstream quark_name;

  quark_name << "quark_construction_n"<< k;
  if(qid==QuarkID::u){quark_name <<"_u.dat";}
  else if(qid==QuarkID::d){quark_name <<"_d.dat";}
  else if(qid==QuarkID::s){quark_name <<"_s.dat";}

  return quark_name.str();
}


void IPSat::TestDump(double T1,double T2){
  std::ofstream density_f;
  std::ostringstream densityname;
  densityname << SETPATH <<"/TestDump_T1_"<< T1<<"_T2_"<<T2<< ".dat" ;
  if(config.get_Verbose()){std::cout<< "Writing dump for " << densityname.str() <<std::endl;}

	GluonParsIPSat parametersG;
	parametersG.sqrts= config.get_collEnergy();
  parametersG.dip= Dip;
  parametersG.T1 = T1;
  parametersG.T2 = T2;

  QuarkParsIPSat parametersQ;
	parametersQ.sqrts= config.get_collEnergy();
	parametersQ.partons = quark_dist;
  parametersQ.dip= Dip;


	gsl_integration_workspace * w= gsl_integration_workspace_alloc ( gen_pars::limit);
	gsl_function F;
	F.function = &IPSatYields::quark_n12;

	double res12u,err12u;
	double res21u,err21u;
	double res12ubar,err12ubar;
	double res21ubar,err21ubar;

  double res12d,err12d;
	double res21d,err21d;
	double res12dbar,err12dbar;
	double res21dbar,err21dbar;


	double resEG=0;
  double resEQ=0;
  double UDens=0;
  double DDens=0;
  // int counter = 0;
  density_f.open(densityname.str());
	for (size_t iy = 0; iy < config.get_NETA(); iy++) {
		double y_t = iy*config.get_dETA() + config.get_ETAMIN();
		parametersG.y =y_t ;
    resEG = IPSatYields::GluonEnergyDens(&parametersG);

    bool is_null12 = check_if_zero_F(y_t,T2);
    bool is_null21 = check_if_zero_F(-y_t,T1);

    ////// ENERGY!!!!!!!
    parametersQ.k=1;
    ///////////////////

    if(is_null12){
      res12u=0;err12u=0;res12ubar=0;err12ubar=0;
      res12d=0;err12u=0;res12dbar=0;err12dbar=0;
    }
    else{
      parametersQ.y = y_t;
      parametersQ.T = T2;
      /////////// UP QUARK
      parametersQ.quark_id = QuarkID::u;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12u,&err12u);

      parametersQ.quark_id = QuarkID::ubar;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12ubar,&err12ubar);

      /////////// DOWN QUARK
      parametersQ.quark_id = QuarkID::d;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12d,&err12d);

      parametersQ.quark_id = QuarkID::dbar;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12dbar,&err12dbar);
    }

    if(is_null21){
      res21u=0;err21u=0;res21ubar=0;err21ubar=0;
      res21d=0;err21d=0;res21dbar=0;err21dbar=0;
    }
    else{
      parametersQ.y = -y_t ;
      parametersQ.T = T1;
      /////////// UP QUARK
      parametersQ.quark_id = QuarkID::u;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21u,&err21u);

      parametersQ.quark_id = QuarkID::ubar;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21ubar,&err21ubar);
      /////////// DOWN QUARK
      parametersQ.quark_id = QuarkID::d;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21d,&err21d);

      parametersQ.quark_id = QuarkID::dbar;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21dbar,&err21dbar);
    }

    resEQ= T1*(res12u+res12ubar+ res12d+res12dbar)+T2*(res21u+res21ubar+ res21d+res21dbar);

    ////// Density!!!!!!!
    parametersQ.k=0;
    ///////////////////

    if(is_null12){
      res12u=0;err12u=0;res12ubar=0;err12ubar=0;
      res12d=0;err12u=0;res12dbar=0;err12dbar=0;
    }
    else{
      parametersQ.y = y_t;
      parametersQ.T = T2;
      /////////// UP QUARK
      parametersQ.quark_id = QuarkID::u;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12u,&err12u);

      parametersQ.quark_id = QuarkID::ubar;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12ubar,&err12ubar);

      /////////// DOWN QUARK
      parametersQ.quark_id = QuarkID::d;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12d,&err12d);

      parametersQ.quark_id = QuarkID::dbar;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12dbar,&err12dbar);
    }

    if(is_null21){
      res21u=0;err21u=0;res21ubar=0;err21ubar=0;
      res21d=0;err21d=0;res21dbar=0;err21dbar=0;
    }
    else{
      parametersQ.y = -y_t ;
      parametersQ.T = T1;
      /////////// UP QUARK
      parametersQ.quark_id = QuarkID::u;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21u,&err21u);

      parametersQ.quark_id = QuarkID::ubar;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21ubar,&err21ubar);
      /////////// DOWN QUARK
      parametersQ.quark_id = QuarkID::d;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21d,&err21d);

      parametersQ.quark_id = QuarkID::dbar;
      F.params = &parametersQ;
      gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21dbar,&err21dbar);
    }

    UDens= T1*(res12u-res12ubar)+T2*(res21u-res21ubar);
    DDens= T1*(res12d-res12dbar)+T2*(res21d-res21dbar);

    density_f<<parametersG.y << "\t"<<resEG<< "\t"<<resEQ<< "\t"<<UDens<< "\t"<<DDens<< "\n";
    if(config.get_Verbose()){std::cout<<parametersG.y << "\t"<<resEG<< "\t"<<resEQ<< "\t"<<UDens<< "\t"<<DDens<< "\n";}
	}
  density_f.close();
  if(config.get_Verbose()){std::cout<<"-> Done"<<std::endl;}

}