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
#include <map>
#include <array>
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
    // static const int VegasInitCalls=1e5;
    bool is_validX(double x, double xmin, double xmax){return (x>xmin) && (x<xmax); }


    // q^2 + k^2 + 2qk cos(phi) cancels catastrophically as phi -> pi with
    // q ~ k: the sum can come out negative and sqrt returns NaN. This form is
    // algebraically identical and stays non-negative by construction.
    inline double P(double q, double k, double phi ){
      double h = cos(0.5*phi);
      return sqrt((q-k)*(q-k) + 4.0*q*k*h*h);
    }

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

    static int GluonEnergyDensIntegrand_non_log(const int *ndim, const double xx[],const int *ncomp, double ff[], void *params) {

          struct GluonParsIPSat * pars = (struct GluonParsIPSat *) params;

          double phi_t=(phi_max-phi_min)*xx[2]+phi_min;
          double JacCuba=pow(IPsat_pars::KTMAX-IPsat_pars::KTMIN,2.);
            
          double  q_t=(IPsat_pars::KTMAX-IPsat_pars::KTMIN)*xx[0]+ IPsat_pars::KTMIN;
          double  k_t=(IPsat_pars::KTMAX-IPsat_pars::KTMIN)*xx[1]+ IPsat_pars::KTMIN;
        
          
          double p_t=  P(q_t,k_t,phi_t); 
          double cutP2 = pow(pars->p_reg,2.) + p_t*p_t  ;
          // double p_t = sqrt(pow(P(q_t,k_t,phi_t),2.) + pow(pcut,2.));

          if(p_t==0){ff[0]=0.0;}
          else{
            double x1=p_t*exp( pars->y)/pars->sqrts;
            double x2=p_t*exp(-pars->y)/pars->sqrts;

              double D1 =(pars->dip)->AdjointDipole_k(x1,q_t,pars->T1);
              double D2 =(pars->dip)->AdjointDipole_k(x2,k_t,pars->T2);
              
              if(D1<0){D1=0.0;}
              if(D2<0){D2=0.0;}
              ff[0]= 2*M_PI * JacCuba* ( gen_pars::pref_glue/ ( 2.0*M_PI) ) * ( p_t * pow( q_t,3.)* pow(k_t,3.)/ cutP2) * D1 * D2 *gen_pars::GeV2_to_fmm2 ;    
          }
          return 0.0 ;
        }

      static int GluonEnergyDensIntegrand(const int *ndim, const double xx[],const int *ncomp, double ff[], void *params) {

          struct GluonParsIPSat * pars = (struct GluonParsIPSat *) params;

          double phi_t=(pars->phi_hi - pars->phi_lo)*xx[2] + pars->phi_lo;
          double UMIN = 0; double UMAX=log(IPsat_pars::KTMAX/IPsat_pars::KTMIN);
          double JacCuba=pow(UMAX-UMIN,2.);
          double u=(UMAX-UMIN)*xx[0]+ UMIN;
          double w=(UMAX-UMIN)*xx[1]+ UMIN;
          double q_t= IPsat_pars::KTMIN*exp(u);
          double k_t= IPsat_pars::KTMIN*exp(w);
      
          double p_t=  P(q_t,k_t,phi_t); 
          double cutP2 = pow(pars->p_reg,2.) + p_t*p_t  ;

          double x1=p_t*exp( pars->y)/pars->sqrts;
          double x2=p_t*exp(-pars->y)/pars->sqrts;

          double D1 =(pars->dip)->AdjointDipole_k(x1,q_t,pars->T1);
          double D2 =(pars->dip)->AdjointDipole_k(x2,k_t,pars->T2);

          if(p_t==0){ff[0]=0.0;}
          else{
            if(D1<0){D1=0.0;}
            if(D2<0){D2=0.0;}
            // The phi jacobian is the width of the interval actually integrated over.
            ff[0] = (pars->phi_hi - pars->phi_lo) * JacCuba* ( gen_pars::pref_glue/ ( 2.0*M_PI) ) * ( p_t * pow( q_t,4.)* pow(k_t,4.)/ cutP2) * D1 * D2 *gen_pars::GeV2_to_fmm2 ;
          }
          return 0.0 ;
        }

    // Same physics as GluonEnergyDensIntegrand, integrated over (q,p) instead
    // of (q,k). Momentum conservation is p = q + k, so trading k for p at fixed
    // q is a translation with unit Jacobian: d2q d2k = d2q d2p exactly. In polar
    // form the measure then carries p dp, which cancels the 1/p from the gluon
    // propagator, so nothing here diverges as p -> 0. That removes the reason
    // for the +-0.99pi azimuthal cut and lets phi span the full circle.
    static int GluonEnergyDensIntegrand_pvar(const int *ndim, const double xx[],const int *ncomp, double ff[], void *params) {

          struct GluonParsIPSat * pars = (struct GluonParsIPSat *) params;

          // q log-mapped over [KTMIN,KTMAX], as in the (q,k) version.
          double UMAX=log(IPsat_pars::KTMAX/IPsat_pars::KTMIN);
          double u=UMAX*xx[0];
          double q_t= IPsat_pars::KTMIN*exp(u);

          // p is the transverse momentum of the produced gluon, an integration
          // variable here rather than the derived |q+k|. Log-mapped for the same
          // reason q is: the integrand peaks near p ~ q, where k = |p-q| goes
          // small and the dipole is largest, and q itself spans five decades, so
          // a linear grid would spend nearly all its samples above 1 GeV.
          //
          // The lower bound is converged, not chosen: dropping it to 1e-6 moves
          // the table by 0.008% on average (0.2% worst). Raising it does matter -
          // 1e-3 costs 0.34% on average and 3.6% at eta=10, 1e-2 costs 5.8% and
          // 35% - so it cannot be tied to KTMIN, and there is nothing to gain by
          // making it configurable.
          double PMIN_p = 1.0e-4;
          double PMAX_p = 2.0*IPsat_pars::KTMAX;
          double VMAX = log(PMAX_p/PMIN_p);
          double v = VMAX*xx[1];
          double p_t = PMIN_p*exp(v);

          // phi is the angle between q and p, and spans the full circle.
          double phi_t = 2.0*M_PI*xx[2] - M_PI;

          double k_t = sqrt(q_t*q_t + p_t*p_t - 2.0*q_t*p_t*cos(phi_t));
          double cutP2 = pow(pars->p_reg,2.) + p_t*p_t;

          double x1=p_t*exp( pars->y)/pars->sqrts;
          double x2=p_t*exp(-pars->y)/pars->sqrts;

          double D1 =(pars->dip)->AdjointDipole_k(x1,q_t,pars->T1);
          double D2 =(pars->dip)->AdjointDipole_k(x2,k_t,pars->T2);

          if(p_t==0){ff[0]=0.0;}
          else{
            if(D1<0){D1=0.0;}
            if(D2<0){D2=0.0;}
            ff[0] = 2*M_PI * UMAX * VMAX * ( gen_pars::pref_glue/ ( 2.0*M_PI) )
                    * ( pow(q_t,4.) * pow(k_t,2.) * (pow(p_t,3.)/cutP2) )
                    * D1 * D2 * gen_pars::GeV2_to_fmm2 ;
          }
          return 0.0 ;
        }

    double GluonEnergyDens_pvar(GluonParsIPSat *  Variables){
      double sqrts_t= Variables->sqrts;
      double y_t= Variables->y;
      double kmin= 0;
      double kmax= sqrts_t*exp(-fabs(y_t));
      double error, result ;
      bool is_valid = (kmin<kmax);
      is_valid = is_valid && ( Variables->T1>0 && Variables->T2>0 );
      if(is_valid){
        Routines::make_cuhre_1C(3, GluonEnergyDensIntegrand_pvar,Variables,result,error);
      }
      else{
        return 0;
      }
      return result;
    }

    double GluonEnergyDens(GluonParsIPSat *  Variables){
      double sqrts_t= Variables->sqrts;
      double y_t= Variables->y;
      //
      // double xmin_dip = (Variables->dip)->get_XMIN();
      // double xmax_dip = 1;

      double kmin= 0;
      double kmax= sqrts_t*exp(-fabs(y_t));

      double error=0;
      double result = 0 ;
      //
      bool is_valid = (kmin<kmax);
      is_valid = is_valid && ( Variables->T1>0 && Variables->T2>0 );
      if ( is_valid && std::fabs(y_t) < std::fabs(Variables->dip->getMaxEta()) ){
          double x1l = Variables->dip->get_KMIN()*exp(y_t)/sqrts_t;
          double x1h = Variables->dip->get_KMAX()*exp(y_t)/sqrts_t;
          double x2l = Variables->dip->get_KMIN()*exp(-y_t)/sqrts_t;
          double x2h = Variables->dip->get_KMAX()*exp(-y_t)/sqrts_t;
          bool is_x1_within_range = (Variables->dip->get_XMIN() <x1h) && (Variables->dip->get_XMAX() > x1l);
          bool is_x2_within_range = (Variables->dip->get_XMIN() <x2h) && (Variables->dip->get_XMAX() > x2l);
          // std::cout<< x1l << "\t"<< x1h << "\t"<< x2l << "\t"<< x2h << "\n";
          if (is_x1_within_range && is_x2_within_range) {Routines::make_cuhre_1C(3, GluonEnergyDensIntegrand,Variables,result,error);} 
          else{return result= 0.0;}
          
      }
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

  // Cuba forks worker processes of its own, defaulting to one per core. The
  // tabulation loop is already threaded over (T1,T2), so leaving that on would
  // fork a set of workers from every thread. Switched off; it is also about
  // 30% faster, since the integrations are small and forking dominates.
  { int ncores = 0, pcores = 10000; cubacores(&ncores, &pcores); }

  config=Config(ConfInput);
  
  p_set= int( config.get_ModelParams(0) );
  xscaling = config.get_ModelParams(1);
  p_reg = config.get_ModelParams(2);
  std::cout<< "[IP-Sat]: Regulator : " << p_reg << std::endl;

  if(config.get_Verbose()> VerboseLevel::None){
    std::cout<< std::endl;
    std::cout<<"[IP-Sat]: Initializing IP-Sat model "<< std::endl;
    std::cout<< std::endl;
  }
  quark_dist = new PDFs(&ConfInput);

  if(config.get_Verbose()> VerboseLevel::None){std::cout<< std::endl;}
  Dip= new Dipole(p_set,config.get_Verbose());
  Dip->set_new_xscaling(xscaling);
  config.set_dT();
  std:: cout<< "The x0_scaling = "<<Dip->get_xscaling()<< std::endl;
  
	//config.set_TMax(gen_pars::TMax);
	//config.set_TMin(gen_pars::TMin);
	//config.set_NT(gen_pars::NT);
	
} 

IPSat::~IPSat(){
  delete Dip;
  delete quark_dist;
}

void IPSat::MakeTable(std::string path_to_set){			
	// MAKES ith SET DIRECTORY. Uncomment when config dump is done
	SETPATH = path_to_set;
  fs::create_directories(SETPATH);
	// Write new config to setpath
	config.set_dump(SETPATH);
  
	if(config.get_Verbose()> VerboseLevel::None){std::cout<<"[IP-Sat]: New config written to "<<SETPATH  << std::endl;}
  if(config.get_Verbose()> VerboseLevel::None){std::cout<<"[IP-Sat]: Tabulating conserved charges in the IP-Sat model framework"<<SETPATH  << std::endl;}
  make_gluon_energy();
	if(config.get_Verbose()> VerboseLevel::None){std::cout<<"\n[IP-Sat]:Gluon Energy written to"<<SETPATH  << std::endl;}

  make_baryon_stopping_all();
	if(config.get_Verbose()> VerboseLevel::None){std::cout<<"\n[IP-Sat]:(All) Quark elements written to"<<SETPATH  << std::endl;}
}


// Where the azimuthal integration is split into a bulk and a tail piece (see
// the table loop below).
static const double PHI_SPLIT = M_PI*0.99;

// Coordinates quantised so a point can be looked up by value rather than by
// its position in the file. Six significant digits is what the tables carry,
// which supports any dT above 1e-4.
static std::array<long,3> gluon_key(double y, double T1, double T2){
	return { llround(y*1e4), llround(T1*1e4), llround(T2*1e4) };
}

void IPSat::make_gluon_energy(){

	std::ofstream density_f;
  std::ostringstream densityname;
  densityname << SETPATH <<"/"<< gluon_energy_table_name ;

	// Points an earlier build already computed. Keyed on coordinates, so the
	// table they came from may be on a different grid - a narrower TMax, a
	// coarser dT, or the same grid left short by a job that ran out of time.
	std::map<std::array<long,3>, double> reused_points;
	if(config.is_reuse()){
		std::ifstream seed(config.get_reuse_table());
		if(seed.is_open()){
			double y,t1,t2,v;
			while(seed >> y >> t1 >> t2 >> v){ reused_points[gluon_key(y,t1,t2)] = v; }
			std::cout<<"[ IPSat ] ReuseTable: "<<reused_points.size()<<" points from "
			         <<config.get_reuse_table()<<std::endl;
		}
		else{
			std::cout<<"[ IPSat ] ReuseTable: nothing at "<<config.get_reuse_table()
			         <<", computing all points"<<std::endl;
		}
	}
  if(config.get_Verbose()==VerboseLevel::Dynamic){ std::cout<< "Total (y,T1,T2)                                                          Local (T1,T2)"<< std::endl;  ;}

	GluonParsIPSat parameters;
	parameters.sqrts= config.get_collEnergy();
  parameters.dip= Dip;
  parameters.p_reg=p_reg;
  
	density_f.open(densityname.str());
	const int NT = config.get_NT();

	// Cuba initialises itself on the first call - reading CUBACORES, setting up
	// its worker state. Do that once here, serially, rather than from whichever
	// thread happens to arrive first.
	{
		GluonParsIPSat warm = parameters;
		warm.y = config.get_ETAMIN();
		warm.T1 = warm.T2 = (config.get_dT() > 0)? config.get_dT() : 1.0;
		warm.phi_lo = -PHI_SPLIT; warm.phi_hi = PHI_SPLIT;
		IPSatYields::GluonEnergyDens(&warm);
	}

	for (int iy = 0; iy < config.get_NETA(); iy++) {
		double y_t = iy*config.get_dETA() + config.get_ETAMIN();
		std::vector<double> plane(long(NT)*NT);

		// Every (T1,T2) point is an independent Cuhre call over a table that is
		// only read, so the rapidity plane parallelises directly. Each thread
		// takes its own copy of the parameter struct. Results go to an array and are
		// written afterwards, so the file does not depend on thread order.
		#pragma omp parallel for collapse(2) schedule(dynamic)
		for (int i1 = 0; i1 < NT; i1++) {
			for (int i2 = 0; i2 < NT; i2++) {
				double T1_t = i1*config.get_dT() + config.get_TMin();
				double T2_t = i2*config.get_dT() + config.get_TMin();

				GluonParsIPSat local = parameters;
				local.y = y_t; local.T1 = T1_t; local.T2 = T2_t;

				double res;
				auto found = reused_points.find(gluon_key(y_t, T1_t, T2_t));
				if(found != reused_points.end()){ res = found->second; }
				else{
					// Integrated in two pieces. The integrand carries a spike of
					// width ~p_reg where phi -> pi and k -> q. Over the full circle
					// that spike is about 1e-5 of the volume, small enough that
					// Cuhre's samples miss it entirely - and a rule that never
					// samples a feature reports convergence rather than error.
					// Split at 0.99*pi the spike is a thirtieth of the tail
					// interval and the first subdivision already resolves it.
					local.phi_lo = -PHI_SPLIT; local.phi_hi = PHI_SPLIT;
					double bulk = IPSatYields::GluonEnergyDens(&local);
					// The integrand is even in phi, so one side of the tail doubled.
					local.phi_lo = PHI_SPLIT; local.phi_hi = M_PI;
					double tail = IPSatYields::GluonEnergyDens(&local);
					res = bulk + 2.0*tail;
				}
				plane[long(i1)*NT + i2] = res;
			}
		}

		for (int i1 = 0; i1 < NT; i1++) {
			double T1_t = i1*config.get_dT() + config.get_TMin();
			for (int i2 = 0; i2 < NT; i2++) {
				double T2_t = i2*config.get_dT() + config.get_TMin();
				density_f<<y_t << "\t"<<T1_t << "\t"<<T2_t << "\t"<<plane[long(i1)*NT + i2]<< "\n";
			}
		}
		if(config.get_Verbose()==VerboseLevel::Dynamic){printProgress(double(iy)/double(config.get_NETA()));}
	}

	if(config.get_Verbose()==VerboseLevel::Dynamic){printProgress(1);}
	Routines::report_cuhre_stats("gluon energy table");
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


	density_f.open(densityname.str());
	for (int iy = 0; iy < config.get_NETA(); iy++) {
		double y_t = iy*config.get_dETA() + config.get_ETAMIN();

		for (int i1 = 0; i1 < config.get_NT(); i1++) {
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
				gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12q,&err12q);
  

				parameters.quark_id = aqid;
        F.params = &parameters;
				gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12aq,&err12aq);

			}

      if(is_null21){
				res21q=0;err21q=0;res21aq=0;err21aq=0;
			}
			else{
				parameters.y = -y_t ;
				parameters.T = T_t;
				parameters.quark_id = qid;
				F.params = &parameters;
				gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21q,&err21q);

				parameters.quark_id = aqid;
        F.params = &parameters;
				gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21aq,&err21aq);
			}

      density_f<<y_t<< "\t"<<T_t<< "\t"<<res12q<< "\t"<<res12aq<< "\t"<<res21q<< "\t"<<res21aq<< "\n";
		}
		density_f<< std::endl;
		if(config.get_Verbose()==VerboseLevel::Dynamic){
			if(remainder(iy,gen_pars::skip)==0){double percentage_done = double(iy)/double(config.get_NETA());printProgress(percentage_done);}
		}

	}
	gsl_integration_workspace_free (w);
 	density_f.close();
	if(config.get_Verbose()==VerboseLevel::Dynamic){printProgress(1);std::cout<<std::endl;}
}

void IPSat::make_baryon_stopping_all(){
	make_baryon_stopping(0, QuarkID::u, QuarkID::ubar);
	make_baryon_stopping(1, QuarkID::u, QuarkID::ubar);
	if(config.get_Verbose()>VerboseLevel::None){std::cout<<"[IP-Sat]: U-Quark elements written \n " << std::endl;}

	make_baryon_stopping(0, QuarkID::d, QuarkID::dbar);
	make_baryon_stopping(1, QuarkID::d, QuarkID::dbar);
	if(config.get_Verbose()>VerboseLevel::None){std::cout<<"[IP-Sat]: D-Quark elements written \n " << std::endl;}

	make_baryon_stopping(0, QuarkID::s, QuarkID::sbar);
	make_baryon_stopping(1, QuarkID::s, QuarkID::sbar);
	if(config.get_Verbose()>VerboseLevel::None){std::cout<<"[IP-Sat]: S-Quark elements written \n " << std::endl;}


}

bool IPSat::check_if_zero_F(double y,double T){

  double x1=gen_pars::PMIN*exp(y)/config.get_collEnergy();

	// No nuclear matter means no charge; any T > 0 carries a small but physical
	// contribution and is computed.
	bool is_null = (T <= 0.0) ;
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
  if(config.get_Verbose()>VerboseLevel::None){std::cout<< "Writing dump for " << densityname.str() <<std::endl;}

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

//Dynamical
  std::cout<<"sqrts="<< parametersG.sqrts<< std::endl;
  std::cout<<"y="<< parametersG.y<< std::endl;
  std::cout<<"qt="<< parametersG.p_reg<< std::endl;


//Geometrical
  std::cout<<"T1="<< parametersG.T1<< std::endl;
  std::cout<<"T2="<< parametersG.T2<< std::endl;


	for (int iy = 0; iy < config.get_NETA(); iy++) {
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
    // if(config.get_Verbose()==VerboseLevel::Dynamic){std::cout<<parametersG.y << "\t"<<resEG<< "\t"<<resEQ<< "\t"<<UDens<< "\t"<<DDens<< "\n";}
	}
  density_f.close();
  if(config.get_Verbose()>VerboseLevel::None){std::cout<<"-> Done"<<std::endl;}

}
