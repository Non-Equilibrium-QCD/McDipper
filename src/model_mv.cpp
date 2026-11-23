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

#include "include/model_mv.h"
#include "include/params_gen.h"


// BESSEL FUNCTIONS //
#ifndef BesselJ0
#define BesselJ0(x) jn(0,x)
#endif

namespace fs = std::filesystem;

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


namespace MVYields{

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

    // static int GluonEnergyDensIntegrand_non_log(const int *ndim, const double xx[],const int *ncomp, double ff[], void *params) {

    //       struct GluonParsMV * pars = (struct GluonParsMV *) params;

    //       double phi_t=(phi_max-phi_min)*xx[2]+phi_min;
    //       double JacCuba=pow(pars->kmax-pars->kmin,2.);
            
    //       double  q_t=(pars->kmax-pars->kmin)*xx[0]+ pars->kmin;
    //       double  k_t=(pars->kmax-pars->kmin)*xx[1]+ pars->kmin;
      
    //       double p_t=  P(q_t,k_t,phi_t); 
    //       double cutP2 = pow(pars->p_reg,2.) + p_t*p_t  ;
          
    //       // double p_t = sqrt(pow(P(q_t,k_t,phi_t),2.) + pow(pcut,2.));
        
    //       if(p_t==0){ff[0]=0.0;}
    //       else{
    //         double x1=p_t*exp( pars->y)/pars->sqrts;
    //         double x2=p_t*exp(-pars->y)/pars->sqrts;
    //         double D1,D2;
              
              
    //           D1 = pars->dip->AdjointDipole_k(x1, q_t, pars->T1);
    //           D2 = pars->dip->AdjointDipole_k(x2, k_t, pars->T2);
              
    //           // double D1 =(pars->dip)->AdjointDipole_k(x1,q_t,pars->T1);
    //           // double D2 =(pars->dip)->AdjointDipole_k(x2,k_t,pars->T2);
              
    //           if(D1<0){D1=0.0;}
    //           if(D2<0){D2=0.0;}
    //           ff[0]= 2*M_PI * JacCuba* ( gen_pars::pref_glue/ ( 2.0*M_PI) ) * ( p_t * pow( q_t,3.)* pow(k_t,3.)/ cutP2) * D1 * D2 *gen_pars::GeV2_to_fmm2 ;    
    //       }
    //       return 0 ;
    //     }

      static int GluonEnergyDensIntegrand(const int *ndim, const double xx[],const int *ncomp, double ff[], void *params) {

          struct GluonParsMV * pars = (struct GluonParsMV *) params;

          double phi_t=(pars->phi_hi - pars->phi_lo)*xx[2] + pars->phi_lo;
          double UMIN = 0; double UMAX=log(pars->kmax/pars->kmin);

          double JacCuba=pow(UMAX-UMIN,2.);
          double u=(UMAX-UMIN)*xx[0]+ UMIN;
          double w=(UMAX-UMIN)*xx[1]+ UMIN;
          double q_t= pars->kmin*exp(u);
          double k_t= pars->kmin*exp(w);

          // std::cout<< "qs ->"<<q_t<< "\t" << k_t<<std::endl;
      
          double p_t=  P(q_t,k_t,phi_t); 
          double cutP2 = pow(pars->p_reg,2.) + p_t*p_t  ;

          double x1=p_t*exp( pars->y)/pars->sqrts;
          double x2=p_t*exp(-pars->y)/pars->sqrts;

          // std::cout<< "xs ->"<<x1<< "\t" << x2 <<std::endl;

          double D1 =(pars->dip)->AdjointDipole_k(x1,q_t,pars->T1);
          double D2 =(pars->dip)->AdjointDipole_k(x2,k_t,pars->T2);
          // std::cout<< "Ds ->"<<D1<< "\t" << D2 <<std::endl;
        // std::cout<< cutP2 <<"\t"<<pars->p_reg <<"\t"<<p_t <<"\t" << JacCuba <<std::endl;
        // std::cout<< "xs ->"<<x1<< "\t" << x2 <<std::endl;
          if(p_t==0){ff[0]=0.0;}
          else{
            // if(!std::isfinite(D1) ){ D1 = 0.0; std::cout<< "NaN!" << std::endl;exit(0);}
            // if(!std::isfinite(D2)){ D2 = 0.0; std::cout<< "NaN!" << std::endl;exit(0);}
            if(D1<0){D1=0.0;}
            if(D2<0){D2=0.0;}
            // std::cerr << D1 <<  "\t" << D2 << "\t" << std::endl
            // The phi jacobian is the width of the interval actually integrated over.
            ff[0] = (pars->phi_hi - pars->phi_lo) * JacCuba* ( gen_pars::pref_glue/ ( 2.0*M_PI) ) * ( p_t * pow( q_t,4.)* pow(k_t,4.)/ cutP2) * D1 * D2 *gen_pars::GeV2_to_fmm2 ;
          }
          return 0 ;
        }

        
    double GluonEnergyDens(GluonParsMV *  Variables){
      // double sqrts_t= Variables->sqrts;
      // double y_t= Variables->y;

      double error, result ;

      bool is_valid = (Variables->kmin<Variables->kmax);
      is_valid = is_valid && ( Variables->T1>0 && Variables->T2>0 );
      // std::cout << Variables->kmin<< "\t"<< Variables->kmax << "\t"<<Variables->T1 << "\t"<<Variables->T2 << std::endl;
      if(is_valid){
        Routines::make_cuhre_1C(3, GluonEnergyDensIntegrand,Variables,result,error);
      }
      else{
        
        return 0;
      }
      return result;
    }


    double quark_n12(double p, void * params){
      QuarkParsMV * pars = ( (struct QuarkParsMV *) params);
      double x1=p*exp( pars->y)/pars->sqrts;
      double x2=p*exp(-pars->y)/pars->sqrts;


      double Dip_t = (pars->dip)->FundamentalDipole_k(x2,p,pars->T);
      double result;
      //Quark contribution
    
      result = pow(p,1+pars->k)* ( (pars->partons)->xfpart(pars->quark_id,x1,p) ) * Dip_t* pow(2*M_PI,-1);
      if(result>0){return result;}
      else{return 0;}
    }

}

MV::MV(){}
MV::MV(Config ConfInput){
  gsl_set_error_handler_off();

  // Cuba forks worker processes of its own, defaulting to one per core. The
  // tabulation loop is already threaded over (T1,T2), so leaving that on would
  // fork a set of workers from every thread. Switched off (as in IP-Sat).
  { int ncores = 0, pcores = 10000; cubacores(&ncores, &pcores); }
  config=Config(ConfInput);

  /// THIS NEEDS CHANGES!
  // HERE GO THE PARAMETERS FOR MV?

  p_reg = config.get_ModelParams(0);
  std::cout<< "[MV]: Regulator : " << p_reg << std::endl;

  if(config.get_Verbose()>VerboseLevel::None){
    std::cout<<"[MV] Initializing MV model  "<< std::endl;
  }
  quark_dist = new PDFs(&ConfInput);

  if(config.get_Verbose()>VerboseLevel::None){std::cout<< std::endl;}
  Dip= new MVDipole(&config);
  config.set_TMax(Dip->getTmax());
  config.set_TMin(Dip->getTmin());
  config.set_NT(Dip->getNT());
  config.set_dT();
  // Dip->dump_test_norm();
} 

MV::~MV(){
  delete Dip;
  delete quark_dist;
}

void MV::MakeTable(std::string path_to_set){			
	// MAKES ith SET DIRECTORY. Uncomment when config dump is done
	SETPATH = path_to_set;
  fs::create_directories(SETPATH);
	// Write new config to setpath
	config.set_dump(SETPATH);
  
	if(config.get_Verbose()>VerboseLevel::None){std::cout<<"[MV]: New config written to "<<SETPATH  << std::endl;}
  if(config.get_Verbose()>VerboseLevel::None){std::cout<<"[MV]: Tabulating conserved charges in the MV model framework"<<SETPATH  << std::endl;}
  make_gluon_energy(); 
	if(config.get_Verbose()>VerboseLevel::None){std::cout<<"\n[MV]: Gluon Energy written to"<<SETPATH  << std::endl;}

  make_baryon_stopping_all();
	if(config.get_Verbose()>VerboseLevel::None){std::cout<<"\n [MV]: (All) Quark elements written to"<<SETPATH  << std::endl;}
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

void MV::make_gluon_energy(){

	std::ofstream density_f;
  std::ostringstream densityname;
  densityname << SETPATH <<"/"<< gluon_energy_table_name ;

	// Points an earlier build already computed. Keyed on coordinates, so the
	// table they came from may be on a different grid or left short by a job
	// that ran out of time.
	std::map<std::array<long,3>, double> reused_points;
	if(config.is_reuse()){
		std::ifstream seed(config.get_reuse_table());
		if(seed.is_open()){
			double y,t1,t2,v;
			while(seed >> y >> t1 >> t2 >> v){ reused_points[gluon_key(y,t1,t2)] = v; }
			std::cout<<"[MV]: ReuseTable: "<<reused_points.size()<<" points from "
			         <<config.get_reuse_table()<<std::endl;
		}
		else{
			std::cout<<"[MV]: ReuseTable: nothing at "<<config.get_reuse_table()
			         <<", computing all points"<<std::endl;
		}
	}
  if(config.get_Verbose()>VerboseLevel::None) std::cout<< "[MV]: Explicit integration for rapidities < " <<Dip->getMaxEta()<<std::endl;

	GluonParsMV parameters;
	parameters.sqrts= config.get_collEnergy();
  parameters.dip= Dip;
  parameters.p_reg=p_reg;
  parameters.kmax=Dip->getKmax();
  parameters.kmin=Dip->getKmin();

  density_f.open(densityname.str());
	const int NT = Dip->getNT();

	// Cuba initialises itself on the first call. Do that once here, serially,
	// rather than from whichever thread happens to arrive first.
	{
		GluonParsMV warm = parameters;
		warm.y = 0.0;
		warm.T1 = warm.T2 = Dip->T_at(NT-1);
		warm.phi_lo = -PHI_SPLIT; warm.phi_hi = PHI_SPLIT;
		MVYields::GluonEnergyDens(&warm);
	}

	for (int iy = 0; iy < config.get_NETA(); iy++) {
		double y_t = iy*config.get_dETA() + config.get_ETAMIN();
		std::vector<double> plane(long(NT)*NT);

		// x-range of the dipole table: same for the whole rapidity plane
		bool in_range = false;
		if ( std::fabs(y_t) < std::fabs(Dip->getMaxEta()) ){
			double x1l = parameters.kmin*exp(y_t)/parameters.sqrts;
			double x1h = parameters.kmax*exp(y_t)/parameters.sqrts;
			double x2l = parameters.kmin*exp(-y_t)/parameters.sqrts;
			double x2h = parameters.kmax*exp(-y_t)/parameters.sqrts;
			bool is_x1_within_range = (Dip->getXmin() <x1h) && (Dip->getXmax() > x1l);
			bool is_x2_within_range = (Dip->getXmin() <x2h) && (Dip->getXmax() > x2l);
			in_range = is_x1_within_range && is_x2_within_range;
		}

		// Every (T1,T2) point is an independent Cuhre call over a table that is
		// only read, so the rapidity plane parallelises directly. Each thread
		// takes its own copy of the parameter struct. Results go to an array and
		// are written afterwards, so the file does not depend on thread order.
		#pragma omp parallel for collapse(2) schedule(dynamic)
		for (int i1 = 0; i1 < NT; i1++) {
			for (int i2 = 0; i2 < NT; i2++) {
				double T1_t = Dip->T_at(i1);
				double T2_t = Dip->T_at(i2);

				GluonParsMV local = parameters;
				local.y = y_t; local.T1 = T1_t; local.T2 = T2_t;

				double res = 0.0;
				auto found = reused_points.find(gluon_key(y_t, T1_t, T2_t));
				if(found != reused_points.end()){ res = found->second; }
				else if(in_range){
					// Integrated in two pieces. The integrand carries a spike of
					// width ~p_reg where phi -> pi and k -> q, which Cuhre misses over
					// the full circle; split at 0.99 pi the first subdivision of the
					// tail resolves it. Makes e(y,T1,T2) = e(-y,T2,T1) exact.
					local.phi_lo = -PHI_SPLIT; local.phi_hi = PHI_SPLIT;
					double bulk = MVYields::GluonEnergyDens(&local);
					// The integrand is even in phi, so one side of the tail doubled.
					local.phi_lo = PHI_SPLIT; local.phi_hi = M_PI;
					double tail = MVYields::GluonEnergyDens(&local);
					res = bulk + 2.0*tail;
				}
				plane[long(i1)*NT + i2] = res;
			}
		}

		for (int i1 = 0; i1 < NT; i1++) {
			for (int i2 = 0; i2 < NT; i2++) {
				density_f<<y_t << "\t"<<Dip->T_at(i1) << "\t"<<Dip->T_at(i2) << "\t"<<plane[long(i1)*NT + i2]<< "\n";
			}
		}
		if(config.get_Verbose()==VerboseLevel::Dynamic){printProgress(double(iy)/double(config.get_NETA()));}
	}

	density_f.close();
	if(config.get_Verbose()==VerboseLevel::Dynamic){printProgress(1);}
	Routines::report_cuhre_stats("gluon energy table");
}

void MV::make_baryon_stopping(int k, QuarkID qid, QuarkID aqid){
  // gsl_set_error_handler_off();
	std::ofstream density_f;
  std::ostringstream densityname;
  densityname << SETPATH <<"/"<< get_quark_construction_file(k, qid) ;
	QuarkParsMV parameters;
	parameters.sqrts= config.get_collEnergy();
	parameters.k=k;
	parameters.partons = quark_dist;
  parameters.dip= Dip;
  parameters.kmax=Dip->getKmax();
  parameters.kmin=Dip->getKmin();

 

	gsl_integration_workspace * w= gsl_integration_workspace_alloc ( gen_pars::limit);
	gsl_function F;
	F.function = &MVYields::quark_n12;//(double p, void * params);

	double res12q,err12q;
	double res21q,err21q;
	double res12aq,err12aq;
	double res21aq,err21aq;

	density_f.open(densityname.str());
	for (int iy = 0; iy < config.get_NETA(); iy++) {
		double y_t = iy*config.get_dETA() + config.get_ETAMIN();

		for (int i1 = 0; i1 < Dip->getNT(); i1++) {
			double T_t = Dip->T_at(i1);

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

void MV::make_baryon_stopping_all(){
	make_baryon_stopping(0, QuarkID::u, QuarkID::ubar);
	make_baryon_stopping(1, QuarkID::u, QuarkID::ubar);
	if(config.get_Verbose()>VerboseLevel::None){std::cout<<"[MV]: U-Quark elements written \n " << std::endl;}

	make_baryon_stopping(0, QuarkID::d, QuarkID::dbar);
	make_baryon_stopping(1, QuarkID::d, QuarkID::dbar);
	if(config.get_Verbose()>VerboseLevel::None){std::cout<<"[MV]: D-Quark elements written \n " << std::endl;}

	make_baryon_stopping(0, QuarkID::s, QuarkID::sbar);
	make_baryon_stopping(1, QuarkID::s, QuarkID::sbar);
	if(config.get_Verbose()>VerboseLevel::None){std::cout<<"[MV]: S-Quark elements written \n " << std::endl;}


}

bool MV::check_if_zero_F(double y,double T){

  double x1=gen_pars::PMIN*exp(y)/config.get_collEnergy();

	// No nuclear matter means no charge; any T > 0 carries a small but physical
	// contribution and is computed.
	bool is_null = (T <= 0.0) ;
  is_null = is_null || (x1 > gen_pars::XdipMax);


	return is_null;
}

void  MV::printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
void MV::printProgress2(double percentage1, double percentage2) {
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

std::string MV::get_quark_construction_file(int k, QuarkID qid){
  std::ofstream quark_name_f;
  std::ostringstream quark_name;

  quark_name << "quark_construction_n"<< k;
  if(qid==QuarkID::u){quark_name <<"_u.dat";}
  else if(qid==QuarkID::d){quark_name <<"_d.dat";}
  else if(qid==QuarkID::s){quark_name <<"_s.dat";}

  return quark_name.str();
}


void MV::TestDump(double T1,double T2){
  std::ofstream density_f;
  std::ostringstream densityname;
  densityname << SETPATH <<"/TestDump_T1_"<< T1<<"_T2_"<<T2<< ".dat" ;
  if(config.get_Verbose()>VerboseLevel::None){std::cout<< "[MVDipole]: Writing dump for " << densityname.str() <<std::endl;}

	GluonParsMV parametersG;
	parametersG.sqrts= config.get_collEnergy();
  parametersG.dip= Dip;
  parametersG.T1 = T1;
  parametersG.T2 = T2;

  QuarkParsMV parametersQ;
	parametersQ.sqrts= config.get_collEnergy();
	parametersQ.partons = quark_dist;
  parametersQ.dip= Dip;


	gsl_integration_workspace * w= gsl_integration_workspace_alloc ( gen_pars::limit);
	gsl_function F;
	F.function = &MVYields::quark_n12;

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
    resEG = MVYields::GluonEnergyDens(&parametersG);

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
    // if(config.get_Verbose()){std::cout<<parametersG.y << "\t"<<resEG<< "\t"<<resEQ<< "\t"<<UDens<< "\t"<<DDens<< "\n";}
	}
  density_f.close();
  if(config.get_Verbose()>VerboseLevel::None){std::cout<<"-> Done"<<std::endl;}

}
