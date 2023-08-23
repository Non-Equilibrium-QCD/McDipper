/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */
#include <iostream>
#include <math.h>
#include <sstream>
#include <gsl/gsl_integration.h>
#include <sys/stat.h>
#include <filesystem>


#include "include/model_gbw.h"
#include "include/params_gen.h"

#define KINMODE 0

namespace fs = std::filesystem;

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

namespace GBW_funcs{
	double Q2F(double x, double Q02, double x0, double lambda, double Ts0){
		if( x>=1 ){return 0;}
		else{return Q02*pow(x0/x,lambda)*(1-x)*Ts0;}
	}

  double Q2A(double x, double Q02, double x0, double lambda, double Ts0){
		return (gen_pars::CA/gen_pars::CF)*Q2F(x,Q02,x0,lambda,Ts0);
	}

	double XF_AVG(double y, double sqrts, double Q02, double x0, double lambda, double Ts0){
		double power_term= (Q02*Ts0*exp(2*y))/(x0*x0*sqrts*sqrts);
		return x0*Q02*Ts0*pow( power_term, 1/(2.+lambda));
	}

	double Q2A_AVG(double y, double sqrts, double Q02, double x0, double lambda, double Ts0){
		double power_term= (gen_pars::CA/gen_pars::CF)*(Q02*Ts0*exp(2*y))/(x0*x0*sqrts*sqrts);
		return (gen_pars::CA/gen_pars::CF)*Q02*Ts0*pow( power_term, -lambda/(2.+lambda));
	}

	double Q2F_AVG(double y, double sqrts, double Q02, double x0, double lambda, double Ts0){
		double power_term= (Q02*Ts0*exp(2*y))/(x0*x0*sqrts*sqrts);
		return Q02*Ts0*pow( power_term, -lambda/(2.+lambda));
	}

  double Dipole(double k, double Q2){
		if(Q2<gen_pars::Q_tolerance){return 0;}
		else{return 4*M_PI *exp(- k*k/Q2 )/(Q2);}
  }
	// gsl_sf_bessel_J0(double x)
	// double gsl_sf_bessel_J1(double x)

  double quark_n12(double p, void * params){
		QuarkPars * pars = ( (struct QuarkPars *) params);
		double Ts0= (pars->Sigma0) * (pars->T);

		#if KINMODE==0
			double x1=p*exp( pars->y)/pars->sqrts;
			double x2=p*exp(-pars->y)/pars->sqrts;
			double Q2= Q2F(x2,pars->Q02,pars->x0,pars->lambda,Ts0);
		#elif KINMODE>0
			double x1= XF_AVG(pars->y, pars->sqrts, pars->Q02,pars->x0,pars->lambda,Ts0);
			double x2= XF_AVG(-pars->y, pars->sqrts, pars->Q02,pars->x0,pars->lambda,Ts0);
			double Q2 = Q2F_AVG(-pars->y, pars->sqrts, pars->Q02,pars->x0,pars->lambda,Ts0);
		#endif

		//Quark contribution
		if( x2>pars->xcut ){return 0;}
		else if( x1>pars->xcut ){return 0;}
		else{return pow(p,1+pars->k)* ( (pars->partons)->xfpart(pars->quark_id,x1,p) ) * Dipole(p,Q2) * pow(2*M_PI,-1); }
  }


	double gluon_density(double p,GluonPars * pars){
		double pref =  gen_pars::dA *pow(M_PI,-3)/(4*gen_pars::alpha_S*gen_pars::NC);
		/* Notice that the p^{-2} factor is not here. It has been added as the energy per mode 
		so to say  dens* p_{jacobian}* p_{first moment} to increase convergence. */
		double T1s0= (pars->Sigma0) * (pars->T1);
		double T2s0= (pars->Sigma0) * (pars->T2);

		#if KINMODE==0
		double x1=p*exp( pars->y)/pars->sqrts;
		double x2=p*exp(-pars->y)/pars->sqrts;
		double Q12= Q2A(x1,pars->Q02,pars->x0,pars->lambda,T1s0);
    	double Q22= Q2A(x2,pars->Q02,pars->x0,pars->lambda,T2s0);
		#elif KINMODE>0
		double Q12= Q2A_AVG( pars->y, pars->sqrts, pars->Q02,pars->x0,pars->lambda,T1s0);
		double Q22= Q2A_AVG(-pars->y, pars->sqrts, pars->Q02,pars->x0,pars->lambda,T2s0);
		#endif

		double QBar2 = Q12+Q22; double Xi2 =pow(Q12-Q22,2);
		double G_f = pow(p,4.)*Q12*Q22 +pow(p,2.)*QBar2*Xi2 +2*Q12*Q22*pow(QBar2,2) ;
		if(QBar2<gen_pars::Q_tolerance){return 0;}
		else if ( x1>pars->xcut){return 0;}
		else if ( x2>pars->xcut){return 0;}
		else{ return pref * exp(- p*p/QBar2 ) * Q12 * Q22 * G_f * pow(QBar2,-5.) * gen_pars::GeV2_to_fmm2;}
  }

	double gluon_energy_integrand(double p, void * params){
		return 2*M_PI*gluon_density(p,((struct GluonPars *) params));
	}

	double gluon_energy_analytical(GluonPars * pars){

		double pref = gen_pars::dA/(16*M_PI*sqrt(M_PI)*gen_pars::NC*gen_pars::alpha_S);

		double T1s0= (pars->Sigma0) * (pars->T1);
		double T2s0= (pars->Sigma0) * (pars->T2);

		double Q12= Q2A_AVG( pars->y, pars->sqrts, pars->Q02,pars->x0,pars->lambda,T1s0);
		double Q22= Q2A_AVG(-pars->y, pars->sqrts, pars->Q02,pars->x0,pars->lambda,T2s0);

		double QBar2 = Q12+Q22;
		double H_f =  2*Q12*Q12 + 7*Q12*Q22 +2*Q22*Q22 ;

		return pref * Q12 * Q22 * pow(QBar2,-5./2.) * H_f * gen_pars::GeV2_to_fmm2;
	}

	

	double gluon_energy(GluonPars * pars)
	{
		double res=0;
		double p_t = gen_pars::PMIN;
		for (size_t ip = 0; ip < gen_pars::NP; ip++) {
			p_t += gen_pars::DP * ip;
			if(ip==0 || ip==gen_pars::NP-1){res += p_t*p_t*gluon_density(p_t, pars) /2. ;}
			else{res += p_t*p_t*gluon_density(p_t, pars) ;}
		}
	  return 2*M_PI*gen_pars::DP * res;
	}


}


///////////////////////////////////////////////////////
//                  GBW MODEL CLASS                  //
///////////////////////////////////////////////////////


GBW::GBW(){}
GBW::GBW(Config ConfInput){
  config=Config(ConfInput);
	quark_dist = new PDFs(&ConfInput);
  if(config.get_Verbose()){std::cout<<"   --> Initializing GBW model  "<< std::endl;}

	config.set_TMax(gen_pars::TMax);
	config.set_TMin(gen_pars::TMin);
	config.set_NT(gen_pars::NT);
	config.set_dT();
}
GBW::~GBW(){}

void GBW::MakeTable(std::string path_to_set){

	// MAKES ith SET DIRECTORY. Uncomment when config dump is done
	SETPATH = path_to_set;
  fs::create_directories(SETPATH);
	// Write new config to setpath
	config.set_dump(SETPATH);
	if(config.get_Verbose()){std::cout<<"New config written to "<<SETPATH  << std::endl;}

	if(config.get_Verbose()){std::cout<<"--> Tabulating conserved charges in the GBW model framework"<<SETPATH  << std::endl;}
	make_gluon_energy();
	if(config.get_Verbose()){std::cout<<"\nGluon Energy written to"<<SETPATH  << std::endl;}
	make_baryon_stopping_all();
	if(config.get_Verbose()){std::cout<<"\n(All) Quark elements written to"<<SETPATH  << std::endl;}

}


void GBW::make_gluon_energy(){
	std::ofstream density_f;
  std::ostringstream densityname;
  densityname << SETPATH <<"/"<< gluon_energy_table_name ;
	//																	sqrts        y       Q02                          x0                        lambda            T1s0   T2s0
	GluonPars parameters;
	parameters.sqrts= config.get_collEnergy();
	parameters.Q02= config.get_ModelParams(0);
	parameters.x0= config.get_ModelParams(1);
	parameters.lambda= config.get_ModelParams(2);
	parameters.xcut= config.get_ModelParams(3);
	parameters.Sigma0= 2*M_PI*config.get_BG();

	#if KINMODE<2
	gsl_integration_workspace * w= gsl_integration_workspace_alloc ( gen_pars::limit);
	gsl_function F;
  F.function = &GBW_funcs::gluon_energy_integrand;
	#endif

	double res,err;
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
				double Q12_t,Q22_t;
				bool is_null;

				if(T1_t==0. || T2_t==0. ){is_null=true;}
				else{is_null = check_if_zero(y_t,T1_t, T2_t,Q12_t,Q22_t);}


				if(is_null){res=0;err=0;}
				else{
					#if KINMODE<2
					F.params = &parameters;
					gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res,&err);
					#elif KINMODE==2
					res=GBW_funcs::gluon_energy_analytical(&parameters);
					#endif
				}
				density_f<<y_t<< "\t"<<T1_t<< "\t"<<T2_t<< "\t"<<res<< "\n";
			}
			// density_f<<"\n";
		}
		density_f<<"\n";
		if(config.get_Verbose()){
			if(remainder(iy,gen_pars::skip)==0){double percentage_done = double(iy)/double(config.get_NETA());printProgress(percentage_done);}
		}

	}
	#if KINMODE<2
	gsl_integration_workspace_free (w);
	#endif
	 //
	density_f.close();
	if(config.get_Verbose()){printProgress(1);}
}

void GBW::make_baryon_stopping(int k, QuarkID qid, QuarkID aqid){
	std::ofstream density_f;
  std::ostringstream densityname;
  densityname << SETPATH <<"/"<< get_quark_construction_file(k, qid) ;
	QuarkPars parameters;
	parameters.sqrts= config.get_collEnergy();
	parameters.Q02= config.get_ModelParams(0);
	parameters.x0= config.get_ModelParams(1);
	parameters.lambda= config.get_ModelParams(2);
	parameters.xcut= config.get_ModelParams(3);
	parameters.Sigma0= 2*M_PI*config.get_BG();
	parameters.k=k;
	parameters.partons = quark_dist;

	gsl_integration_workspace * w= gsl_integration_workspace_alloc ( gen_pars::limit);
	gsl_function F;
	F.function = &GBW_funcs::quark_n12;//(double p, void * params);

	double res12q,err12q;
	double res21q,err21q;
	double res12aq,err12aq;
	double res21aq,err21aq;

	density_f.open(densityname.str());
	for (size_t iy = 0; iy < config.get_NETA(); iy++) {
		double y_t = iy*config.get_dETA() + config.get_ETAMIN();

		for (size_t i1 = 0; i1 < config.get_NT(); i1++) {
			double T_t = i1*config.get_dT() + config.get_TMin();

			if(T_t==0.){
				res12q=0;err12q=0;
				res21q=0;err21q=0;
				res12aq=0;err12aq=0;
				res21aq=0;err21aq=0;
			}
			else{
				double x1min=gen_pars::PMIN *exp(y_t)/ config.get_collEnergy();
				double x2min=gen_pars::PMIN *exp(-y_t)/ config.get_collEnergy();
				parameters.y = y_t ;
				parameters.T = T_t;
				parameters.quark_id = qid;
				F.params = &parameters;

				gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12q,&err12q);
				parameters.y = -y_t ;
				F.params = &parameters;

				gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21q,&err21q);


				parameters.y = y_t ;
				F.params = &parameters;
				parameters.quark_id = aqid;

				gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res12aq,&err12aq);
				parameters.y = -y_t ;
				F.params = &parameters;

				gsl_integration_qag(&F,gen_pars::PMIN, gen_pars::PMAX,gen_pars::epsabs, gen_pars::epsrel, gen_pars::limit, gen_pars::routine, w, &res21aq,&err21aq);
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
  
void GBW::make_baryon_stopping_all(){
	make_baryon_stopping(0, QuarkID::u, QuarkID::ubar);
	make_baryon_stopping(1, QuarkID::u, QuarkID::ubar);
	if(config.get_Verbose()){std::cout<<"\n U-Quark elements written " << std::endl;}

	make_baryon_stopping(0, QuarkID::d, QuarkID::dbar);
	make_baryon_stopping(1, QuarkID::d, QuarkID::dbar);
	if(config.get_Verbose()){std::cout<<"\n D-Quark elements written " << std::endl;}

	make_baryon_stopping(0, QuarkID::s, QuarkID::sbar);
	make_baryon_stopping(1, QuarkID::s, QuarkID::sbar);
	if(config.get_Verbose()){std::cout<<"\n S-Quark elements written " << std::endl;}
}

bool GBW::check_if_zero(double y,double T1, double T2, double &Q12, double &Q22){

	double x1=gen_pars::PMIN*exp( y)/config.get_collEnergy();
	double x2=gen_pars::PMIN*exp(-y)/config.get_collEnergy();

	bool is_null = (x1 > gen_pars::XdipMax);
	is_null = is_null || (x2 > gen_pars::XdipMax);

	return is_null;
}


bool GBW::check_if_zero_F(double y,double T1, double T2, double &Q12, double &Q22){

	double x1=gen_pars::PMIN*exp( y)/config.get_collEnergy();
	double x2=gen_pars::PMIN*exp(-y)/config.get_collEnergy();
	Q12= GBW_funcs::Q2F(x1,config.get_ModelParams(0),config.get_ModelParams(1),config.get_ModelParams(2),T1);
	Q22= GBW_funcs::Q2F(x2,config.get_ModelParams(0),config.get_ModelParams(1),config.get_ModelParams(2),T2);

	bool is_null = (Q12 < gen_pars::Q_tolerance);
	is_null = is_null ||  (Q22 < gen_pars::Q_tolerance);
	return is_null;
}

void  GBW::printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

std::string GBW::get_quark_construction_file(int k, QuarkID qid){
  std::ofstream quark_name_f;
  std::ostringstream quark_name;

  quark_name << "quark_construction_n"<< k;
  if(qid==QuarkID::u){quark_name <<"_u.dat";}
  else if(qid==QuarkID::d){quark_name <<"_d.dat";}
  else if(qid==QuarkID::s){quark_name <<"_s.dat";}

  return quark_name.str();
}
