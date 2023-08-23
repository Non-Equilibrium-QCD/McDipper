/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */

#include <iostream>
#include <math.h>
#include <random>
#include <filesystem>
#include <sstream>

#include "include/event.h"
#include "include/params_gen.h"
// #include "include/params_gen.h"

#define OPTICAL 1
namespace fs = std::filesystem;

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

// Constructor and Destructor

Event::Event(Config ConfInput){

	config=Config(ConfInput);
	N1.A=config.get_A(1);N1.Z=config.get_Z(1);N1.mode=config.get_NuclearMode(1);
	N2.A=config.get_A(2);N2.Z=config.get_Z(2);N2.mode=config.get_NuclearMode(2);

	PrimalSeed=config.get_seed();

	b1=new double[2];
	b2=new double[2];

	sqrtsNN=config.get_collEnergy();

	NX=config.get_NX();
	NY=config.get_NY();
	MakeGrid();

	ChargeMaker = new Charges(config);

	Initialize_output();

	cell_trans_volume=config.get_dX()*config.get_dY();

}

Event::~Event(){
	delete[] b1;
	delete[] b2;
	delete[] T1p_ptr;
	delete[] T1n_ptr;
	delete[] T2p_ptr;
	delete[] T2n_ptr;
}

// Structural Functions
void Event::MakeGrid(){
	T1p_ptr = new double[config.get_NX() * config.get_NY()];
	T1n_ptr = new double[config.get_NX() * config.get_NY()];
	T2p_ptr = new double[config.get_NX() * config.get_NY()];
	T2n_ptr = new double[config.get_NX() * config.get_NY()];
}

// Interfaces of T1 and T2
/* These functions allows the user to use the Thickness pointers safely and easily. */

const double& Event::T1p (int64_t nx, int64_t ny) const {return (T1p_ptr)[ NY*nx + ny];}
double& Event:: T1p (int64_t nx, int64_t ny){return (T1p_ptr)[ NY*nx + ny];}

const double& Event::T1n (int64_t nx, int64_t ny) const {return (T1n_ptr)[ NY*nx + ny];}
double& Event:: T1n (int64_t nx, int64_t ny){return (T1n_ptr)[ NY*nx + ny];}

const double& Event::T2p (int64_t nx, int64_t ny) const {return (T2p_ptr)[ NY*nx + ny];}
double& Event:: T2p (int64_t nx, int64_t ny){return (T2p_ptr)[ NY*nx + ny];}

const double& Event::T2n (int64_t nx, int64_t ny) const {return (T2n_ptr)[ NY*nx + ny];}
double& Event:: T2n (int64_t nx, int64_t ny){return (T2n_ptr)[ NY*nx + ny];}

double Event::T1(int64_t nx, int64_t ny){return  T1p (nx,ny)+T1n (nx,ny); }
double Event::T2(int64_t nx, int64_t ny){return  T2p (nx,ny)+T2n (nx,ny); }

// Evaluation
// Functions

void Event::NewEvent(int EventID_in){
	EventID=EventID_in;
	EventSeed =std::rand(); //  CHECK THIS !

	srand48(EventSeed);// Dump this in the config?
	// Create Nuclei
	bool is_event_valid=false;
	if(config.get_Verbose()){
		std::cout<< "|---------------------------------- New Event: "<< EventID+1<<"/"<< config.get_NEvents() << " -------------------------------------|\n";
		std::cout<< "    Event ID: " <<EventID << "       Event Seed: " << EventSeed << "\n"<<std::endl;}

	while(!is_event_valid){

		Nucleus A1(N1.A,N1.Z,N1.mode);
		Nucleus A2(N2.A,N2.Z,N2.mode);
		//Sample b
		if(config.get_ImpactMode()== ImpSample::Fixed){get_impact_from_value(config.get_ImpactValue(),0);}
		else if(config.get_ImpactMode()== ImpSample::dbSampled){sample_db_impact(config.get_MinImpact(), config.get_MaxImpact());}
		else if(config.get_ImpactMode()== ImpSample::bdbSampled){sample_bdb_impact(config.get_MinImpact(), config.get_MaxImpact());}
		else{
			std::cerr<<"Error! Impact Sampling mode invalid"<<std::endl;
			std::cerr<<"       Use - 0 for fixed value b"<<std::endl;
			std::cerr<<"           - 1 for uniform b-sampling (P(b) = K )"<<std::endl;
			std::cerr<<"           - 2 for quadratic b-sampling (P(b) = K*b )"<<std::endl;
			exit(EXIT_FAILURE);
		}

		A1.shift_nucleus_by_impact(b1[0],b1[1]);
		A2.shift_nucleus_by_impact(b2[0],b2[1]);

		#if OPTICAL==0
		if(config.get_Verbose()){std::cout<< "    Running for impact parameter ( bx = " << b1[0]-b2[0] << " , by = " << b1[1]-b2[1] << ")"<<std::endl ;}
		is_event_valid=true;
		for (size_t ix = 0; ix < config.get_NX(); ix++) {
			for (size_t iy = 0; iy < config.get_NY(); iy++) {
				double x_t= get_x(ix);
				double y_t= get_y(iy);
				T1p (ix,iy)= A1.get_Z()*A1.nuclear_thickness_optical(x_t-b1[0], y_t-b1[1]);
				T1n (ix,iy)= (A1.get_A()-A1.get_Z() )*A1.nuclear_thickness_optical(x_t-b1[0], y_t-b1[1]);
				T2p (ix,iy)= A2.get_Z()*A2.nuclear_thickness_optical(x_t-b2[0], y_t-b2[1]);
				T2n (ix,iy)= (A2.get_A()-A2.get_Z() )*A2.nuclear_thickness_optical(x_t-b2[0], y_t-b2[1]);
			}
		}
		MakeGlobalQuantities();

		#elif OPTICAL==1


		CheckParticipants(&A1,&A2);

		is_event_valid = (A1.get_number_of_participants()>0) && (A2.get_number_of_participants()>0);
		if(is_event_valid){
			if(config.get_Verbose()){
				std::cout<< "    Running for impact parameter ( bx = " << b1[0]-b2[0] << " , by = " << b1[1]-b2[1] << ")"<<std::endl ;
				std::cout<< "    Number of Participants: " <<A1.get_number_of_participants() + A2.get_number_of_participants() <<  "\n";
				std::cout<< "    Number of Collisions: " <<A1.get_NColl()+A2.get_NColl() <<  "\n";
			}

			dump_nucleon_pos(&A1,&A2);

			if(config.get_Verbose()){std::cout<<"    Nucleon Positions written out\n";}

			for (size_t ix = 0; ix < config.get_NX(); ix++) {
				for (size_t iy = 0; iy < config.get_NY(); iy++) {
					double x_t= get_x(ix);
					double y_t= get_y(iy);
					T1p (ix,iy)= A1.GetThickness_p(x_t,y_t, config.get_BG());
					T1n (ix,iy)= A1.GetThickness_n(x_t,y_t, config.get_BG());
					T2p (ix,iy)= A2.GetThickness_p(x_t,y_t, config.get_BG());
					T2n (ix,iy)= A2.GetThickness_n(x_t,y_t, config.get_BG());
				}
			}
			print_glauber_data_to_file(&A1,&A2);

			MakeGlobalQuantities();
			bool is_charge_output= config.is_format("EMoments") || config.is_format("NMoments") || config.is_format("Charges");
			if (is_charge_output){MakeChargeOutput();}
		}
		#endif
	}

	if(config.get_Verbose()){
		std::cout<<"    Thickness Functions computed!\n";
		std::cout<< "\n|------------------------------------- End Event ----------------------------------------|\n";
	}
}


void Event::MakeEventByEvent(){

	/* This function iterates over the number of listed EVENTS */
	if( ev0 < config.get_NEvents()){
		std::cerr<< "[ Warning::Event ]: Event generation starting from event # "<< ev0 << std::endl;
		for (size_t ev = ev0; ev < config.get_NEvents(); ev++) {NewEvent(get_ID(ev,config.get_NEvents()));}
	}
	else{std::cerr<< "[ Error::Event ]: Initial event ID="<< ev0 << " larger than Nevents in config file. Exiting."<<std::endl;}
}


void Event::get_impact_from_value(double AbsB, double ThetaB){
	b1[0] =  AbsB*cos(ThetaB)/2.;
	b1[1] =  AbsB*sin(ThetaB)/2.;
	b2[0] = -AbsB*cos(ThetaB)/2.;
	b2[1] = -AbsB*sin(ThetaB)/2.;
}

void Event::sample_db_impact(double bmin, double bmax){
	double ThetaB = 0; // Without loss of generality, we
	double AbsB = (bmax-bmin)*uni_ev_rn() + bmin;
	get_impact_from_value(AbsB,ThetaB);
}

void Event::sample_bdb_impact(double bmin, double bmax){
	double ThetaB = 0; // Without loss of generality, we
	double AbsB = sqrt( ( pow(bmax,2) - pow(bmin,2) )*uni_ev_rn() + pow(bmin,2) );
	get_impact_from_value(AbsB,ThetaB);
}

double Event::get_x(double ix){return ix*config.get_dX() + config.get_XMIN();}
double Event::get_y(double iy){return iy*config.get_dY() + config.get_YMIN();}
double Event::get_eta(double ieta){return ieta*config.get_dETA() + config.get_ETAMIN();}

////////////////////

void Event::Initialize_output(){

	if(config.is_output_named()){

		std::ostringstream dirname_t;
		dirname_t << config.get_out_path()<<"/"<< config.get_run_name();
		std::string dirpath = dirname_t.str();
		OUTPATH = dirpath;
		if(config.get_Verbose()){std::cout<< OUTPATH<< std::endl;}
		ev0=0;
		if ( !(fs::create_directories(dirpath))){
			std::cerr<< "[ Warning::Event ]: Output directory ("<< OUTPATH<< ") already exists!"<< std::endl ;
			for (size_t iev = 0; iev < config.get_NEvents(); iev++)
			{
				std::ostringstream filename_t;
				filename_t << config.get_out_path()<<"/" <<config.get_run_name() << "/global_"<< iev<< ".dat"	;
				std::string filepath = filename_t.str();

				if(fs::exists(filepath) && iev==config.get_NEvents()-1){ev0=config.get_NEvents();}
				if(!fs::exists(filepath) && iev>0){ev0=iev-1;break;}//set the initial event to the last event
				
			}
			
		}

	}
	else{
		int dircount = 0;
		bool direx = false;

		while (!direx) {
			std::ostringstream dirname_t;
			dirname_t << config.get_out_path()<<"/"<< "Run" << dircount;
			std::string dirpath = dirname_t.str();
			direx = fs::create_directories(dirpath);
			dircount++;
		}
		std::ostringstream dirname;
		dirname << config.get_out_path()<<"/"<< "Run"<< dircount-1;
		ev0=0; //set the initial event to 0
		OUTPATH = dirname.str();
		if(config.get_Verbose()){std::cout<< OUTPATH<< std::endl;}
		
	}
	fs::create_directories(OUTPATH);
	config.dump(OUTPATH);

}

void Event::dump_nucleon_pos(Nucleus *A1,Nucleus *A2){

	std::ofstream pos_f;
  std::ostringstream posname;
  posname << OUTPATH << "/Event_"<<EventID<<"_Nuclei.dat";
  pos_f.open(posname.str());

	double x_t,y_t,z_t;
  for (size_t i = 0; i < A1->get_A(); i++) {
		A1->get_position_nucleon(i, x_t,y_t,z_t);
		pos_f<< i <<"\t"<< x_t<<"\t"<< y_t<<"\t"<< z_t<<"\t"<< A1->get_ParticipantStatus(i)<<std::endl;
	}
		pos_f<<std::endl;
	for (size_t i = 0; i < A2->get_A(); i++) {
		A2->get_position_nucleon(i, x_t,y_t,z_t);
		pos_f<< i <<"\t"<< x_t<<"\t"<< y_t<<"\t"<< z_t<<"\t"<< A2->get_ParticipantStatus(i)<<std::endl;
	}
pos_f.close();
}

void Event::MakeGlobalQuantities(){


	std::ofstream global_f;
	std::ostringstream global_name;
	global_name << OUTPATH <<"/global_"<<EventID<< ".dat" ;
	global_f.open(global_name.str(), std::ios_base::app);

	double eg_t,eq_t,nu_t,nd_t,ns_t;
	double t1p_t,t1n_t,t2p_t,t2n_t;
	double t1_t,t2_t;

	double total_energy_event_q=0;
	double total_energy_event_g=0;

	double total_u_event=0;
	double total_d_event=0;
	double total_s_event=0;
	double total_energy_event,total_q_event,total_B_event;

	x_cm.assign(config.get_NETA(),0.0);
	y_cm.assign(config.get_NETA(),0.0);
	std::vector<double> x_cm_unnorm(config.get_NETA(),0.0);
	std::vector<double> y_cm_unnorm(config.get_NETA(),0.0);

	std::vector<double> egtau0_2o3_midrap(config.get_NETA(),0.0);
	std::vector<double> eqtau0_2o3_midrap(config.get_NETA(),0.0);
	std::vector<double> etottau0_2o3_midrap(config.get_NETA(),0.0);

	int izero = int( -config.get_ETAMIN()/config.get_dETA());

	dEdeta.assign(config.get_NETA(),0.0);
	dEgdeta.assign(config.get_NETA(),0.0);
	dEqdeta.assign(config.get_NETA(),0.0);
	dnudeta.assign(config.get_NETA(),0.0);
	dnddeta.assign(config.get_NETA(),0.0);
	dnsdeta.assign(config.get_NETA(),0.0);

	double x_cm_global_unnorm=0.;
	double y_cm_global_unnorm=0.;
 	x_cm_global=0.;
	y_cm_global=0.;

	if(config.get_Verbose()){std::cout<< "\n Making global observables and observables for Event " << EventID <<std::endl;}

	for (size_t ieta = 0; ieta < config.get_NETA(); ieta++) {
		double eta_t = ieta*config.get_dETA() + config.get_ETAMIN();
		for (size_t ix = 0; ix < config.get_NX(); ix++) {
			for (size_t iy = 0; iy < config.get_NY(); iy++) {

				double x_t= get_x(ix);
				double y_t= get_y(iy);

				t1p_t=T1p(ix,iy);
				t1n_t=T1n(ix,iy);
				t2p_t=T2p(ix,iy);
				t2n_t=T2n(ix,iy);
				t1_t = t1p_t+t1n_t;
				t2_t = t2p_t+t2n_t;

				eg_t=ChargeMaker->gluon_energy(eta_t, t1_t, t2_t) * config.get_KFactor() ;
				eq_t=ChargeMaker->quark_energy(eta_t, t1_t, t2_t) ;
				nu_t=ChargeMaker->u_density(eta_t, t1p_t, t1n_t, t2p_t, t2n_t) ;
				nd_t=ChargeMaker->d_density(eta_t, t1p_t, t1n_t, t2p_t, t2n_t) ;
				ns_t=ChargeMaker->s_density(eta_t, t1_t, t2_t) ;

				x_cm_unnorm[ieta] += cell_trans_volume * x_t * (eg_t+eq_t);
				y_cm_unnorm[ieta] += cell_trans_volume * y_t * (eg_t+eq_t);


				dEgdeta[ieta]+= cell_trans_volume*(eg_t);
				dEqdeta[ieta]+= cell_trans_volume*(eq_t);
				dEdeta[ieta]+= cell_trans_volume*(eg_t+eq_t);


				dnudeta[ieta]+= cell_trans_volume*nu_t;
				dnddeta[ieta]+= cell_trans_volume*nd_t;
				dnsdeta[ieta]+= cell_trans_volume*ns_t;

				total_energy_event_g += config.get_dETA() * cell_trans_volume * eg_t;
				total_energy_event_q += config.get_dETA() * cell_trans_volume * eq_t;

				egtau0_2o3_midrap[ieta] += cell_trans_volume * gen_pars::GeV2_to_fmm2 * pow(gen_pars::fmm2_to_GeV2*eg_t,2/3.);//cell_trans_volume* gen_pars::GeV2_to_fmm2*pow( gen_pars::fmm2_to_GeV2 * (eg_t + eq_t), 2./3.);
				eqtau0_2o3_midrap[ieta] += cell_trans_volume * gen_pars::GeV2_to_fmm2 * pow(gen_pars::fmm2_to_GeV2*eq_t,2/3.);
				etottau0_2o3_midrap[ieta] += cell_trans_volume * gen_pars::GeV2_to_fmm2 * pow( gen_pars::fmm2_to_GeV2*(eg_t+eq_t),2/3.);


			}
		}

		if(dEdeta[ieta]==0){
			x_cm[ieta]=0;
			y_cm[ieta]=0;
		}
		else{
			x_cm[ieta] += x_cm_unnorm[ieta]/dEdeta[ieta];
			y_cm[ieta] += y_cm_unnorm[ieta]/dEdeta[ieta];
		}

  

		x_cm_global_unnorm += config.get_dETA() *x_cm_unnorm[ieta] ;
		y_cm_global_unnorm += config.get_dETA() *y_cm_unnorm[ieta] ;

		total_u_event += config.get_dETA() * dnudeta[ieta];
		total_d_event += config.get_dETA() * dnddeta[ieta];
		total_s_event += config.get_dETA() * dnsdeta[ieta];

		if(config.get_Verbose()){
			double percentage_done=double(ieta)/double(config.get_NETA());
			printProgress(percentage_done);
		}
	}

	if(config.get_Verbose()){printProgress(1);}

	total_energy_event=total_energy_event_q+total_energy_event_g;

	if(total_energy_event==0){
		x_cm_global=0;
		y_cm_global=0;
	}
	else{
		x_cm_global=x_cm_global_unnorm/total_energy_event;
		y_cm_global=y_cm_global_unnorm/total_energy_event;
	}

	total_q_event= (2*total_u_event-total_d_event)/3;
	total_B_event= (total_u_event+total_d_event)/3;

	global_f<<  "E_g\t" << total_energy_event_g<<std::endl;
	global_f<<  "E_q\t" << total_energy_event_q<<std::endl;
	global_f<<  "E_tot\t" << total_energy_event<<std::endl;

	global_f<<  "dEg_deta\t" << dEgdeta[izero]<<std::endl;
	global_f<<  "dEq_deta\t" << dEqdeta[izero]<<std::endl;
	global_f<<  "dE_deta\t" << dEdeta[izero]<<std::endl;

	global_f<<  "Int(egtau)2/3\t" << egtau0_2o3_midrap[izero] <<std::endl;
	global_f<<  "Int(eqtau)2/3\t" << eqtau0_2o3_midrap[izero]<<std::endl;
	global_f<<  "Int(egtau+eqtau)2/3\t" << etottau0_2o3_midrap[izero]<<std::endl;

	global_f<<  "N_u\t" << total_u_event<<std::endl;
	global_f<<  "N_d\t" << total_d_event<<std::endl;
	global_f<<  "N_s\t" << total_s_event<<std::endl;

	global_f<<  "dNu_deta(eta=0)\t" << dnudeta[izero]<<std::endl;
	global_f<<  "dNd_deta(eta=0)\t" << dnddeta[izero]<<std::endl;
	global_f<<  "dNs_deta(eta=0)\t" << dnsdeta[izero]<<std::endl;

	global_f<<  "Q\t" << total_q_event<<std::endl;
	global_f<<  "B\t" << total_B_event<<std::endl;

	global_f<<  "dQ_deta(eta=0)\t" << (2*dnudeta[izero]-dnddeta[izero] )/3.<<std::endl;
	global_f<<  "dB_deta(eta=0)\t" << (dnudeta[izero]+dnddeta[izero] )/3.<<std::endl;

	global_f<<  "x_cm_global\t" << x_cm_global<<std::endl;
	global_f<<  "y_cm_global\t" << y_cm_global<<std::endl;

	global_f<<  "x_cm_mid\t" << x_cm[izero]<<std::endl;
	global_f<<  "y_cm_mid\t" << y_cm[izero]<<std::endl;


	if(config.get_Verbose()){
		std::cout << "\n\n    E_g = " << total_energy_event_g << "\tE_q = " << total_energy_event_q << "\tE_tot = " << total_energy_event<<std::endl ;
		std::cout<<  "    dE_deta = " << dEdeta[izero]<<  "\t with \tdEg_deta = " << dEgdeta[izero]<<  "    and    dEq_deta = " << dEqdeta[izero]<<std::endl;
		std::cout << "    N_u = " << total_u_event<< "\tN_d = " << total_d_event << "\tN_s = " << total_s_event <<std::endl ;
		std::cout<<  "    dN_i_deta(eta=0)     u: " << dnudeta[izero];
		std::cout<<  "\t d: " << dnddeta[izero]<< "\t s: " << dnsdeta[izero]<<std::endl;
		std::cout << "    Q = " << total_q_event<< "\tB = " << total_B_event <<std::endl ;
		std::cout << "\nGlobal CoM coordinates:  ( x = " << x_cm_global << " ,y = " << y_cm_global<< ")"<<std::endl ;

	}

	global_f.close();

}

void Event::MakeChargeOutput(){
	std::ofstream charges_f; std::ostringstream charges_name;
	std::ofstream e_moments_f; std::ostringstream e_moments_name;
	std::ofstream nu_moments_f; std::ostringstream nu_moments_name;
	std::ofstream nd_moments_f; std::ostringstream nd_moments_name;

	if(config.is_format("Charges")){
		charges_name << OUTPATH <<"/Event_"<<EventID<< "_charges.dat" ;
		charges_f.open(charges_name.str());
	}
	
	if(config.is_format("EMoments")){
		e_moments_name << OUTPATH <<"/EnergyMoments_"<<EventID<< ".dat" ;
		e_moments_f.open(e_moments_name.str());
	}

	if(config.is_format("NMoments")){
		nu_moments_name << OUTPATH <<"/NuMoments_"<<EventID<< ".dat" ;
		nu_moments_f.open(nu_moments_name.str());

		nd_moments_name << OUTPATH <<"/NdMoments_"<<EventID<< ".dat" ;
		nd_moments_f.open(nd_moments_name.str());
	}

	double eg_t,eq_t,nu_t,nd_t,ns_t;
	double t1p_t,t1n_t,t2p_t,t2n_t;
	double t1_t,t2_t;


	if(config.get_Verbose()){std::cout<< "\nWriting out charges and moments for Event " << EventID <<std::endl;}

	for (size_t ieta = 0; ieta < config.get_NETA(); ieta++) {
		double eta_t = ieta*config.get_dETA() + config.get_ETAMIN();

		double lattice_sum_dint23dy = 0.;

		std::vector<double> lattice_sum_eg(n_max+1, 0.0);
		std::vector<double> lattice_sum_eq(n_max+1, 0.0);

		std::vector<double> lattice_sum_eg_cos(n_max+1, 0.0);
		std::vector<double> lattice_sum_eq_cos(n_max+1, 0.0);

		std::vector<double> lattice_sum_eg_sin(n_max+1, 0.0);
		std::vector<double> lattice_sum_eq_sin(n_max+1, 0.0);

		std::vector<double> lattice_sum_nu(n_max+1, 0.0);
		std::vector<double> lattice_sum_nd(n_max+1, 0.0);

		std::vector<double> lattice_sum_nu_cos(n_max+1, 0.0);
		std::vector<double> lattice_sum_nd_cos(n_max+1, 0.0);

		std::vector<double> lattice_sum_nu_sin(n_max+1, 0.0);
		std::vector<double> lattice_sum_nd_sin(n_max+1, 0.0);


		for (size_t ix = 0; ix < config.get_NX(); ix++) {
			for (size_t iy = 0; iy < config.get_NY(); iy++) {

				double x_t= get_x(ix)-x_cm_global;
				double y_t= get_y(iy)-y_cm_global;

				t1p_t=T1p(ix,iy);
				t1n_t=T1n(ix,iy);
				t2p_t=T2p(ix,iy);
				t2n_t=T2n(ix,iy);
				t1_t = t1p_t+t1n_t;
				t2_t = t2p_t+t2n_t;


				eg_t=ChargeMaker->gluon_energy(eta_t, t1_t, t2_t) * config.get_KFactor()  ;
				eq_t=ChargeMaker->quark_energy(eta_t, t1_t, t2_t) ;
				nu_t=ChargeMaker->u_density(eta_t, t1p_t, t1n_t, t2p_t, t2n_t) ;
				nd_t=ChargeMaker->d_density(eta_t, t1p_t, t1n_t, t2p_t, t2n_t) ;
				ns_t=ChargeMaker->s_density(eta_t, t1_t, t2_t) ;

				double r_t = sqrt(x_t*x_t + y_t*y_t);
				double phi_t = atan2(y_t,x_t);
				if(config.is_format("EMoments")){
					lattice_sum_dint23dy +=  cell_trans_volume* gen_pars::GeV2_to_fmm2*pow( gen_pars::fmm2_to_GeV2 * (eg_t + eq_t), 2./3.);
				}	

				for (int j = 0; j <= n_max; j++)
				{
					double jj = (double) j;
					if(config.is_format("EMoments")){
						lattice_sum_eg[j] += cell_trans_volume*pow(r_t,jj) * eg_t;
						lattice_sum_eq[j] += cell_trans_volume*pow(r_t,jj) * eq_t;
						//
						lattice_sum_eg_cos[j] += cell_trans_volume*pow(r_t,jj) * cos(jj*phi_t) * eg_t ;
						lattice_sum_eq_cos[j] += cell_trans_volume*pow(r_t,jj) * cos(jj*phi_t) * eq_t;

						lattice_sum_eg_sin[j] += cell_trans_volume*pow(r_t,jj) * sin(jj*phi_t) * eg_t;
						lattice_sum_eq_sin[j] += cell_trans_volume*pow(r_t,jj) * sin(jj*phi_t) * eq_t;
					}
					if(config.is_format("NMoments")){
						lattice_sum_nu[j] += cell_trans_volume*pow(r_t,jj) * nu_t;
						lattice_sum_nd[j] += cell_trans_volume*pow(r_t,jj) * nd_t;
						//
						lattice_sum_nu_cos[j] += cell_trans_volume*pow(r_t,jj) * cos(jj*phi_t) * nu_t ;
						lattice_sum_nd_cos[j] += cell_trans_volume*pow(r_t,jj) * cos(jj*phi_t) * nd_t;

						lattice_sum_nu_sin[j] += cell_trans_volume*pow(r_t,jj) * sin(jj*phi_t) * nu_t;
						lattice_sum_nd_sin[j] += cell_trans_volume*pow(r_t,jj) * sin(jj*phi_t) * nd_t;
					}
				}

				if(config.is_format("Charges")){
					charges_f<< eta_t <<"\t"<< x_t <<"\t"<< y_t <<"\t"<< eg_t  <<"\t"<< eq_t<<"\t"<< nu_t <<"\t"<< nd_t<<"\t"<< ns_t <<std::endl;
				}
			}
		}
		if(config.is_format("EMoments")){e_moments_f << eta_t << "\t"<< lattice_sum_dint23dy ;}
		if(config.is_format("NMoments")){
			nu_moments_f << eta_t;
			nd_moments_f << eta_t;
		}	

		for (size_t j = 0; j <= n_max; j++) {
			if(config.is_format("EMoments")){
				e_moments_f << "\t" << lattice_sum_eg[j]<< "\t" <<lattice_sum_eq[j];
				e_moments_f << "\t" << lattice_sum_eg_cos[j]<< "\t" <<lattice_sum_eq_cos[j];
				e_moments_f << "\t" << lattice_sum_eg_sin[j]<< "\t" <<lattice_sum_eq_sin[j];
			}
			if(config.is_format("NMoments")){
				nu_moments_f << "\t" << lattice_sum_nu[j]<< "\t" <<lattice_sum_nu_cos[j]<< "\t" << lattice_sum_nu_sin[j];
				nd_moments_f << "\t" << lattice_sum_nd[j]<< "\t" <<lattice_sum_nd_cos[j]<< "\t" << lattice_sum_nd_sin[j];
			}
		}
		if(config.is_format("EMoments")){e_moments_f << std::endl;}
		if(config.is_format("NMoments")){
			nu_moments_f << std::endl;
			nd_moments_f << std::endl;
		}

		if(config.get_Verbose()){
			double percentage_done=double(ieta)/double(config.get_NETA());
			printProgress(percentage_done);
		}
	}

	if(config.is_format("EMoments")){e_moments_f.close();}
	if(config.is_format("EMoments")){
		nu_moments_f.close();
		nd_moments_f.close();
	}

	if(config.is_format("Charges")){charges_f.close();}

	if(config.get_Verbose()){printProgress(1);}

}

void Event::MakeChargeOutput_Transverse(double eta){

	std::ofstream charges_f;
	std::ostringstream chargesname;
	chargesname << OUTPATH <<"/Event_"<<EventID<< "_Charges_eta_"<< eta << ".dat" ;
	charges_f.open(chargesname.str());

	double eg_t,eq_t,nu_t,nd_t,ns_t;
	double t1p_t,t1n_t,t2p_t,t2n_t;
	double t1_t,t2_t;

	for (size_t ix = 0; ix < config.get_NX(); ix++) {
		for (size_t iy = 0; iy < config.get_NY(); iy++) {

			double x_t= get_x(ix);
			double y_t= get_y(iy);

			t1p_t=T1p(ix,iy);
			t1n_t=T1n(ix,iy);
			t2p_t=T2p(ix,iy);
			t2n_t=T2n(ix,iy);
			t1_t = t1p_t+t1n_t;
			t2_t = t2p_t+t2n_t;

			eg_t=ChargeMaker->gluon_energy(eta, t1_t, t2_t)  ;
			eq_t=ChargeMaker->quark_energy(eta, t1_t, t2_t) ;
			nu_t=ChargeMaker->u_density(eta, t1p_t, t1n_t, t2p_t, t2n_t) ;
			nd_t=ChargeMaker->d_density(eta, t1p_t, t1n_t, t2p_t, t2n_t) ;
			ns_t=ChargeMaker->s_density(eta, t1_t, t2_t) ;

			charges_f<< x_t <<"\t"<< y_t <<"\t"<< eg_t  <<"\t"<< eq_t<<"\t"<< nu_t <<"\t"<< nd_t<<"\t"<< ns_t <<std::endl;
		}
		charges_f<<std::endl;
	}

	charges_f.close();
}

void Event::MakeThicknessOutput(){

		std::ofstream charges_f;
	  std::ostringstream chargesname;
	  chargesname << OUTPATH <<"/Event_"<<EventID<< "_Thickness.dat" ;
	  charges_f.open(chargesname.str());

		for (size_t ix = 0; ix < config.get_NX(); ix++) {
			for (size_t iy = 0; iy < config.get_NY(); iy++) {
				double x_t= get_x(ix);
				double y_t= get_y(iy);
				charges_f<< x_t <<"\t"<< y_t <<"\t"<< T1p(ix,iy) <<"\t"<< T1n(ix,iy) <<"\t"<< T2p(ix,iy) <<"\t"<< T2n(ix,iy) <<std::endl;
			}
			charges_f <<std::endl;
		}
		charges_f.close();
}

double Event::compute_internucleon_distance(Nucleus *N1,int ind1,Nucleus *N2,int ind2){
	double x1,y1,x2,y2;
	N1->get_trans_position_nucleon(ind1,x1,y1);
	N2->get_trans_position_nucleon(ind2,x2,y2);

	return sqrt( pow(x1-x2,2)+pow(y1-y2,2) );   /////////CHECK THIS!!!!!!!!!!
}

double Event::ComputeIndividualOverlap_Gaussian(Nucleus *N1,int ind1,Nucleus *N2,int ind2){
		double d_12=Event::compute_internucleon_distance(N1,ind1,N2,ind2);
    return exp(-d_12*d_12/(4.0*config.get_BG()))/(4.0*M_PI*config.get_BG());
}

double Event::ComputeIndividualOverlap_Exponential (Nucleus *N1,int ind1,Nucleus *N2,int ind2){
		double d_12=Event::compute_internucleon_distance(N1,ind1,N2,ind2);
    return d_12; //exp(-dsqr/(4.0*config.get_BG()))/(4.0*M_PI*config.get_BG());
}

double Event::ComputeInteractionProbability(double CollisionOverlap){
    // DETERMINE NUCLEON-NUCLEON INELASTIC CROSS-SECTION //
		double sigma0=2*M_PI*config.get_BG();
		double sigma_inelastic= sigma_inelastic_fm2(config.get_collEnergy());
		double sigmaNN =sigma0 * FifthOrderPoly(sigma_inelastic/sigma0);

		if(!is_sigma_in_range(sigma_inelastic/sigma0)){
			// WARNING MESSAGE //
			if(config.get_Verbose() && is_sigma_output ){
				std::cerr << std::endl;
				std::cerr << "[[ WARNING! ]] -- Solution for SIGMA(SIGMA_INELASTIC) not trusted in the specified range! Re-check assumptions and values! " << std::endl;is_sigma_output=false;
				std::cerr << std::endl;
			}
		}

		if(config.get_Verbose() && is_sigma_output ){
			std::cerr << std::endl;
			std::cerr << "#  SIGMA_INEL=" << sigma_inelastic * fm2_to_mb << " mb,  SIGMA_NN=" << sigmaNN * fm2_to_mb << " mb,     CENTER OF MASS ENERGY " << config.get_collEnergy() << " GEV" << std::endl;is_sigma_output=false;
			std::cerr << std::endl;
		}
    // COMPUTE INTERACTION PROBABILITY //
    return 1.0-exp(-sigmaNN*CollisionOverlap);
}


bool Event::check_interaction_status(Nucleus *N1,int ind1,Nucleus *N2,int ind2){
	// CHECK THAT intERACTION PROBABILITY IS SUFFICIENTLY LARGE TO BE RELEVANT //
	bool is_interacting=false;
	if(config.get_GlauberAcceptance() == GlauberMode::Standard){
		double d_12=Event::compute_internucleon_distance(N1,ind1,N2,ind2);
		is_interacting= d_12 <=  (sqrt( sigma_inelastic_fm2(config.get_collEnergy())/M_PI) ) ;
	}
	else{
			double NNOverlap=0.;
			//Now we compute the overlap
			if(config.get_GlauberAcceptance() == GlauberMode::Gaussian){NNOverlap=ComputeIndividualOverlap_Gaussian(N1,ind1,N2,ind2);}
			else if(config.get_GlauberAcceptance() == GlauberMode::Exponential){NNOverlap=ComputeIndividualOverlap_Exponential(N1,ind1,N2,ind2);}
			//and the interaction probability coming from the overlap
			double NNInteractionProbability=ComputeInteractionProbability(NNOverlap);
			// Sampling to see if interaction happens
			is_interacting=uni_nu_rn()<NNInteractionProbability;
	}
	return is_interacting;
}

void Event::CheckParticipants(Nucleus* N1,Nucleus *N2){

    // NUCLEON-NUCLEON INTERACTION STATUS //
    int **InteractionStatus= new int*[N1->get_A()];
    for(int s1=0;s1<N1->get_A();s1++){InteractionStatus[s1]=new int[N2->get_A()];}

    // COMPUTE NUCLEON-NUCLEON intERACTION PROBABILITIES AND SAMPLE intERACTION //
    for(int s1=0;s1<N1->get_A();s1++){
        for(int s2=0;s2<N2->get_A();s2++){
            if(check_interaction_status(N1,s1,N2,s2)){InteractionStatus[s1][s2]=1;}
            else{InteractionStatus[s1][s2]=0;}
        }
    }

    // CHECK FOR EACH NUCLEON IN NUCLEUS ONE THE NUMBER OF COLLISIONS AND PARTICIPANT STATUS //
    N1->set_NPart(0); N1->set_NColl(0);
    for(int s1=0;s1<N1->get_A();s1++){
        // NUMBER OF INDIVIDUAL COLLISIONS //
        int NIndColl=0;
        for(int s2=0;s2<N2->get_A();s2++){
            if(InteractionStatus[s1][s2]==1){NIndColl++;}
        }
        // UPDATE NUMBER OF COLLISIONS //
        N1->set_CollisionNumber(s1,NIndColl); N1->update_NColl(NIndColl);
        // UPDATE PARTICIPANT STATUS //
        if(NIndColl>0){N1->set_ParticipantStatus(s1,1); N1->add_Participant();}
        else{N1->set_ParticipantStatus(s1,0);}
    }

    // CHECK FOR EACH NUCLEON IN NUCLEUS TWO THE NUMBER OF COLLISIONS AND PARTICIPANT STATUS //
    N2->set_NPart(0); N2->set_NColl(0);

    for(int s2=0;s2<N2->get_A();s2++){
        // NUMBER OF INDIVIDUAL COLLISIONS //
        int NIndColl=0;
        for(int s1=0;s1<N1->get_A();s1++){
            if(InteractionStatus[s1][s2]==1){NIndColl++;}
        }

        // UPDATE NUMBER OF COLLISIONS //
        N2->set_CollisionNumber(s2,NIndColl); N2->update_NColl(NIndColl);
        // UPDATE PARTICIPANT STATUS //
        if(NIndColl>0){N2->set_ParticipantStatus(s2,1); N2->add_Participant();}
        else{N2->set_ParticipantStatus(s2,0);}

    }

    // CLEAN-UP //
    for(int s1=0;s1<N1->get_A();s1++){delete InteractionStatus[s1];}

    delete[] InteractionStatus;
}


void  Event::printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

void Event::print_glauber_data_to_file(Nucleus * N1,Nucleus * N2){
	/* This function prints the MC-Glauber*/
	std::ofstream global_f;
	std::ostringstream global_name;
	global_name << OUTPATH <<"/global_"<<EventID<< ".dat" ;
	global_f.open(global_name.str());

	global_f<<  "N_part_1\t"<< N1->get_number_of_participants() << std::endl;
	global_f<<  "N_part_2\t"<< N2->get_number_of_participants() << std::endl;

	global_f<<  "N_coll_1\t"<< N1->get_NColl() << std::endl;
	global_f<<  "N_coll_2\t"<< N2->get_NColl() << std::endl;

	global_f<<  "bx\t"<< abs(b1[0]-b2[0]) << std::endl;
	global_f<<  "by\t"<< abs(b1[1]-b2[1])  << std::endl;


	global_f.close();
}
