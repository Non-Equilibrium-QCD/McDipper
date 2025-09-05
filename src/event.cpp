/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */

#include <iostream>
#include <math.h>
#include <random>
#include <filesystem>
#include <sstream>
#include "include/random.h"
#include "include/event.h"
#include "include/params_gen.h"

#define OPTICAL 0
#define AVG_EVENT 0
namespace fs = std::filesystem;

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

// Constructor and Destructor

Event::Event(Config ConfInput){
	config=Config(ConfInput);
	
	N1.A=config.get_A(1);N1.Z=config.get_Z(1);N1.mode=config.get_NuclearMode(1);
    N1.is_hotspots_fluct=config.is_hotspots_fluct();
    N1.Nq=config.get_Nq();
    N1.Bq=config.get_Bq();
    N1.Br=config.get_Br();
    
    N1.is_thick_fluct=config.is_thick_fluct();
    N1.sigma=config.get_sigma();
    N1.fluct_mode =config.get_fluct_mode();

	if(N1.mode==3){
		N1.inputFile=config.get_NucleusInput(1);
		N1.IsospinSpecified=config.get_IsospinDefinition(1);
		N1.NConf=config.get_NConf(1);
		N1.is_weights=config.get_weight(1);
	}


	N2.A=config.get_A(2);N2.Z=config.get_Z(2);N2.mode=config.get_NuclearMode(2);
    N2.is_hotspots_fluct=config.is_hotspots_fluct();
    N2.Nq=config.get_Nq();
    N2.Bq=config.get_Bq();
    N2.Br=config.get_Br();
    
    N2.is_thick_fluct=config.is_thick_fluct();
    N2.sigma=config.get_sigma();
    N2.fluct_mode =config.get_fluct_mode();

	if(N2.mode==3){
		N2.inputFile=config.get_NucleusInput(2);
		N2.IsospinSpecified=config.get_IsospinDefinition(2);
		N2.NConf=config.get_NConf(2);
		N2.is_weights=config.get_weight(2);
	}

	PrimalSeed=config.get_seed();
	srand48(PrimalSeed);
    engine.seed(static_cast<Engine::result_type>(PrimalSeed));

	b1=new double[2];
	b2=new double[2];

	sqrtsNN=config.get_collEnergy();

	NX=config.get_NX();
	NY=config.get_NY();
	NETA= config.get_NETA();
	MakeGrid();
	if(config.print_avg_event()>0){InitializeAverageEvent();}

	ChargeMaker = new Charges(config);
	if (!config.is_output_suppressed()) {
		Initialize_output();
	}

	cell_trans_volume=config.get_dX()*config.get_dY();
}

Event::~Event(){
	delete[] b1;
	delete[] b2;
	delete[] T1p_ptr;
	delete[] T1n_ptr;
	delete[] T2p_ptr;
	delete[] T2n_ptr;

	if(config.print_avg_event()>0){
		delete[] EgAvg_ptr;
		delete[] EqAvg_ptr;
		delete[] nuAvg_ptr;
		delete[] ndAvg_ptr;
		delete[] nsAvg_ptr;
	}
	delete ChargeMaker;
}

// Structural Functions
void Event::MakeGrid(){

	T1p_ptr = new double[config.get_NX() * config.get_NY()];
	T1n_ptr = new double[config.get_NX() * config.get_NY()];
	T2p_ptr = new double[config.get_NX() * config.get_NY()];
	T2n_ptr = new double[config.get_NX() * config.get_NY()];
	
	if(config.print_avg_event()>0){
		EgAvg_ptr = new double[config.get_NETA() * config.get_NX() * config.get_NY()];
		EqAvg_ptr = new double[config.get_NETA() * config.get_NX() * config.get_NY()];
		nuAvg_ptr = new double[config.get_NETA() * config.get_NX() * config.get_NY()];
		ndAvg_ptr = new double[config.get_NETA() * config.get_NX() * config.get_NY()];
		nsAvg_ptr = new double[config.get_NETA() * config.get_NX() * config.get_NY()];
	}
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


////////// average event interface 

const double& Event::EgAvg (int64_t ieta, int64_t nx, int64_t ny) const {return (EgAvg_ptr)[ NETA*(NY*nx + ny)+ieta];}
double& Event:: EgAvg (int64_t ieta, int64_t nx, int64_t ny){return (EgAvg_ptr)[ NETA*(NY*nx + ny)+ieta];}

const double& Event::EqAvg (int64_t ieta, int64_t nx, int64_t ny) const {return (EqAvg_ptr)[ NETA*(NY*nx + ny)+ieta];}
double& Event:: EqAvg (int64_t ieta, int64_t nx, int64_t ny){return (EqAvg_ptr)[ NETA*(NY*nx + ny)+ieta];}

const double& Event::nuAvg (int64_t ieta, int64_t nx, int64_t ny) const {return (nuAvg_ptr)[ NETA*(NY*nx + ny)+ieta];}
double& Event:: nuAvg (int64_t ieta, int64_t nx, int64_t ny){return (nuAvg_ptr)[ NETA*(NY*nx + ny)+ieta];}

const double& Event::ndAvg (int64_t ieta, int64_t nx, int64_t ny) const {return (ndAvg_ptr)[ NETA*(NY*nx + ny)+ieta];}
double& Event:: ndAvg (int64_t ieta, int64_t nx, int64_t ny){return (ndAvg_ptr)[ NETA*(NY*nx + ny)+ieta];}

const double& Event::nsAvg (int64_t ieta, int64_t nx, int64_t ny) const {return (nsAvg_ptr)[ NETA*(NY*nx + ny)+ieta];}
double& Event:: nsAvg (int64_t ieta, int64_t nx, int64_t ny){return (nsAvg_ptr)[ NETA*(NY*nx + ny)+ieta];}

////////////////////////////////

// Evaluation
// Functions
		
void Event::NewEvent(int EventID_in, Nucleus *A1, Nucleus *A2){
	EventID=EventID_in;

	// Create Nuclei
	bool is_event_valid=false;
	if(config.get_Verbose()){
		std::cout<< "|---------------------------------- New Event: "<< EventID+1<<"/"<< config.get_NEvents() << " -------------------------------------|\n";
		std::cout<< "                       Event ID: " <<EventID <<std::endl;
	}

	while(!is_event_valid){
		A1->refresh_positions();
		A2->refresh_positions();
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

		A1->shift_nucleus_by_impact(b1[0],b1[1]);
		A2->shift_nucleus_by_impact(b2[0],b2[1]);

		#if OPTICAL==1
		if(config.get_Verbose()){std::cout<< "    Running for impact parameter ( bx = " << b1[0]-b2[0] << " , by = " << b1[1]-b2[1] << ")"<<std::endl ;}
		is_event_valid=true;
		for (int ix = 0; ix < config.get_NX(); ix++) {
			for (int iy = 0; iy < config.get_NY(); iy++) {
				double x_t= get_x(ix);
				double y_t= get_y(iy);
				T1p (ix,iy)= A1->get_Z()*A1->nuclear_thickness_optical(x_t-b1[0], y_t-b1[1]);
				T1n (ix,iy)= (A1->get_A()-A1->get_Z() )*A1->nuclear_thickness_optical(x_t-b1[0], y_t-b1[1]);
				T2p (ix,iy)= A2->get_Z()*A2->nuclear_thickness_optical(x_t-b2[0], y_t-b2[1]);
				T2n (ix,iy)= (A2->get_A()-A2->get_Z() )*A2->nuclear_thickness_optical(x_t-b2[0], y_t-b2[1]);
			}
		}
		MakeGlobalQuantities();
		bool is_charge_output= config.is_format("EMoments") || config.is_format("NMoments") || config.is_format("Charges");
		if (is_charge_output){MakeChargeOutput();}

		#elif OPTICAL==0
		CheckParticipants(A1,A2);

		is_event_valid = (A1->get_number_of_participants()>0) && (A2->get_number_of_participants()>0);
		if(is_event_valid){
			if(config.get_Verbose()){
				std::cout<< "    Running for impact parameter ( bx = " << b1[0]-b2[0] << " , by = " << b1[1]-b2[1] << ")"<<std::endl ;
				std::cout<< "    Number of Participants: " <<A1->get_number_of_participants() + A2->get_number_of_participants() <<  "\n";
				std::cout<< "    Number of Collisions: " <<A1->get_NColl()+A2->get_NColl() <<  "\n";
			}

			if (config.is_format("Nuclei") ){print_nucleon_pos(A1,A2,EventID);}

			if(config.get_Verbose()){std::cout<<"    Nucleon Positions written out\n";}
            if (config.is_thick_fluct()){
              A1->Thickness_fluct();
              A2->Thickness_fluct();
            }
			for (int ix = 0; ix < config.get_NX(); ix++) {
				for (int iy = 0; iy < config.get_NY(); iy++) {
					double x_t= get_x(ix);
					double y_t= get_y(iy);
					T1p (ix,iy)= A1->GetThickness_p(x_t,y_t, config.get_BG());
					T1n (ix,iy)= A1->GetThickness_n(x_t,y_t, config.get_BG());
					T2p (ix,iy)= A2->GetThickness_p(x_t,y_t, config.get_BG());
					T2n (ix,iy)= A2->GetThickness_n(x_t,y_t, config.get_BG());
				}
			}
			print_glauber_data_to_file(A1,A2);
			//Print weights for the NLEFT cases.
			if(config.get_weight(1) || config.get_weight(2)){print_weights(A1,A2,EventID);}

			MakeGlobalQuantities();
			bool is_charge_output= config.is_format("EMoments") || config.is_format("NMoments") || config.is_format("Charges");
			if (is_charge_output){
				if(config.is_boost_invariant() ) {MakeChargeOutputMidrapidity();}
				else{MakeChargeOutput();}
			}
			
		}
		#endif
	}

	if(config.get_Verbose()){
		std::cout<<"    Thickness Functions computed!\n";
		std::cout<< "\n|------------------------------------- End Event ----------------------------------------|\n";
	}
}

void Event::MakeBlock(int EventID){
		std::ofstream block_f;
	std::ostringstream block_name;
	block_name << OUTPATH <<"/block.txt" ;
	block_f.open(block_name.str());
	block_f<< EventID;
	block_f.close();
}

void Event::RemoveBlock(){
	std::ostringstream block_name;
	block_name << OUTPATH <<"/block.txt" ;
	std::filesystem::remove(block_name.str());
}

void Event::MakeEventByEvent(){

	/* Create the nuclei for the evaluation */
	Nucleus A1=Nucleus(N1);
	Nucleus A2=Nucleus(N2);

	
	if( ev0 < config.get_NEvents()){
		std::cerr<< "[ Warning::Event ]: Event generation starting from event # "<< ev0 << std::endl;
		for (int ev = ev0; ev < config.get_NEvents(); ev++) {
			NewEvent(get_ID(ev,config.get_NEvents()),&A1,&A2);
			MakeBlock(ev);
		}
	}
	else{std::cerr<< "[ Error::Event ]: Initial event ID="<< ev0 << " larger than Nevents in config file. Exiting."<<std::endl;}
	
	if(config.print_avg_event()>0){
		MakeGlobalQuantities_AverageEvent();
		MakeChargeOutput_AverageEvent();
	}
	// RemoveBlock();
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
void Event::MakeOutputFiles(){
	/// Make Nuclei File
	std::ofstream pos_f;
	std::ostringstream posname;
	posname << OUTPATH << "/Nuclei.dat";
	pos_f.open(posname.str());
	pos_f.close();	

	// Make Global Data
	std::ofstream global_f;
	std::ostringstream global_name;
	global_name << OUTPATH <<"/Global.dat" ;
	global_f.open(global_name.str());
	global_f.close();

	std::ofstream e_moments_f; std::ostringstream e_moments_name;
	std::ofstream nu_moments_f; std::ostringstream nu_moments_name;
	std::ofstream nd_moments_f; std::ostringstream nd_moments_name;

	if(config.get_weight(1) || config.get_weight(2)){
		std::ofstream weights_f;
		std::ostringstream weights_name;
		weights_name << OUTPATH <<"/Weights.dat" ;
		weights_f.open(weights_name.str());
		weights_f <<"# EventID \t W(1) \t W(2) \t W(Event) " << std::endl;
		weights_f.close();
	}

	if(config.is_format("EMoments")){
		e_moments_name << OUTPATH <<"/EnergyMoments.dat" ;
		e_moments_f.open(e_moments_name.str());
		e_moments_f.close();
	}

	if(config.is_format("NMoments")){
		nu_moments_name << OUTPATH <<"/NuMoments.dat" ;
		nu_moments_f.open(nu_moments_name.str());
		nu_moments_f.close();

		nd_moments_name << OUTPATH <<"/NdMoments.dat" ;
		nd_moments_f.open(nd_moments_name.str());
		nd_moments_f.close();
	}

}



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
			std::ostringstream block_name;
			block_name << OUTPATH <<"/block.txt" ;
	
			std::ifstream f(block_name.str());
			if ( !f.is_open() ) {std::cerr << "Error opening the file!";}

			std::string s;
			getline(f, s);
			ev0=std::stoi(s) +1;
			f.close();

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

void Event::print_nucleon_pos(Nucleus *A1,Nucleus *A2, int EventID){

	std::ofstream pos_f;
	std::ostringstream posname;
	posname << OUTPATH << "/Nuclei.dat";
	pos_f.open(posname.str(), std::ios_base::app);
	pos_f<<"# Start Event "<< EventID << std::endl;
	double x_t,y_t,z_t;
	for (int i = 0; i < A1->get_A(); i++) {
		A1->get_position_nucleon(i, x_t,y_t,z_t);
		pos_f<< 1 <<"\t"<< i <<"\t"<< x_t<<"\t"<< y_t<<"\t"<< z_t<<"\t"<< A1->get_ParticipantStatus(i)<<"\t"<< A1->get_nucleon_type(i) <<std::endl;
	}
	for (int i = 0; i < A2->get_A(); i++) {
		A2->get_position_nucleon(i, x_t,y_t,z_t);
		pos_f<< 2 <<"\t"<< i <<"\t"<< x_t<<"\t"<< y_t<<"\t"<< z_t<<"\t"<< A2->get_ParticipantStatus(i)<<"\t"<< A2->get_nucleon_type(i)<<std::endl;
	}
	pos_f<<"# End Event "<< EventID << std::endl;
	pos_f.close();	

}

void Event::MakeGlobalQuantities(){

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

	for (int ieta = 0; ieta < config.get_NETA(); ieta++) {
		double eta_t = ieta*config.get_dETA() + config.get_ETAMIN();
		for (int ix = 0; ix < config.get_NX(); ix++) {
			for (int iy = 0; iy < config.get_NY(); iy++) {

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

				if(config.print_avg_event()>0)
				{
					EgAvg (ieta, ix, iy)+=eg_t/config.get_NEvents();
					EqAvg (ieta, ix, iy)+=eq_t/config.get_NEvents();
					nuAvg (ieta, ix, iy)+=nu_t/config.get_NEvents();
					ndAvg (ieta, ix, iy)+=nd_t/config.get_NEvents();
					nsAvg (ieta, ix, iy)+=ns_t/config.get_NEvents();
				}

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

	std::ofstream global_f;
	std::ostringstream global_name;
	global_name << OUTPATH <<"/Global.dat" ;
	global_f.open(global_name.str(), std::ios_base::app);
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
	global_f<<"# End Event "<< EventID << std::endl;
	global_f.close();


	if(config.get_Verbose()){
		std::cout << "\n\n    E_g = " << total_energy_event_g << "\tE_q = " << total_energy_event_q << "\tE_tot = " << total_energy_event<<std::endl ;
		std::cout<<  "    dE_deta = " << dEdeta[izero]<<  "\t with \tdEg_deta = " << dEgdeta[izero]<<  "    and    dEq_deta = " << dEqdeta[izero]<<std::endl;
		std::cout << "    N_u = " << total_u_event<< "\tN_d = " << total_d_event << "\tN_s = " << total_s_event <<std::endl ;
		std::cout<<  "    dN_i_deta(eta=0)     u: " << dnudeta[izero];
		std::cout<<  "\t d: " << dnddeta[izero]<< "\t s: " << dnsdeta[izero]<<std::endl;
		std::cout << "    Q = " << total_q_event<< "\tB = " << total_B_event <<std::endl ;
		std::cout << "\nGlobal CoM coordinates:  ( x = " << x_cm_global << " ,y = " << y_cm_global<< ")"<<std::endl ;
	}

	

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
		e_moments_name << OUTPATH <<"/EnergyMoments.dat" ;
		e_moments_f.open(e_moments_name.str(), std::ios_base::app);
		e_moments_f<<"# Start Event "<< EventID << std::endl;
	}

	if(config.is_format("NMoments")){
		nu_moments_name << OUTPATH <<"/NuMoments.dat" ;
		nu_moments_f.open(nu_moments_name.str(), std::ios_base::app);
		nu_moments_f<<"# Start Event "<< EventID << std::endl;

		nd_moments_name << OUTPATH <<"/NdMoments.dat" ;
		nd_moments_f.open(nd_moments_name.str(), std::ios_base::app);
		nd_moments_f<<"# Start Event "<< EventID << std::endl;
	}

	double eg_t,eq_t,nu_t,nd_t,ns_t;
	double t1p_t,t1n_t,t2p_t,t2n_t;
	double t1_t,t2_t;


	if(config.get_Verbose()){std::cout<< "\nWriting out charges and moments for Event " << EventID <<std::endl;}

	for (int ieta = 0; ieta < config.get_NETA(); ieta++) {
		double eta_t = ieta*config.get_dETA() + config.get_ETAMIN();

		double lattice_sum_dint23dy = 0.;
		double lattice_sum_dint43dy = 0.;

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


		for (int ix = 0; ix < config.get_NX(); ix++) {
			for (int iy = 0; iy < config.get_NY(); iy++) {

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
					lattice_sum_dint43dy +=  cell_trans_volume* pow(gen_pars::GeV2_to_fmm2 , 2.0)*pow( gen_pars::fmm2_to_GeV2 * (eg_t + eq_t), 4./3.);  // in fm^{-4}
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
					charges_f<< eta_t <<"\t"<< get_x(ix) <<"\t"<< get_y(iy) <<"\t"<< eg_t  <<"\t"<< eq_t<<"\t"<< nu_t <<"\t"<< nd_t<<"\t"<< ns_t <<std::endl;
				}
			}
		}
		if(config.is_format("EMoments")){e_moments_f << eta_t << "\t"<< lattice_sum_dint23dy<< "\t"<< lattice_sum_dint43dy ;}
		if(config.is_format("NMoments")){
			nu_moments_f << eta_t;
			nd_moments_f << eta_t;
		}	

		for (int j = 0; j <= n_max; j++) {
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

	if(config.is_format("EMoments")){
		e_moments_f<<"# End Event "<< EventID << std::endl;
		e_moments_f.close();}
	if(config.is_format("EMoments")){
		nu_moments_f<<"# End Event "<< EventID << std::endl;
		nd_moments_f<<"# End Event "<< EventID << std::endl;
		nu_moments_f.close();
		nd_moments_f.close();
	}

	if(config.is_format("Charges")){charges_f.close();}

	if(config.get_Verbose()){printProgress(1);}

}


void Event::MakeChargeOutputMidrapidity(){
	std::ofstream charges_f; std::ostringstream charges_name;
	std::ofstream e_moments_f; std::ostringstream e_moments_name;
	std::ofstream nu_moments_f; std::ostringstream nu_moments_name;
	std::ofstream nd_moments_f; std::ostringstream nd_moments_name;

	if(config.is_format("Charges")){
		charges_name << OUTPATH <<"/Event_"<<EventID<< "_charges.dat" ;
		charges_f.open(charges_name.str());
	}
	
	if(config.is_format("EMoments")){
		e_moments_name << OUTPATH <<"/EnergyMoments.dat" ;
		if(EventID==0){ 
			e_moments_f.open(e_moments_name.str());
			e_moments_f<< "# Moment Output at y=0"<< std::endl;
		}
		else{e_moments_f.open(e_moments_name.str(),std::ofstream::app);}
	}

	if(config.is_format("NMoments")){
		nu_moments_name << OUTPATH <<"/NuMoments.dat" ;
		if(EventID==0){ 
			nu_moments_f.open(nu_moments_name.str());
			nu_moments_f<< "# Moment Output at y=0"<< std::endl;
		}
		else{nu_moments_f.open(nu_moments_name.str(),std::ofstream::app);}

		nd_moments_name << OUTPATH <<"/NdMoments.dat" ;
		if(EventID==0){ 
			nd_moments_f.open(nd_moments_name.str());
			nd_moments_f<< "# Moment Output at y=0"<< std::endl;
		}
		else{nd_moments_f.open(nd_moments_name.str(),std::ofstream::app);}
	}

	double eg_t,eq_t,nu_t,nd_t,ns_t;
	double t1p_t,t1n_t,t2p_t,t2n_t;
	double t1_t,t2_t;

	double eta_t =0.0;


	if(config.get_Verbose()){std::cout<< "\nWriting out -midrapidity- charges and moments for Event " << EventID <<std::endl;}
	double lattice_sum_dint23dy = 0.;
	double lattice_sum_dint43dy = 0.;

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


	for (int ix = 0; ix < config.get_NX(); ix++) {
		for (int iy = 0; iy < config.get_NY(); iy++) {

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
				lattice_sum_dint43dy +=  cell_trans_volume* pow(gen_pars::GeV2_to_fmm2 , 2.0) * pow( gen_pars::fmm2_to_GeV2 * (eg_t + eq_t), 4./3.);
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
	if(config.is_format("EMoments")){e_moments_f << EventID << "\t"<< lattice_sum_dint23dy<< "\t"<< lattice_sum_dint43dy ;}
	if(config.is_format("NMoments")){
		nu_moments_f << EventID;
		nd_moments_f << EventID;
	}	

	for (int j = 0; j <= n_max; j++) {
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

	for (int ix = 0; ix < config.get_NX(); ix++) {
		for (int iy = 0; iy < config.get_NY(); iy++) {

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

		for (int ix = 0; ix < config.get_NX(); ix++) {
			for (int iy = 0; iy < config.get_NY(); iy++) {
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
	if (!config.is_hotspots_fluct()){
        double d_12=Event::compute_internucleon_distance(N1,ind1,N2,ind2);
        return exp(-d_12*d_12/(4.0*config.get_BG())) / (4.0*M_PI*config.get_BG());
      }
    else{
        double NNOverlap = 0.0;
	    double x1,y1,x2,y2,bsq_hotspot; double Bq=config.get_Bq();
        for (int i=0; i<N1->get_Nq(); i++){
            for (int j=0; j<N2->get_Nq(); j++){
	            N1->get_trans_position_hotspot(ind1,x1,y1,i);
	            N2->get_trans_position_hotspot(ind2,x2,y2,j);
                bsq_hotspot=pow(x1-x2,2)+pow(y1-y2,2);

                NNOverlap+=exp(-bsq_hotspot/(4.0*Bq));//suppose A and B have the same hotspots width.
              }  
          }
        NNOverlap /=  4.0*M_PI*Bq*N1->get_Nq()*N2->get_Nq();
        return NNOverlap;
    }
}

double Event::ComputeIndividualOverlap_Exponential (Nucleus *N1,int ind1,Nucleus *N2,int ind2){
    //left for furture
		double d_12=Event::compute_internucleon_distance(N1,ind1,N2,ind2);
    return d_12; //exp(-dsqr/(4.0*config.get_BG()))/(4.0*M_PI*config.get_BG());
}

double Event::ComputeInteractionProbability(double CollisionOverlap){
    // DETERMINE NUCLEON-NUCLEON INELASTIC CROSS-SECTION //
		double sigma0=4*M_PI*config.get_BG();
		double sigma_inelastic= sigma_inelastic_fm2(config.get_collEnergy());
        double sigmaNN=0.0;

        if (!config.is_hotspots_fluct())
		{ sigmaNN =sigma0 * FitInverseInelXSec(sigma_inelastic/sigma0);}
        else if (config.get_Nq()==3)
        { sigmaNN =FitInverse_parton_sigma(sigma_inelastic); }
        else 
        {std::cerr << "Fitted cross-section for hotspots only apply for nq=3 case." << std::endl;
         std::cerr << "The parton cross section should be fit inside the programe." <<std::endl;}

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
	// CHECK THAT INTERACTION PROBABILITY IS SUFFICIENTLY LARGE TO BE RELEVANT //
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
	global_name << OUTPATH <<"/Global.dat" ;
	global_f.open(global_name.str(), std::ios_base::app);
	global_f<<"# Start Event "<< EventID << std::endl;

	global_f<<  "N_part_1\t"<< N1->get_number_of_participants() << std::endl;
	global_f<<  "N_part_2\t"<< N2->get_number_of_participants() << std::endl;

	global_f<<  "N_coll_1\t"<< N1->get_NColl() << std::endl;
	global_f<<  "N_coll_2\t"<< N2->get_NColl() << std::endl;

	global_f<<  "bx\t"<< abs(b1[0]-b2[0]) << std::endl;
	global_f<<  "by\t"<< abs(b1[1]-b2[1])  << std::endl;


	global_f.close();
}

// CREATE A NUCLEUS OBJECT - THIS IS FOR THE USE AS AN EXTERNAL LIBRARY //
Nucleus Event::CreateNucleusObject(NucStruct N){
	Nucleus Nucl = Nucleus(N);
	return Nucl;
}

// CREATE A DENSITY PROFILE FOR AN EXTERNALLY DEFINED GRID - THIS IS FOR THE USE AS AN EXTERNAL LIBRARY //
/*
- 'ExtGrid' is a struct that contains NX_EXT, NY_EXT, NZ_EXT, XMIN_EXT, XMAX_EXT, YMIN_EXT, YMAX_EXT, ETAMIN_EXT, ETAMAX_EXT.
  The external grid can have equal dimensions as the McDipper grid speciefied in the config, but it can also be
  smaller / coarser than it.
- 'density' is a pointer to an array of size NX_EXT*NY_EXT*NZ_EXT, which is filled with the density values by the function.
- 'mode' is the parameter for the type of the density returned by the function:
	- 0 = energy density
	- 1 = baryon density
	- 2 = charge density
	- 3 = strangeness density
	- any other value will cause an error and the code exits
- All densities are returned in units of [fm^-3]
*/
void Event::EventDensityCustomGrid(ExternalGrid ExtGrid, double *density, int mode){
	double T1p_tmp;
	double T1n_tmp;
	double T2p_tmp;
	double T2n_tmp;

	// Create Nuclei
	bool is_event_valid=false;
	while(!is_event_valid){

		Nucleus A1(N1);
		Nucleus A2(N2);

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

		CheckParticipants(&A1,&A2);
		
		is_event_valid = (A1.get_number_of_participants()>0) && (A2.get_number_of_participants()>0);
		if(is_event_valid){
			double dETA_EXT = (ExtGrid.ETAMAX_EXT-ExtGrid.ETAMIN_EXT)/(ExtGrid.NETA_EXT-1);
			double dX_EXT = (ExtGrid.XMAX_EXT-ExtGrid.XMIN_EXT)/(ExtGrid.NX_EXT-1);
			double dY_EXT = (ExtGrid.YMAX_EXT-ExtGrid.YMIN_EXT)/(ExtGrid.NY_EXT-1);
			int super_index;
			double eg_t,eq_t,nu_t,nd_t,ns_t;

			if (config.is_thick_fluct()){
				A1.Thickness_fluct();
				A2.Thickness_fluct();
			}

			double BG = config.get_BG();
			double KFactor = config.get_KFactor();
			for (int ieta = 0; ieta < ExtGrid.NETA_EXT; ieta++) {
        		double eta = ExtGrid.ETAMIN_EXT + ieta * dETA_EXT;

		        for (int ix = 0; ix < ExtGrid.NX_EXT; ix++) {
        		    double x = ExtGrid.XMIN_EXT + ix * dX_EXT;

            		for (int iy = 0; iy < ExtGrid.NY_EXT; iy++) {
                		double y = ExtGrid.YMIN_EXT + iy * dY_EXT;
						
						T1p_tmp = A1.GetThickness_p(x,y,BG);
						T1n_tmp = A1.GetThickness_n(x,y,BG);
						T2p_tmp = A2.GetThickness_p(x,y,BG);
						T2n_tmp = A2.GetThickness_n(x,y,BG);

						super_index = ix+ExtGrid.NX_EXT*iy+ExtGrid.NX_EXT*ExtGrid.NY_EXT*ieta;
						if(mode == 0){ // energy density
							eg_t = ChargeMaker->gluon_energy(eta, T1p_tmp+T1n_tmp, T2p_tmp+T2n_tmp) * KFactor;
							eq_t = ChargeMaker->quark_energy(eta, T1p_tmp+T1n_tmp, T2p_tmp+T2n_tmp);
							density[super_index] = (eg_t + eq_t) * gen_pars::GeV_to_fmm1; // final unit [fm^-3]
						} else if(mode == 1){ // baryon density
							nu_t = ChargeMaker->u_density(eta, T1p_tmp, T1n_tmp, T2p_tmp, T2n_tmp);
							nd_t = ChargeMaker->d_density(eta, T1p_tmp, T1n_tmp, T2p_tmp, T2n_tmp);
							density[super_index] = ((nu_t+nd_t)/3.) * gen_pars::GeV_to_fmm1; 
						} else if(mode == 2){ // charge density
							nu_t = ChargeMaker->u_density(eta, T1p_tmp, T1n_tmp, T2p_tmp, T2n_tmp);
							nd_t = ChargeMaker->d_density(eta, T1p_tmp, T1n_tmp, T2p_tmp, T2n_tmp);
							density[super_index] = ((2.*nu_t-nd_t)/3.) * gen_pars::GeV_to_fmm1;
						} else if(mode == 3){ // strangeness density
							ns_t = ChargeMaker->s_density(eta, T1p_tmp+T1n_tmp, T2p_tmp+T2n_tmp);
							density[super_index] = ns_t * gen_pars::GeV_to_fmm1;
						} else {
							std::cout << "Requested mode in 'EventDensityCustomGrid' not available." << std::endl;
							exit(1);
						}
					}
				}
			}
		}
	}
}

void Event::EventDensityCustomGridQuarkAndGluonContribution(ExternalGrid ExtGrid, double *densityGluons, double *densityQuarks){
	double T1p_tmp;
	double T1n_tmp;
	double T2p_tmp;
	double T2n_tmp;

	// Create Nuclei
	bool is_event_valid=false;
	while(!is_event_valid){

		Nucleus A1(N1);
		Nucleus A2(N2);

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

		CheckParticipants(&A1,&A2);
		
		is_event_valid = (A1.get_number_of_participants()>0) && (A2.get_number_of_participants()>0);
		if(is_event_valid){
			double dETA_EXT = (ExtGrid.ETAMAX_EXT-ExtGrid.ETAMIN_EXT)/(ExtGrid.NETA_EXT-1);
			double dX_EXT = (ExtGrid.XMAX_EXT-ExtGrid.XMIN_EXT)/(ExtGrid.NX_EXT-1);
			double dY_EXT = (ExtGrid.YMAX_EXT-ExtGrid.YMIN_EXT)/(ExtGrid.NY_EXT-1);
			int super_index;
			double eg_t,eq_t,nu_t,nd_t,ns_t;

			if (config.is_thick_fluct()){
				A1.Thickness_fluct();
				A2.Thickness_fluct();
			}

			double BG = config.get_BG();
			double KFactor = config.get_KFactor();
			std::cout << "KFactor: " << KFactor << std::endl;
			for (int ieta = 0; ieta < ExtGrid.NETA_EXT; ieta++) {
				double eta = ExtGrid.ETAMIN_EXT + ieta * dETA_EXT;

				for (int ix = 0; ix < ExtGrid.NX_EXT; ix++) {
					double x = ExtGrid.XMIN_EXT + ix * dX_EXT;

					for (int iy = 0; iy < ExtGrid.NY_EXT; iy++) {
						double y = ExtGrid.YMIN_EXT + iy * dY_EXT;
						
						T1p_tmp = A1.GetThickness_p(x,y,BG);
						T1n_tmp = A1.GetThickness_n(x,y,BG);
						T2p_tmp = A2.GetThickness_p(x,y,BG);
						T2n_tmp = A2.GetThickness_n(x,y,BG);

						super_index = ix+ExtGrid.NX_EXT*iy+ExtGrid.NX_EXT*ExtGrid.NY_EXT*ieta;
						// energy density contributions for quarks and gluons
						eg_t = ChargeMaker->gluon_energy(eta, T1p_tmp+T1n_tmp, T2p_tmp+T2n_tmp) * KFactor;
						eq_t = ChargeMaker->quark_energy(eta, T1p_tmp+T1n_tmp, T2p_tmp+T2n_tmp);
						densityGluons[super_index] = eg_t * gen_pars::GeV_to_fmm1;
						densityQuarks[super_index] = eq_t * gen_pars::GeV_to_fmm1;
					}
				}
			}
		}
	}
}

void Event::print_weights(Nucleus * N1,Nucleus * N2, int EventID){
	/* This function prints the weights in a compact way */
	std::ofstream weights_f;
	std::ostringstream weights_name;
	weights_name << OUTPATH <<"/Weights.dat" ;
	weights_f.open(weights_name.str(), std::ios_base::app);
	weights_f << EventID<<"\t"<<N1->get_weight()<<"\t"<< N2->get_weight() <<"\t"<< N1->get_weight()*N2->get_weight() << std::endl;
	weights_f.close();
}

void Event::InitializeAverageEvent(){
	for (int ieta = 0; ieta < config.get_NETA(); ieta++) {
		for (int ix = 0; ix < config.get_NX(); ix++) {
			for (int iy = 0; iy < config.get_NY(); iy++) {
				EgAvg (ieta, ix, iy)=0.0;
				EqAvg (ieta, ix, iy)=0.0;
				nuAvg (ieta, ix, iy)=0.0;
				ndAvg (ieta, ix, iy)=0.0;
				nsAvg (ieta, ix, iy)=0.0;
			}
		}
	}
}


void Event::MakeGlobalQuantities_AverageEvent(){

	std::ofstream global_f;
	std::ostringstream global_name;
	global_name << OUTPATH <<"/AverageEvent_global.dat" ;
	global_f.open(global_name.str(), std::ios_base::app);

	double eg_t,eq_t,nu_t,nd_t,ns_t;

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

	if(config.get_Verbose()){std::cout<< "\n Making global observables and observables for Average Event "<<std::endl;}

	for (int ieta = 0; ieta < config.get_NETA(); ieta++) {
		// double eta_t = ieta*config.get_dETA() + config.get_ETAMIN();
		for (int ix = 0; ix < config.get_NX(); ix++) {
			for (int iy = 0; iy < config.get_NY(); iy++) {

				double x_t= get_x(ix);
				double y_t= get_y(iy);


				eg_t= EgAvg (ieta, ix, iy);
				eq_t= EqAvg (ieta, ix, iy);
				nu_t= nuAvg (ieta, ix, iy);
				nd_t= ndAvg (ieta, ix, iy);
				ns_t= nsAvg (ieta, ix, iy);

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


void Event::MakeChargeOutput_AverageEvent(){
	std::ofstream charges_f; std::ostringstream charges_name;
	std::ofstream e_moments_f; std::ostringstream e_moments_name;
	std::ofstream nu_moments_f; std::ostringstream nu_moments_name;
	std::ofstream nd_moments_f; std::ostringstream nd_moments_name;

	if(config.print_avg_event()>1){
		charges_name << OUTPATH <<"/AverageEvent_Event_charges.dat" ;
		charges_f.open(charges_name.str());
	}
	
	if(config.print_avg_event()>0){
		e_moments_name << OUTPATH <<"/AverageEvent_EnergyMoments.dat" ;
		e_moments_f.open(e_moments_name.str());
	}

	if(config.print_avg_event()>0){
		nu_moments_name << OUTPATH <<"/AverageEvent_NuMoments.dat" ;
		nu_moments_f.open(nu_moments_name.str());

		nd_moments_name << OUTPATH <<"/AverageEvent_NdMoments.dat" ;
		nd_moments_f.open(nd_moments_name.str());
	}

	double eg_t,eq_t,nu_t,nd_t,ns_t;


	if(config.get_Verbose()){std::cout<< "\nWriting out charges and moments for Average Event "<<std::endl;}

	for (int ieta = 0; ieta < config.get_NETA(); ieta++) {
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


		for (int ix = 0; ix < config.get_NX(); ix++) {
			for (int iy = 0; iy < config.get_NY(); iy++) {

				double x_t= get_x(ix)-x_cm_global;
				double y_t= get_y(iy)-y_cm_global;

				eg_t= EgAvg (ieta, ix, iy);
				eq_t= EqAvg (ieta, ix, iy);
				nu_t= nuAvg (ieta, ix, iy);
				nd_t= ndAvg (ieta, ix, iy);
				ns_t= nsAvg (ieta, ix, iy);


				double r_t = sqrt(x_t*x_t + y_t*y_t);
				double phi_t = atan2(y_t,x_t);
				if(config.print_avg_event()>0){
					lattice_sum_dint23dy +=  cell_trans_volume* gen_pars::GeV2_to_fmm2*pow( gen_pars::fmm2_to_GeV2 * (eg_t + eq_t), 2./3.);
				}	

				for (int j = 0; j <= n_max; j++)
				{
					double jj = (double) j;
					if(config.print_avg_event()>0){
						lattice_sum_eg[j] += cell_trans_volume*pow(r_t,jj) * eg_t;
						lattice_sum_eq[j] += cell_trans_volume*pow(r_t,jj) * eq_t;
						//
						lattice_sum_eg_cos[j] += cell_trans_volume*pow(r_t,jj) * cos(jj*phi_t) * eg_t ;
						lattice_sum_eq_cos[j] += cell_trans_volume*pow(r_t,jj) * cos(jj*phi_t) * eq_t;

						lattice_sum_eg_sin[j] += cell_trans_volume*pow(r_t,jj) * sin(jj*phi_t) * eg_t;
						lattice_sum_eq_sin[j] += cell_trans_volume*pow(r_t,jj) * sin(jj*phi_t) * eq_t;
					}
					if(config.print_avg_event()>0){
						lattice_sum_nu[j] += cell_trans_volume*pow(r_t,jj) * nu_t;
						lattice_sum_nd[j] += cell_trans_volume*pow(r_t,jj) * nd_t;
						//
						lattice_sum_nu_cos[j] += cell_trans_volume*pow(r_t,jj) * cos(jj*phi_t) * nu_t ;
						lattice_sum_nd_cos[j] += cell_trans_volume*pow(r_t,jj) * cos(jj*phi_t) * nd_t;

						lattice_sum_nu_sin[j] += cell_trans_volume*pow(r_t,jj) * sin(jj*phi_t) * nu_t;
						lattice_sum_nd_sin[j] += cell_trans_volume*pow(r_t,jj) * sin(jj*phi_t) * nd_t;
					}
				}

				if(config.print_avg_event()>1){
					charges_f<< eta_t <<"\t"<< x_t <<"\t"<< y_t <<"\t"<< eg_t  <<"\t"<< eq_t<<"\t"<< nu_t <<"\t"<< nd_t<<"\t"<< ns_t <<std::endl;
				}
			}
		}
		if(config.print_avg_event()>0){e_moments_f << eta_t << "\t"<< lattice_sum_dint23dy ;}
		if(config.print_avg_event()>0){
			nu_moments_f << eta_t;
			nd_moments_f << eta_t;
		}	

		for (int j = 0; j <= n_max; j++) {
			if(config.print_avg_event()>0){
				e_moments_f << "\t" << lattice_sum_eg[j]<< "\t" <<lattice_sum_eq[j];
				e_moments_f << "\t" << lattice_sum_eg_cos[j]<< "\t" <<lattice_sum_eq_cos[j];
				e_moments_f << "\t" << lattice_sum_eg_sin[j]<< "\t" <<lattice_sum_eq_sin[j];
			}
			if(config.print_avg_event()>0){
				nu_moments_f << "\t" << lattice_sum_nu[j]<< "\t" <<lattice_sum_nu_cos[j]<< "\t" << lattice_sum_nu_sin[j];
				nd_moments_f << "\t" << lattice_sum_nd[j]<< "\t" <<lattice_sum_nd_cos[j]<< "\t" << lattice_sum_nd_sin[j];
			}
		}
		if(config.print_avg_event()>0){e_moments_f << std::endl;}
		if(config.print_avg_event()>0){
			nu_moments_f << std::endl;
			nd_moments_f << std::endl;
		}

		if(config.get_Verbose()){
			double percentage_done=double(ieta)/double(config.get_NETA());
			printProgress(percentage_done);
		}
	}

	if(config.print_avg_event()>0){e_moments_f.close();}
	if(config.print_avg_event()>0){
		nu_moments_f.close();
		nd_moments_f.close();
	}

	if(config.print_avg_event()>1){charges_f.close();}

	if(config.get_Verbose()){printProgress(1);}

}


