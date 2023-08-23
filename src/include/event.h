/* Copyright (c) 2023 Oscar Garcia-Montero
 * All rights reserved. */
#ifndef EVENT_H
#define EVENT_H
#include <iostream>

#include "config.h"
#include "charges.h"
#include "nucleus.h"

// Constructor and destructor

enum class WeightMode : int { epsGluon = 0, epsQuark = 1, epsTot = 2};

class Event{
	public:
		Event(Config ConfInput);
		virtual ~Event();

		void NewEvent(int EventID);
		void MakeGrid();
		void get_impact_from_value(double AbsB, double ThetaB);
		void sample_db_impact(double bmin, double bmax);
		void sample_bdb_impact(double bmin, double bmax);

		void MakeEventByEvent();

		double uni_nu_rn(){return drand48();}

		// Interfaces of T1 and T2
		const double& T1p (int64_t nx, int64_t ny) const ;
		double& T1p (int64_t nx, int64_t ny);

		const double& T1n (int64_t nx, int64_t ny) const ;
		double& T1n (int64_t nx, int64_t ny);

		const double& T2p (int64_t nx, int64_t ny) const;
		double& T2p (int64_t nx, int64_t ny);

		const double& T2n (int64_t nx, int64_t ny) const;
		double& T2n (int64_t nx, int64_t ny);

		double T1(int64_t nx, int64_t ny);
		double T2(int64_t nx, int64_t ny);

		// OUTPUT
		void Initialize_output();
		void dump_nucleon_pos(Nucleus * A1,Nucleus * A2);
		void MakeChargeOutput();
		void MakeGlobalQuantities();
		void MakeThicknessOutput();
		void MakeChargeOutput_Transverse(double eta);
		void Make_Event_Output();

		

		////   TOOLS
		double uni_ev_rn(){return drand48();}

		double get_x(double ix);
		double get_y(double iy);
		double get_eta(double ieta);
		int get_ID(int ev,int Nruns){return Nruns*id + ev; }
		int get_nmax(){return n_max;}
		void set_nmax(double n_max_new){n_max=n_max_new;}

		//Spectator Selection

		double ComputeIndividualOverlap_Gaussian(Nucleus *N1,int in1,Nucleus *N2,int in2);
		double ComputeIndividualOverlap_Exponential(Nucleus *N1,int in1,Nucleus *N2,int in2);
		double ComputeInteractionProbability(double CollisionOverlap);
		void CheckParticipants(Nucleus* N1,Nucleus *N2);

		// TOOLS!

		double sigma_inelastic_fm2(double sqrts){return mb_to_fm2*(25.2 + 33.8/pow(sqrts,1.1) + 45.2/pow(sqrts,0.9) + 0.05*log(sqrts) + 0.56*pow(log(sqrts),2));}
		double FifthOrderPoly(double x){return a5 * x + b5 * pow(x,2) + c5 * pow(x,3) + d5 * pow(x,4) + f5 * pow(x,5) + g5 ; }



	private:
		//int Nthreads=1;
		int id=0;
		int n_max=4;
		int ev0;

		NucStruct N1;
		NucStruct N2;

		int NX;
		int NY;
		double cell_trans_volume;

		int EventID;
		int PrimalSeed;
		int EventSeed;

		double sqrtsNN;

		double *b1;
		double *b2;

		double *T1p_ptr;
		double *T1n_ptr;
		double *T2p_ptr;
		double *T2n_ptr;

		std::vector<double> x_cm;
		std::vector<double> y_cm;
		std::vector<double> dEgdeta;
		std::vector<double> dEqdeta;
		std::vector<double> dEdeta;

		std::vector<double> dInt23deta;

		std::vector<double> dnudeta;
		std::vector<double> dnddeta;
		std::vector<double> dnsdeta;
		double x_cm_global,y_cm_global;

		Config config;

		std::string OUTPATH;
		bool is_sigma_output = true;

		Charges * ChargeMaker;

		void printProgress(double percentage);

		//This are parameters for the MC-Glauber Sampling. 
		double a5=0.339446;
		double b5=0.00776864;
		double c5=0.00819359;
		double d5=-0.000767424;
		double f5=0.0000886261;
		double g5=-0.00637744;

		double mb_to_fm2 = 0.1;
	  	double fm2_to_mb = 1/mb_to_fm2;

		double norm_sigma_inelastic_min= 0.5;
		double norm_sigma_inelastic_max= 10.5;

		/// Acceptance Functions//
		bool is_sigma_in_range(double sigma){return (sigma > norm_sigma_inelastic_min && sigma < norm_sigma_inelastic_max ); }
		double compute_internucleon_distance(Nucleus *N1,int ind1,Nucleus *N2,int ind2);
		bool check_interaction_status(Nucleus *N1,int ind1,Nucleus *N2,int ind2);

		// Output
		void print_glauber_data_to_file(Nucleus * N1,Nucleus * N2);
};

#endif /* event_h */ 