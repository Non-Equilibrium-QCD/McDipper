/* Copyright (c) 2023 Oscar Garcia-Montero
 * All rights reserved. */
#ifndef EVENT_H
#define EVENT_H
#include <iostream>

#include "config.h"
#include "charges.h"
#include "nucleus.h"
#include "h5_out.h"

// Constructor and destructor

enum class WeightMode : int { epsGluon = 0, epsQuark = 1, epsTot = 2};

class Event{
	public:
		Event(Config ConfInput);
		virtual ~Event();

		void NewEvent(int EventID, Nucleus * A1, Nucleus *A2);
		void MakeGrid();
		void get_impact_from_value(double AbsB, double ThetaB);
		void sample_db_impact(double bmin, double bmax);
		void sample_bdb_impact(double bmin, double bmax);

		void MakeEventByEvent();
		void InitializeAverageEvent();

		
		

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

		//AVG EVENT 

		const double& EgAvg (int64_t neta, int64_t nx, int64_t ny) const;
		double& EgAvg (int64_t neta, int64_t nx, int64_t ny);

		const double& EqAvg (int64_t neta, int64_t nx, int64_t ny) const;
		double& EqAvg (int64_t neta, int64_t nx, int64_t ny);

		const double& nuAvg (int64_t neta, int64_t nx, int64_t ny) const;
		double& nuAvg (int64_t neta, int64_t nx, int64_t ny);

		const double& ndAvg (int64_t neta, int64_t nx, int64_t ny) const;
		double& ndAvg (int64_t neta, int64_t nx, int64_t ny);

		const double& nsAvg (int64_t neta, int64_t nx, int64_t ny) const;
		double& nsAvg (int64_t neta, int64_t nx, int64_t ny);


		// OUTPUT
		void Initialize_output();
		void print_nucleon_pos(Nucleus * A1,Nucleus * A2, int EventID);
		void print_weights(Nucleus * N1,Nucleus * N2, int EventID);
		void MakeChargeOutput();
		void MakeGlobalQuantities();
		void MakeGlobalQuantities_AverageEvent();
		void MakeChargeOutput_AverageEvent();
		void MakeThicknessOutput();
		void MakeChargeOutput_Transverse(double eta);
		void Make_Event_Output();
		void MakeBlock(int EventID);
		void RemoveBlock();

		void MakeGlobalQuantities_Midrapidity();
		void MakeChargeOutputMidrapidity();

		void MakeOutputFiles();

		
		

		////   TOOLS
		double uni_ev_rn();
		int uni_ev_int();

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
		double solve_parton_sigma();
		void CheckParticipants(Nucleus* N1,Nucleus *N2);

		// TOOLS!

		double sigma_inelastic_fm2(double sqrts){return mb_to_fm2*(25.2 + 33.8/pow(sqrts,1.1) + 45.2/pow(sqrts,0.9) + 0.05*log(sqrts) + 0.56*pow(log(sqrts),2));}
		//double sigma_inelastic_fm2(double sqrts){return 3.1253+0.1280 * pow(log(sqrts)-2.0412, 1.8231);}

		// This function is the inverse to the dimensionless solution to the matching (sigmaInel/4piBG = F[sig_g/(4pi BG)])
		double FitInverseInelXSec(double x){return  0.561459*exp(x) + 0.00615644*x - 0.0533597 ; } // Gives F^{-1}(sigmaInel/4piBG)
    
		// Replaced by Event::solve_parton_sigma(), which does the matching at run
		// time for any Nq, BG and Bq. Kept for reference: each line is a fit to a
		// set of TRENTo runs, Nq=3 only, at the one (w,v) its comment names.
		//double FitInverse_parton_sigma(double x)
		//{return exp(0.000307231898*pow(x,5)-0.0094106092*pow(x,4)+0.110155447*pow(x,3)-0.556793005*pow(x,2)+1.95110038*x-1.30207);}//w=0.4fm,v=0.2fm
        	//{return exp(0.000154975735*pow(x,5)-0.0050539279*pow(x,4)+0.063947562*pow(x,3)-0.380394993*pow(x,2)+1.60700817*x-1.12212);}//w=0.4fm,v=0.3fm
		//{return exp(0.000668144874*pow(x,5)-0.0179029924*pow(x,4)+0.178578559*pow(x,3)-0.621370047*pow(x,2)+2.05178913*x-1.24023);}//w=0.4fm,v=0.116fm
		//{return exp(0.000114104509*pow(x,5)-0.0037398475*pow(x,4)+0.049565563*pow(x,3)-0.333671377*pow(x,2)+1.45320013*x-1.06351);}//w=0.6fm,v=0.3fm
		//{return exp(0.000111855876*pow(x,5)-0.0037261861*pow(x,4)+0.049413354*pow(x,3)-0.336687571*pow(x,2)+1.42577757*x-1.06848);}//w=0.8fm,v=0.4fm
		//{return exp(0.000303478666*pow(x,5)-0.0092639543*pow(x,4)+0.107327778*pow(x,3)-0.529384823*pow(x,2)+1.89562230*x-1.25462);}//w=0.3795,v=0.2
		//{return exp(0.000121814035*pow(x,5)-0.0041679746*pow(x,4)+0.056055270*pow(x,3)-0.357869035*pow(x,2)+1.55038338*x-1.10125);}//w=0.4405,v=0.3
		//{return exp(0.000110121625*pow(x,5)-0.0036891554*pow(x,4)+0.049334386*pow(x,3)-0.331558266*pow(x,2)+1.45530546*x-1.06012);}//w=0.5138,v=0.4
		//{return exp(0.000119036963*pow(x,5)-0.0039188703*pow(x,4)+0.051408643*pow(x,3)-0.344204653*pow(x,2)+1.45851547*x-1.07931);}//w=0.595,v=0.5
		//{return exp(0.000117366864*pow(x,5)-0.0038643018*pow(x,4)+0.050640694*pow(x,3)-0.340879880*pow(x,2)+1.43311119*x-1.06750);}//w=0.6812,v=0.6
	private:
		//int Nthreads=1;
		int id=0;
		int n_max=4;
		int ev0;

		NucStruct N1;
		NucStruct N2;

		int NX;
		int NY;
		int NETA;
		double cell_trans_volume;

		int EventID;
		int PrimalSeed;

		double sqrtsNN;

		double *b1;
		double *b2;

		double *T1p_ptr;
		double *T1n_ptr;
		double *T2p_ptr;
		double *T2n_ptr;

		double * EgAvg_ptr;
		double * EqAvg_ptr;
		double * nuAvg_ptr;
		double * ndAvg_ptr;
		double * nsAvg_ptr;

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
		double sigma_parton = 0.0;	// solved once, in the constructor

		Charges * ChargeMaker;

		// One file per event, opened before the first output function runs and
		// closed after the last. Null unless Output:H5 is set.
		H5Event * h5out;

		void printProgress(double percentage);



		double mb_to_fm2 = 0.1;
	  	double fm2_to_mb = 1/mb_to_fm2;

		double norm_sigma_inelastic_min= 1.0;
		double norm_sigma_inelastic_max= 20.0;

		/// Acceptance Functions//
		bool is_sigma_in_range(double sigma){return (sigma > norm_sigma_inelastic_min && sigma < norm_sigma_inelastic_max ); }
		double compute_internucleon_distance(Nucleus *N1,int ind1,Nucleus *N2,int ind2);
		bool check_interaction_status(Nucleus *N1,int ind1,Nucleus *N2,int ind2);

		// Output
		void print_glauber_data_to_file(Nucleus * N1,Nucleus * N2);
};

#endif /* event_h */ 
