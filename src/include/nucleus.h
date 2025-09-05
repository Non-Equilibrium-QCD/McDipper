/* Copyright (c) 2023 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef NUCLEUS_H
#define NUCLEUS_H
#include <iostream>
#include <random>
#include "random.h"

#include "config.h"
enum class Nucleon : int { proton = 1, neutron=0};

class Nucleus{
	/* This is the nuleus class. It creates a nucleus as an object. It contains knowledge of a single realization of the nucleus, 
	according to Wood-Saxon's sampling, and it can compute the thickness for such realization. The nucleus class does not contain 
	any memory allocated for a "Thickness Grid", that is in the event class. */
	public:
		Nucleus(NucStruct NucIn);
		virtual ~Nucleus();
		// 

		double nuclear_density(double x,double y,double z);
		double nuclear_thickness_optical(double x,double y);
		double random_position();
		void sample_single_position(double * x_t);
        void sample_hotspots(double * x_t);
		void set_nucleon_positions();
		void rotate_nucleus();
        void Thickness_fluct();
		void refresh_positions();
		void shift_nucleus_by_impact(double bx,double by);

		//Nuclear structure functions 	
		void import_nuclear_configurations();
		double get_weight(){return weight;}
		int ConfIndex(int ie,int n,int ix);
		//Retrievers
		const double& Configuration(int64_t ie, int64_t n, int64_t ix) const;
		double& Configuration (int64_t ie, int64_t n, int64_t ix);
	
		// Thickness functions
		double NucleonThickness(double x,double y,double x0,double y0,int n,double BG);
		double GetThickness(double xt,double yt,double BG);
		double GetThickness_p(double xt,double yt,double BG);
		double GetThickness_n(double xt,double yt,double BG);

		// Nucleon positions
		void get_position_nucleon(int n, double &x, double &y,  double &z);
		void get_trans_position_nucleon(int n, double &x, double &y);
		void get_trans_position_hotspot(int n, double &x, double &y, int hotspoti);

		// Tools
		double uni_nu_rn(){return drand48();}
		int uni_nu_int(){return lrand48();}

        std::gamma_distribution<double> gamma_dist;
        std::lognormal_distribution<double> lognorm_dist;
        std::normal_distribution<double> hotspots_posi_dist;


		int get_A(){return A;}
		int get_Z(){return Z;}

		// Retrieving functions 
		int get_NColl(){return NumberOfCollisions;}
		int get_NPart(){return NumberOfParticipants;}
		int get_ParticipantStatus(int n){return ParticipantStatus[n];}
		int get_CollisionNumber(int n){return CollisionNumber[n];}
		int get_nucleon_type(int n){
			if(NucleonType[n]==Nucleon::neutron){return 0;}
			else if(NucleonType[n]==Nucleon::proton){return 1;}
			else{return 0;}
		}

		void set_ParticipantStatus(int n, int ind){ParticipantStatus[n]=ind;}
		void set_CollisionNumber(int n, int Ncoll_new){CollisionNumber[n]=Ncoll_new;}

		void set_NPart(int NPart_new){NumberOfParticipants=NPart_new;}
		void set_NColl(int NColl_new){NumberOfCollisions=NColl_new;}
		void update_NColl(int NColl_update){NumberOfCollisions+= NColl_update;}
		void add_Participant(){NumberOfParticipants++;}

		int get_number_of_participants(){return NumberOfParticipants;}
		int get_Nq(){return Nq;}

	private:
		int A;  // Atomic Number
		int Z; // Charge Number
 		int mode; //mode 0-> Spherical, 1-> deformed
		std::string modeStr;
		std::string InputName;
		bool IsIsospinSpecified;
		// bool IsSpinSpecified;
        bool is_hotspots_fluct=false;
        bool is_thick_fluct=false;

		double * Configurations_ptr;
		double * Config_weights;
		int NConf;
		bool is_weights=false;

        int Nq=0;
        double Bq=0.0;
        double Br=0.0;
        
        double sigma=0.0;
        std::string fluct_mode="";

		double * NucPars;

		double ** r;// r is a vector which contains each nucleon's information (x,y,z, hotspot_1x, hotspot_1y, hotspot_2x, hotspot_2y, ...)
		double * rBar;

        double ** w;

		Nucleon * NucleonType;
		int * ParticipantStatus;
		int * CollisionNumber;

    	int NumberOfCollisions;
    	int NumberOfParticipants;

		double Y02pref;
		double Y04pref;

		double RSampFactor;

		void rotate_X_axis(double *r, double theta);
		void rotate_Y_axis(double *r, double theta);
		void rotate_Z_axis(double *r, double theta);

		void shuffle(int *array, size_t n);

		double weight= 1.0; // Statistical weight of the nuclear configuration, DEF=1, but can change due to, i.e. NLEFT, configs.
	};

#endif /* nucleus */
 
