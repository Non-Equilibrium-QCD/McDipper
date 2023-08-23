/* Copyright (c) 2023 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef NUCLEUS_H
#define NUCLEUS_H
#include <iostream>

enum class Nucleon : int { proton = 0, neutron=1};

class Nucleus{
	/* This is the nuleus class. It creates a nucleus as an object. It contains knowledge of a single realization of the nucleus, 
	according to Wood-Saxon's sampling, and it can compute the thickness for such realization. The nucleus class does not contain 
	any memory allocated for a "Thickness Grid", that is in the event class. */
	public:
		Nucleus(int AtomicNumber,int ChargeNumber,int NucleusType);
		virtual ~Nucleus();
		// 

		double nuclear_density(double x,double y,double z);
		double nuclear_thickness_optical(double x,double y);
		double random_position();
		void sample_single_position(double * x_t);
		void set_nucleon_positions();
		void rotate_nucleus();
		void shift_nucleus_by_impact(double bx,double by);

		// Thickness functions
		double NucleonThickness(double x,double y,double x0,double y0,double BG);
		double GetThickness(double xt,double yt,double BG);
		double GetThickness_p(double xt,double yt,double BG);
		double GetThickness_n(double xt,double yt,double BG);

		// Nucleon positions
		void get_position_nucleon(int n, double &x, double &y,  double &z);
		void get_trans_position_nucleon(int n, double &x, double &y);

		// Tools
		double uni_nu_rn(){return drand48();}

		int get_A(){return A;}
		int get_Z(){return Z;}

		// Retrieving functions 
		int get_NColl(){return NumberOfCollisions;}
		int get_NPart(){return NumberOfParticipants;}
		int get_ParticipantStatus(int n){return ParticipantStatus[n];}
		int get_CollisionNumber(int n){return CollisionNumber[n];}

		void set_ParticipantStatus(int n, int ind){ParticipantStatus[n]=ind;}
		void set_CollisionNumber(int n, int Ncoll_new){CollisionNumber[n]=Ncoll_new;}

		void set_NPart(int NPart_new){NumberOfParticipants=NPart_new;}
		void set_NColl(int NColl_new){NumberOfCollisions=NColl_new;}
		void update_NColl(int NColl_update){NumberOfCollisions+= NColl_update;}
		void add_Participant(){NumberOfParticipants++;}

		int get_number_of_participants(){return NumberOfParticipants;}

	private:
		int A;  // Atomic Number
		int Z; // Charge Number
 		int mode; //mode 0-> Spherical, 1-> deformed
		std::string modeStr;


		double * NucPars;

		double ** r;
		double * rBar;

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

};

#endif /* nucleus */
 