/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */

#include <iostream>
#include <math.h>

#include "include/nuclear/H2.cpp"
#include "include/nuclear/He3.cpp"

#include "include/nuclear_data.h"
#include "include/nucleus.h"

// Constructor and destructor

Nucleus::Nucleus(int AtomicNumber,int ChargeNumber,int NucleusType): A(AtomicNumber),Z(ChargeNumber),mode(NucleusType){
		// Chose parameter structure according to
		if(A==1){modeStr="Proton";}
		else if(A==2 && Z==1){modeStr="Deuteron";}
		else if(A==3 && Z==2){modeStr="He3";}
		else if(A>3){
			switch (mode) {
				case 0:
					NucPars=new double[2];
					modeStr="Spherical";break;
				case 1:
					NucPars=new double[4];
					modeStr="Deformed";break;
				default:
					std::cerr<<"Error: Nucleus type not implemented. Exiting.";exit(EXIT_FAILURE);
			}
			Y02pref = sqrt(5/(16*M_PI));
			Y04pref = (3/16)/sqrt(M_PI);
			RSampFactor=4.0;
	    // set parameters //
	     NuclearData::getNPars(A,Z,mode,NucPars);
		 }
		 else{std::cerr<<"Error: Nucleus type not implemented. Exiting.";exit(EXIT_FAILURE); }


		r=new double*[A];
		rBar=new double[3];
		ParticipantStatus=new int[A];
		CollisionNumber=new int[A];
		NucleonType = new Nucleon[A];

		for(int n=0;n<A;n++){r[n]=new double[3];}

    // SET NUCLEON POSITIONS //
    set_nucleon_positions();
	// ROTATE RANDOMLY (EULER ANGLES)	
	rotate_nucleus();
}

Nucleus::~Nucleus(){
  for(int n=0;n<A;n++){
      delete[] r[n];
  }
  delete[] r;
	delete[] rBar;
	delete[] NucleonType;
	delete[] ParticipantStatus;
	delete[] CollisionNumber;

  // delete participants;
  // delete collisionNum;
}

// Functions

double Nucleus::nuclear_density(double x,double y, double z){
	double r_t = sqrt(x*x + y*y +z*z );
	double R_t;
	if(mode==0){R_t=NucPars[0];}
	if(mode==1){
		double Cos_th = z/r_t;
		double Y02 = Y02pref *(3*pow(Cos_th,2) -1);
		double Y04 = Y04pref *(35*pow(Cos_th,4)-30*pow(Cos_th,2) + 3);
		R_t=NucPars[0]*(1 + NucPars[2]*Y02 + NucPars[3]*Y04 ) ;
	}
	return (1+exp(-R_t/NucPars[1]))/(1+exp( (r_t - R_t)/NucPars[1]));

}

double Nucleus::nuclear_thickness_optical(double x,double y){

	if(mode!=0){std::cerr<<" [ Error ]: Nucleus type not implemented for optical thickness! Exiting.";exit(EXIT_FAILURE);}
	else{
		double rho0,CC;
		if(A==208){rho0=0.000771107;}
		if(A==197){rho0=0.000859619;}

		double dz= 0.01*NucPars[0];
		double ZMaxAbs=3*NucPars[0];
		int NZ= int(2*ZMaxAbs/dz);

		double thickness_unnorm=0;
		for (size_t iz = 0; iz < NZ; iz++) {
			double zz = dz*iz - ZMaxAbs;
			if(iz==0 or iz==NZ-1){CC=1/2.;}
			else{CC=1.;}
			thickness_unnorm += dz*CC*nuclear_density(x, y, zz);
		}

		return rho0 * thickness_unnorm;
	}
	return 0;

}


double Nucleus::random_position(){
	return RSampFactor*NucPars[0]*(uni_nu_rn()-0.5);
}


void Nucleus::sample_single_position(double * x_t){
	int Accept=0;
  while(Accept==0){
	  // Sample positions //
	  x_t[0]=random_position();
	  x_t[1]=random_position();
	  x_t[2]=random_position();
    // get rejection probability //
    double Probability=nuclear_density(x_t[0],x_t[1],x_t[2]);
    if(uni_nu_rn()<Probability){Accept=1;} //Check rejection criterion
	}
}

void Nucleus::set_nucleon_positions(){
	if(A==1){r[0][0]=0.0; r[0][1]=0.0; r[0][2]=0.0;NucleonType[0]=Nucleon::proton;}
	else if(A==2){
		H2::Init();H2::GetNucleonPositions(r[0],r[1]);
		NucleonType[0]=Nucleon::proton;
		NucleonType[1]=Nucleon::neutron;
	}
	else if(A==3){
		He3::Init();He3::GetNucleonPositions(r[0],r[1],r[2]);
		NucleonType[0]=Nucleon::proton;
		NucleonType[1]=Nucleon::proton;
		NucleonType[2]=Nucleon::neutron;
	}
	else{
		for(int n=0;n<A;n++){
			sample_single_position(r[n]);
			if(n<=Z){NucleonType[n]=Nucleon::proton;}
			else{NucleonType[n]=Nucleon::neutron;}
		}
	}  ///Sample Glauber Positions

	// Locate Center of Mass
	rBar[0]=0.0; rBar[1]=0.0; rBar[2]=0.0;
	for(int n=0;n<A;n++){
		for(int i=0;i<3;i++){
			rBar[i]+=r[n][i]/double(A);
		}
	}
	// Shift by the center of mass
	for(int n=0;n<A;n++){
		for(int i=0;i<3;i++){
			r[n][i]-=rBar[i];
		}
	}

	rBar[0]=0.0; rBar[1]=0.0; rBar[2]=0.0; //Clean the Center of Mass vector
}

void Nucleus::rotate_nucleus(){
	double thetaX=2*M_PI*uni_nu_rn();
	double thetaY=M_PI*(uni_nu_rn()-0.5);
	double thetaZ=2*M_PI*uni_nu_rn();
	for(int n=0;n<A;n++){
		rotate_X_axis(r[n],thetaX);
		rotate_Y_axis(r[n],thetaY);
		rotate_Z_axis(r[n],thetaZ);
	}
}

double Nucleus::NucleonThickness(double x,double y,double x0,double y0,double BG){
    double r2 = pow(x-x0,2) + pow(y-y0,2);
    return exp(-0.5*r2/BG)/(2.0*M_PI*BG);
}

double Nucleus::GetThickness(double xt,double yt,double BG){
    double TValue=0.0;
		// SUM ALL NUCLEONS //
    for(int n=0;n<A;n++){
			if (ParticipantStatus[n]==1){TValue+=NucleonThickness(xt,yt,r[n][0],r[n][1],BG);}
		}
    return TValue;
}

double Nucleus::GetThickness_p(double xt,double yt, double BG){
    double TValue=0.0;
		// SUM ALL NUCLEONS //
    for(int n=0;n<A;n++){
			if (ParticipantStatus[n]==1 && NucleonType[n]==Nucleon::proton){TValue+=NucleonThickness(xt,yt,r[n][0],r[n][1],BG);}
		}
    return TValue;
}

double Nucleus::GetThickness_n(double xt,double yt,double BG){
    double TValue=0.0;
		// SUM ALL NUCLEONS //
    for(int n=0;n<A;n++){
			if (ParticipantStatus[n]==1 && NucleonType[n]==Nucleon::neutron){TValue+=NucleonThickness(xt,yt,r[n][0],r[n][1],BG);}
		}
    return TValue;
}



void Nucleus::get_position_nucleon(int n, double &x, double &y,  double &z){x=r[n][0];y=r[n][1];z=r[n][2];}
void Nucleus::get_trans_position_nucleon(int n, double &x, double &y){x=r[n][0];y=r[n][1];}


void Nucleus::rotate_X_axis(double *r, double theta){
	double rp[3]={r[0],r[1],r[2]};
	r[0]=rp[0];
	r[1]=rp[1]*cos(theta)-rp[2]*sin(theta);
	r[2]=rp[1]*sin(theta)+rp[2]*cos(theta);
}

void Nucleus::rotate_Y_axis(double *r, double theta){
	double rp[3]={r[0],r[1],r[2]};
	r[0]=rp[0]*cos(theta)+rp[2]*sin(theta);
	r[1]=rp[1];
	r[2]=rp[2]*cos(theta)-rp[0]*sin(theta);
}

void Nucleus::rotate_Z_axis(double *r, double theta){
	double rp[3]={r[0],r[1],r[2]};
	r[0]=rp[0]*cos(theta)-rp[1]*sin(theta);
	r[1]=rp[0]*sin(theta)+rp[1]*cos(theta);
	r[2]=rp[2];
};

void Nucleus::shift_nucleus_by_impact(double bx,double by){
	for(int n=0;n<A;n++){
		r[n][0] -= bx;
		r[n][1] -= by;
	}
}


  