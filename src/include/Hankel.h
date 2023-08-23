/* Copyright (c) 2023 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef HANKEL_H
#define HANKEL_H
#include <iostream>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_dht.h>
#include <vector>

namespace HankelParameters{
  const double MinFactor=0.05;
  const double MaxFactor=10;
}

class Hankel{
public:
  Hankel(){
    // Initializes Default
    tf = gsl_dht_new(NPoints,TransformOrder,MaxRadius);
    create_arrays();
    output_creation_message();
  };

  Hankel(int TransformOrder_new){
    // Initializes custom order with default points + Domain
    TransformOrder= TransformOrder_new;
    tf = gsl_dht_new(NPoints,TransformOrder,MaxRadius);
    create_arrays();
    output_creation_message();
  };

  Hankel(int TransformOrder_new, bool Verbose){
    // Initializes custom order with default points + Domain
    HankelVerbose=Verbose;
    TransformOrder= TransformOrder_new;
    tf = gsl_dht_new(NPoints,TransformOrder,MaxRadius);
    create_arrays();
    output_creation_message();
  };

  Hankel(double NPoints_new,double MaxRadius_new, int TransformOrder_new, bool Verbose){
    // Rehashes the default parameters
    HankelVerbose=Verbose;
    TransformOrder= TransformOrder_new;
    NPoints= NPoints_new;
    MaxRadius =MaxRadius_new;
    tf = gsl_dht_new( NPoints,TransformOrder,MaxRadius );
    create_arrays();
    output_creation_message();
  }

  virtual ~Hankel(){
    gsl_dht_free(tf);}

	//Initialization Functions

	void ReInit(double NPoints_new,double MaxRadius_new){
    // Order Can't change!
    NPoints= NPoints_new;
    MaxRadius =MaxRadius_new;
    create_arrays();
    gsl_dht_free(tf);
    tf = gsl_dht_new( NPoints,TransformOrder,MaxRadius );
  }

	void Initialize_Domain_w_CharScale(double CharacteristicDomain){
    // Order does not change!
    CharScale = CharacteristicDomain;
    MaxRadius = MaxFactor*CharacteristicDomain;
    MinRadius = MinFactor*CharacteristicDomain;

    // TransformOrder= TransformOrder_new;
    int NPoints_new= ceil( (2 *TransformOrder*(DomainRatio-1)+3*DomainRatio + 1)/4. );
    // if(NPoints!=NPoints_new){create_arrays();}
    NPoints=NPoints_new;
    create_arrays();
    tf = gsl_dht_new( NPoints, TransformOrder, MaxRadius );
  }


  double get_max_radius(){return MaxRadius;}
  double get_min_radius(){return MinRadius;}
  double get_char_scale(){return CharScale;}
  int get_Npoints(){return NPoints;}

  void output_state(){
    std::cout<< "Hankel Object of order:  "<<TransformOrder << std::endl;
    std::cout<< "Minimum Radius:  "<< MinRadius << std::endl;
    std::cout<< "Maximum Radius:  "<< MaxRadius << std::endl;
    std::cout<< "Number of Points:  "<< NPoints << std::endl;
    std::cout<< std::endl;
  }

	//Sample
  double get_X(int i){return gsl_dht_x_sample(tf,i);} // Returns x_n = ( j_{\nu,n+1} / j_{\nu,M}) X.
  double get_K(int i){return gsl_dht_k_sample(tf,i);} // Returns k_n = ( j_{\nu,n+1} / X )

  void set_FX(int i,double fx){Fx.at(i)=fx;} // Returns Fx_n
  void set_FK(int i,double fk){Fk.at(i)=fk;}

  double get_FX(int i){return Fx.at(i);} // Returns Fx_n
  double get_FK(int i){return Fk.at(i);} // Returns Fk_n


	void get_min_Zero(){
		if(TransformOrder==0){MinZero=gsl_sf_bessel_zero_J0(1);}
		else if(TransformOrder==1){MinZero=gsl_sf_bessel_zero_J1(1);}
		else{MinZero=gsl_sf_bessel_zero_Jnu(TransformOrder, 1);}
	}

	// Transform
  void Transform(){
    double Fx_t[NPoints];
    double Fk_t[NPoints];
    for (size_t i = 0; i < NPoints; i++) {Fx_t[i]=Fx.at(i);}
		gsl_dht_apply(tf,Fx_t,Fk_t);
    for (size_t i = 0; i < NPoints; i++) {Fk.at(i)=Fk_t[i];}
  }

  //Auxiliary
  double get_min_factor(){return MinFactor;}
  double get_max_factor(){return MaxFactor;}
  // double get_min_factor(){return MinFactor;}


private:
  bool HankelVerbose=false;
  gsl_dht * tf;
  int TransformOrder=0;
  int NPoints= 50;
  double MaxRadius = 50;
  double MinRadius = 0.01;
  double CharScale;

	double MinFactor=HankelParameters::MinFactor;
  double MaxFactor=HankelParameters::MaxFactor;
  double DomainRatio=MaxFactor/MinFactor;

	int MinZero;

  std::vector<double> Fx;
  std::vector<double> Fk;

  void output_creation_message(){
    if(HankelVerbose){
      std::cout<< "Created Hankel Object for nu = " << TransformOrder <<", Rmax = " << MaxRadius <<" and NPoints = " << NPoints << std::endl;
    }
  }

  void create_arrays(){
    Fx.resize(NPoints);
    Fk.resize(NPoints);
  }



};

#endif /* Hankel */
 