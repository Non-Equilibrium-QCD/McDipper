/* Copyright (c) 2023 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef CHARGES_H
#define CHARGES_H
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "config.h"
#include "model_gbw.h"
#include "model_ipsat.h"


class Charges{
  /*  This is the Charges class. In case parameter set has not been precomputed, it calls the models.  
  When computed, it interpolates the energy/charge tabulations necessary to produce the events */
   
  public:
    Charges();
    Charges(Config ConfInput);
    virtual ~Charges();

    //Fill-in Functions
    // Grid Related Functions

    void MakeGrid();

    // Model Dependent Functions
    bool import_charges();
		bool call_models();

    bool check_sets();
    bool check_config();

    //  Model Independent Functions
    bool read_in_energy_gluons();
    bool read_in_nK_quark(int k, QuarkID qid);

    // INTERPOLATION FUNCTIONS

    double gluon_energy(double eta, double T1, double T2);
    double quark_energy(double eta, double T1, double T2);
    double u_density(double eta, double T1p,double T1n, double T2p, double T2n);
    double d_density(double eta, double T1p,double T1n, double T2p, double T2n);
    double s_density(double eta, double T1, double T2);

    // Auxiliary functions

    void reset_warnings(){T_warning=true;eta_warning=true;}
    void set_T_warning_off(){T_warning=false;}
    void set_eta_warning_off(){eta_warning=false;}


    // TEST (OUTPUT) FUNCTIONS 
    void dump_charges(double T1,double T2);
    void dump_charges_eta(double eta );
    void set_output_tests();

    /////// OBSERVABLES!!
 

  private:
    Config config;
    Config config_set;

    bool T_warning=true;
    bool eta_warning=true;

    int NX,NY,NETA;
    int NQComp=2*2;

    std::string path_to_tabs;
    std::string path_to_set;

    // Tools
    void check_for_tabs_folders();
    bool IsPathExist(const std::string &s);
    std::string get_set_name(int n);
    std::string get_quark_construction_file(int k, QuarkID qid);

    void get_name_and_value(std::string testline,std::string &name, std::string &value);
    double double_tolerance=0.001;

    std::string gluon_energy_table_name = "gluon_energy_density.dat";

    // Interpolators
    gsl_spline2d ** e_g_spl ;

    gsl_spline2d ** N0u_spl ;  //2*2
    gsl_spline2d ** N0d_spl ;
    gsl_spline2d ** N0s_spl ;

    gsl_spline2d ** N1u_spl ;  //2*2
    gsl_spline2d ** N1d_spl ;
    gsl_spline2d ** N1s_spl ;

    gsl_interp_accel ** xaccEG ;
    gsl_interp_accel ** yaccEG ;

    gsl_interp_accel ** xaccN0u ;
    gsl_interp_accel ** yaccN0u ;

    gsl_interp_accel ** xaccN0d ;
    gsl_interp_accel ** yaccN0d ;

    gsl_interp_accel ** xaccN0s ;
    gsl_interp_accel ** yaccN0s ;

    gsl_interp_accel ** xaccN1u ;
    gsl_interp_accel ** yaccN1u ;

    gsl_interp_accel ** xaccN1d ;
    gsl_interp_accel ** yaccN1d ;

    gsl_interp_accel ** xaccN1s ;
    gsl_interp_accel ** yaccN1s ;

    std::vector<double> T_vector;


    bool initialized=false;

};

#endif  /* charges */
