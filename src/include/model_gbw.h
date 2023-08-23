/* Copyright (c) 2023 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef GBW_H
#define GBW_H
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

#include "config.h"
#include "pdfs_lhapdf.h"


struct GluonPars{
//Dynamical
  double sqrts;
  double y;
//GBW
  double Q02;
  double x0;
  double lambda;
  double xcut;
//Geometrical
  double T1;
  double T2;
  //general
  double Sigma0;
};

struct QuarkPars{
//Dynamical
  double sqrts;
  double y;
//GBW
  double Q02;
  double x0;
  double lambda;
  double xcut;
//Geometrical
  double T;
  //general
  double Sigma0;
//Quarks
  PDFs * partons;
  QuarkID quark_id;
  int k;
};

class GBW{
	public:
		GBW();
		GBW(Config ConfInput);
		virtual ~GBW();

		void MakeTable(std::string path_to_set);
		void make_gluon_energy();
		// void make_quark_energy();
		void make_baryon_stopping(int k, QuarkID qid, QuarkID aqid);
    void make_baryon_stopping_all();


    ////// GET FUNCTIONS FROM NAMESPACE
    double get_Q2(double x, double Ts0);
    double get_Q2F(double x, double Ts0);
    double get_Dipole(double k, double Q2);

    bool check_if_zero(double y,double T1, double T2,double &Q12, double &Q22);
    bool check_if_zero_F(double y,double T1, double T2,double &Q12, double &Q22);


	private:
		Config config;
		PDFs * quark_dist;
		std::string path_to_set;

		double dT;

		bool IsPathExist(const std::string &s);
    std::string get_quark_construction_file(int k, QuarkID qid);

		std::string SETPATH;

		std::string gluon_energy_table_name = "gluon_energy_density.dat";
		std::string quark_energy_table_name = "quark_energy_density.dat";
		std::string charge_density_table_name = "charge_density.dat";
		std::string baryon_density_table_name = "baryon_density.dat";

    void printProgress(double percentage);
    int skip=3;

};
#endif /* gbw */
 