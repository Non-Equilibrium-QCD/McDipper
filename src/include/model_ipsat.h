/* Copyright (c) 2023 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef IPSAT_H
#define IPSAT_H
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

#include "config.h"
#include "pdfs_lhapdf.h"
#include "dipole_ipsat.h"


struct GluonParsIPSat{
//Dynamical
  double sqrts;
  double y;
  double qt;
  double kt;

//Geometrical
  double T1;
  double T2;
  Dipole * dip;
};

struct QuarkParsIPSat{
//Dynamical
  double sqrts;
  double y;
//Geometrical
  double T;
//Quarks
  PDFs *partons;
  QuarkID quark_id;
  int k;
  Dipole * dip;
};

class IPSat{
	public:
		IPSat();
		IPSat(Config ConfInput);
		virtual ~IPSat();

		void MakeTable(std::string path_to_set);
    void make_gluon_energy();

    void make_baryon_stopping(int k, QuarkID qid, QuarkID aqid);
    void make_baryon_stopping_all();

    bool check_if_zero(double y,double sqrts,double kmin,double kmax, double xmin, double xmax);
    bool check_if_zero_F(double y,double T);



	private:
		Config config;

		PDFs * quark_dist;

    Dipole * Dip;
		std::string path_to_set;
    int p_set;
    double xscaling;

    double S0;
		double dT;
    double TMax;
    double TMin;

		bool IsPathExist(const std::string &s);
    std::string get_quark_construction_file(int k, QuarkID qid);

		std::string SETPATH;

    std::string gluon_energy_table_name = "gluon_energy_density.dat";
    std::string quark_energy_table_name = "quark_energy_density.dat";
    std::string charge_density_table_name = "charge_density.dat";
    std::string baryon_density_table_name = "baryon_density.dat";


    void printProgress(double percentage);
    void printProgress2(double percentage1,double percentage2);

    int skip=3;

    void TestDump(double T1,double T2);



};
#endif /* ipsat */
 