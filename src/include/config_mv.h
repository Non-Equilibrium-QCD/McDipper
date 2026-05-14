/* Copyright (c) 2025 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef CONFIG_MV_H
#define CONFIG_MV_H

#include <iostream>
#include <fstream>
#include <string>

#include "params_gen.h"

class MVConfig{
	/* This is the MV-model config reader class. It reads in input from the */
	public:
		MVConfig();
		MVConfig(std::string PATH);
		virtual ~MVConfig();
		bool operator==(const MVConfig& other) const { return compare(other); }
		bool operator!=(const MVConfig& other) const { return !compare(other); }

		 /*  **** General functions ****   */ 
		void process_parameters();
		// General
		void process_general_parameters(std::string testline);
		//Grid
		void process_grid_parameters(std::string testline);
		void make_momentum_parameters(); 
		// Model Parameters
		void process_MV_parameters(std::string testline);
		// Evolution
		void process_BK_parameters(std::string testline);
		// processing tools 
		int count_indent(std::string testline);
		std::string get_header(std::string testline);
		//Comparison

		bool compare(const MVConfig& other) const;
		//Output 
		void terminal_setup_output();
		void write_config(std::string path_to_write);	

		// ------- Getters -------
		// General
		int    getNC()        const { return NC; }
		int    getNF()        const { return NF; }
		double getLambdaQCD() const { return LambdaQCD; }

		// MV
		double getQs02()      const { return Qs02; }
		double getGamma()     const { return gamma; }
		double getLambdaMV()  const { return LambdaMV; }
		double getSigma0()    const { return sigma0; }
		double getEc()        const { return ec; }

		double getQs02_in_fm()      const { return Qs02*gen_pars::GeV2_to_fmm2; }
		double getLambdaMV_in_fm()  const { return LambdaMV*gen_pars::GeV_to_fmm1; }


		// BK
		double getLambdaBK()  const { return LambdaBK; }
		double getMu2()       const { return mu2; }
		double getC2()        const { return C2; }
		double getZeta()      const { return zeta; }
		double getMp()        const { return mp; }
		double getNGC()       const { return nGC; }
		const std::string& getKtype() const { return Ktype; }

		// Grid sizes
		int    getNr() const { return Nr; }
		int    getNT() const { return NT; }
		int    getNY() const { return NY; }
		int    getNK() const { return NK; }

		// r-grid / q-grid
		double getRmin() const { return rmin; }
		double getRmax() const { return rmax; }
		double getQmin() const { return qmin; }
		double getQmax() const { return qmax; }
		double getDq()   const { return dq; }

		// Y-grid
		double getX0()   const { return x0; }
		double getYmin() const { return Ymin; }
		double getYmax() const { return Ymax; }
		double getdY()   const { return dY; }
		


		// T-grid
		double getTmin() const { return Tmin; }
		double getTmax() const { return Tmax; }
		double getdT()   const { return dT; }

		// k/u grid
		double getKmax() const { return kmax; }
		double getKmin() const { return kmin; }
		double getUmax() const { return umax; }
		double getUmin() const { return umin; }
		double getDu()   const { return du; }


		//transformesrs
		double getq(double r){return log(r/rmin);}
		double getu(double k){return log(k/kmin);}
		double getY(double x){return log(x0/x);}
		double getx(double Y){return x0*exp(-Y);}

		double q_at(int i){return qmin+i*dq; }
		double u_at(int i){return umin+i*du; }
		double T_at(int i){return Tmin+i*dT; }
		double Y_at(int i){return Ymin+i*dY; }
		double r_at(int i){return rmin*exp(qmin+i*dq);}
		double k_at(int i){return kmin*exp(umin+i*du);}
		double x_at(int i){return x0*exp(-(Ymin+i*dY));}
	
		
		double getXmin() const { return xmin;}

	private:
		std::string subheader;
		std::string path;
		//Config Variables 
		// General
		int NC, NF;
		double LambdaQCD; 
		// MV parameters
		double Qs02,gamma,LambdaMV,sigma0,ec;
		// BK-Parameters
		double LambdaBK,mu2,C2,zeta,mp,nGC;
		std::string Ktype; 
		// Grid
		int Nr, NT, NY, NK;
		double rmin,rmax,qmin,qmax,dq;
		double x0, Ymin, Ymax, dY, xmin;
		double Tmin, Tmax, dT;
		double kmax,kmin, umax,umin,du;

		void get_name_and_value(std::string testline,std::string &name, std::string &value);
		// void make_momentum_parameters(bool in_run);
		// void quadratic_extrapolation(double x1,double f1,double x2,double f2,double x3,double f3,double &a,double &b,double &c);

		//checkers
		bool check_dipole_sets();
		// bool check_config();
		// void read_in_config();
		// bool is_in_range(double x, double x1,double x2);

		bool reldiff(double a, double b, double relEpsilon = 1e-9)const {
			double diff = std::fabs(a - b);
			double largest = std::max(std::fabs(a), std::fabs(b));
			return diff <= largest * relEpsilon;
		}

		

 };

#endif /* config_mv */
