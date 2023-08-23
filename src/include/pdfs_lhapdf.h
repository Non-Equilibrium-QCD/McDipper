/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef PDFS_H
#define PDFS_H
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/Paths.h"
#include "config.h"

// #include "params_pdfs_lhapdf.h"


class PDFs{
	public:
	PDFs();
	// PDFs(int iparton, bool Verbose);
	PDFs(Config * ConfInput);
	virtual ~PDFs();
	// --------------------------------------- GET LIMITS--------------------------------------- //
	double get_xMin(){return xMin;}
	double get_xMax(){return xMax;}
	double get_QMin(){return QMin;}
	double get_QMax(){return QMax;}

	// --------------------------------------- SET LIMITS--------------------------------------- //
	void set_xMin(double xMin_t){xMin=xMin_t;}
	void set_xMax(double xMax_t){xMax=xMax_t;}
	void set_QMin(double QMin_t){QMin=QMin_t;}
	void set_QMax(double QMax_t){QMax=QMax_t;}

	// --------------------------------------- EVALUATION OF PDFS --------------------------------------- //
	double xfpart(QuarkID ID,double xx, double QQ);

	// --------------------------------------- OTHER INFOS --------------------------------------- //

	int get_NF(){return NF;}
	private:
	bool PDF_Verbose;
	LHAPDF::PDF* PDF_ptr;
	std::string PDFName;
	int Mode;
	double QMin,QMax;
	double xMin,xMax;
	std::vector< int > flavour_vector;
	int NF;

	std::string get_parton_name(int qid){
		std::ostringstream parton_name;
		if(qid<0){parton_name << "anti-";}
		switch (abs(qid)) {
			case 1:
				parton_name << "down";
				break;
			case 2:
				parton_name << "up";
				break;
			case 3:
				parton_name << "strange";
				break;
			case 4:
				parton_name << "charm";
				break;
			case 5:
				parton_name << "bottom";
				break;
			case 6:
				parton_name << "top";
				break;
			case 0:
				parton_name << "gluon";
				break;
			case 21:
				parton_name << "gluon";
				break;
			case 9:
				parton_name << "gluon";
				break;
			default:
				std::cout << "Parton identifier, "<< qid <<" is not valid\n";
				exit(1);
				break;
		}
		return parton_name.str();
	}

};

#endif /* pdfs_h */
 