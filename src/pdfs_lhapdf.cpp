/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */

#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/Paths.h"

#include "include/pdfs_lhapdf.h"

PDFs::PDFs(){}
PDFs::PDFs(Config * ConfInput){
	PDFName=ConfInput->get_PDFSet();
	Mode=ConfInput->get_ForcedPositive();
	PDF_Verbose=false;

	LHAPDF::setVerbosity(0);
	PDF_ptr=LHAPDF::mkPDF(PDFName);
	PDF_ptr->setForcePositive(Mode);


	set_xMin(PDF_ptr->xMin());
	set_xMax(PDF_ptr->xMax());
	set_QMin(PDF_ptr->qMin());
	set_QMax(PDF_ptr->qMax());

	flavour_vector= PDF_ptr->flavors();
	NF=flavour_vector.size();
	if(PDF_Verbose){
		std::cout<< PDFName << " distributions loaded for"<< std::endl ;
		for (size_t i = 0; i < NF; i++) {
			std::cout<< "  " << get_parton_name(flavour_vector[i]);
		}
		std::cout << std::endl;
	}

}
PDFs::~PDFs(){}

	// --------------------------------------- EVALUATION OF PDFS --------------------------------------- //

double PDFs::xfpart(QuarkID ID,double xx, double QQ){

	if(xx <= xMin || xx >= xMax){return 0;}
	else if(QQ >= QMax ){return 0;}
	else if(QQ <= QMin){return PDF_ptr->xfxQ ((int) ID, xx, QMin);}
	else{return PDF_ptr->xfxQ ((int) ID, xx, QQ);}
}
 