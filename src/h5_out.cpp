/* Copyright (c) 2026 */

#include <iostream>
#include <cstdlib>
#include "include/h5_out.h"

H5Event::H5Event(){}

H5Event::H5Event(std::string path){
	file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if(file<0){
		std::cerr<<"[ Error::H5 ]: Cannot create "<<path<<std::endl;
		exit(EXIT_FAILURE);
	}
	charges = -1;
	charges_space = -1;
	charges_row = 0;
}

H5Event::~H5Event(){
	if(charges_space>=0){H5Sclose(charges_space);}
	if(charges>=0){H5Dclose(charges);}
	H5Fclose(file);
}

void H5Event::open_charges(long nrows, long chunk){
	hsize_t dims[2] = {(hsize_t) nrows, 5};
	hsize_t cdims[2] = {(hsize_t) chunk, 5};

	charges_space = H5Screate_simple(2, dims, NULL);

	// Shuffle groups the bytes of equal significance together before gzip,
	// which is what makes float32 field data compress: 253 MB -> 210 MB on a
	// central PbPb event, and faster to write than without it.
	hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(plist, 2, cdims);
	H5Pset_shuffle(plist);
	H5Pset_deflate(plist, 5);

	charges = H5Dcreate2(file, "charges", H5T_IEEE_F32LE, charges_space,
	                     H5P_DEFAULT, plist, H5P_DEFAULT);
	H5Pclose(plist);
}

void H5Event::write_charges(double * rows, long n){
	hsize_t offset[2] = {(hsize_t) charges_row, 0};
	hsize_t count[2] = {(hsize_t) n, 5};

	H5Sselect_hyperslab(charges_space, H5S_SELECT_SET, offset, NULL, count, NULL);
	hid_t mem = H5Screate_simple(2, count, NULL);
	H5Dwrite(charges, H5T_NATIVE_DOUBLE, mem, charges_space, H5P_DEFAULT, rows);
	H5Sclose(mem);

	charges_row += n;
}

void H5Event::write_table(std::string name, double * data, long nrows, long ncols){
	hsize_t dims[2] = {(hsize_t) nrows, (hsize_t) ncols};
	hid_t space = H5Screate_simple(2, dims, NULL);

	hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(plist, 2, dims);
	H5Pset_deflate(plist, 5);

	hid_t set = H5Dcreate2(file, name.c_str(), H5T_IEEE_F32LE, space,
	                       H5P_DEFAULT, plist, H5P_DEFAULT);
	H5Dwrite(set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	H5Dclose(set);
	H5Pclose(plist);
	H5Sclose(space);
}

void H5Event::attr(std::string name, double value){
	hid_t space = H5Screate(H5S_SCALAR);
	hid_t a = H5Acreate2(file, name.c_str(), H5T_IEEE_F64LE, space,
	                     H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(a, H5T_NATIVE_DOUBLE, &value);
	H5Aclose(a);
	H5Sclose(space);
}

void H5Event::attr(std::string name, double * values, long n){
	hsize_t dims[1] = {(hsize_t) n};
	hid_t space = H5Screate_simple(1, dims, NULL);
	hid_t a = H5Acreate2(file, name.c_str(), H5T_IEEE_F64LE, space,
	                     H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(a, H5T_NATIVE_DOUBLE, values);
	H5Aclose(a);
	H5Sclose(space);
}

void H5Event::attr(std::string name, std::string value){
	hid_t type = H5Tcopy(H5T_C_S1);
	H5Tset_size(type, value.size()+1);
	hid_t space = H5Screate(H5S_SCALAR);
	hid_t a = H5Acreate2(file, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(a, type, value.c_str());
	H5Aclose(a);
	H5Sclose(space);
	H5Tclose(type);
}
