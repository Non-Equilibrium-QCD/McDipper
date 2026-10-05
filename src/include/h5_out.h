/* Copyright (c) 2026 */

#ifndef H5_OUT_H
#define H5_OUT_H
#include <string>
#include <hdf5.h>

class H5Event{
	public:
		H5Event();
		H5Event(std::string path);
		virtual ~H5Event();

		// charges is NETA*NX*NY x 5 - the densities eg, eq, nu, nd, ns. The
		// coordinates are not stored: rows run eta, then x, then y, all
		// ascending, so the grid follows from NETA/NX/NY and the ranges in the
		// config the run writes out. Created once, filled one slice at a time.
		void open_charges(long nrows, long chunk);
		void write_charges(double * rows, long n);
		void write_table(std::string name, double * data, long nrows, long ncols);

		void attr(std::string name, double value);
		void attr(std::string name, double * values, long n);
		void attr(std::string name, std::string value);

	private:
		hid_t file;
		hid_t charges;
		hid_t charges_space;
		long charges_row;
};

#endif /* h5_out.h */
