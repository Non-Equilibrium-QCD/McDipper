/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef ROUTINES_H
#define ROUTINES_H
#include <iostream>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "../external/cuba.h"

namespace Routines {

    // Constants
    extern const double chispar;
    extern const int calls_miser_min;

    extern const int INCREASE;
    extern const int loopcut;
    extern const int vegas_iter_cut;
    extern const double fEpsRel;
    extern const double fEpsAbs;

    // Cuba
    extern const double EpsRelCuba;
    extern const double EpsAbsCuba;
    extern const int nVec;
    extern const int flags;
    extern const int minEval;
    extern const int maxEval;
    extern void* spin;
    extern void* statefile;
    extern const int key;

    extern const int nnew;
    extern const int nmin;
    extern const double flatness;
    extern const int seed;

    // Function declarations only
    void make_vegas_gen(gsl_rng *r_ptr, gsl_monte_function f2b,
                        double * xl, double * xu, int DIM,
                        int mc_calls_def, double &result, double &error);

    void make_cuhre_1C(int DIM, integrand_t integrand, void * parameters,
                       double &result, double &error);

    void make_suave_1C(int DIM, integrand_t integrand, void * parameters,
                       double &result, double &error);

    // Cuhre returns fail>0 when the requested accuracy was not reached within
    // maxEval. Keep a tally of those and report it once per table.
    inline long   & cuhre_calls()    { static long   n = 0;   return n; }
    inline long   & cuhre_fails()    { static long   n = 0;   return n; }
    inline double & cuhre_worst_rel(){ static double e = 0.0; return e; }
    // make_cuhre_1C replaces a NaN result by zero. Count these so they are
    // reported together with the failed calls.
    inline long   & cuhre_nans()     { static long   n = 0;   return n; }

    void report_cuhre_stats(const char * where);

}


#endif /* routines.h */
 