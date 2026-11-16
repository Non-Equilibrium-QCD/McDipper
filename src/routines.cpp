/* Copyright (c) 2025 Oscar Garcia-Montero
 * All rights reserved. */
#include "include/routines.h"

namespace Routines{

    const double chispar = 0.5;
    const int calls_miser_min    = pow(10,3);

    // GSL VEGAS  args
    const int INCREASE =10;
    const int loopcut = 6;
    const int vegas_iter_cut= 5;
    const double fEpsRel = 1e-2;
    const double fEpsAbs = 1e-10;

    //Cuhre Args
    const double EpsRelCuba=1e-4;
    const double EpsAbsCuba=1e-15;
    const int nVec =1;
    const int flags = 0;
    const int minEval= 100;
    // Cuhre stops as soon as EpsRelCuba is met, so this cap only costs time on
    // the hard points. With 5e6 the gluon table is accurate to about 0.4% (the
    // azimuthal tail to about 1.5%); with 5e4 entries carried up to 23% error.
    const int maxEval= 5000000;
    void * spin =NULL;
    void * statefile =NULL;
    const int key=13;
    
    const int nnew =1000;
    const int nmin =2;
    const double flatness= 25.;
    const int seed=0;

    

    void make_vegas_gen( gsl_rng *r_ptr, gsl_monte_function f2b, double * xl, double * xu, int DIM, int mc_calls_def, double &result, double &error)
    {
        bool MCcheck = true;
        int counter = 0;

        int MC_CALLS=mc_calls_def;
        do {
            // ====================== THIS IS THE ROUTINE =============================//
            gsl_monte_vegas_state *mcstate = gsl_monte_vegas_alloc(DIM);
    		
            gsl_monte_vegas_integrate (&f2b, xl, xu, DIM, MC_CALLS, r_ptr, mcstate, &result, &error);
            int whilecounter=0;
            do{
                gsl_monte_vegas_integrate (&f2b, xl, xu, DIM, MC_CALLS/2, r_ptr, mcstate, &result, &error);
                if (whilecounter>=vegas_iter_cut){break;}
                whilecounter++;
            }
            while (fabs (gsl_monte_vegas_chisq (mcstate) - 1.0) > chispar);
            gsl_monte_vegas_free (mcstate);

            // ====================== THIS IS THE ROUTINE =============================//

            if(error/result < fEpsRel){
                MCcheck = false;
    						break;
            }
            else if(result==0){
              MCcheck = false;
              break;
            }

            if(counter>=loopcut && result!=0){
    					std::cout<< "Result not achieved at tryout "<< counter<< " res = "<<result <<" +-"<< error <<std::endl;
    					break;}
            counter++;
        } while(MCcheck);
        if(MCcheck && result!=0){std::cerr << "\033[1mRoutines.h\033[0m : \033[1;31merror\033[0m : integration failure! Relative Error not achieved after " << loopcut<<" tryouts. Result = " << result << "+_" << error << std::endl;}
      }


    void make_cuhre_1C(int DIM, integrand_t integrand,void * parameters,double &result, double &error){

      int nregions, neval, fail;
      double integral_arr[1], error_arr[1], prob_arr[1];

      Cuhre(DIM,1, integrand, parameters, nVec,EpsRelCuba, EpsAbsCuba, flags, minEval, maxEval, key,"", spin, &nregions, &neval, &fail, integral_arr, error_arr, prob_arr);
      result= integral_arr[0];
      error= error_arr[0];
      // The tabulation loop calls this from several threads at once.
      #pragma omp atomic
      cuhre_calls()++;
      if(fail > 0){
        #pragma omp atomic
        cuhre_fails()++;
        double rel = (result!=0.0)? std::fabs(error/result) : 0.0;
        // Read, compare, write - not something atomic can cover.
        #pragma omp critical(cuhre_worst)
        { if(rel > cuhre_worst_rel()){ cuhre_worst_rel() = rel; } }
      }
      if(result!=result){
        #pragma omp atomic
        cuhre_nans()++;
        result=0.;error=0.0;
      }
    }

    void report_cuhre_stats(const char * where){
      if(cuhre_fails() > 0){
        std::cerr << "\n[ Cuhre ] " << where << ": " << cuhre_fails() << " of "
                  << cuhre_calls() << " integrations did not reach EpsRel="
                  << EpsRelCuba << " within maxEval=" << maxEval
                  << " (worst relative error " << cuhre_worst_rel() << ")" << std::endl;
      }
      else{
        std::cerr << "\n[ Cuhre ] " << where << ": all " << cuhre_calls()
                  << " integrations converged." << std::endl;
      }
      if(cuhre_nans() > 0){
        std::cerr << "[ Cuhre ] " << where << ": " << cuhre_nans() << " of "
                  << cuhre_calls() << " integrations returned NaN and were"
                  << " written as zero." << std::endl;
      }
    }

    void make_suave_1C(int DIM, integrand_t integrand,void * parameters,double &result, double &error){

      int nregions, neval, fail;
      double integral_arr[1], error_arr[1], prob_arr[1];

      Suave(DIM,1,integrand,parameters, nVec,EpsRelCuba, EpsAbsCuba,flags, seed, minEval, maxEval, nnew,nmin, flatness, "", spin,&nregions,&neval, &fail,integral_arr, error_arr, prob_arr);

      result= integral_arr[0];
      error= error_arr[0];
      if(result!=result){
        result=0.;error=0.0;}
    }
  }