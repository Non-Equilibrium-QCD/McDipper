/* Copyright (c) 2023 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef CONFIG_H
#define CONFIG_H
#include <iostream>
#include <fstream>
#include <string>

struct Impact{ double b1[2];double b2[2];};
enum class ImpSample : int { Fixed = 0, dbSampled = 1, bdbSampled = 2 };
enum class Model : int { GBW = 0, IPSat = 1, MV = 2, GBWSimp = 3};
enum class QuarkID : int { u = 2, d = 1, s = 3,ubar = -2, dbar = -1, sbar = -3, g=21 };
enum class GlauberMode : int { Standard = 0, Gaussian=1, Exponential=2};

struct NucStruct{ 
    int A; int Z;int mode;
    std::string inputFile;
    bool IsospinSpecified; int NConf; 
    bool is_weights;
    bool is_thick_fluct; 
    bool is_hotspots_fluct;
    int Nq; double Bq; double Br;
    std::string fluct_mode; double sigma;
};

class Config{
  /*  This is the config class. It reads in configuration files necessary for the EbE HICs initial conditions */ 
  public:
    Config();
    Config(std::string configfile_path);
    Config(const Config& OldConf);
    virtual ~Config();

    /*  **** General functions ****   */ 
    void process_parameters();
    //Version
    void check_version(std::string testline);
    //Logging
    void process_logging_parameters(std::string testline);
    // General
    void process_general_parameters(std::string testline);
    //Grid
    void process_grid_parameters(std::string testline);
    void fill_up_grid_params(); 
    // Model Parameters
    void process_model_parameters(std::string testline);
    void process_gbw_parameters(std::string testline);
    void process_ipsat_parameters(std::string testline);
    void process_mv_parameters(std::string testline);
    //Output
    void process_output_parameters(std::string testline);
    void write_config_file(std::string target_folder);
    void terminal_setup_output();
    void dump(std::string outpath);
    void set_dump(std::string outpath);

    //Output
    void process_thickness_parameters(std::string testline);

    // Operational Instructions

    bool compare_grid_parameters (Config *OldConf, double tolerance);
    bool compare_Thickness_parameters (Config *OldConf, double tolerance);
    bool compare_model_parameters (Config *OldConf, double tolerance);
    bool compare_PDF_parameters (Config *OldConf, double tolerance);

    //  RETRIEVING FUNCTIONS
    bool get_Verbose(){return Verbose;};
    // bool& getVerbose(){return Verbose;};
    //
    Model get_Model(){return cModel;};
    int get_A(int i){int A_t=-1;if(i==1){A_t=A1;}else if(i==2){A_t=A2;}return A_t;};
    int get_Z(int i){int Z_t=-1;if(i==1){Z_t=Z1;}else if(i==2){Z_t=Z2;}return Z_t;};
    int get_NuclearMode(int i){int mode_t=-1; if(i==1){mode_t=mode1;}else if(i==2){mode_t=mode2;}return mode_t;};
    std::string get_NucleusInput(int i){
      std::string input=""; 
      if(i==1 && mode1==3){input=Nucleus1ListFile;}
      else if(i==2 && mode2==3){input=Nucleus2ListFile;}
      return input;
    };
    bool get_IsospinDefinition(int i){
      bool IsoDef_t=false;
      if(i==1 && mode1==3){IsoDef_t= N1IsospinSpec;}
      else if(i==2 && mode2 ==3){IsoDef_t= N2IsospinSpec;}
      return IsoDef_t;
    }

    int get_NConf(int i){int NConf_t=0;if(i==1){NConf_t=NConf1;}else if(i==2){NConf_t=NConf2;}return NConf_t;};


    bool get_weight(int i){
      bool DefWeight_t=false;
      if(i==1 && mode1==3){DefWeight_t= N1_weights;}
      else if(i==2 && mode2 ==3){DefWeight_t= N2_weights;}
      return DefWeight_t;}

    bool is_hotspots_fluct(){return hotspots_fluct;}
    int get_Nq() {return Nq;}
    double get_Bq() {return Bq;}
    double get_Br() {
      //The sampling width of hotspots. Br=(BG-Bq)/(1-1/Nq)
      if (Nq==1){return 0.0;}
      double Br = (BG-Bq) / (1.0-1.0/Nq);
      return Br;
    }

    std::string get_PDFSet(){return cPDFSetStr;}
    int get_ForcedPositive(){return cForcedMode;}

    double get_collEnergy(){return sqrtsNN;}
    bool compare_collEnergy(double CollEnergy, double tolerance){ return ( NDiff(sqrtsNN,CollEnergy) < tolerance ) ;}

    ImpSample get_ImpactMode(){return ImpactMode;}
    double get_ImpactValue(){return ImpactValue;}
    double get_MaxImpact(){return bMax;}
    double get_MinImpact(){return bMin;}
    double get_KFactor(){return KFactor;}
    int get_NEvents(){return NEvents;}
    GlauberMode get_GlauberAcceptance(){return GMode;}
    
    // setter functions if the code is used as a library
    void set_ImpactMode(ImpSample val){ImpactMode = val;}
    void set_ImpactValue(double val){ImpactValue = val;}
    void set_MaxImpact(double val){bMax = val;}
    void set_MinImpact(double val){bMin = val;}
    void set_KFactor(double val){KFactor = val;}
    void set_Seed(int val){seed = val;}

    //Grid!

    int get_NX(){return NX;}
    int get_NY(){return NY;}
    int get_NETA(){return NETA;}
    double get_XMIN(){return XMIN;}
    double get_XMAX(){return XMAX;}
    double get_YMIN(){return YMIN;}
    double get_YMAX(){return YMAX;}
    double get_ETAMIN(){return ETAMIN;}
    double get_ETAMAX(){return ETAMAX;}

    void set_NX(int val){NX=val;}
    void set_NY(int val){NY=val;}
    void set_NETA(int val){NETA=val;}
    void set_XMIN(double val){XMIN=val;}
    void set_XMAX(double val){XMAX=val;}
    void set_YMIN(double val){YMIN=val;}
    void set_YMAX(double val){YMAX=val;}
    void set_ETAMIN(double val){ETAMIN=val;}
    void set_ETAMAX(double val){ETAMAX=val;}

    double get_dX(){return dX;}
    double get_dY(){return dY;}
    double get_dETA(){return dETA;}

    double get_BG(){return BG;}

    double get_ModelParams(int i){return ModelPars[i];}

    int get_seed(){return seed;}
    void set_seed();

    Model get_model(){return cModel;}
    int get_mpars_number(){return NModelParams;}

    //Thickess

    double get_TMax(){return TMax;}
    double get_TMin(){return TMin;}
    int get_NT(){return NT;}
    double get_dT(){return dT;}
    bool is_thick_fluct(){return thick_fluct;}
    double get_sigma(){return sigma;}
    std::string get_fluct_mode() {return fluct_mode;}

    void set_TMax(double TMax_new){TMax=TMax_new;}
    void set_TMin(double TMin_new){TMin=TMin_new;}
    void set_NT(int NT_new){NT=NT_new;}
    void set_dT(){dT=(TMax-TMin)/(NT-1);}

    //Output
    std::string get_config_path(){return path_to_configfile;}
    std::string get_out_path(){return path_to_output;}
    std::string get_run_name(){return run_name;}
    std::string get_version(){return version;}
    bool is_output_named(){return is_run_renamed;}
    bool is_format(std::string test_format){
      for (int i = 0; i < n_formats; i++){
        if(format[i]==test_format){return true;}
      }
      return false;
    }
    int print_avg_event(){return print_avg;}
    bool is_boost_invariant(){return boost_invariant;}
    void set_suppress_output(bool flag){ suppress_output = flag; }
    bool is_output_suppressed() const { return suppress_output; }

  private:

    std::string path_to_configfile;
    //Version
    std::string version;
    //Logging
    bool Verbose=false;
    // General
    Model cModel;
    bool hotspots_fluct=false;
    int Nq = 0;
    double Bq = 0.0;

    bool thick_fluct=false;
    double sigma=-1.0;
    std::string fluct_mode="Uniform";
    
    std::string cModelStr;
    int A1,A2;
    int Z1,Z2;
    int mode1, mode2;
    std::string Nucleus1ListFile="";
    std::string Nucleus2ListFile="";
    bool N1IsospinSpec=false;
    bool N2IsospinSpec=false;
    int NConf1 = 0;
    int NConf2 = 0;
    bool N1_weights=false;
    bool N2_weights=false;
    
    std::string mode1name,mode2name;
    double sqrtsNN;
    GlauberMode GMode=GlauberMode::Standard;
    std::string GModeStr="Standard";

    // int ImpactMode=3;
    ImpSample ImpactMode=ImpSample::bdbSampled;
    double ImpactValue;
    double bMin, bMax;
    double KFactor=1.;
    int seed=-1;
    int NEvents=1;


    //PDFs
    std::string cPDFSetStr;
    int cForcedMode=0;
    //Grid
    int NX,NY,NETA;
    double XMAX,XMIN; /* fm */
    double YMAX,YMIN; /* fm */
    double ETAMAX,ETAMIN;
    double dX,dY,dETA; /* fm */
    double BG; /// In fm^2
    // // Model Parameters
    double ModelPars[10];
    int NModelParams;
    // Output
    int n_formats;
    std::string path_to_output;
    std::string run_name;
    bool is_run_renamed=false;
    std::string * format;
    int print_avg=0;
    bool boost_invariant=false;
    bool suppress_output=false;

    // Thickness_Parameters
    double TMax=-1;
    double TMin=-1;
    int NT=-1;
    double dT;
    
    //Tools
    int count_indent(std::string testline);
    std::string get_header(std::string testline);
    void get_name_and_value(std::string testline,std::string &name, std::string &value);
    bool make_bool(std::string boolstring);

    std::string bool_string(bool testbool){if(testbool){return "True";}else{return "False";}}

    void remove_char(std::string &str_rem,char rem_char);
    void retrieve_range(std::string rangestring, double &xmin, double &xmax);

    double NDiff(double a, double b)                       
    {                                                    
        if (a==0.0 || b==0.0) {return fabs(a-b);}        
        else {return fabs( (a-b)/a );}                   
    }

    //Temps
    std::string subheader;
    int cModel_int;
};

#endif  // config //
