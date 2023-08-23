/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */
 #include <iostream>
 #include <math.h>
 #include <sstream>
 #include <sys/stat.h>
 #include <filesystem>


#include "include/charges.h"
#include "include/model_gbw.h"

namespace fs = std::filesystem;

Charges::Charges(){}
Charges::Charges(Config ConfInput){
  config=Config(ConfInput);

  NETA=config.get_NETA();

  check_for_tabs_folders();
  bool is_imported=import_charges(); // Import charges from computed table.

  if(!is_imported){
    call_models(); // Calls the models if there are no matching set of tables
    is_imported=import_charges();// Imports the newly constructed tables
    if(is_imported){initialized = true;}
  }
  if(!is_imported){std::cerr<<"Error after importing newly constructed charges!"<<std::endl;exit(EXIT_FAILURE);}
  set_output_tests();
}

Charges::~Charges(){
  for (size_t i = 0; i < NETA; i++) {
    gsl_spline2d_free(e_g_spl[i]);
    gsl_interp_accel_free (xaccEG[i]);
    gsl_interp_accel_free (yaccEG[i]);
	}

  for (size_t i = 0; i < NQComp; i++) {
    gsl_spline2d_free(N0u_spl[i]);
    gsl_spline2d_free(N0d_spl[i]);
    gsl_spline2d_free(N0s_spl[i]);

    gsl_spline2d_free(N1u_spl[i]);
    gsl_spline2d_free(N1d_spl[i]);
    gsl_spline2d_free(N1s_spl[i]);

    gsl_interp_accel_free (xaccN0u[i]);
    gsl_interp_accel_free (yaccN0u[i]);
    gsl_interp_accel_free (xaccN1u[i]);
    gsl_interp_accel_free (yaccN1u[i]);

    gsl_interp_accel_free (xaccN0d[i]);
    gsl_interp_accel_free (yaccN0d[i]);
    gsl_interp_accel_free (xaccN1d[i]);
    gsl_interp_accel_free (yaccN1d[i]);

    gsl_interp_accel_free (xaccN0s[i]);
    gsl_interp_accel_free (yaccN0s[i]);
    gsl_interp_accel_free (xaccN1s[i]);
    gsl_interp_accel_free (yaccN1s[i]);
  }
}

//Fill-in Functions

void Charges::MakeGrid(){

  NX=config_set.get_NT();NY=config_set.get_NT();

  e_g_spl = new gsl_spline2d*[NETA] ;

  xaccEG = new gsl_interp_accel*[NETA] ;
  yaccEG = new gsl_interp_accel*[NETA] ;

  N0u_spl = new gsl_spline2d*[NQComp] ;
  N1u_spl = new gsl_spline2d*[NQComp] ;

  N0d_spl = new gsl_spline2d*[NQComp] ;
  N1d_spl = new gsl_spline2d*[NQComp] ;

  N0s_spl = new gsl_spline2d*[NQComp] ;
  N1s_spl = new gsl_spline2d*[NQComp] ;

  xaccN0u = new gsl_interp_accel*[NQComp] ;
  yaccN0u = new gsl_interp_accel*[NQComp] ;

  xaccN0d = new gsl_interp_accel*[NQComp] ;
  yaccN0d = new gsl_interp_accel*[NQComp] ;

  xaccN0s = new gsl_interp_accel*[NQComp] ;
  yaccN0s = new gsl_interp_accel*[NQComp] ;

  xaccN1u = new gsl_interp_accel*[NQComp] ;
  yaccN1u = new gsl_interp_accel*[NQComp] ;

  xaccN1d = new gsl_interp_accel*[NQComp] ;
  yaccN1d = new gsl_interp_accel*[NQComp] ;

  xaccN1s = new gsl_interp_accel*[NQComp] ;
  yaccN1s = new gsl_interp_accel*[NQComp] ;

  if(config.get_Verbose()){std::cout<<"[ Charges ] Grid created for Charges."<<std::endl;}

}
 
bool Charges::import_charges(){

  bool import_is_success=true;
  bool is_table_set=check_sets();
  // If set is found, import
  if(is_table_set){
    MakeGrid();
    import_is_success = import_is_success && read_in_energy_gluons();
    if(config.get_Verbose()){std::cout<< "[ Charges ] Gluon energy density    -->  Imported "<< std::endl;}
    import_is_success = import_is_success && read_in_nK_quark(0, QuarkID::u );
    import_is_success = import_is_success && read_in_nK_quark(1, QuarkID::u );
    if(config.get_Verbose()){std::cout<< "[ Charges ] U-Quark constructions   -->  Imported "<< std::endl;}
    import_is_success = import_is_success && read_in_nK_quark(0, QuarkID::d );
    import_is_success = import_is_success && read_in_nK_quark(1, QuarkID::d );
    if(config.get_Verbose()){std::cout<< "[ Charges ] D-Quark constructions   -->  Imported "<< std::endl;}
    import_is_success = import_is_success && read_in_nK_quark(0, QuarkID::s );
    import_is_success = import_is_success && read_in_nK_quark(1, QuarkID::s );
    if(config.get_Verbose()){std::cout<< "[ Charges ] S-Quark constructions   -->  Imported "<< std::endl;}
    import_is_success=true;
  }
  else{import_is_success=false;}

  return import_is_success;
}

bool Charges::call_models(){

  bool call_is_success=false;
  if(config.get_model()== Model::GBW ){
    GBW gbw(config);
    gbw.MakeTable(path_to_set);
  }
  if(config.get_model()== Model::IPSat ){
    IPSat ipsat(config);
    ipsat.MakeTable(path_to_set);
  }
  if(config.get_model()== Model::MV ){}
  return call_is_success;
}

bool Charges::check_sets(){
  bool is_computed =false;
  int set=0;// counter for the sets (which are numbered as Setn, n={0,...,N})
  path_to_set=get_set_name(set);
  bool is_path=IsPathExist(path_to_set);

  while (is_path) {
    is_computed = check_config();// Checks configfiles match
    if(is_computed){break;}
    else{
      set++;
      path_to_set=get_set_name(set);// goes to next Set and
      is_path=IsPathExist(path_to_set); // checks if exists

    }
  }

  return is_computed;
}

bool Charges::check_config(){
  bool is_config=true;
  // std::string line,name_t,value_t;

  std::ostringstream path_to_config;
  path_to_config << path_to_set<< "/config.yaml" ;

  Config * set_conf= new Config(path_to_config.str());

  /// Grid Comparison
  is_config = is_config && config.compare_grid_parameters(set_conf,double_tolerance);
  is_config = is_config && config.compare_model_parameters(set_conf,double_tolerance);
  is_config = is_config && config.compare_collEnergy(set_conf->get_collEnergy(), double_tolerance);
  is_config = is_config && config.compare_PDF_parameters(set_conf, double_tolerance);

  if(is_config){config_set = Config(path_to_config.str());}
  return is_config;
}


//////////////////////////////
////////// TOOLS /////////////
//////////////////////////////

bool Charges::IsPathExist(const std::string &s){
	struct stat buffer;
	return (stat (s.c_str(), &buffer) == 0);
}

void Charges::check_for_tabs_folders(){
  std::ostringstream tabspath;
  fs::path PWD = fs::current_path();
  tabspath << PWD.string() << "/tabs" ;
  if(!IsPathExist(tabspath.str())){ fs::create_directories(tabspath.str()); }

  if(config.get_model()== Model::GBW ){
    std::ostringstream gbwpath;
    gbwpath << PWD.string()<< "/tabs/GBW" ;
    path_to_tabs = gbwpath.str();
    if(!IsPathExist(path_to_tabs)){ fs::create_directories(path_to_tabs );}
  }
  if(config.get_model()== Model::IPSat ){
    std::ostringstream ipsatpath;
    ipsatpath << PWD.string()<< "/tabs/IPSat" ;
    path_to_tabs = ipsatpath.str();
    if(!IsPathExist(path_to_tabs)){ fs::create_directories( path_to_tabs);}
  }
  if(config.get_model()== Model::MV ){
    std::ostringstream mvpath;
    mvpath << PWD.string()<< "/tabs/MV" ;
    path_to_tabs = mvpath.str();
    if(!IsPathExist(path_to_tabs)){ fs::create_directories(path_to_tabs );}
  }

}

std::string Charges::get_set_name(int n){
  std::ostringstream temppath;
  temppath << path_to_tabs<< "/Set"<< n;
  return temppath.str();
}


void Charges::get_name_and_value(std::string testline,std::string &name, std::string &value){
  std::string line_temp= testline;
  size_t pos = line_temp.find(":");
  name = line_temp.substr(0,pos);
  value = line_temp.substr(pos+1, line_temp.length()-pos );
  name.erase(std::remove(name.begin(),name.end(), ' '),name.end());
  value.erase(std::remove(value.begin(),value.end(), ' '),value.end());
}

////// Read in routines

bool Charges::read_in_energy_gluons(){

  bool is_read=true;
  std::ostringstream tablename;
  tablename << path_to_set <<"/"<< gluon_energy_table_name ;

  FILE *table = fopen(tablename.str().c_str(),"r");
  double y_t,T1_t,T2_t,edensG_t;

  double T1Grid[NX*NY];
  double T2Grid[NX*NY];
  double T1array[NX];
  double T2array[NY];
  double etau0_g[NX*NY];

  for (int iY=0; iY < NETA; iY++)
  {
    for (int i1=0; i1 < NX; i1++)
    {
      for (int i2=0; i2 < NY; i2++)
      {
        if (fscanf(table,"%lf %lf %lf %lf", &y_t,&T1_t,&T2_t,&edensG_t) != 4){printf("Error reading distribution table!\n");exit(1);}
        T1Grid[NX*i2+i1]=T1_t;
        T2Grid[NX*i2+i1]=T2_t;
        etau0_g[NX*i2+i1]=edensG_t;
      }
    }
    for (size_t i1 = 0; i1 < NX; i1++) {T1array[i1]=T1Grid[i1]; }
    for (size_t i2 = 0; i2 < NY; i2++) {T2array[i2]=T2Grid[NX*i2]; }

    size_t nx = sizeof(T1array) / sizeof(T1array[0]);
    size_t ny = sizeof(T2array) / sizeof(T2array[0]);

    const gsl_interp2d_type *T= gsl_interp2d_bilinear;
    e_g_spl[iY] = gsl_spline2d_alloc(T, nx, ny);
    xaccEG[iY] = gsl_interp_accel_alloc();
    yaccEG[iY] = gsl_interp_accel_alloc();
    gsl_spline2d_init(e_g_spl[iY], T1array, T2array, &etau0_g[0], nx, ny);
  }
  fclose(table);
  // test_f.close();
  return is_read;
}

bool Charges::read_in_nK_quark(int k, QuarkID qid){


  bool is_read=true;
  std::ostringstream tablename;
  tablename << path_to_set <<"/"<< get_quark_construction_file(k,qid) ;

  FILE *table = fopen(tablename.str().c_str(),"r");
  double y_t,T_t,nk_12_q,nk_21_q,nk_12_qbar,nk_21_qbar;

  double Yarray[NETA];
  double Tarray[NX];



  double Nk_t[NQComp][NX*NETA];

  for (int iY=0; iY < NETA; iY++)
  {
    for (int i1=0; i1 < NX; i1++)
    {
        if (fscanf(table,"%lf %lf %lf %lf %lf %lf", &y_t,&T_t,&nk_12_q,&nk_12_qbar,&nk_21_q,&nk_21_qbar) != 6){printf("Error reading distribution table!\n");exit(1);}
        if(iY==0){Tarray[i1]=T_t;};
        if(i1==0){Yarray[iY]=y_t;};

        Nk_t[0][NETA*i1+iY]=nk_12_q;
        Nk_t[1][NETA*i1+iY]=nk_12_qbar;
        Nk_t[2][NETA*i1+iY]=nk_21_q;
        Nk_t[3][NETA*i1+iY]=nk_21_qbar;
    }
  }
  fclose(table);

  size_t neta = sizeof(Yarray) / sizeof(Yarray[0]);
  size_t nx = sizeof(Tarray) / sizeof(Tarray[0]);

  for (size_t ik = 0; ik < NQComp; ik++) {
    const gsl_interp2d_type *T= gsl_interp2d_bilinear;
    if(k==0){
      if(qid==QuarkID::u){
        N0u_spl[ik] = gsl_spline2d_alloc(T, neta, nx);
        xaccN0u[ik] = gsl_interp_accel_alloc();yaccN0u[ik] = gsl_interp_accel_alloc();
        gsl_spline2d_init(N0u_spl[ik], Yarray, Tarray, &Nk_t[ik][0],neta, nx);
      }else if(qid==QuarkID::d){
        N0d_spl[ik] = gsl_spline2d_alloc(T, neta, nx);
        xaccN0d[ik] = gsl_interp_accel_alloc();yaccN0d[ik] = gsl_interp_accel_alloc();
        gsl_spline2d_init(N0d_spl[ik], Yarray, Tarray, &Nk_t[ik][0],neta, nx);
      }else if(qid==QuarkID::s){
        N0s_spl[ik] = gsl_spline2d_alloc(T, neta, nx);
        xaccN0s[ik] = gsl_interp_accel_alloc();yaccN0s[ik] = gsl_interp_accel_alloc();
        gsl_spline2d_init(N0s_spl[ik], Yarray, Tarray, &Nk_t[ik][0],neta, nx);
      }
    }else if(k==1){
      if(qid==QuarkID::u){
        N1u_spl[ik] = gsl_spline2d_alloc(T, neta, nx);
        xaccN1u[ik] = gsl_interp_accel_alloc();yaccN1u[ik] = gsl_interp_accel_alloc();
        gsl_spline2d_init(N1u_spl[ik], Yarray, Tarray, &Nk_t[ik][0],neta, nx);
      }else if(qid==QuarkID::d){
        N1d_spl[ik] = gsl_spline2d_alloc(T, neta, nx);
        xaccN1d[ik] = gsl_interp_accel_alloc();yaccN1d[ik] = gsl_interp_accel_alloc();
        gsl_spline2d_init(N1d_spl[ik], Yarray, Tarray, &Nk_t[ik][0],neta, nx);
      }else if(qid==QuarkID::s){
        N1s_spl[ik] = gsl_spline2d_alloc(T, neta, nx);
        xaccN1s[ik] = gsl_interp_accel_alloc();yaccN1s[ik] = gsl_interp_accel_alloc();
        gsl_spline2d_init(N1s_spl[ik], Yarray, Tarray, &Nk_t[ik][0],neta, nx);
      }
    }



  }
  return is_read;
}


////// Evaluate


double Charges::gluon_energy(double eta, double T1, double T2){
  if ( T1 < config_set.get_TMin() ){return 0.0;}
  if ( T1 > config_set.get_TMax() ){
    if(T_warning){
      std::cerr<<"[ Charges::Warning ] Accessed T(1) larger than range. Increase T range to avoid missing nuclear matter.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  if ( T2 < config_set.get_TMin() ){return 0.0;}
  if ( T2 > config_set.get_TMax() ){
    if(T_warning){
      std::cerr<<"[ Charges::Warning ] Accessed T(2) larger than range. Increase T range to avoid missing nuclear matter.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  if ( eta < config.get_ETAMIN() || eta > config.get_ETAMAX() ){
    if(eta_warning){
      std::cerr<<"[ Charges::Warning ] Accessed eta not in range. Increase eta range to include forward/backward rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  else if ( eta == config.get_ETAMIN() ){return gsl_spline2d_eval(e_g_spl[0],T1,T2,xaccEG[0], yaccEG[0]);}
  else if ( eta == config.get_ETAMAX() ){return gsl_spline2d_eval(e_g_spl[config.get_NETA()-1],T1,T2,xaccEG[config.get_NETA()-1], yaccEG[config.get_NETA()-1]);}
  else{
    int iETA;
    double ETA1,ETA2;
    double tmp1,tmp2;

    iETA = (int) ( (eta - config.get_ETAMIN() )/ config.get_dETA() );
    ETA1= iETA * config.get_dETA() + config.get_ETAMIN() ;
    ETA2= (iETA+1) * config.get_dETA() + config.get_ETAMIN() ;

    tmp1 = gsl_spline2d_eval(e_g_spl[iETA],T1,T2,xaccEG[iETA], yaccEG[iETA]);
    tmp2 = gsl_spline2d_eval(e_g_spl[iETA+1],T1,T2,xaccEG[iETA+1], yaccEG[iETA+1]);

    return ( tmp1 + (tmp2-tmp1) * (eta-ETA1)/(ETA2-ETA1) ) ; //Linear interp in ETA
  }
}

double Charges::quark_energy(double eta, double T1, double T2){
  if ( T1 < config_set.get_TMin() ){return 0.0;}
  if ( T1 > config_set.get_TMax() ){
    if(T_warning){
      std::cerr<<"[ Charges::Warning ] Accessed T(1) larger than range. Increase T range to avoid missing rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  if ( T2 < config_set.get_TMin() ){return 0.0;}
  if ( T2 > config_set.get_TMax() ){
    if(T_warning){
      std::cerr<<"[ Charges::Warning ] Accessed T(2) larger than range. Increase T range to avoid missing rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  if ( eta < config.get_ETAMIN() || eta > config.get_ETAMAX() ){
    if(eta_warning){
      std::cerr<<"[ Charges::Warning ] Accessed eta not in range. Increase eta range to include forward/backward rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  else{

    double TEMP=0;
    // U QUARK
    TEMP += T1*gsl_spline2d_eval(N1u_spl[0],eta,T2,xaccN1u[0], yaccN1u[0]); // 12 u
    TEMP += T1*gsl_spline2d_eval(N1u_spl[1],eta,T2,xaccN1u[1], yaccN1u[1]); // 12 ubar
    TEMP += T2*gsl_spline2d_eval(N1u_spl[2],eta,T1,xaccN1u[2], yaccN1u[2]); // 21 u
    TEMP += T2*gsl_spline2d_eval(N1u_spl[3],eta,T1,xaccN1u[3], yaccN1u[3]); // 21 ubar
    // D QUARK
    TEMP += T1*gsl_spline2d_eval(N1d_spl[0],eta,T2,xaccN1d[0], yaccN1d[0]); // 12 d
    TEMP += T1*gsl_spline2d_eval(N1d_spl[1],eta,T2,xaccN1d[1], yaccN1d[1]); // 12 dbar
    TEMP += T2*gsl_spline2d_eval(N1d_spl[2],eta,T1,xaccN1d[2], yaccN1d[2]); // 21 d
    TEMP += T2*gsl_spline2d_eval(N1d_spl[3],eta,T1,xaccN1d[3], yaccN1d[3]); // 21 dbar
    // S QUARK
    TEMP += T1*gsl_spline2d_eval(N1s_spl[0],eta,T2,xaccN1s[0], yaccN1s[0]); // 12 s
    TEMP += T1*gsl_spline2d_eval(N1s_spl[1],eta,T2,xaccN1s[1], yaccN1s[1]); // 12 sbar
    TEMP += T2*gsl_spline2d_eval(N1s_spl[2],eta,T1,xaccN1s[2], yaccN1s[2]); // 21 s
    TEMP += T2*gsl_spline2d_eval(N1s_spl[3],eta,T1,xaccN1s[3], yaccN1s[3]); // 21 sbar

    return TEMP ;
  }
}

double Charges::u_density(double eta, double T1p,double T1n, double T2p, double T2n){
  // TO conserve charges, we apply isospin symmetry
  // meaning ->  p:u -> n:d. Therefore, for the u density
  // n_u = n_u(Tp)+n_d(Tn)
  double T1 = T1p + T1n;
  double T2 = T2p + T2n;

  if ( T1 < config_set.get_TMin() ){return 0.0;}
  if ( T1 > config_set.get_TMax() ){
    if(T_warning){
      std::cerr<<"[ Charges::Warning ] Accessed T(1) larger than range. Increase T range to avoid missing rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  if ( T2 < config_set.get_TMin() ){return 0.0;}
  if ( T2 > config_set.get_TMax() ){
    if(T_warning){
      std::cerr<<"[ Charges::Warning ] Accessed T(2) larger than range. Increase T range to avoid missing rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  if ( eta < config.get_ETAMIN() || eta > config.get_ETAMAX() ){
    if(eta_warning){
      std::cerr<<"[ Charges::Warning ] Accessed eta not in range. Increase eta range to include forward/backward rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  else{

    double TEMP=0;
    // U QUARK - p
    TEMP += T1p*gsl_spline2d_eval(N0u_spl[0],eta,T2,xaccN0u[0], yaccN0u[0]); // 12 u
    TEMP -= T1p*gsl_spline2d_eval(N0u_spl[1],eta,T2,xaccN0u[1], yaccN0u[1]); // 12 ubar
    TEMP += T2p*gsl_spline2d_eval(N0u_spl[2],eta,T1,xaccN0u[2], yaccN0u[2]); // 21 u
    TEMP -= T2p*gsl_spline2d_eval(N0u_spl[3],eta,T1,xaccN0u[3], yaccN0u[3]); // 21 ubar

    // D QUARK - n
    TEMP += T1n*gsl_spline2d_eval(N0d_spl[0],eta,T2,xaccN0d[0], yaccN0d[0]); // 12 u
    TEMP -= T1n*gsl_spline2d_eval(N0d_spl[1],eta,T2,xaccN0d[1], yaccN0d[1]); // 12 ubar
    TEMP += T2n*gsl_spline2d_eval(N0d_spl[2],eta,T1,xaccN0d[2], yaccN0d[2]); // 21 u
    TEMP -= T2n*gsl_spline2d_eval(N0d_spl[3],eta,T1,xaccN0d[3], yaccN0d[3]); // 21 ubar

    return TEMP ;
  }
}

double Charges::d_density(double eta, double T1p,double T1n, double T2p, double T2n){
  // TO conserve charges, we apply isospin symmetry
  // meaning ->  p:d -> n:u. Therefore, for the d density
  // n_d = n_d(Tp)+n_u(Tn)
  double T1 = T1p + T1n;
  double T2 = T2p + T2n;

  if ( T1 < config_set.get_TMin() ){return 0.0;}
  if ( T1 > config_set.get_TMax() ){
    if(T_warning){
      std::cerr<<"[ Charges::Warning ] Accessed T(1) larger than range. Increase T range to avoid missing rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  if ( T2 < config_set.get_TMin() ){return 0.0;}
  if ( T2 > config_set.get_TMax() ){
    if(T_warning){
      std::cerr<<"[ Charges::Warning ] Accessed T(2) larger than range. Increase T range to avoid missing rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  if ( eta < config.get_ETAMIN() || eta > config.get_ETAMAX() ){
    if(eta_warning){
      std::cerr<<"[ Charges::Warning ] Accessed eta not in range. Increase eta range to include forward/backward rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  else{

    double TEMP=0;
    // U QUARK - n
    TEMP += T1n*gsl_spline2d_eval(N0u_spl[0],eta,T2,xaccN0u[0], yaccN0u[0]); // 12 u
    TEMP -= T1n*gsl_spline2d_eval(N0u_spl[1],eta,T2,xaccN0u[1], yaccN0u[1]); // 12 ubar
    TEMP += T2n*gsl_spline2d_eval(N0u_spl[2],eta,T1,xaccN0u[2], yaccN0u[2]); // 21 u
    TEMP -= T2n*gsl_spline2d_eval(N0u_spl[3],eta,T1,xaccN0u[3], yaccN0u[3]); // 21 ubar

    // D QUARK - p
    TEMP += T1p*gsl_spline2d_eval(N0d_spl[0],eta,T2,xaccN0d[0], yaccN0d[0]); // 12 u
    TEMP -= T1p*gsl_spline2d_eval(N0d_spl[1],eta,T2,xaccN0d[1], yaccN0d[1]); // 12 ubar
    TEMP += T2p*gsl_spline2d_eval(N0d_spl[2],eta,T1,xaccN0d[2], yaccN0d[2]); // 21 u
    TEMP -= T2p*gsl_spline2d_eval(N0d_spl[3],eta,T1,xaccN0d[3], yaccN0d[3]); // 21 ubar

    return TEMP ;
  }
}

double Charges::s_density(double eta, double T1, double T2){
  if ( T1 < config_set.get_TMin() ){return 0.0;}
  if ( T1 > config_set.get_TMax() ){
    if(T_warning){
      std::cerr<<"[ Charges::Warning ] Accessed T(1) larger than range. Increase T range to avoid missing rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  if ( T2 < config_set.get_TMin() ){return 0.0;}
  if ( T2 > config_set.get_TMax() ){
    if(T_warning){
      std::cerr<<"[ Charges::Warning ] Accessed T(2) larger than range. Increase T range to avoid missing rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  if ( eta < config.get_ETAMIN() || eta > config.get_ETAMAX() ){
    if(eta_warning){
      std::cerr<<"[ Charges::Warning ] Accessed eta not in range. Increase eta range to include forward/backward rapidities.";
      set_eta_warning_off();
    }
    return 0.0;
  }
  else{
    double TEMP=0;
    // S QUARK
    TEMP += T1*gsl_spline2d_eval(N0s_spl[0],eta,T2,xaccN0s[0], yaccN0s[0]); // 12 s
    TEMP -= T1*gsl_spline2d_eval(N0s_spl[1],eta,T2,xaccN0s[1], yaccN0s[1]); // 12 sbar
    TEMP += T2*gsl_spline2d_eval(N0s_spl[2],eta,T1,xaccN0s[2], yaccN0s[2]); // 21 s
    TEMP -= T2*gsl_spline2d_eval(N0s_spl[3],eta,T1,xaccN0s[3], yaccN0s[3]); // 21 sbar

    return TEMP ;
  }
}

std::string Charges::get_quark_construction_file(int k, QuarkID qid){
  std::ofstream quark_name_f;
  std::ostringstream quark_name;

  quark_name << "quark_construction_n"<< k;
  if(qid==QuarkID::u){quark_name <<"_u.dat";}
  else if(qid==QuarkID::d){quark_name <<"_d.dat";}
  else if(qid==QuarkID::s){quark_name <<"_s.dat";}

  return quark_name.str();
}


void Charges::dump_charges(double T1,double T2){
  std::ofstream density_f;
  std::ostringstream densityname;
  densityname << path_to_set <<"/charges_dump_T1_"<< T1<< "_T2_"<< T2<<".dat";
	density_f.open(densityname.str());

	for (size_t iy = 0; iy < config.get_NETA(); iy++) {
		double y_t = iy*config.get_dETA() + config.get_ETAMIN();
		// density_f<< std::endl;
    density_f<<y_t<< "\t"<<gluon_energy(y_t, T1,T2);
    density_f<< "\t"<<quark_energy(y_t, T1,T2);
    density_f<< "\t"<<u_density(y_t, T1,0, T2, 0);
    density_f<< "\t"<<d_density(y_t, T1,0, T2, 0);
    density_f<< "\t"<<s_density(y_t, T1,T2)<< std::endl;
	}
 	density_f.close();
}


void Charges::dump_charges_eta(double eta ){
  std::ofstream density_f;
  std::ostringstream densityname;
  densityname << path_to_set <<"/charges_dump_eta_"<< eta <<".dat";
	density_f.open(densityname.str());

  double Tmax_t=8.;
  double Tmin_t=0;
  double dT_t= (Tmax_t-Tmin_t)/(double(NX) - 1.);
  for (size_t i1 = 0; i1 < NX; i1++) {
    double T1= i1*dT_t + Tmin_t;
    for (size_t i2 = 0; i2 < NY; i2++) {
      double T2= i2*dT_t + Tmin_t;
      density_f<<T1<< "\t"<<T2<< "\t"<<gluon_energy(eta, T1,T2);
      density_f<< "\t"<<quark_energy(eta, T1,T2);
      density_f<< "\t"<<u_density(eta, T1,0,T2,0);
      density_f<< "\t"<<d_density(eta, T1,0,T2,0);
      density_f<< "\t"<<s_density(eta, T1,T2)<< std::endl;
    }
  }

 	density_f.close();
}

void Charges::set_output_tests(){
  dump_charges(1,1);
  dump_charges(1,4);
  dump_charges(4,4);
  dump_charges(2,2);

  dump_charges_eta(0);
  dump_charges_eta(-7);
  dump_charges_eta(-3);
  dump_charges_eta(-1);
  dump_charges_eta(3);
  dump_charges_eta(7);
  dump_charges_eta(1);
  
  if(config.get_Verbose()){std::cout<< "[ Charges ] Charge test files written out" << std::endl;}

}
