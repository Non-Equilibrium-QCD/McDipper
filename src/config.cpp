/* Copyright (c) 2022 Oscar Garcia-Montero
 * All rights reserved. */
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <time.h>       /* time */
#include <random>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>


#include "include/config.h"

Config::Config(){}
Config::Config(std::string configfile_path){
  path_to_configfile = configfile_path;
  NModelParams=0;
  process_parameters();
  set_seed();
  if(Verbose){terminal_setup_output();}
}

Config::~Config(){}

void Config::process_parameters(){
  std::string line;
  std::string header;

  std::ifstream configfile(path_to_configfile);
  if (configfile.is_open())
  {
    while ( getline (configfile,line) )
    {
      int indent_level = count_indent(line);
      if(line.length()>0){
        if(indent_level == 0){
          if (line[0]=='V'){check_version(line);}
          else{header=get_header(line);}
        }
        else if(indent_level>0){
          // subheader=get_header(line);
          if(header=="Logging") {process_logging_parameters(line);}
          else if(header=="General") {
            process_general_parameters(line);}
          else if(header=="Grid") {process_grid_parameters(line);}
          else if(header=="Model_Parameters") {process_model_parameters(line);}
          else if(header=="Output") {process_output_parameters(line);}
          else if(header=="Thickness") {process_thickness_parameters(line);}
          else{std::cerr<< "Config Error! Indentation in configuration is not correct! "<<std::endl; exit(EXIT_FAILURE);}
        }
        else if(indent_level%4!=0){std::cerr<< "Config Error! Indentation in configuration is not correct! "<<std::endl; exit(EXIT_FAILURE);}
      }
      else{continue;}
    }
  }
  else{std::cerr<<"Config Error! Failed to open config-file!    Config path:   "<< path_to_configfile <<std::endl; exit(EXIT_FAILURE); }
  // Extra processing
  fill_up_grid_params();
  configfile.close();
}

///////////////////// Version Check  ///////////////////////

void Config::check_version(std::string testline){
  std::string line_temp= testline;
  int pos = line_temp.find(":");
  version = line_temp.substr(pos+1, line_temp.length()-pos );
  version.erase(std::remove(version.begin(), version.end(), ' '), version.end());
  /// Later on, checks version and if not compatible, exits. For now this is just a placeholder.
}

////////////////////////  Logging  ///////////////////////
////////////////////////  Logging  ///////////////////////

void Config::process_logging_parameters(std::string testline){
  std::string subheader_t,value_t;
  get_name_and_value(testline,subheader_t, value_t);
  Verbose = make_bool(value_t);
}
////////////////////////  General  ///////////////////////
////////////////////////  General  ///////////////////////
void Config::process_general_parameters(std::string testline){
  std::string name_t,value_t;
  get_name_and_value(testline,name_t,value_t);
  if(value_t.length()>0){
    //Process Parameters
    if(name_t=="SqrtsNN"){sqrtsNN=std::stod(value_t);}
    if(name_t=="Model"){
      cModel_int=std::stoi(value_t);
      switch(cModel_int){
        case 0: cModel=Model::GBW;cModelStr="GBW";break;
        case 1: cModel=Model::IPSat;cModelStr="IP-Sat";break;
        case 2: cModel=Model::MV;cModelStr="MV";break;
        case 3: cModel=Model::GBWSimp;cModelStr="GBW-Simplified";break;
      }
    }
    if(subheader=="Fluctuations"){
      if (name_t=="Thickness_fluct"){thick_fluct=make_bool(value_t);}
      if (name_t=="Fluct_mode"){fluct_mode=value_t;}
      if (name_t=="Sigma"){sigma=std::stod(value_t);}
      if (name_t=="Hotspots_fluct"){hotspots_fluct=make_bool(value_t);}
      if (name_t=="Bq"){Bq=std::stod(value_t);}
      if (name_t=="Nq"){Nq=std::stoi(value_t);}
    }
    if(name_t=="A"){
      if(subheader=="Nucleus1"){A1=std::stoi(value_t);}
      if(subheader=="Nucleus2"){A2=std::stoi(value_t);}
    }
    if(name_t=="Z"){
      if(subheader=="Nucleus1"){Z1=std::stoi(value_t);}
      if(subheader=="Nucleus2"){Z2=std::stoi(value_t);}
    }
    if(name_t=="mode"){
      if(subheader=="Nucleus1"){
        mode1=std::stoi(value_t);
        if(mode1==0){mode1name="Spherical";}
        else if(mode1==1){mode1name="Deformed";}
        else if(mode1==2){mode1name="Non-Spherical";}
        else if(mode1==3){mode1name="Input-Sampling";}
      }
      if(subheader=="Nucleus2"){
        mode2=std::stoi(value_t);
        if(mode2==0){mode2name="Spherical";}
        else if(mode2==1){mode2name="Deformed";}
        else if(mode2==2){mode2name="Non-Spherical";}
        else if(mode2==3){mode2name="Input-Sampling";}
        }
    }
    if(name_t== "InputFile"){
      if(subheader=="Nucleus1"){Nucleus1ListFile=value_t;}
      if(subheader=="Nucleus2"){Nucleus2ListFile=value_t;}
    }
    if(name_t== "IsospinSpecified"){
      if(subheader=="Nucleus1"){N1IsospinSpec=make_bool(value_t);}
      if(subheader=="Nucleus2"){N2IsospinSpec=make_bool(value_t);}
    }
    if(name_t== "Configurations"){
      if(subheader=="Nucleus1"){NConf1=std::stoi(value_t);}
      if(subheader=="Nucleus2"){NConf2=std::stoi(value_t);}
    }
    if(name_t== "Weights"){
      if(subheader=="Nucleus1"){N1_weights=make_bool(value_t);}
      if(subheader=="Nucleus2"){N2_weights=make_bool(value_t);}
      
    }

    if(name_t=="GlauberAcceptance"){
      if(value_t=="Standard"){GMode=GlauberMode::Standard;GModeStr="Standard";}
      else if(value_t=="Gaussian"){GMode=GlauberMode::Gaussian;GModeStr="Gaussian";}
      else if(value_t=="Exponential"){GMode=GlauberMode::Exponential;GModeStr="Exponential";}
      else {std::cerr<< "Config Error! Glauber Acceptance mode is not valid! "<<std::endl; exit(EXIT_FAILURE);}
    }
    if(subheader=="Impact" && name_t=="Value"){ ImpactValue = std::stod(value_t); ImpactMode=ImpSample::Fixed;}
    if(subheader=="Impact" && name_t=="Range"){retrieve_range(value_t,bMin, bMax);}
    if(subheader=="Impact" && name_t=="Sampling"){
      if(value_t=="Quadratic"){ImpactMode=ImpSample::bdbSampled;}
      else if(value_t=="Uniform"){ImpactMode=ImpSample::dbSampled;}
      else{std::cerr<< "Config Error! Impact sampling mode is not valid! "<<std::endl; exit(EXIT_FAILURE);}
    }
    if(name_t=="Seed"){seed=std::stoi(value_t);}
    if(subheader=="PDFs" && name_t=="PDFSet"){cPDFSetStr=value_t;}
    if(subheader=="PDFs" && name_t=="ForcePositive"){cForcedMode=std::stoi(value_t);}
    if(name_t=="K-Factor"){KFactor=std::stod(value_t);}
    if(name_t=="Events"){NEvents=std::stoi(value_t);}
  }
  else{subheader=name_t;}
}

////////////////////////  Grid  ///////////////////////
////////////////////////  Grid  ///////////////////////

void Config::process_grid_parameters(std::string testline){
  std::string name_t,value_t;
  get_name_and_value(testline,name_t,value_t);
  if(name_t=="NX"){NX=std::stoi(value_t);}
  if(name_t=="NY"){NY=std::stoi(value_t);}
  if(name_t=="NETA"){NETA=std::stoi(value_t);}
  if(name_t=="X_RANGE"){retrieve_range(value_t,XMIN,XMAX);}
  if(name_t=="Y_RANGE"){retrieve_range(value_t,YMIN,YMAX);}
  if(name_t=="ETA_RANGE"){retrieve_range(value_t,ETAMIN,ETAMAX);}
  if(name_t=="BG"){BG=std::stod(value_t);} // In fm2
}

void Config::fill_up_grid_params(){
  dX=(XMAX-XMIN)/(NX-1.);
  dY=(YMAX-YMIN)/(NY-1.);
  dETA=(ETAMAX-ETAMIN)/(NETA-1.);
};

//////////////////  Model Parameters  //////////////////
//////////////////  Model Parameters  //////////////////
void Config::process_model_parameters(std::string testline){
  if(cModel==Model::GBW || cModel==Model::GBWSimp ){process_gbw_parameters(testline);}
  if(cModel==Model::IPSat){process_ipsat_parameters(testline);}
  if(cModel==Model::MV){process_mv_parameters(testline);}
}
void Config::process_gbw_parameters(std::string testline){
  std::string name_t,value_t;
  get_name_and_value(testline,name_t,value_t);
  if(name_t=="Q02"){ModelPars[0]=std::stod(value_t);NModelParams++;}
  if(name_t=="x0"){ModelPars[1]=std::stod(value_t);NModelParams++;}
  if(name_t=="lambda"){ModelPars[2]=std::stod(value_t);NModelParams++;}
  if(name_t=="XCut"){ModelPars[3]=std::stod(value_t);NModelParams++;}
}

void Config::process_ipsat_parameters(std::string testline){
  std::string name_t,value_t;
  get_name_and_value(testline,name_t,value_t);
  if(name_t=="Set"){ModelPars[0]=std::stod(value_t);NModelParams++;}
  if(name_t=="XScaling"){
    ModelPars[1]=std::stod(value_t);
    // IPsat_pars::x0_scaling =ModelPars[1];  
    NModelParams++;}
  if(name_t=="P_reg"){ModelPars[2]=std::stod(value_t);NModelParams++;}
}
void Config::process_mv_parameters(std::string testline){}


///////////////////// Thickness //////////////////////
///////////////////// Thickness  /////////////////////

void Config::process_thickness_parameters(std::string testline){
  std::string name_t,value_t;
  get_name_and_value(testline,name_t,value_t);
  if(name_t=="TMax"){TMax=std::stod(value_t);}
  if(name_t=="TMin"){TMin=std::stod(value_t);}
  if(name_t=="NT"){NT=std::stoi(value_t);}
}


/////////////////////  Output  //////////////////////
/////////////////////  Output   /////////////////////
void Config::process_output_parameters(std::string testline){
  std::string name_t, value_t,token;
  std::string delimiter= ",";
  get_name_and_value(testline,name_t,value_t);
  if(name_t == "path_to_output"){
    if (value_t.length() > 0)
    {
      std::string::iterator it = value_t.end() - 1;
      if (*it == '/'){value_t.erase(it);} 
    }
    path_to_output= value_t;
  }
  if(name_t == "run_name"){run_name=value_t; is_run_renamed=true;}
  if(name_t == "Format"){
    n_formats = 1 + std::count(value_t.begin(), value_t.end(), ',');
    format=new std::string[n_formats];
    remove_char(value_t,'[');remove_char(value_t,']');
    int pos = 0;
    for (int i = 0; i < n_formats-1; i++) {
      pos = value_t.find(delimiter);
      token = value_t.substr(0, pos);
      value_t.erase(0, pos + delimiter.length());
      remove_char(token,'"');
      format[i]=token;
    }
    remove_char(value_t,'"');
    format[n_formats-1] = value_t;
  }
  if(name_t == "PrintAvg"){
    if(value_t == "False"){print_avg=0;}
    else if(value_t == "ObservablesOnly"){print_avg=1;}
    else if(value_t == "True"){print_avg=2;}
    else{std::cerr<<"Error: Average-Event output mode not implemented! Exiting.";exit(EXIT_FAILURE);}
  }
  if(name_t == "BoostInvariant"){ boost_invariant= make_bool(value_t);}
  if(name_t == "suppress_output"){ suppress_output= make_bool(value_t);}
  // 
 
}
 
void Config::write_config_file(std::string target_folder){}

void Config::terminal_setup_output(){
  std::cout<< "|----------------------------------------------------------------------------------------|"<<std::endl;
  std::cout<< "|----------------------------------------------------------------------------------------|"<<std::endl;
  std::cout<< "                             Running McDipper. Version " << version <<std::endl;
  std::cout<< "|----------------------------------------------------------------------------------------|"<<std::endl;
  std::cout<< "|------------------------------- General Parameters -------------------------------------|\n";
  std::cout<< "       Running model: " << cModelStr << " for sqrt(S)= "<<sqrtsNN << "\n";
  if(A1>3 && A2>3){std::cout<< "       Nucleus 1  ->  A="<<A1<< " Z="<< Z1 <<" ("<< mode1name  << ")        Nucleus 2  ->  A="<<A2<< " Z=" << Z2 <<" ("<< mode2name  << ")\n";}
  if(A1<=3 && A2>3){std::cout<< "       Nucleus 1  ->  A="<<A1<< " Z="<< Z1 <<"        Nucleus 2  ->  A="<<A2<< " Z=" << Z2 <<" ("<< mode2name  << ")\n";}
  if(A1>3 && A2<=3){std::cout<< "       Nucleus 1  ->  A="<<A1<< " Z="<< Z1<<" ("<< mode1name  << ")        Nucleus 2  ->  A="<<A2<< " Z=" << Z2 <<" \n";}
  if(A1<=3 && A2<=3){std::cout<< "       Nucleus 1  ->  A="<<A1<< " Z="<< Z1<<"           Nucleus 2  ->  A="<<A2<< " Z=" << Z2 <<"  \n";}
  std::cout<< "       Nucleon-Nucleon selection mode : "  << GModeStr << "\n";

  if(ImpactMode==ImpSample::Fixed){std::cout<< "       Impact Parameter Mode: 'Value'  for  b = "<<ImpactValue <<" fm  \n";}
  if(ImpactMode==ImpSample::bdbSampled){std::cout<< "       Impact Parameter Mode: 'Range'  for  b = [ "<<bMin<<" , "<< bMax <<" ] fm, Sampling: Quadratic  \n";}
  if(ImpactMode==ImpSample::dbSampled){std::cout<< "       Impact Parameter Mode: 'Range'  for  b = [ "<<bMin<<" , "<< bMax <<" ] fm, Sampling: Uniform \n";}
  std::cout<< "       Number of Events: "<< NEvents <<" , using a K-factor = "<< KFactor <<" for the gluon energy\n";
  std::cout<< "       Running seed: "<< seed <<"\n";
  if (thick_fluct){
    std::cout<< "               Thickness fluctuation: True\n";
    std::cout<< "               Thickness fluctuation distribution: " << fluct_mode << "\n";
    std::cout<< "               Thickness fluctuation width parameter: " << sigma << "\n"; 
  }
  std::cout<< "|--------------------------------- Grid Parameters --------------------------------------|\n";
  std::cout<< "                       X=["<<XMIN<<","<< XMAX<< "] ,   NX = "<<NX<<" ,   dX = "<< dX << " fm \n";
  std::cout<< "                       Y=["<<YMIN<<","<< YMAX<< "] ,   NY = "<<NY<<" ,   dY = "<< dY << " fm \n";
  std::cout<< "                     ETA=["<<ETAMIN<<","<< ETAMAX<< "] , NETA = "<<NETA<<" , dETA = "<< dETA << "\n";
  std::cout<< "                     Nucleonic smearing parameter B_G = "<< BG <<"fm^2 \n";
  if (hotspots_fluct){
  std::cout<< "                 Hotspots  fluctuation: True\n";
  std::cout<< "                 Hotspots  smearing parameter Bq = "<< Bq << "fm^2 \n";
  std::cout<< "                 Hotspots  number = " << Nq << "\n";
  }
  std::cout<< "|--------------------------------- Model Parameters -------------------------------------|\n";
  if(cModel==Model::GBW){
  std::cout<< "               Q02="<<ModelPars[0]<<" GeV^2,    x0="<<ModelPars[1]<<",   lambda="<<ModelPars[2]<<",   XCut="<<ModelPars[3]<<" \n";
  }
  if(cModel==Model::IPSat){
  std::cout<< "               ParameterSet -> "<<ModelPars[0]<<",    x0_scaling="<<ModelPars[1]<<"\n";
  }
  std::cout<< "|-------------------------------------- Output ------------------------------------------|\n";
  std::cout<< "  Writing output to : " << path_to_output << " \n";
  std::cout<< "  Boost invariant output : " ;
  if(boost_invariant){ std::cout<<"True\n";}
  else{ std::cout<<"False\n";}
  std::cout<< "  Chosen output Formats : ";
  for (int i = 0; i < n_formats; i++) { if(i==n_formats-1){std::cout<< format[i] << "\n";}else{std::cout<< format[i] << ", ";}}
  std::cout<< "  Average Event output: ";
  if(print_avg==0){std::cout<< "None \n";}
  else if(print_avg==1){std::cout<< "Only Observables \n";}
  else if(print_avg==2){std::cout<< "All \n";}
  else{std::cerr<<"Error: Average-Event output mode not implemented! Exiting.";exit(EXIT_FAILURE);}
  std::cout<< "|----------------------------------------------------------------------------------------|\n";
}

void Config::dump(std::string OUTPATH){
  std::ofstream config_f;
  std::ostringstream configname;
  configname << OUTPATH << "/config.yaml";
  config_f.open(configname.str());

  config_f << "Version:"<<version<<"\n";
  config_f << "\n";
  config_f << "Logging:\n";
  config_f << "    Verbose: "<< bool_string(Verbose)<<"\n";
  config_f << "\n";
  config_f << "General:\n";
  config_f << "    SqrtsNN: "<< sqrtsNN<<"\n";
  config_f << "    Nucleus1:\n";
  config_f << "        A: "<< A1<<"\n";
  config_f << "        Z: "<< Z1<<"\n";
  if(A1>3){
    config_f << "        mode: "<< mode1<<"\n";
    if(mode1==3){
      config_f << "        InputFile: "<< Nucleus1ListFile<<"\n";
      config_f << "        IsospinSpecified: "<< N1IsospinSpec<<"\n";
      config_f << "        Configurations: "<< NConf1<<"\n";
    } 
  }
  config_f << "    Nucleus2:\n";
  config_f << "        A: "<< A2<<"\n";
  config_f << "        Z: "<< Z2<<"\n";
  if(A2>3){
    config_f << "        mode: "<< mode2<<"\n";
    if(mode2==3){
      config_f << "        InputFile: "<< Nucleus2ListFile<<"\n";
      config_f << "        IsospinSpecified: "<< N2IsospinSpec<<"\n";
      config_f << "        Configurations: "<< NConf2<<"\n";
    }
  }

  config_f << " GlauberAcceptance: " << GModeStr << "\n";
  config_f << "    Events: "<< NEvents << "\n";
  if(cModel==Model::GBW){config_f << "    Model: 0\n";}
  if(cModel==Model::IPSat){config_f << "    Model: 1\n";}
  if(cModel==Model::MV){config_f << "    Model: 2\n";}
  if(cModel==Model::GBWSimp){config_f << "    Model: 3\n";}
  config_f << "    Impact:\n";
  if(ImpactMode==ImpSample::Fixed){config_f << "        Value: "<< ImpactValue<<"\n";}
  if(ImpactMode!=ImpSample::Fixed){
    config_f << "        Range:  ["<< bMin<<","<< bMax<<"]\n";
    if(ImpactMode==ImpSample::bdbSampled){config_f << "        Sampling:  Quadratic\n";}
    if(ImpactMode==ImpSample::dbSampled){config_f << "        Sampling:  Uniform\n";}
  }
  config_f << "    Seed: "<< seed<<"\n";
  config_f << "    PDFs: "<<"\n";
  config_f << "        PDFSet: "<< cPDFSetStr<<"\n";
  config_f << "        ForcePositive: "<< cForcedMode<<"\n";
  config_f << "    K-Factor: "<< KFactor << "\n";
  config_f << "    Fluctuations:\n";
  config_f << "        Thickness_fluct: "<< thick_fluct << "\n";
  config_f << "        Fluct_mode: "     << fluct_mode << "\n";
  config_f << "        Sigma: "          << sigma << "\n";
  config_f << "        Hotspots_fluct: " << hotspots_fluct << "\n";
  config_f << "        Bq: "             << Bq << "\n";
  config_f << "        Nq: "             << Nq << "\n";
  config_f << "Grid:\n";
  config_f << "    NX: "<< NX<<"\n";
  config_f << "    NY: "<< NY<<"\n";
  config_f << "    NETA: "<< NETA<<"\n";
  config_f << "    X_RANGE: ["<< XMIN<<","<< XMAX <<"]\n";
  config_f << "    Y_RANGE: ["<< YMIN<<","<< YMAX <<"]\n";
  config_f << "    ETA_RANGE: ["<< ETAMIN<<","<< ETAMAX <<"]\n";
  config_f << "    BG: "<< BG<<"\n";
  config_f << "Model_Parameters:\n";
  if(cModel==Model::GBW){
    config_f << "    Q02: "<<ModelPars[0]<<  "\n";
    config_f << "    x0: "<<ModelPars[1]<<  "\n";
    config_f << "    lambda: "<<ModelPars[2]<<  "\n";
    config_f << "    XCut: "<<ModelPars[3]<<  "\n";
  }
  if(cModel==Model::IPSat){
    config_f << "    Set: "<<ModelPars[0]<<  "\n";
    config_f << "    XScaling: "<<ModelPars[1]<<  "\n";
    config_f << "    P_reg: "<<ModelPars[2]<<  "\n";

  }
  config_f << "\n";
  config_f << "Output:\n";
  config_f << "    path_to_output: "<< path_to_output<<"\n";
  config_f << "    Format:  [";
  for (int i = 0; i < n_formats; i++) { if(i==n_formats-1){config_f << "\"" << format[i] << "\"]\n";}else{config_f << "\"" << format[i] << "\", ";}}
  config_f << "PrintAvg: ";
  if(print_avg==0){config_f<< "False\n";}
  else if(print_avg==1){config_f<< "ObservablesOnly\n";}
  else if(print_avg==2){config_f<< "True\n";}
  config_f<< "  BoostInvariant: " ;
  if(boost_invariant){ config_f<<"True\n";}
  else{ config_f<<"False\n";}
  config_f << "\n";
config_f.close();
}

void Config::set_dump(std::string OUTPATH){
  // This writes out the configuration needed to compute the CHARGE TABLES
  std::ofstream config_f;
  std::ostringstream configname;
  configname << OUTPATH << "/config.yaml";
  config_f.open(configname.str());

  config_f << "Version: "<<version<<"\n";
  config_f << "\n";
  config_f << "General:\n";
  config_f << "    SqrtsNN: "<< sqrtsNN<<"\n";
  if(cModel==Model::GBW){config_f << "    Model: 0\n";}
  if(cModel==Model::IPSat){config_f << "    Model: 1\n";}
  if(cModel==Model::MV){config_f << "    Model: 2\n";}
  if(cModel==Model::GBWSimp){config_f << "    Model: 3\n";}
  config_f << "    PDFs: "<<"\n";
  config_f << "        PDFSet: "<< cPDFSetStr<<"\n";
  config_f << "        ForcePositive: "<< cForcedMode<<"\n";
  config_f << "\n";
  config_f << "Grid:\n";
  config_f << "    NX: "<< NX<<"\n";
  config_f << "    NY: "<< NY<<"\n";
  config_f << "    NETA: "<< NETA<<"\n";
  config_f << "    X_RANGE: ["<< XMIN<<","<< XMAX <<"]\n";
  config_f << "    Y_RANGE: ["<< YMIN<<","<< YMAX <<"]\n";
  config_f << "    ETA_RANGE: ["<< ETAMIN<<","<< ETAMAX <<"]\n";
  config_f << "    BG: "<< BG<<"\n";
  config_f << "\n";
  config_f << "Model_Parameters:\n";
  if(cModel==Model::GBW){
    config_f << "    Q02: "<<ModelPars[0]<<  "\n";
    config_f << "    x0: "<<ModelPars[1]<<  "\n";
    config_f << "    lambda: "<<ModelPars[2]<<  "\n";
    config_f << "    XCut: "<<ModelPars[3]<<  "\n";
  }
  else if(cModel==Model::IPSat){
    config_f << "    Set: "<<ModelPars[0]<<  "\n";
    config_f << "    XScaling: "<<ModelPars[1]<<  "\n";
    config_f << "    P_reg: "<<ModelPars[2]<<  "\n";
  }
  config_f << "\n";
  config_f << "Thickness:\n";
  config_f << "    TMax: "<< TMax<<"\n";
  config_f << "    TMin: "<< TMin<<"\n";
  config_f << "    NT: "<< NT<<"\n";
  config_f << "\n";
config_f.close();
}

/////////////////////   Tools   /////////////////////
/////////////////////   Tools   /////////////////////
/////////////////////   Tools   /////////////////////

int Config::count_indent(std::string testline){
  char cset[] = " ";
  return strspn (testline.c_str(),cset);
}

std::string Config::get_header(std::string testline){
  std::string line_temp= testline;
  line_temp.erase(std::remove(line_temp.begin(), line_temp.end(), ' '), line_temp.end());
  line_temp.erase(std::remove(line_temp.begin(), line_temp.end(), ':'), line_temp.end());
  return line_temp;
}

void Config::get_name_and_value(std::string testline,std::string &name, std::string &value){
  std::string line_temp= testline;
  int pos = line_temp.find(":");
  name = line_temp.substr(0,pos);
  value = line_temp.substr(pos+1, line_temp.length()-pos );
  name.erase(std::remove(name.begin(),name.end(), ' '),name.end());
  value.erase(std::remove(value.begin(),value.end(), ' '),value.end());
}

bool Config::make_bool(std::string boolstr){
  bool testbool=false;
  if(boolstr=="True" || boolstr=="true"){testbool=true;}
  else if(boolstr=="False" || boolstr=="false"){testbool=false;}
  else{std::cerr<< "Config Error! Assignment of boolean variable is not correct ! "<<std::endl; exit(EXIT_FAILURE);}
  return testbool;
}

void Config::retrieve_range(std::string rangestring, double &xmin, double &xmax){
    std::string tempstr=rangestring;
    remove_char(tempstr,'[');remove_char(tempstr,']');
    int pos = tempstr.find(",");
    std::string min = tempstr.substr(0,pos);
    std::string max = tempstr.substr(pos+1, tempstr.length()-pos );
    min.erase(std::remove(min.begin(),min.end(), ' '),min.end());
    max.erase(std::remove(max.begin(),max.end(), ' '),max.end());
    xmin=std::stod(min);xmax=std::stod(max);
}


void Config::remove_char(std::string &str_rem,char rem_char){
  str_rem.erase(std::remove(str_rem.begin(),str_rem.end(), rem_char),str_rem.end());
}

void Config::set_seed(){
    if(seed < 0){
      int time_seed = time(NULL);
      std::srand(time_seed); seed=std::rand();
    }
  }

Config::Config(const Config& OldConf){
   Verbose=OldConf.Verbose;
   path_to_configfile=OldConf.path_to_configfile;
   //Version
   version=OldConf.version;
   //General
   cModel=OldConf.cModel;
   cModelStr=OldConf.cModelStr;
   A1=OldConf.A1;A2=OldConf.A2;
   Z1=OldConf.Z1;Z2=OldConf.Z2;
   mode1=OldConf.mode1;mode2=OldConf.mode2;
   mode1name=OldConf.mode1name;mode2name=OldConf.mode2name;
   Nucleus1ListFile=OldConf.Nucleus1ListFile; 
   Nucleus2ListFile=OldConf.Nucleus2ListFile;
   N1IsospinSpec=OldConf.N1IsospinSpec;
   N2IsospinSpec=OldConf.N2IsospinSpec;
   NConf1=OldConf.NConf1;
   NConf2=OldConf.NConf2;
   N1_weights=OldConf.N1_weights;
   N2_weights=OldConf.N2_weights;

   GMode=OldConf.GMode;
   GModeStr=OldConf.GModeStr;
  
   sqrtsNN=OldConf.sqrtsNN;
   ImpactMode=OldConf.ImpactMode;
   ImpactValue=OldConf.ImpactValue;
   bMin=OldConf.bMin;bMax=OldConf.bMax;
   seed=OldConf.seed;
   cPDFSetStr = OldConf.cPDFSetStr;
   cForcedMode = OldConf.cForcedMode;
   KFactor = OldConf.KFactor;
   thick_fluct=OldConf.thick_fluct;
   fluct_mode =OldConf.fluct_mode;
   sigma      =OldConf.sigma;
   hotspots_fluct=OldConf.hotspots_fluct;
   Bq=OldConf.Bq;
   Nq=OldConf.Nq;
   NEvents = OldConf.NEvents;
   //Grid
   NX=OldConf.NX;NY=OldConf.NY;NETA=OldConf.NETA;
   XMAX=OldConf.XMAX;XMIN=OldConf.XMIN;
   YMAX=OldConf.YMAX;YMIN=OldConf.YMIN;
   ETAMAX=OldConf.ETAMAX;ETAMIN=OldConf.ETAMIN;
   dX=OldConf.dX;dY=OldConf.dY;dETA=OldConf.dETA; /* fm */
   BG = OldConf.BG; /// In fm^2
   // Model Parameters
   for (int i = 0; i < 10; i++) {ModelPars[i]=OldConf.ModelPars[i];}
   NModelParams=OldConf.NModelParams;
   //OutPut params
   n_formats= OldConf.n_formats;
   path_to_output=OldConf.path_to_output;
   run_name=OldConf.run_name;
   is_run_renamed=OldConf.is_run_renamed;
   format=new std::string[n_formats];
   print_avg=OldConf.print_avg;
   boost_invariant=OldConf.boost_invariant;
   suppress_output=OldConf.suppress_output;
   for (int i = 0; i < n_formats; i++) {format[i]=OldConf.format[i];}
   //Thickness_Parameters
   TMax= OldConf.TMax;
   TMin= OldConf.TMin;
   NT= OldConf.NT;
}

bool Config::compare_grid_parameters (Config *OldConf, double tolerance){
  bool is_equal= ( NETA==OldConf->get_NETA()) ;
  is_equal= is_equal && ( NDiff(ETAMIN,OldConf->get_ETAMIN())<tolerance ) ;
  is_equal= is_equal && ( NDiff(ETAMAX,OldConf->get_ETAMAX())<tolerance ) ;
  return is_equal;
}

bool Config::compare_Thickness_parameters (Config *OldConf, double tolerance){
    bool is_equal= ( NT==OldConf->get_NT());
    is_equal= is_equal && ( NDiff(TMin,OldConf->get_TMin())<tolerance ) ;
    is_equal= is_equal && ( NDiff(TMax,OldConf->get_TMax())<tolerance ) ;
    return is_equal;
}

bool Config::compare_PDF_parameters (Config *OldConf, double tolerance){
  bool is_equal= ( cPDFSetStr==OldConf->get_PDFSet()) ;
  is_equal= is_equal && ( cForcedMode==OldConf->get_ForcedPositive() ) ;
  return is_equal;
}

bool Config::compare_model_parameters (Config * OldConf, double tolerance){
  bool is_equal=false;
  if(cModel == OldConf->get_model()){
    is_equal=true;
    for (int i = 0; i < get_mpars_number(); i++) {
      if(ModelPars[i]==0 && OldConf->get_ModelParams(i) == 0){continue;}
      else{
        double diff = NDiff(ModelPars[i],OldConf->get_ModelParams(i));
        is_equal=is_equal && (diff<tolerance );
      }
    }

  }
  return is_equal;
}
