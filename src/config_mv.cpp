/* Copyright (c) 2025 Oscar Garcia-Montero
 * All rights reserved. */
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <random>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>

#include  "include/config_mv.h"
#include  "include/params_gen.h"

MVConfig::MVConfig(){}
MVConfig::MVConfig(std::string PATH){
  path = PATH;
  process_parameters();
  
}
MVConfig::~MVConfig(){}

void MVConfig::get_name_and_value(std::string testline,std::string &name, std::string &value){
  std::string line_temp= testline;
  int pos = line_temp.find(":");
  name = line_temp.substr(0,pos);
  value = line_temp.substr(pos+1, line_temp.length()-pos );
  value = value.substr(0, value.find('#'));   // drop trailing "# comment"
  name.erase(std::remove(name.begin(),name.end(), ' '),name.end());
  value.erase(std::remove(value.begin(),value.end(), ' '),value.end());
}

//////////////////////////////
////////// CONFIG ////////////
//////////////////////////////

int MVConfig::count_indent(std::string testline){
  char cset[] = " ";
  return strspn (testline.c_str(),cset);
}

std::string MVConfig::get_header(std::string testline){
  std::string line_temp= testline;
  line_temp.erase(std::remove(line_temp.begin(), line_temp.end(), ' '), line_temp.end());
  line_temp.erase(std::remove(line_temp.begin(), line_temp.end(), ':'), line_temp.end());
  return line_temp;
}


void MVConfig::process_parameters(){
  std::string line;
  std::string header;
  std::ostringstream configpath;
  configpath << path<< "/config.yaml";

  std::ifstream configfile(configpath.str());
  if (configfile.is_open())
  {
    while ( getline (configfile,line) )
    {
      int indent_level = count_indent(line);
      if(line.length()>0){
        if(indent_level == 0){
          header=get_header(line);
        }
        else if(indent_level>0){
          // subheader=get_header(line);
          if(header=="General") {process_general_parameters(line);}
          else if(header=="Grid") {process_grid_parameters(line);}
          else if(header=="MV-Parameters") {process_MV_parameters(line);}
          else if(header=="BK-Parameters") {process_BK_parameters(line);}
          else if(header=="McDipper-Settings") {process_McDipper_settings(line);}
          else{std::cerr<< "Config Error! Indentation in configuration is not correct! "<<std::endl; exit(EXIT_FAILURE);}
        }
        else if(indent_level%4!=0){std::cerr<< "Config Error! Indentation in configuration is not correct! "<<std::endl; exit(EXIT_FAILURE);}
      }
      else{continue;}
    }
  }
  else{std::cerr<<"Config Error! Failed to open config-file!    Config path:   "<< configpath.str() <<std::endl; exit(EXIT_FAILURE); }
  // Extra processing
  xmin=x0*exp(-Ymax);
  make_momentum_parameters();
  configfile.close();
}


void MVConfig::process_general_parameters(std::string testline){
	std::string name_t,value_t;
	get_name_and_value(testline,name_t,value_t);
	if(name_t=="Nc"){NC=std::stoi(value_t);}
	if(name_t=="Nf"){NF=std::stoi(value_t);}
	if(name_t=="LambdaQCD"){LambdaQCD =std::stod(value_t);}
	if(name_t=="LambdaEff"){LambdaEff =std::stod(value_t);}
}

void MVConfig::process_MV_parameters(std::string testline){
	std::string name_t,value_t;
	get_name_and_value(testline,name_t,value_t);
	if(name_t=="IC"){IC=value_t;}
	if(name_t=="Qs02"){ Qs02=std::stod(value_t);}
	if(name_t=="gamma"){gamma=std::stod(value_t);}
	if(name_t=="sigma0"){sigma0=std::stod(value_t);}
	if(name_t=="ec"){ec =std::stod(value_t);}
	if(name_t=="Lambda"){LambdaMV =std::stod(value_t);}
	if(name_t=="kappa"){kappa=std::stod(value_t);}
	if(name_t=="m_IC"){m_IC=std::stod(value_t);}
}


void MVConfig::process_BK_parameters(std::string testline){
	std::string name_t,value_t;
	get_name_and_value(testline,name_t,value_t);
	if(name_t=="mu2"){mu2=std::stod(value_t);}
	if(name_t=="C2"){C2=std::stod(value_t);}
	if(name_t=="zeta"){zeta=std::stod(value_t);}
	if(name_t=="mp"){mp =std::stod(value_t);}
	if(name_t=="Lambda"){LambdaBK =std::stod(value_t);}
	if(name_t=="nGC"){nGC=std::stoi(value_t);}
	if(name_t=="Ktype"){Ktype =value_t;}
	if(name_t=="Adjoint"){Adjoint=value_t;}
	if(name_t=="Coupling"){Coupling=value_t;}
	if(name_t=="alpha_freeze"){alpha_freeze=std::stod(value_t);}
}

void MVConfig::process_McDipper_settings(std::string testline){
	std::string name_t,value_t;
	get_name_and_value(testline,name_t,value_t);
	if(name_t=="LargeX"){LargeXMode=std::stoi(value_t);}
	if(name_t=="Npoints"){NpointsLX=std::stoi(value_t);}
	if(name_t=="xmax"){xmaxLX=std::stod(value_t);}
	if(name_t=="HankelMode"){HankelMode=std::stoi(value_t);}
	if(name_t=="LargeXPars"){
		LargeXPars.clear();
		std::stringstream ss(value_t); std::string item;
		while(std::getline(ss,item,',')){ if(item.length()>0){LargeXPars.push_back(std::stod(item));} }
	}
}

void MVConfig::set_McDipper_settings(int LargeX_t, int Npoints_t, double xmax_t, const std::vector<double>& pars_t, int HankelMode_t){
	LargeXMode=LargeX_t; NpointsLX=Npoints_t; xmaxLX=xmax_t; LargeXPars=pars_t; HankelMode=HankelMode_t;
}

void MVConfig::process_grid_parameters(std::string testline){
  std::string name_t,value_t;
  get_name_and_value(testline,name_t,value_t);
  if(value_t.length()>0){
    //Process Parameters
    if(subheader=="r-grid"){
		if (name_t=="Nr"){Nr=std::stoi(value_t);}
		if (name_t=="rmin"){rmin=std::stod(value_t);}
		if (name_t=="rmax"){rmax=std::stod(value_t);}
    }
	if(subheader=="T-grid"){
		if (name_t=="NT"){NT=std::stoi(value_t);}
		if (name_t=="Tmin"){Tmin=std::stod(value_t);}
		if (name_t=="Tmax"){Tmax=std::stod(value_t);}
		if (name_t=="dT"){dT=std::stod(value_t);}
    }
	if(subheader=="Y-grid"){
		if (name_t=="NY"){NY=std::stoi(value_t);}
		if (name_t=="Ymin"){Ymin=std::stod(value_t);}
		if (name_t=="Ymax"){Ymax=std::stod(value_t);}
		if (name_t=="x0"){x0=std::stod(value_t);}
		if (name_t=="dY"){dY=std::stod(value_t);}
    }
  }
  else{subheader=name_t;}
}

void MVConfig::make_momentum_parameters(){ 
	qmin = 0;
	qmax =log(rmax/rmin);
	dq = log(rmax/rmin)/(Nr-1.);

	NK=Nr;
	kmin= 1./(rmax*gen_pars::fm_to_GeVm1);
	kmax= 0.25/(rmin*gen_pars::fm_to_GeVm1);
	umin =0;
	umax= log(kmax/kmin);
	du = (umax-umin)/(NK-1.);
} 

void MVConfig::terminal_setup_output(){
    std::cout<< "|--------------------------------- MV-like Dipole Parameters -------------------------------------|\n";
    std::cerr << "[General]:\n";
    std::cerr <<"    Nc = " << NC<< ", Nf = " << NF<< "\n";
    std::cerr <<"    LambdaQCD: " << LambdaQCD << " GeV,\tLambdaEff: " << LambdaEff << " GeV\n"; // GeV
    std::cerr << "[MV-Parameters]:\n";
    std::cerr << "    IC = " << IC << "\n";
    if(IC=="MV-Regulated"){
    std::cerr << "    kappa = " << kappa << ",\tm_IC = " << m_IC << " GeV\n";
    } else {
    std::cerr << "    Qs02= " << Qs02 << " GeV^2,\tgamma = " << gamma<< ",\tLambdaMV = " << LambdaMV<< " GeV\n"; 
    std::cerr << "    sigma0 = " << sigma0 << " fm^2 ,\tec = " << ec<< "\n";
    }
    std::cerr << "[BK-Parameters]:\n";
    std::cerr << "    LambdaBK = " << LambdaBK <<" GeV,\tmu2 = " << mu2<<" GeV^2,\tC2 = " << C2<<"\n";
    std::cerr << "    zeta = " << zeta<<",\tmp = " << mp <<" GeV,\tnGC = " << nGC << "\n";
    std::cerr << "    Ktype = " << Ktype << ",\tCoupling = " << Coupling;
    if(Coupling=="Frozen"){std::cerr << " (alpha_freeze = " << alpha_freeze << ")";}
    std::cerr << ",\tAdjoint = " << Adjoint << "\n";
    std::cerr << "[Grid]:\n";
    std::cerr << "r-grid:\n";
    std::cerr << "    Nr = " << Nr<<", rmin = " << rmin<<" fm, rmax = " << rmax<<" fm\n";
    std::cerr << "        \t qmin = " << qmin<<", qmax = " << qmax<<", dq =" << dq<<" \n";
    std::cerr << "Y-grid:\n";
    std::cerr << "    x0 = " << x0<<"\n";
    std::cerr << "    NY = " << NY<<", Ymin = " << Ymin<<", Ymax = " << Ymax<<" => dY = " << dY<<"\n";    
    std::cerr << "T-grid:\n";
    std::cerr << "    NT = " << NT<<", Tmin = " << Tmin<<" fm^{-2}, Tmax = " << Tmax<<" fm^{-2} => dT = " << dT<<" fm^{-2}\n"; 
    if(HankelMode>=0){
    std::cerr << "[McDipper-Settings]:\n";
    std::cerr << "    LargeX = " << LargeXMode << ",\tNpoints = " << NpointsLX << ",\txmax = " << xmaxLX << ",\tLargeXPars =";
    for(double v: LargeXPars){std::cerr << " " << v;}
    std::cerr << ",\tHankelMode = " << HankelMode << "\n";
    }
    std::cerr << "[k-Parameters]:\n"; 
    std::cerr << "    NK = " << NK<<", kmin = " << kmin<<" GeV, kmax = " << kmax<<" GeV\n";
    std::cerr << "       \t umin = " << umin<<", umax = " << umax<<", du = " << du<<" \n";
    std::cout<< "|-------------------------------------------------------------------------------------------------|\n";
} 

void MVConfig::write_config(std::string path_to_write){
  std::ofstream config_f;
  std::ostringstream config_f_name;
  config_f_name << path_to_write<< "/config.yaml";
  config_f.open (config_f_name.str());
    config_f << std::setprecision(12);
    config_f << "General:\n";
    config_f << "    Nc:        " << NC<< "\n";
    config_f << "    Nf:        " << NF<< "\n";
    config_f << "    LambdaQCD: " << LambdaQCD << " # GeV \n"; // GeV
    config_f << "    LambdaEff: " << LambdaEff << " # GeV \n";
    config_f << std::endl;
    config_f << "MV-Parameters:\n";
    config_f << "    IC:     " << IC << "\n";
    if(IC=="MV-Regulated"){
    config_f << "    kappa:  " << kappa << "\n";
    config_f << "    m_IC:   " << m_IC << " # GeV \n";
    } else {
    config_f << "    Qs02:   " << Qs02  << " # GeV^2\n"; // GeV
    config_f << "    gamma:  " << gamma<< "\n";
    config_f << "    Lambda: " << LambdaMV << " # GeV \n";// GeV
    config_f << "    sigma0: " << sigma0 << "\n";
    config_f << "    ec:     " << ec<< "\n";
    }
    config_f << std::endl;
    config_f << "BK-Parameters:\n";
    config_f << "    Lambda: " << LambdaBK <<" # GeV \n";
    config_f << "    mu2:    " << mu2 <<" # GeV^2\n";
    config_f << "    C2:     " << C2<<"\n";
    config_f << "    zeta:   " << zeta<<"\n";
    config_f << "    mp:     " << mp <<" # GeV \n";
    config_f << "    dY:     " << dY << "\n";
    config_f << "    nGC:    " << nGC << "\n";
    config_f << "    Ktype:  " << Ktype << "\n";
    config_f << "    Coupling: " << Coupling << "\n";
    if(Coupling=="Frozen"){config_f << "    alpha_freeze: " << alpha_freeze << "\n";}
    config_f << "    Adjoint: " << Adjoint << "\n";
    config_f << "Grid:\n";
    config_f << "    r-grid:\n";
    config_f << "        Nr:   " << Nr<<"\n";
    config_f << "        rmin: " << rmin<<" # in fm \n";
    config_f << "        rmax: " << rmax<<"\n";
    config_f << "        dr: Evolved in log-scale\n";
    config_f << "    Y-grid:\n";
    config_f << "        x0:   " << x0 <<"\n";
    config_f << "        NY:   " << NY<<"\n";
    config_f << "        Ymin: " << Ymin<<"\n";
    config_f << "        Ymax: " << Ymax<<"\n";
    config_f << "        dY:   " << dY<<"\n";
    config_f << "    T-grid:\n";
    config_f << "        NT:   " << NT<<"\n";
    config_f << "        Tmin: " << Tmin<<" # in fm^{-2}\n";
    config_f << "        Tmax: " << Tmax<<" # in fm^{-2}\n";
    config_f << "        dT:   " << dT<<"\n";
    config_f << "    k-grid:\n";
    config_f << "        Nr:   " << NK<<"\n";
    config_f << "        kmin: " << kmin<<" # in GeV \n";
    config_f << "        kmax: " << kmax<<"\n";
    config_f << "        dk: Transformed in log-scale\n";
    config_f << std::endl;
    config_f << "McDipper-Settings:\n";
    config_f << "    LargeX:     " << LargeXMode << " # 0 Frozen, 1 Simple, 2 TMD, 3 PDFMatch\n";
    config_f << "    Npoints:    " << NpointsLX << "\n";
    config_f << "    xmax:       " << xmaxLX << "\n";
    config_f << "    LargeXPars: ";
    for(size_t i=0;i<LargeXPars.size();i++){config_f << (i>0 ? ", " : "") << LargeXPars[i];}
    config_f << " # Simple: beta\n";
    config_f << "    HankelMode: " << HankelMode << " # 0 DHT, 1 J0, 2 J1\n";
    config_f.close(); 
}

bool MVConfig::compare(const MVConfig& other) const {
    bool equal = true;
    double relEPS=1e-3;
    std::vector<std::string> diffs;
    // General
    if (NC != other.NC) { diffs.emplace_back("NC"); equal = false; }
    if (NF != other.NF) { diffs.emplace_back("NF"); equal = false; }
    if (!reldiff(LambdaQCD, other.LambdaQCD,relEPS)) { diffs.emplace_back("LambdaQCD"); equal = false; }

    // MV parameters
    if (IC != other.IC) { diffs.emplace_back("IC"); equal = false; }
    if (!reldiff(kappa, other.kappa,relEPS)) { diffs.emplace_back("kappa"); equal = false; }
    if (!reldiff(m_IC, other.m_IC,relEPS)) { diffs.emplace_back("m_IC"); equal = false; }
    if (!reldiff(Qs02, other.Qs02,relEPS)) { diffs.emplace_back("Qs02"); equal = false; }
    if (!reldiff(gamma, other.gamma,relEPS)) { diffs.emplace_back("gamma"); equal = false; }
    if (!reldiff(LambdaMV, other.LambdaMV,relEPS)) { diffs.emplace_back("LambdaMV"); equal = false; }
    if (!reldiff(sigma0, other.sigma0,relEPS)) { diffs.emplace_back("sigma0"); equal = false; }
    if (!reldiff(ec, other.ec,relEPS)) { diffs.emplace_back("ec"); equal = false; }

    // BK parameters
    if (!reldiff(LambdaBK, other.LambdaBK,relEPS)) { diffs.emplace_back("LambdaBK"); equal = false; }
    if (!reldiff(mu2, other.mu2,relEPS)) { diffs.emplace_back("mu2"); equal = false; }
    if (!reldiff(C2, other.C2,relEPS)) { diffs.emplace_back("C2"); equal = false; }
    if (!reldiff(zeta, other.zeta,relEPS)) { diffs.emplace_back("zeta"); equal = false; }
    if (!reldiff(mp, other.mp,relEPS)) { diffs.emplace_back("mp"); equal = false; }
    if (!reldiff(nGC, other.nGC,relEPS)) { diffs.emplace_back("nGC"); equal = false; }
    if (Ktype != other.Ktype) { diffs.emplace_back("Ktype"); equal = false; }
    if (Adjoint != other.Adjoint) { diffs.emplace_back("Adjoint"); equal = false; }
    if (Coupling != other.Coupling) { diffs.emplace_back("Coupling"); equal = false; }
    if (!reldiff(alpha_freeze, other.alpha_freeze,relEPS)) { diffs.emplace_back("alpha_freeze"); equal = false; }

    // Grid
    if (Nr != other.Nr) { diffs.emplace_back("Nr"); equal = false; }
    if (NT != other.NT) { diffs.emplace_back("NT"); equal = false; }
    if (NY != other.NY) { diffs.emplace_back("NY"); equal = false; }
    if (NK != other.NK) { diffs.emplace_back("NK"); equal = false; }

    if (!reldiff(rmin, other.rmin,relEPS)) { diffs.emplace_back("rmin"); equal = false; }
    if (!reldiff(rmax, other.rmax,relEPS)) { diffs.emplace_back("rmax"); equal = false; }
    if (!reldiff(qmin, other.qmin,relEPS)) { diffs.emplace_back("qmin"); equal = false; }
    if (!reldiff(qmax, other.qmax,relEPS)) { diffs.emplace_back("qmax"); equal = false; }
    if (!reldiff(dq, other.dq,relEPS)) { diffs.emplace_back("dq"); equal = false; }

    if (!reldiff(x0, other.x0,relEPS)) { diffs.emplace_back("x0"); equal = false; }
    if (!reldiff(Ymin, other.Ymin,relEPS)) { diffs.emplace_back("Ymin"); equal = false; }
    if (!reldiff(Ymax, other.Ymax,relEPS)) { diffs.emplace_back("Ymax"); equal = false; }
    if (!reldiff(dY, other.dY,relEPS)) { diffs.emplace_back("dY"); equal = false; }

    if (!reldiff(Tmin, other.Tmin,relEPS)) { diffs.emplace_back("Tmin"); equal = false; }
    if (!reldiff(Tmax, other.Tmax,relEPS)) { diffs.emplace_back("Tmax"); equal = false; }
    if (!reldiff(dT, other.dT,relEPS)) { diffs.emplace_back("dT"); equal = false; }

    if (!reldiff(kmax, other.kmax,relEPS)) { diffs.emplace_back("kmax"); equal = false; }
    if (!reldiff(kmin, other.kmin,relEPS)) { diffs.emplace_back("kmin"); equal = false; }
    if (!reldiff(umax, other.umax,relEPS)) { diffs.emplace_back("umax"); equal = false; }
    if (!reldiff(umin, other.umin,relEPS)) { diffs.emplace_back("umin"); equal = false; }
    // if (!reldiff(du, other.du,relEPS)) { diffs.emplace_back("du"); equal = false; }

    // McDipper-Settings
    if (LargeXMode != other.LargeXMode) { diffs.emplace_back("LargeX"); equal = false; }
    if (NpointsLX != other.NpointsLX) { diffs.emplace_back("Npoints"); equal = false; }
    if (!reldiff(xmaxLX, other.xmaxLX,relEPS)) { diffs.emplace_back("xmax"); equal = false; }
    if (HankelMode != other.HankelMode) { diffs.emplace_back("HankelMode"); equal = false; }
    if (LargeXPars.size() != other.LargeXPars.size()) { diffs.emplace_back("LargeXPars"); equal = false; }
    else{
      for (size_t i = 0; i < LargeXPars.size(); i++) {
        if (!reldiff(LargeXPars[i], other.LargeXPars[i],relEPS)) { diffs.emplace_back("LargeXPars"); equal = false; break; }
      }
    }

   if (!equal) {std::cerr<< "Different variables:\n";
    for (const auto& name : diffs) {
        std::cerr << " - " << name ;
    }
    std::cerr<< std::endl;
    }
    return equal;
}





