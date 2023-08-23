
#ifndef NUCLEAR_DATA_H
#define NUCLEAR_DATA_H

namespace  NuclearData { ///TO DO, get more elements!

    // Nuclear Parameters from 2011.14909
    // Format:  R(fm) a(fm) beta2 beta4
    // ( multipole parameters are dimensionless)
    // (deviation from sphericity still to be done)
    // double Silicon[5]={3.34, 0.580,-0.233,0.,0.};

    //Spherical Nuclei (denoted by sph subindex)
    double O_sph[2]={2.608, 0.5130}; // A=16, Z=8
    double Ar_sph[2]={3.53, 0.542};  // A=40
    double Cu_sph[2]={4.2, 0.5130};  // A=63
    double Xe_sph[2]={5.42, 0.57};   // A=129
    double Au_sph[2]={6.38, 0.535};  // A=197
    double Pb_sph[2]={6.62, 0.546};  // A=208

    // Spherically Deviated Nuclei (denoted by dev subindex) //TODO
    // Needs new parametrization

    // Deformed Nuclei (denoted by def subindex)
    double Al_def[4]={5.36, 0.59,-0.448, -0.239};   // A=27
    double Cu_def[4]={4.2, 0.596, 0.162, -0.006};   // A=63
    double Xe_def[4]={5.36, 0.59,0.161,-0.003};   // A=129  parameters from  Z. Physik 270 (1974) 113.
    double Au_def[4]={6.38, 0.535,-0.131,-0.031}; // A=197
    double U_def[4]={6.670, 0.44,0.280,0.093};    // A=238
    // double ATest[4]={6.670,0.44,-1.,0.093}; // 512 (Fake Test Nucleus)

    void getNPars(int A,int Z,  int mode,double * pars){
      // A is the nuclear number, mode corresponds to the type of Nucleus
      // mode 0 -> spherical, mode 1 -> deformed, mode 2-> sphericity dev. (NOT IMPLEMENTED!)
      switch (Z) {
        case 8: // Oxygen 
          if(A==16){
            if(mode==0){for (size_t i = 0; i < 2; i++) { pars[i]=O_sph[i];} }
            else{std::cerr<< "Error! Non-spherical parametrization not implemented for O-16!" ;exit(EXIT_FAILURE);}
            }
          break;
        case 13:
          if(A==27){
            if(mode==1){for (size_t i = 0; i < 4; i++) { pars[i]=Al_def[i];} }
            else{std::cerr<< "Error! Non-deformed parametrization not implemented for Al-27!" ;exit(EXIT_FAILURE);}
          }
          break;
        case 18:
          if(A==40){
            if(mode==0){for (size_t i = 0; i < 2; i++) { pars[i]=Ar_sph[i];} }
            else{std::cerr<< "Error! Non-spherical parametrization not implemented for Ar-40!" ;exit(EXIT_FAILURE);}
          }
          break;
        case 29:
          if(A==63){
            if(mode==0){for (size_t i = 0; i < 2; i++) { pars[i]=Cu_sph[i];} }
            else if(mode==1){for (size_t i = 0; i < 4; i++) { pars[i]=Cu_def[i];} }
            else{std::cerr<< "Error! Selected parametrization not implemented for Cu-63!" ;exit(EXIT_FAILURE);}
          }
          break;
        case 54:
          if(A==129){
            if(mode==0){for (size_t i = 0; i < 2; i++) { pars[i]=Xe_sph[i];} }
            else if(mode==1){for (size_t i = 0; i < 4; i++) { pars[i]=Xe_def[i];} }
            else{std::cerr<< "Error! Selected parametrization not implemented for Xe-129!" ;exit(EXIT_FAILURE);}
          }
          break;
        case 79:
          if(A==197){
            if(mode==0){for (size_t i = 0; i < 2; i++) { pars[i]=Au_sph[i];} }
            else if(mode==1){for (size_t i = 0; i < 4; i++) { pars[i]=Au_def[i];} }
            else{std::cerr<< "Error! Selected parametrization not implemented for Au-197!" ;exit(EXIT_FAILURE);}
          }
          break;
        case 82:
          if(A==208){
            if(mode==0){for (size_t i = 0; i < 2; i++) { pars[i]=Pb_sph[i];} }
            else{std::cerr<< "Error! Non-spherical parametrization not implemented for Pb-208!" ;exit(EXIT_FAILURE);}
          }          
          break;
        case 92:
          if(A==238){
            if(mode==1){for (size_t i = 0; i < 4; i++) { pars[i]=U_def[i];} }
            else{std::cerr<< "Error! Non-deformed parametrization not implemented for U-238!" ;exit(EXIT_FAILURE);}
          }
          break;
        default:
          std::cerr<< "Error! Nucleus not implemented!" ;
          exit(EXIT_FAILURE);

      }
    }

  }

  #endif /* nuclear_data */
 