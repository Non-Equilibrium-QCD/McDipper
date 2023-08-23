namespace H2{

int NumberOfNucleons=2;
int NumberOfSamples=300;

double *xPos;

int Index(int s,int i,int d){
return d+3*(i+NumberOfNucleons*s);
}

void Init(){

xPos=new double[3*NumberOfSamples*NumberOfNucleons];

xPos[Index( 0 ,0,0)]= 0.284959 ;  xPos[Index( 0 ,0,1)]= 3.8088 ;  xPos[Index( 0 ,0,2)]= -0.213055 ;
xPos[Index( 0 ,1,0)]= -0.284959 ;  xPos[Index( 0 ,1,1)]= -3.8088 ;  xPos[Index( 0 ,1,2)]= 0.213055 ;

xPos[Index( 1 ,0,0)]= 0.0821029 ;  xPos[Index( 1 ,0,1)]= -0.919647 ;  xPos[Index( 1 ,0,2)]= -0.461124 ;
xPos[Index( 1 ,1,0)]= -0.0821029 ;  xPos[Index( 1 ,1,1)]= 0.919647 ;  xPos[Index( 1 ,1,2)]= 0.461124 ;

xPos[Index( 2 ,0,0)]= -0.995285 ;  xPos[Index( 2 ,0,1)]= -1.49994 ;  xPos[Index( 2 ,0,2)]= 1.83565 ;
xPos[Index( 2 ,1,0)]= 0.995285 ;  xPos[Index( 2 ,1,1)]= 1.49994 ;  xPos[Index( 2 ,1,2)]= -1.83565 ;

xPos[Index( 3 ,0,0)]= 0.813212 ;  xPos[Index( 3 ,0,1)]= -0.989215 ;  xPos[Index( 3 ,0,2)]= -0.762924 ;
xPos[Index( 3 ,1,0)]= -0.813212 ;  xPos[Index( 3 ,1,1)]= 0.989215 ;  xPos[Index( 3 ,1,2)]= 0.762924 ;

xPos[Index( 4 ,0,0)]= -0.00350246 ;  xPos[Index( 4 ,0,1)]= 0.245365 ;  xPos[Index( 4 ,0,2)]= 0.64326 ;
xPos[Index( 4 ,1,0)]= 0.00350246 ;  xPos[Index( 4 ,1,1)]= -0.245365 ;  xPos[Index( 4 ,1,2)]= -0.64326 ;

xPos[Index( 5 ,0,0)]= -0.558697 ;  xPos[Index( 5 ,0,1)]= -0.355127 ;  xPos[Index( 5 ,0,2)]= 1.12605 ;
xPos[Index( 5 ,1,0)]= 0.558697 ;  xPos[Index( 5 ,1,1)]= 0.355127 ;  xPos[Index( 5 ,1,2)]= -1.12605 ;

xPos[Index( 6 ,0,0)]= 0.537212 ;  xPos[Index( 6 ,0,1)]= -0.0559806 ;  xPos[Index( 6 ,0,2)]= -1.45247 ;
xPos[Index( 6 ,1,0)]= -0.537212 ;  xPos[Index( 6 ,1,1)]= 0.0559806 ;  xPos[Index( 6 ,1,2)]= 1.45247 ;

xPos[Index( 7 ,0,0)]= 0.180559 ;  xPos[Index( 7 ,0,1)]= -1.02056 ;  xPos[Index( 7 ,0,2)]= -0.969269 ;
xPos[Index( 7 ,1,0)]= -0.180559 ;  xPos[Index( 7 ,1,1)]= 1.02056 ;  xPos[Index( 7 ,1,2)]= 0.969269 ;

xPos[Index( 8 ,0,0)]= -0.119586 ;  xPos[Index( 8 ,0,1)]= -1.91054 ;  xPos[Index( 8 ,0,2)]= 0.311179 ;
xPos[Index( 8 ,1,0)]= 0.119586 ;  xPos[Index( 8 ,1,1)]= 1.91054 ;  xPos[Index( 8 ,1,2)]= -0.311179 ;

xPos[Index( 9 ,0,0)]= -0.254927 ;  xPos[Index( 9 ,0,1)]= -1.53615 ;  xPos[Index( 9 ,0,2)]= -0.0965124 ;
xPos[Index( 9 ,1,0)]= 0.254927 ;  xPos[Index( 9 ,1,1)]= 1.53615 ;  xPos[Index( 9 ,1,2)]= 0.0965124 ;

xPos[Index( 10 ,0,0)]= -0.588479 ;  xPos[Index( 10 ,0,1)]= 0.768619 ;  xPos[Index( 10 ,0,2)]= -0.146882 ;
xPos[Index( 10 ,1,0)]= 0.588479 ;  xPos[Index( 10 ,1,1)]= -0.768619 ;  xPos[Index( 10 ,1,2)]= 0.146882 ;

xPos[Index( 11 ,0,0)]= 1.7728 ;  xPos[Index( 11 ,0,1)]= 0.611907 ;  xPos[Index( 11 ,0,2)]= 1.15434 ;
xPos[Index( 11 ,1,0)]= -1.7728 ;  xPos[Index( 11 ,1,1)]= -0.611907 ;  xPos[Index( 11 ,1,2)]= -1.15434 ;

xPos[Index( 12 ,0,0)]= -0.0216875 ;  xPos[Index( 12 ,0,1)]= 0.384991 ;  xPos[Index( 12 ,0,2)]= -1.15852 ;
xPos[Index( 12 ,1,0)]= 0.0216875 ;  xPos[Index( 12 ,1,1)]= -0.384991 ;  xPos[Index( 12 ,1,2)]= 1.15852 ;

xPos[Index( 13 ,0,0)]= 0.106948 ;  xPos[Index( 13 ,0,1)]= 0.904691 ;  xPos[Index( 13 ,0,2)]= -0.33013 ;
xPos[Index( 13 ,1,0)]= -0.106948 ;  xPos[Index( 13 ,1,1)]= -0.904691 ;  xPos[Index( 13 ,1,2)]= 0.33013 ;

xPos[Index( 14 ,0,0)]= 1.08897 ;  xPos[Index( 14 ,0,1)]= -1.17781 ;  xPos[Index( 14 ,0,2)]= -0.0363303 ;
xPos[Index( 14 ,1,0)]= -1.08897 ;  xPos[Index( 14 ,1,1)]= 1.17781 ;  xPos[Index( 14 ,1,2)]= 0.0363303 ;

xPos[Index( 15 ,0,0)]= -0.0513882 ;  xPos[Index( 15 ,0,1)]= -1.45462 ;  xPos[Index( 15 ,0,2)]= 0.327864 ;
xPos[Index( 15 ,1,0)]= 0.0513882 ;  xPos[Index( 15 ,1,1)]= 1.45462 ;  xPos[Index( 15 ,1,2)]= -0.327864 ;

xPos[Index( 16 ,0,0)]= 0.677903 ;  xPos[Index( 16 ,0,1)]= 0.636809 ;  xPos[Index( 16 ,0,2)]= 1.12366 ;
xPos[Index( 16 ,1,0)]= -0.677903 ;  xPos[Index( 16 ,1,1)]= -0.636809 ;  xPos[Index( 16 ,1,2)]= -1.12366 ;

xPos[Index( 17 ,0,0)]= -1.22361 ;  xPos[Index( 17 ,0,1)]= -0.350834 ;  xPos[Index( 17 ,0,2)]= 0.727709 ;
xPos[Index( 17 ,1,0)]= 1.22361 ;  xPos[Index( 17 ,1,1)]= 0.350834 ;  xPos[Index( 17 ,1,2)]= -0.727709 ;

xPos[Index( 18 ,0,0)]= -2.0259 ;  xPos[Index( 18 ,0,1)]= -0.763525 ;  xPos[Index( 18 ,0,2)]= 1.11937 ;
xPos[Index( 18 ,1,0)]= 2.0259 ;  xPos[Index( 18 ,1,1)]= 0.763525 ;  xPos[Index( 18 ,1,2)]= -1.11937 ;

xPos[Index( 19 ,0,0)]= -1.72149 ;  xPos[Index( 19 ,0,1)]= 0.274465 ;  xPos[Index( 19 ,0,2)]= 0.479987 ;
xPos[Index( 19 ,1,0)]= 1.72149 ;  xPos[Index( 19 ,1,1)]= -0.274465 ;  xPos[Index( 19 ,1,2)]= -0.479987 ;

xPos[Index( 20 ,0,0)]= 0.314667 ;  xPos[Index( 20 ,0,1)]= -0.197382 ;  xPos[Index( 20 ,0,2)]= 0.499684 ;
xPos[Index( 20 ,1,0)]= -0.314667 ;  xPos[Index( 20 ,1,1)]= 0.197382 ;  xPos[Index( 20 ,1,2)]= -0.499684 ;

xPos[Index( 21 ,0,0)]= -1.34138 ;  xPos[Index( 21 ,0,1)]= -0.147979 ;  xPos[Index( 21 ,0,2)]= -0.73972 ;
xPos[Index( 21 ,1,0)]= 1.34138 ;  xPos[Index( 21 ,1,1)]= 0.147979 ;  xPos[Index( 21 ,1,2)]= 0.73972 ;

xPos[Index( 22 ,0,0)]= 0.510871 ;  xPos[Index( 22 ,0,1)]= -0.169312 ;  xPos[Index( 22 ,0,2)]= 1.00629 ;
xPos[Index( 22 ,1,0)]= -0.510871 ;  xPos[Index( 22 ,1,1)]= 0.169312 ;  xPos[Index( 22 ,1,2)]= -1.00629 ;

xPos[Index( 23 ,0,0)]= -1.65536 ;  xPos[Index( 23 ,0,1)]= 0.898053 ;  xPos[Index( 23 ,0,2)]= 0.0904812 ;
xPos[Index( 23 ,1,0)]= 1.65536 ;  xPos[Index( 23 ,1,1)]= -0.898053 ;  xPos[Index( 23 ,1,2)]= -0.0904812 ;

xPos[Index( 24 ,0,0)]= -0.579595 ;  xPos[Index( 24 ,0,1)]= -0.613429 ;  xPos[Index( 24 ,0,2)]= 2.65964 ;
xPos[Index( 24 ,1,0)]= 0.579595 ;  xPos[Index( 24 ,1,1)]= 0.613429 ;  xPos[Index( 24 ,1,2)]= -2.65964 ;

xPos[Index( 25 ,0,0)]= 1.47806 ;  xPos[Index( 25 ,0,1)]= 1.48947 ;  xPos[Index( 25 ,0,2)]= -0.0935041 ;
xPos[Index( 25 ,1,0)]= -1.47806 ;  xPos[Index( 25 ,1,1)]= -1.48947 ;  xPos[Index( 25 ,1,2)]= 0.0935041 ;

xPos[Index( 26 ,0,0)]= -0.296233 ;  xPos[Index( 26 ,0,1)]= 0.532012 ;  xPos[Index( 26 ,0,2)]= -0.0373773 ;
xPos[Index( 26 ,1,0)]= 0.296233 ;  xPos[Index( 26 ,1,1)]= -0.532012 ;  xPos[Index( 26 ,1,2)]= 0.0373773 ;

xPos[Index( 27 ,0,0)]= 0.794997 ;  xPos[Index( 27 ,0,1)]= -0.549669 ;  xPos[Index( 27 ,0,2)]= -0.591306 ;
xPos[Index( 27 ,1,0)]= -0.794997 ;  xPos[Index( 27 ,1,1)]= 0.549669 ;  xPos[Index( 27 ,1,2)]= 0.591306 ;

xPos[Index( 28 ,0,0)]= 0.76119 ;  xPos[Index( 28 ,0,1)]= 2.62427 ;  xPos[Index( 28 ,0,2)]= -0.660497 ;
xPos[Index( 28 ,1,0)]= -0.76119 ;  xPos[Index( 28 ,1,1)]= -2.62427 ;  xPos[Index( 28 ,1,2)]= 0.660497 ;

xPos[Index( 29 ,0,0)]= 0.597405 ;  xPos[Index( 29 ,0,1)]= -0.55899 ;  xPos[Index( 29 ,0,2)]= -0.382729 ;
xPos[Index( 29 ,1,0)]= -0.597405 ;  xPos[Index( 29 ,1,1)]= 0.55899 ;  xPos[Index( 29 ,1,2)]= 0.382729 ;

xPos[Index( 30 ,0,0)]= -0.0223945 ;  xPos[Index( 30 ,0,1)]= -0.473286 ;  xPos[Index( 30 ,0,2)]= 1.23679 ;
xPos[Index( 30 ,1,0)]= 0.0223945 ;  xPos[Index( 30 ,1,1)]= 0.473286 ;  xPos[Index( 30 ,1,2)]= -1.23679 ;

xPos[Index( 31 ,0,0)]= -0.52068 ;  xPos[Index( 31 ,0,1)]= -0.203374 ;  xPos[Index( 31 ,0,2)]= -0.0474003 ;
xPos[Index( 31 ,1,0)]= 0.52068 ;  xPos[Index( 31 ,1,1)]= 0.203374 ;  xPos[Index( 31 ,1,2)]= 0.0474003 ;

xPos[Index( 32 ,0,0)]= -0.639296 ;  xPos[Index( 32 ,0,1)]= 1.20445 ;  xPos[Index( 32 ,0,2)]= 0.429774 ;
xPos[Index( 32 ,1,0)]= 0.639296 ;  xPos[Index( 32 ,1,1)]= -1.20445 ;  xPos[Index( 32 ,1,2)]= -0.429774 ;

xPos[Index( 33 ,0,0)]= 0.187993 ;  xPos[Index( 33 ,0,1)]= -0.882847 ;  xPos[Index( 33 ,0,2)]= 0.186962 ;
xPos[Index( 33 ,1,0)]= -0.187993 ;  xPos[Index( 33 ,1,1)]= 0.882847 ;  xPos[Index( 33 ,1,2)]= -0.186962 ;

xPos[Index( 34 ,0,0)]= 0.556676 ;  xPos[Index( 34 ,0,1)]= -0.00638006 ;  xPos[Index( 34 ,0,2)]= 0.199829 ;
xPos[Index( 34 ,1,0)]= -0.556676 ;  xPos[Index( 34 ,1,1)]= 0.00638006 ;  xPos[Index( 34 ,1,2)]= -0.199829 ;

xPos[Index( 35 ,0,0)]= 0.563306 ;  xPos[Index( 35 ,0,1)]= 0.0792969 ;  xPos[Index( 35 ,0,2)]= 0.323359 ;
xPos[Index( 35 ,1,0)]= -0.563306 ;  xPos[Index( 35 ,1,1)]= -0.0792969 ;  xPos[Index( 35 ,1,2)]= -0.323359 ;

xPos[Index( 36 ,0,0)]= -0.120886 ;  xPos[Index( 36 ,0,1)]= 0.590625 ;  xPos[Index( 36 ,0,2)]= 1.21924 ;
xPos[Index( 36 ,1,0)]= 0.120886 ;  xPos[Index( 36 ,1,1)]= -0.590625 ;  xPos[Index( 36 ,1,2)]= -1.21924 ;

xPos[Index( 37 ,0,0)]= -0.139468 ;  xPos[Index( 37 ,0,1)]= 0.0120337 ;  xPos[Index( 37 ,0,2)]= 1.06795 ;
xPos[Index( 37 ,1,0)]= 0.139468 ;  xPos[Index( 37 ,1,1)]= -0.0120337 ;  xPos[Index( 37 ,1,2)]= -1.06795 ;

xPos[Index( 38 ,0,0)]= 1.30372 ;  xPos[Index( 38 ,0,1)]= -0.23841 ;  xPos[Index( 38 ,0,2)]= -0.173698 ;
xPos[Index( 38 ,1,0)]= -1.30372 ;  xPos[Index( 38 ,1,1)]= 0.23841 ;  xPos[Index( 38 ,1,2)]= 0.173698 ;

xPos[Index( 39 ,0,0)]= -0.602664 ;  xPos[Index( 39 ,0,1)]= -1.58745 ;  xPos[Index( 39 ,0,2)]= -0.511846 ;
xPos[Index( 39 ,1,0)]= 0.602664 ;  xPos[Index( 39 ,1,1)]= 1.58745 ;  xPos[Index( 39 ,1,2)]= 0.511846 ;

xPos[Index( 40 ,0,0)]= 0.745889 ;  xPos[Index( 40 ,0,1)]= -1.34572 ;  xPos[Index( 40 ,0,2)]= 0.447314 ;
xPos[Index( 40 ,1,0)]= -0.745889 ;  xPos[Index( 40 ,1,1)]= 1.34572 ;  xPos[Index( 40 ,1,2)]= -0.447314 ;

xPos[Index( 41 ,0,0)]= -1.22214 ;  xPos[Index( 41 ,0,1)]= -0.207401 ;  xPos[Index( 41 ,0,2)]= -0.736717 ;
xPos[Index( 41 ,1,0)]= 1.22214 ;  xPos[Index( 41 ,1,1)]= 0.207401 ;  xPos[Index( 41 ,1,2)]= 0.736717 ;

xPos[Index( 42 ,0,0)]= 0.912006 ;  xPos[Index( 42 ,0,1)]= 0.620062 ;  xPos[Index( 42 ,0,2)]= -1.39723 ;
xPos[Index( 42 ,1,0)]= -0.912006 ;  xPos[Index( 42 ,1,1)]= -0.620062 ;  xPos[Index( 42 ,1,2)]= 1.39723 ;

xPos[Index( 43 ,0,0)]= -0.0502409 ;  xPos[Index( 43 ,0,1)]= -1.48962 ;  xPos[Index( 43 ,0,2)]= -0.654058 ;
xPos[Index( 43 ,1,0)]= 0.0502409 ;  xPos[Index( 43 ,1,1)]= 1.48962 ;  xPos[Index( 43 ,1,2)]= 0.654058 ;

xPos[Index( 44 ,0,0)]= 0.213393 ;  xPos[Index( 44 ,0,1)]= -1.4747 ;  xPos[Index( 44 ,0,2)]= 0.69472 ;
xPos[Index( 44 ,1,0)]= -0.213393 ;  xPos[Index( 44 ,1,1)]= 1.4747 ;  xPos[Index( 44 ,1,2)]= -0.69472 ;

xPos[Index( 45 ,0,0)]= -0.921806 ;  xPos[Index( 45 ,0,1)]= -0.816478 ;  xPos[Index( 45 ,0,2)]= -0.959842 ;
xPos[Index( 45 ,1,0)]= 0.921806 ;  xPos[Index( 45 ,1,1)]= 0.816478 ;  xPos[Index( 45 ,1,2)]= 0.959842 ;

xPos[Index( 46 ,0,0)]= 0.200827 ;  xPos[Index( 46 ,0,1)]= 0.0135871 ;  xPos[Index( 46 ,0,2)]= 1.5861 ;
xPos[Index( 46 ,1,0)]= -0.200827 ;  xPos[Index( 46 ,1,1)]= -0.0135871 ;  xPos[Index( 46 ,1,2)]= -1.5861 ;

xPos[Index( 47 ,0,0)]= -0.256865 ;  xPos[Index( 47 ,0,1)]= 0.721754 ;  xPos[Index( 47 ,0,2)]= -0.817415 ;
xPos[Index( 47 ,1,0)]= 0.256865 ;  xPos[Index( 47 ,1,1)]= -0.721754 ;  xPos[Index( 47 ,1,2)]= 0.817415 ;

xPos[Index( 48 ,0,0)]= -0.329594 ;  xPos[Index( 48 ,0,1)]= 0.656078 ;  xPos[Index( 48 ,0,2)]= -0.726994 ;
xPos[Index( 48 ,1,0)]= 0.329594 ;  xPos[Index( 48 ,1,1)]= -0.656078 ;  xPos[Index( 48 ,1,2)]= 0.726994 ;

xPos[Index( 49 ,0,0)]= 0.219162 ;  xPos[Index( 49 ,0,1)]= -1.35735 ;  xPos[Index( 49 ,0,2)]= 0.45384 ;
xPos[Index( 49 ,1,0)]= -0.219162 ;  xPos[Index( 49 ,1,1)]= 1.35735 ;  xPos[Index( 49 ,1,2)]= -0.45384 ;

xPos[Index( 50 ,0,0)]= -0.243331 ;  xPos[Index( 50 ,0,1)]= 0.026992 ;  xPos[Index( 50 ,0,2)]= -0.900626 ;
xPos[Index( 50 ,1,0)]= 0.243331 ;  xPos[Index( 50 ,1,1)]= -0.026992 ;  xPos[Index( 50 ,1,2)]= 0.900626 ;

xPos[Index( 51 ,0,0)]= 0.0958838 ;  xPos[Index( 51 ,0,1)]= 0.248912 ;  xPos[Index( 51 ,0,2)]= 0.681637 ;
xPos[Index( 51 ,1,0)]= -0.0958838 ;  xPos[Index( 51 ,1,1)]= -0.248912 ;  xPos[Index( 51 ,1,2)]= -0.681637 ;

xPos[Index( 52 ,0,0)]= 0.0923201 ;  xPos[Index( 52 ,0,1)]= -0.0273433 ;  xPos[Index( 52 ,0,2)]= 0.0273337 ;
xPos[Index( 52 ,1,0)]= -0.0923201 ;  xPos[Index( 52 ,1,1)]= 0.0273433 ;  xPos[Index( 52 ,1,2)]= -0.0273337 ;

xPos[Index( 53 ,0,0)]= -0.0587248 ;  xPos[Index( 53 ,0,1)]= 0.380068 ;  xPos[Index( 53 ,0,2)]= 0.252059 ;
xPos[Index( 53 ,1,0)]= 0.0587248 ;  xPos[Index( 53 ,1,1)]= -0.380068 ;  xPos[Index( 53 ,1,2)]= -0.252059 ;

xPos[Index( 54 ,0,0)]= 0.623514 ;  xPos[Index( 54 ,0,1)]= 0.562993 ;  xPos[Index( 54 ,0,2)]= 0.026217 ;
xPos[Index( 54 ,1,0)]= -0.623514 ;  xPos[Index( 54 ,1,1)]= -0.562993 ;  xPos[Index( 54 ,1,2)]= -0.026217 ;

xPos[Index( 55 ,0,0)]= -0.300407 ;  xPos[Index( 55 ,0,1)]= 0.589541 ;  xPos[Index( 55 ,0,2)]= 0.349212 ;
xPos[Index( 55 ,1,0)]= 0.300407 ;  xPos[Index( 55 ,1,1)]= -0.589541 ;  xPos[Index( 55 ,1,2)]= -0.349212 ;

xPos[Index( 56 ,0,0)]= -1.58212 ;  xPos[Index( 56 ,0,1)]= -0.339552 ;  xPos[Index( 56 ,0,2)]= 1.93865 ;
xPos[Index( 56 ,1,0)]= 1.58212 ;  xPos[Index( 56 ,1,1)]= 0.339552 ;  xPos[Index( 56 ,1,2)]= -1.93865 ;

xPos[Index( 57 ,0,0)]= 0.350936 ;  xPos[Index( 57 ,0,1)]= 0.00175317 ;  xPos[Index( 57 ,0,2)]= -0.425785 ;
xPos[Index( 57 ,1,0)]= -0.350936 ;  xPos[Index( 57 ,1,1)]= -0.00175317 ;  xPos[Index( 57 ,1,2)]= 0.425785 ;

xPos[Index( 58 ,0,0)]= 0.820896 ;  xPos[Index( 58 ,0,1)]= 0.582331 ;  xPos[Index( 58 ,0,2)]= 0.42267 ;
xPos[Index( 58 ,1,0)]= -0.820896 ;  xPos[Index( 58 ,1,1)]= -0.582331 ;  xPos[Index( 58 ,1,2)]= -0.42267 ;

xPos[Index( 59 ,0,0)]= 1.33022 ;  xPos[Index( 59 ,0,1)]= -1.08085 ;  xPos[Index( 59 ,0,2)]= -1.84108 ;
xPos[Index( 59 ,1,0)]= -1.33022 ;  xPos[Index( 59 ,1,1)]= 1.08085 ;  xPos[Index( 59 ,1,2)]= 1.84108 ;

xPos[Index( 60 ,0,0)]= 0.629628 ;  xPos[Index( 60 ,0,1)]= 1.05532 ;  xPos[Index( 60 ,0,2)]= -0.234912 ;
xPos[Index( 60 ,1,0)]= -0.629628 ;  xPos[Index( 60 ,1,1)]= -1.05532 ;  xPos[Index( 60 ,1,2)]= 0.234912 ;

xPos[Index( 61 ,0,0)]= 0.328252 ;  xPos[Index( 61 ,0,1)]= 1.09065 ;  xPos[Index( 61 ,0,2)]= -0.0541522 ;
xPos[Index( 61 ,1,0)]= -0.328252 ;  xPos[Index( 61 ,1,1)]= -1.09065 ;  xPos[Index( 61 ,1,2)]= 0.0541522 ;

xPos[Index( 62 ,0,0)]= -0.116297 ;  xPos[Index( 62 ,0,1)]= -0.850493 ;  xPos[Index( 62 ,0,2)]= -0.0851858 ;
xPos[Index( 62 ,1,0)]= 0.116297 ;  xPos[Index( 62 ,1,1)]= 0.850493 ;  xPos[Index( 62 ,1,2)]= 0.0851858 ;

xPos[Index( 63 ,0,0)]= -0.78202 ;  xPos[Index( 63 ,0,1)]= 2.01367 ;  xPos[Index( 63 ,0,2)]= -0.291797 ;
xPos[Index( 63 ,1,0)]= 0.78202 ;  xPos[Index( 63 ,1,1)]= -2.01367 ;  xPos[Index( 63 ,1,2)]= 0.291797 ;

xPos[Index( 64 ,0,0)]= 1.48094 ;  xPos[Index( 64 ,0,1)]= 0.313099 ;  xPos[Index( 64 ,0,2)]= 0.177203 ;
xPos[Index( 64 ,1,0)]= -1.48094 ;  xPos[Index( 64 ,1,1)]= -0.313099 ;  xPos[Index( 64 ,1,2)]= -0.177203 ;

xPos[Index( 65 ,0,0)]= 0.0268278 ;  xPos[Index( 65 ,0,1)]= -0.342364 ;  xPos[Index( 65 ,0,2)]= -0.295339 ;
xPos[Index( 65 ,1,0)]= -0.0268278 ;  xPos[Index( 65 ,1,1)]= 0.342364 ;  xPos[Index( 65 ,1,2)]= 0.295339 ;

xPos[Index( 66 ,0,0)]= 0.98258 ;  xPos[Index( 66 ,0,1)]= 0.148609 ;  xPos[Index( 66 ,0,2)]= 0.346044 ;
xPos[Index( 66 ,1,0)]= -0.98258 ;  xPos[Index( 66 ,1,1)]= -0.148609 ;  xPos[Index( 66 ,1,2)]= -0.346044 ;

xPos[Index( 67 ,0,0)]= 0.749975 ;  xPos[Index( 67 ,0,1)]= 0.0452311 ;  xPos[Index( 67 ,0,2)]= -0.118104 ;
xPos[Index( 67 ,1,0)]= -0.749975 ;  xPos[Index( 67 ,1,1)]= -0.0452311 ;  xPos[Index( 67 ,1,2)]= 0.118104 ;

xPos[Index( 68 ,0,0)]= -0.462317 ;  xPos[Index( 68 ,0,1)]= 0.60775 ;  xPos[Index( 68 ,0,2)]= -0.728648 ;
xPos[Index( 68 ,1,0)]= 0.462317 ;  xPos[Index( 68 ,1,1)]= -0.60775 ;  xPos[Index( 68 ,1,2)]= 0.728648 ;

xPos[Index( 69 ,0,0)]= -0.0921222 ;  xPos[Index( 69 ,0,1)]= 0.464356 ;  xPos[Index( 69 ,0,2)]= -0.0381018 ;
xPos[Index( 69 ,1,0)]= 0.0921222 ;  xPos[Index( 69 ,1,1)]= -0.464356 ;  xPos[Index( 69 ,1,2)]= 0.0381018 ;

xPos[Index( 70 ,0,0)]= -0.450977 ;  xPos[Index( 70 ,0,1)]= -0.494756 ;  xPos[Index( 70 ,0,2)]= 0.194223 ;
xPos[Index( 70 ,1,0)]= 0.450977 ;  xPos[Index( 70 ,1,1)]= 0.494756 ;  xPos[Index( 70 ,1,2)]= -0.194223 ;

xPos[Index( 71 ,0,0)]= -0.76624 ;  xPos[Index( 71 ,0,1)]= -0.539475 ;  xPos[Index( 71 ,0,2)]= -1.1294 ;
xPos[Index( 71 ,1,0)]= 0.76624 ;  xPos[Index( 71 ,1,1)]= 0.539475 ;  xPos[Index( 71 ,1,2)]= 1.1294 ;

xPos[Index( 72 ,0,0)]= -0.176049 ;  xPos[Index( 72 ,0,1)]= 1.20794 ;  xPos[Index( 72 ,0,2)]= -1.9421 ;
xPos[Index( 72 ,1,0)]= 0.176049 ;  xPos[Index( 72 ,1,1)]= -1.20794 ;  xPos[Index( 72 ,1,2)]= 1.9421 ;

xPos[Index( 73 ,0,0)]= -0.627749 ;  xPos[Index( 73 ,0,1)]= -1.32281 ;  xPos[Index( 73 ,0,2)]= -0.590148 ;
xPos[Index( 73 ,1,0)]= 0.627749 ;  xPos[Index( 73 ,1,1)]= 1.32281 ;  xPos[Index( 73 ,1,2)]= 0.590148 ;

xPos[Index( 74 ,0,0)]= -0.0513461 ;  xPos[Index( 74 ,0,1)]= 0.441765 ;  xPos[Index( 74 ,0,2)]= 1.45987 ;
xPos[Index( 74 ,1,0)]= 0.0513461 ;  xPos[Index( 74 ,1,1)]= -0.441765 ;  xPos[Index( 74 ,1,2)]= -1.45987 ;

xPos[Index( 75 ,0,0)]= -0.494321 ;  xPos[Index( 75 ,0,1)]= -0.822041 ;  xPos[Index( 75 ,0,2)]= -0.447285 ;
xPos[Index( 75 ,1,0)]= 0.494321 ;  xPos[Index( 75 ,1,1)]= 0.822041 ;  xPos[Index( 75 ,1,2)]= 0.447285 ;

xPos[Index( 76 ,0,0)]= 0.546237 ;  xPos[Index( 76 ,0,1)]= 0.739223 ;  xPos[Index( 76 ,0,2)]= -0.245953 ;
xPos[Index( 76 ,1,0)]= -0.546237 ;  xPos[Index( 76 ,1,1)]= -0.739223 ;  xPos[Index( 76 ,1,2)]= 0.245953 ;

xPos[Index( 77 ,0,0)]= 0.182618 ;  xPos[Index( 77 ,0,1)]= -1.5907 ;  xPos[Index( 77 ,0,2)]= -0.581011 ;
xPos[Index( 77 ,1,0)]= -0.182618 ;  xPos[Index( 77 ,1,1)]= 1.5907 ;  xPos[Index( 77 ,1,2)]= 0.581011 ;

xPos[Index( 78 ,0,0)]= 0.926892 ;  xPos[Index( 78 ,0,1)]= 1.69769 ;  xPos[Index( 78 ,0,2)]= -0.550845 ;
xPos[Index( 78 ,1,0)]= -0.926892 ;  xPos[Index( 78 ,1,1)]= -1.69769 ;  xPos[Index( 78 ,1,2)]= 0.550845 ;

xPos[Index( 79 ,0,0)]= -0.0829944 ;  xPos[Index( 79 ,0,1)]= -0.44859 ;  xPos[Index( 79 ,0,2)]= -0.163362 ;
xPos[Index( 79 ,1,0)]= 0.0829944 ;  xPos[Index( 79 ,1,1)]= 0.44859 ;  xPos[Index( 79 ,1,2)]= 0.163362 ;

xPos[Index( 80 ,0,0)]= -0.43475 ;  xPos[Index( 80 ,0,1)]= 0.22518 ;  xPos[Index( 80 ,0,2)]= -0.532497 ;
xPos[Index( 80 ,1,0)]= 0.43475 ;  xPos[Index( 80 ,1,1)]= -0.22518 ;  xPos[Index( 80 ,1,2)]= 0.532497 ;

xPos[Index( 81 ,0,0)]= -2.02045 ;  xPos[Index( 81 ,0,1)]= 0.93484 ;  xPos[Index( 81 ,0,2)]= 0.581666 ;
xPos[Index( 81 ,1,0)]= 2.02045 ;  xPos[Index( 81 ,1,1)]= -0.93484 ;  xPos[Index( 81 ,1,2)]= -0.581666 ;

xPos[Index( 82 ,0,0)]= -0.454221 ;  xPos[Index( 82 ,0,1)]= 1.80082 ;  xPos[Index( 82 ,0,2)]= -0.580749 ;
xPos[Index( 82 ,1,0)]= 0.454221 ;  xPos[Index( 82 ,1,1)]= -1.80082 ;  xPos[Index( 82 ,1,2)]= 0.580749 ;

xPos[Index( 83 ,0,0)]= 0.0656841 ;  xPos[Index( 83 ,0,1)]= 0.989576 ;  xPos[Index( 83 ,0,2)]= -0.638361 ;
xPos[Index( 83 ,1,0)]= -0.0656841 ;  xPos[Index( 83 ,1,1)]= -0.989576 ;  xPos[Index( 83 ,1,2)]= 0.638361 ;

xPos[Index( 84 ,0,0)]= -0.31396 ;  xPos[Index( 84 ,0,1)]= -0.427186 ;  xPos[Index( 84 ,0,2)]= 0.666653 ;
xPos[Index( 84 ,1,0)]= 0.31396 ;  xPos[Index( 84 ,1,1)]= 0.427186 ;  xPos[Index( 84 ,1,2)]= -0.666653 ;

xPos[Index( 85 ,0,0)]= 0.632491 ;  xPos[Index( 85 ,0,1)]= 1.22383 ;  xPos[Index( 85 ,0,2)]= -1.53904 ;
xPos[Index( 85 ,1,0)]= -0.632491 ;  xPos[Index( 85 ,1,1)]= -1.22383 ;  xPos[Index( 85 ,1,2)]= 1.53904 ;

xPos[Index( 86 ,0,0)]= -0.299993 ;  xPos[Index( 86 ,0,1)]= -0.214711 ;  xPos[Index( 86 ,0,2)]= -0.31192 ;
xPos[Index( 86 ,1,0)]= 0.299993 ;  xPos[Index( 86 ,1,1)]= 0.214711 ;  xPos[Index( 86 ,1,2)]= 0.31192 ;

xPos[Index( 87 ,0,0)]= -0.656514 ;  xPos[Index( 87 ,0,1)]= -0.392781 ;  xPos[Index( 87 ,0,2)]= 0.0374674 ;
xPos[Index( 87 ,1,0)]= 0.656514 ;  xPos[Index( 87 ,1,1)]= 0.392781 ;  xPos[Index( 87 ,1,2)]= -0.0374674 ;

xPos[Index( 88 ,0,0)]= -0.142159 ;  xPos[Index( 88 ,0,1)]= 0.697954 ;  xPos[Index( 88 ,0,2)]= -1.45954 ;
xPos[Index( 88 ,1,0)]= 0.142159 ;  xPos[Index( 88 ,1,1)]= -0.697954 ;  xPos[Index( 88 ,1,2)]= 1.45954 ;

xPos[Index( 89 ,0,0)]= -1.01068 ;  xPos[Index( 89 ,0,1)]= -0.354703 ;  xPos[Index( 89 ,0,2)]= 0.443815 ;
xPos[Index( 89 ,1,0)]= 1.01068 ;  xPos[Index( 89 ,1,1)]= 0.354703 ;  xPos[Index( 89 ,1,2)]= -0.443815 ;

xPos[Index( 90 ,0,0)]= 0.982134 ;  xPos[Index( 90 ,0,1)]= -0.34248 ;  xPos[Index( 90 ,0,2)]= -1.07409 ;
xPos[Index( 90 ,1,0)]= -0.982134 ;  xPos[Index( 90 ,1,1)]= 0.34248 ;  xPos[Index( 90 ,1,2)]= 1.07409 ;

xPos[Index( 91 ,0,0)]= 0.113364 ;  xPos[Index( 91 ,0,1)]= -0.83525 ;  xPos[Index( 91 ,0,2)]= 1.91865 ;
xPos[Index( 91 ,1,0)]= -0.113364 ;  xPos[Index( 91 ,1,1)]= 0.83525 ;  xPos[Index( 91 ,1,2)]= -1.91865 ;

xPos[Index( 92 ,0,0)]= -1.07662 ;  xPos[Index( 92 ,0,1)]= -0.660303 ;  xPos[Index( 92 ,0,2)]= 0.296062 ;
xPos[Index( 92 ,1,0)]= 1.07662 ;  xPos[Index( 92 ,1,1)]= 0.660303 ;  xPos[Index( 92 ,1,2)]= -0.296062 ;

xPos[Index( 93 ,0,0)]= 0.266191 ;  xPos[Index( 93 ,0,1)]= 0.321094 ;  xPos[Index( 93 ,0,2)]= 0.939163 ;
xPos[Index( 93 ,1,0)]= -0.266191 ;  xPos[Index( 93 ,1,1)]= -0.321094 ;  xPos[Index( 93 ,1,2)]= -0.939163 ;

xPos[Index( 94 ,0,0)]= 0.597757 ;  xPos[Index( 94 ,0,1)]= 0.403252 ;  xPos[Index( 94 ,0,2)]= 0.20717 ;
xPos[Index( 94 ,1,0)]= -0.597757 ;  xPos[Index( 94 ,1,1)]= -0.403252 ;  xPos[Index( 94 ,1,2)]= -0.20717 ;

xPos[Index( 95 ,0,0)]= 0.546199 ;  xPos[Index( 95 ,0,1)]= -0.614284 ;  xPos[Index( 95 ,0,2)]= -0.622442 ;
xPos[Index( 95 ,1,0)]= -0.546199 ;  xPos[Index( 95 ,1,1)]= 0.614284 ;  xPos[Index( 95 ,1,2)]= 0.622442 ;

xPos[Index( 96 ,0,0)]= 0.559434 ;  xPos[Index( 96 ,0,1)]= -0.615724 ;  xPos[Index( 96 ,0,2)]= 0.650919 ;
xPos[Index( 96 ,1,0)]= -0.559434 ;  xPos[Index( 96 ,1,1)]= 0.615724 ;  xPos[Index( 96 ,1,2)]= -0.650919 ;

xPos[Index( 97 ,0,0)]= 0.414285 ;  xPos[Index( 97 ,0,1)]= 0.211526 ;  xPos[Index( 97 ,0,2)]= -0.975547 ;
xPos[Index( 97 ,1,0)]= -0.414285 ;  xPos[Index( 97 ,1,1)]= -0.211526 ;  xPos[Index( 97 ,1,2)]= 0.975547 ;

xPos[Index( 98 ,0,0)]= 1.37762 ;  xPos[Index( 98 ,0,1)]= 1.17456 ;  xPos[Index( 98 ,0,2)]= -0.352817 ;
xPos[Index( 98 ,1,0)]= -1.37762 ;  xPos[Index( 98 ,1,1)]= -1.17456 ;  xPos[Index( 98 ,1,2)]= 0.352817 ;

xPos[Index( 99 ,0,0)]= -0.485994 ;  xPos[Index( 99 ,0,1)]= -0.263722 ;  xPos[Index( 99 ,0,2)]= -0.991692 ;
xPos[Index( 99 ,1,0)]= 0.485994 ;  xPos[Index( 99 ,1,1)]= 0.263722 ;  xPos[Index( 99 ,1,2)]= 0.991692 ;

xPos[Index( 100 ,0,0)]= -0.0648751 ;  xPos[Index( 100 ,0,1)]= 0.558885 ;  xPos[Index( 100 ,0,2)]= 0.73511 ;
xPos[Index( 100 ,1,0)]= 0.0648751 ;  xPos[Index( 100 ,1,1)]= -0.558885 ;  xPos[Index( 100 ,1,2)]= -0.73511 ;

xPos[Index( 101 ,0,0)]= -0.380593 ;  xPos[Index( 101 ,0,1)]= 0.169672 ;  xPos[Index( 101 ,0,2)]= -0.923179 ;
xPos[Index( 101 ,1,0)]= 0.380593 ;  xPos[Index( 101 ,1,1)]= -0.169672 ;  xPos[Index( 101 ,1,2)]= 0.923179 ;

xPos[Index( 102 ,0,0)]= -0.691077 ;  xPos[Index( 102 ,0,1)]= -1.20428 ;  xPos[Index( 102 ,0,2)]= -0.0542309 ;
xPos[Index( 102 ,1,0)]= 0.691077 ;  xPos[Index( 102 ,1,1)]= 1.20428 ;  xPos[Index( 102 ,1,2)]= 0.0542309 ;

xPos[Index( 103 ,0,0)]= 0.723559 ;  xPos[Index( 103 ,0,1)]= 0.807037 ;  xPos[Index( 103 ,0,2)]= -0.169394 ;
xPos[Index( 103 ,1,0)]= -0.723559 ;  xPos[Index( 103 ,1,1)]= -0.807037 ;  xPos[Index( 103 ,1,2)]= 0.169394 ;

xPos[Index( 104 ,0,0)]= -0.353194 ;  xPos[Index( 104 ,0,1)]= 0.446116 ;  xPos[Index( 104 ,0,2)]= 0.43473 ;
xPos[Index( 104 ,1,0)]= 0.353194 ;  xPos[Index( 104 ,1,1)]= -0.446116 ;  xPos[Index( 104 ,1,2)]= -0.43473 ;

xPos[Index( 105 ,0,0)]= -0.935439 ;  xPos[Index( 105 ,0,1)]= 0.0846375 ;  xPos[Index( 105 ,0,2)]= 0.236976 ;
xPos[Index( 105 ,1,0)]= 0.935439 ;  xPos[Index( 105 ,1,1)]= -0.0846375 ;  xPos[Index( 105 ,1,2)]= -0.236976 ;

xPos[Index( 106 ,0,0)]= -0.674447 ;  xPos[Index( 106 ,0,1)]= -0.0699463 ;  xPos[Index( 106 ,0,2)]= -0.715557 ;
xPos[Index( 106 ,1,0)]= 0.674447 ;  xPos[Index( 106 ,1,1)]= 0.0699463 ;  xPos[Index( 106 ,1,2)]= 0.715557 ;

xPos[Index( 107 ,0,0)]= 0.246195 ;  xPos[Index( 107 ,0,1)]= 0.551032 ;  xPos[Index( 107 ,0,2)]= 1.95453 ;
xPos[Index( 107 ,1,0)]= -0.246195 ;  xPos[Index( 107 ,1,1)]= -0.551032 ;  xPos[Index( 107 ,1,2)]= -1.95453 ;

xPos[Index( 108 ,0,0)]= 0.616339 ;  xPos[Index( 108 ,0,1)]= -0.0208583 ;  xPos[Index( 108 ,0,2)]= 0.0428958 ;
xPos[Index( 108 ,1,0)]= -0.616339 ;  xPos[Index( 108 ,1,1)]= 0.0208583 ;  xPos[Index( 108 ,1,2)]= -0.0428958 ;

xPos[Index( 109 ,0,0)]= 0.659109 ;  xPos[Index( 109 ,0,1)]= 0.0777134 ;  xPos[Index( 109 ,0,2)]= 0.486317 ;
xPos[Index( 109 ,1,0)]= -0.659109 ;  xPos[Index( 109 ,1,1)]= -0.0777134 ;  xPos[Index( 109 ,1,2)]= -0.486317 ;

xPos[Index( 110 ,0,0)]= -0.0802103 ;  xPos[Index( 110 ,0,1)]= -0.837556 ;  xPos[Index( 110 ,0,2)]= 0.117072 ;
xPos[Index( 110 ,1,0)]= 0.0802103 ;  xPos[Index( 110 ,1,1)]= 0.837556 ;  xPos[Index( 110 ,1,2)]= -0.117072 ;

xPos[Index( 111 ,0,0)]= 0.14576 ;  xPos[Index( 111 ,0,1)]= -0.394205 ;  xPos[Index( 111 ,0,2)]= -0.589472 ;
xPos[Index( 111 ,1,0)]= -0.14576 ;  xPos[Index( 111 ,1,1)]= 0.394205 ;  xPos[Index( 111 ,1,2)]= 0.589472 ;

xPos[Index( 112 ,0,0)]= -0.616533 ;  xPos[Index( 112 ,0,1)]= 0.174206 ;  xPos[Index( 112 ,0,2)]= -0.037556 ;
xPos[Index( 112 ,1,0)]= 0.616533 ;  xPos[Index( 112 ,1,1)]= -0.174206 ;  xPos[Index( 112 ,1,2)]= 0.037556 ;

xPos[Index( 113 ,0,0)]= 0.241262 ;  xPos[Index( 113 ,0,1)]= -0.349889 ;  xPos[Index( 113 ,0,2)]= 0.325662 ;
xPos[Index( 113 ,1,0)]= -0.241262 ;  xPos[Index( 113 ,1,1)]= 0.349889 ;  xPos[Index( 113 ,1,2)]= -0.325662 ;

xPos[Index( 114 ,0,0)]= -0.429288 ;  xPos[Index( 114 ,0,1)]= 0.16188 ;  xPos[Index( 114 ,0,2)]= 0.883866 ;
xPos[Index( 114 ,1,0)]= 0.429288 ;  xPos[Index( 114 ,1,1)]= -0.16188 ;  xPos[Index( 114 ,1,2)]= -0.883866 ;

xPos[Index( 115 ,0,0)]= 0.776383 ;  xPos[Index( 115 ,0,1)]= 0.230963 ;  xPos[Index( 115 ,0,2)]= -1.19386 ;
xPos[Index( 115 ,1,0)]= -0.776383 ;  xPos[Index( 115 ,1,1)]= -0.230963 ;  xPos[Index( 115 ,1,2)]= 1.19386 ;

xPos[Index( 116 ,0,0)]= -0.590654 ;  xPos[Index( 116 ,0,1)]= 1.36057 ;  xPos[Index( 116 ,0,2)]= 0.811453 ;
xPos[Index( 116 ,1,0)]= 0.590654 ;  xPos[Index( 116 ,1,1)]= -1.36057 ;  xPos[Index( 116 ,1,2)]= -0.811453 ;

xPos[Index( 117 ,0,0)]= 0.224944 ;  xPos[Index( 117 ,0,1)]= -0.75038 ;  xPos[Index( 117 ,0,2)]= -1.02035 ;
xPos[Index( 117 ,1,0)]= -0.224944 ;  xPos[Index( 117 ,1,1)]= 0.75038 ;  xPos[Index( 117 ,1,2)]= 1.02035 ;

xPos[Index( 118 ,0,0)]= -0.543347 ;  xPos[Index( 118 ,0,1)]= -0.00499034 ;  xPos[Index( 118 ,0,2)]= -0.450391 ;
xPos[Index( 118 ,1,0)]= 0.543347 ;  xPos[Index( 118 ,1,1)]= 0.00499034 ;  xPos[Index( 118 ,1,2)]= 0.450391 ;

xPos[Index( 119 ,0,0)]= -0.847965 ;  xPos[Index( 119 ,0,1)]= -1.30189 ;  xPos[Index( 119 ,0,2)]= 0.560485 ;
xPos[Index( 119 ,1,0)]= 0.847965 ;  xPos[Index( 119 ,1,1)]= 1.30189 ;  xPos[Index( 119 ,1,2)]= -0.560485 ;

xPos[Index( 120 ,0,0)]= 0.286725 ;  xPos[Index( 120 ,0,1)]= -0.0945088 ;  xPos[Index( 120 ,0,2)]= -0.828817 ;
xPos[Index( 120 ,1,0)]= -0.286725 ;  xPos[Index( 120 ,1,1)]= 0.0945088 ;  xPos[Index( 120 ,1,2)]= 0.828817 ;

xPos[Index( 121 ,0,0)]= 0.904756 ;  xPos[Index( 121 ,0,1)]= 0.693947 ;  xPos[Index( 121 ,0,2)]= 0.442764 ;
xPos[Index( 121 ,1,0)]= -0.904756 ;  xPos[Index( 121 ,1,1)]= -0.693947 ;  xPos[Index( 121 ,1,2)]= -0.442764 ;

xPos[Index( 122 ,0,0)]= -1.00002 ;  xPos[Index( 122 ,0,1)]= -1.00647 ;  xPos[Index( 122 ,0,2)]= -0.124055 ;
xPos[Index( 122 ,1,0)]= 1.00002 ;  xPos[Index( 122 ,1,1)]= 1.00647 ;  xPos[Index( 122 ,1,2)]= 0.124055 ;

xPos[Index( 123 ,0,0)]= 0.280381 ;  xPos[Index( 123 ,0,1)]= 0.137264 ;  xPos[Index( 123 ,0,2)]= -0.935339 ;
xPos[Index( 123 ,1,0)]= -0.280381 ;  xPos[Index( 123 ,1,1)]= -0.137264 ;  xPos[Index( 123 ,1,2)]= 0.935339 ;

xPos[Index( 124 ,0,0)]= 1.13942 ;  xPos[Index( 124 ,0,1)]= 0.715675 ;  xPos[Index( 124 ,0,2)]= -0.241377 ;
xPos[Index( 124 ,1,0)]= -1.13942 ;  xPos[Index( 124 ,1,1)]= -0.715675 ;  xPos[Index( 124 ,1,2)]= 0.241377 ;

xPos[Index( 125 ,0,0)]= -1.02758 ;  xPos[Index( 125 ,0,1)]= 0.803924 ;  xPos[Index( 125 ,0,2)]= -0.623545 ;
xPos[Index( 125 ,1,0)]= 1.02758 ;  xPos[Index( 125 ,1,1)]= -0.803924 ;  xPos[Index( 125 ,1,2)]= 0.623545 ;

xPos[Index( 126 ,0,0)]= -0.698397 ;  xPos[Index( 126 ,0,1)]= -0.436969 ;  xPos[Index( 126 ,0,2)]= -0.379779 ;
xPos[Index( 126 ,1,0)]= 0.698397 ;  xPos[Index( 126 ,1,1)]= 0.436969 ;  xPos[Index( 126 ,1,2)]= 0.379779 ;

xPos[Index( 127 ,0,0)]= -1.62607 ;  xPos[Index( 127 ,0,1)]= -1.18098 ;  xPos[Index( 127 ,0,2)]= 0.113314 ;
xPos[Index( 127 ,1,0)]= 1.62607 ;  xPos[Index( 127 ,1,1)]= 1.18098 ;  xPos[Index( 127 ,1,2)]= -0.113314 ;

xPos[Index( 128 ,0,0)]= 0.626782 ;  xPos[Index( 128 ,0,1)]= -1.56422 ;  xPos[Index( 128 ,0,2)]= 1.29163 ;
xPos[Index( 128 ,1,0)]= -0.626782 ;  xPos[Index( 128 ,1,1)]= 1.56422 ;  xPos[Index( 128 ,1,2)]= -1.29163 ;

xPos[Index( 129 ,0,0)]= 0.0916085 ;  xPos[Index( 129 ,0,1)]= 0.242811 ;  xPos[Index( 129 ,0,2)]= -0.289402 ;
xPos[Index( 129 ,1,0)]= -0.0916085 ;  xPos[Index( 129 ,1,1)]= -0.242811 ;  xPos[Index( 129 ,1,2)]= 0.289402 ;

xPos[Index( 130 ,0,0)]= 0.411882 ;  xPos[Index( 130 ,0,1)]= 0.429862 ;  xPos[Index( 130 ,0,2)]= -0.0959314 ;
xPos[Index( 130 ,1,0)]= -0.411882 ;  xPos[Index( 130 ,1,1)]= -0.429862 ;  xPos[Index( 130 ,1,2)]= 0.0959314 ;

xPos[Index( 131 ,0,0)]= -0.223265 ;  xPos[Index( 131 ,0,1)]= -0.0808424 ;  xPos[Index( 131 ,0,2)]= 0.3096 ;
xPos[Index( 131 ,1,0)]= 0.223265 ;  xPos[Index( 131 ,1,1)]= 0.0808424 ;  xPos[Index( 131 ,1,2)]= -0.3096 ;

xPos[Index( 132 ,0,0)]= 0.188764 ;  xPos[Index( 132 ,0,1)]= -0.316927 ;  xPos[Index( 132 ,0,2)]= 0.618108 ;
xPos[Index( 132 ,1,0)]= -0.188764 ;  xPos[Index( 132 ,1,1)]= 0.316927 ;  xPos[Index( 132 ,1,2)]= -0.618108 ;

xPos[Index( 133 ,0,0)]= 0.927329 ;  xPos[Index( 133 ,0,1)]= 0.069507 ;  xPos[Index( 133 ,0,2)]= -0.846125 ;
xPos[Index( 133 ,1,0)]= -0.927329 ;  xPos[Index( 133 ,1,1)]= -0.069507 ;  xPos[Index( 133 ,1,2)]= 0.846125 ;

xPos[Index( 134 ,0,0)]= -0.239547 ;  xPos[Index( 134 ,0,1)]= -1.4106 ;  xPos[Index( 134 ,0,2)]= 1.03836 ;
xPos[Index( 134 ,1,0)]= 0.239547 ;  xPos[Index( 134 ,1,1)]= 1.4106 ;  xPos[Index( 134 ,1,2)]= -1.03836 ;

xPos[Index( 135 ,0,0)]= 1.30851 ;  xPos[Index( 135 ,0,1)]= -0.0591983 ;  xPos[Index( 135 ,0,2)]= -2.06177 ;
xPos[Index( 135 ,1,0)]= -1.30851 ;  xPos[Index( 135 ,1,1)]= 0.0591983 ;  xPos[Index( 135 ,1,2)]= 2.06177 ;

xPos[Index( 136 ,0,0)]= 0.390295 ;  xPos[Index( 136 ,0,1)]= 0.365763 ;  xPos[Index( 136 ,0,2)]= 1.06244 ;
xPos[Index( 136 ,1,0)]= -0.390295 ;  xPos[Index( 136 ,1,1)]= -0.365763 ;  xPos[Index( 136 ,1,2)]= -1.06244 ;

xPos[Index( 137 ,0,0)]= 1.65678 ;  xPos[Index( 137 ,0,1)]= -0.159712 ;  xPos[Index( 137 ,0,2)]= -2.45234 ;
xPos[Index( 137 ,1,0)]= -1.65678 ;  xPos[Index( 137 ,1,1)]= 0.159712 ;  xPos[Index( 137 ,1,2)]= 2.45234 ;

xPos[Index( 138 ,0,0)]= 0.0801618 ;  xPos[Index( 138 ,0,1)]= 0.536841 ;  xPos[Index( 138 ,0,2)]= -1.0002 ;
xPos[Index( 138 ,1,0)]= -0.0801618 ;  xPos[Index( 138 ,1,1)]= -0.536841 ;  xPos[Index( 138 ,1,2)]= 1.0002 ;

xPos[Index( 139 ,0,0)]= -0.368529 ;  xPos[Index( 139 ,0,1)]= 0.618203 ;  xPos[Index( 139 ,0,2)]= -0.382807 ;
xPos[Index( 139 ,1,0)]= 0.368529 ;  xPos[Index( 139 ,1,1)]= -0.618203 ;  xPos[Index( 139 ,1,2)]= 0.382807 ;

xPos[Index( 140 ,0,0)]= -0.281971 ;  xPos[Index( 140 ,0,1)]= -0.0959094 ;  xPos[Index( 140 ,0,2)]= 0.30254 ;
xPos[Index( 140 ,1,0)]= 0.281971 ;  xPos[Index( 140 ,1,1)]= 0.0959094 ;  xPos[Index( 140 ,1,2)]= -0.30254 ;

xPos[Index( 141 ,0,0)]= 0.48177 ;  xPos[Index( 141 ,0,1)]= 0.69277 ;  xPos[Index( 141 ,0,2)]= 0.73907 ;
xPos[Index( 141 ,1,0)]= -0.48177 ;  xPos[Index( 141 ,1,1)]= -0.69277 ;  xPos[Index( 141 ,1,2)]= -0.73907 ;

xPos[Index( 142 ,0,0)]= -1.13262 ;  xPos[Index( 142 ,0,1)]= -0.435083 ;  xPos[Index( 142 ,0,2)]= -0.523314 ;
xPos[Index( 142 ,1,0)]= 1.13262 ;  xPos[Index( 142 ,1,1)]= 0.435083 ;  xPos[Index( 142 ,1,2)]= 0.523314 ;

xPos[Index( 143 ,0,0)]= -0.0480551 ;  xPos[Index( 143 ,0,1)]= 0.216169 ;  xPos[Index( 143 ,0,2)]= -0.714194 ;
xPos[Index( 143 ,1,0)]= 0.0480551 ;  xPos[Index( 143 ,1,1)]= -0.216169 ;  xPos[Index( 143 ,1,2)]= 0.714194 ;

xPos[Index( 144 ,0,0)]= 0.0717889 ;  xPos[Index( 144 ,0,1)]= 0.242582 ;  xPos[Index( 144 ,0,2)]= -0.0080524 ;
xPos[Index( 144 ,1,0)]= -0.0717889 ;  xPos[Index( 144 ,1,1)]= -0.242582 ;  xPos[Index( 144 ,1,2)]= 0.0080524 ;

xPos[Index( 145 ,0,0)]= -0.536269 ;  xPos[Index( 145 ,0,1)]= 0.115084 ;  xPos[Index( 145 ,0,2)]= 0.0933168 ;
xPos[Index( 145 ,1,0)]= 0.536269 ;  xPos[Index( 145 ,1,1)]= -0.115084 ;  xPos[Index( 145 ,1,2)]= -0.0933168 ;

xPos[Index( 146 ,0,0)]= -0.891798 ;  xPos[Index( 146 ,0,1)]= -1.10243 ;  xPos[Index( 146 ,0,2)]= 0.014153 ;
xPos[Index( 146 ,1,0)]= 0.891798 ;  xPos[Index( 146 ,1,1)]= 1.10243 ;  xPos[Index( 146 ,1,2)]= -0.014153 ;

xPos[Index( 147 ,0,0)]= 1.41352 ;  xPos[Index( 147 ,0,1)]= 0.334181 ;  xPos[Index( 147 ,0,2)]= 1.05089 ;
xPos[Index( 147 ,1,0)]= -1.41352 ;  xPos[Index( 147 ,1,1)]= -0.334181 ;  xPos[Index( 147 ,1,2)]= -1.05089 ;

xPos[Index( 148 ,0,0)]= -0.820338 ;  xPos[Index( 148 ,0,1)]= 0.00338418 ;  xPos[Index( 148 ,0,2)]= 0.879796 ;
xPos[Index( 148 ,1,0)]= 0.820338 ;  xPos[Index( 148 ,1,1)]= -0.00338418 ;  xPos[Index( 148 ,1,2)]= -0.879796 ;

xPos[Index( 149 ,0,0)]= 0.328 ;  xPos[Index( 149 ,0,1)]= 1.7887 ;  xPos[Index( 149 ,0,2)]= -0.31815 ;
xPos[Index( 149 ,1,0)]= -0.328 ;  xPos[Index( 149 ,1,1)]= -1.7887 ;  xPos[Index( 149 ,1,2)]= 0.31815 ;

xPos[Index( 150 ,0,0)]= -0.0736092 ;  xPos[Index( 150 ,0,1)]= -0.459078 ;  xPos[Index( 150 ,0,2)]= -0.645065 ;
xPos[Index( 150 ,1,0)]= 0.0736092 ;  xPos[Index( 150 ,1,1)]= 0.459078 ;  xPos[Index( 150 ,1,2)]= 0.645065 ;

xPos[Index( 151 ,0,0)]= -0.295782 ;  xPos[Index( 151 ,0,1)]= 0.855101 ;  xPos[Index( 151 ,0,2)]= -0.109701 ;
xPos[Index( 151 ,1,0)]= 0.295782 ;  xPos[Index( 151 ,1,1)]= -0.855101 ;  xPos[Index( 151 ,1,2)]= 0.109701 ;

xPos[Index( 152 ,0,0)]= 1.46014 ;  xPos[Index( 152 ,0,1)]= 0.19623 ;  xPos[Index( 152 ,0,2)]= 1.56914 ;
xPos[Index( 152 ,1,0)]= -1.46014 ;  xPos[Index( 152 ,1,1)]= -0.19623 ;  xPos[Index( 152 ,1,2)]= -1.56914 ;

xPos[Index( 153 ,0,0)]= -0.342913 ;  xPos[Index( 153 ,0,1)]= -0.473827 ;  xPos[Index( 153 ,0,2)]= -1.3995 ;
xPos[Index( 153 ,1,0)]= 0.342913 ;  xPos[Index( 153 ,1,1)]= 0.473827 ;  xPos[Index( 153 ,1,2)]= 1.3995 ;

xPos[Index( 154 ,0,0)]= 0.168011 ;  xPos[Index( 154 ,0,1)]= -0.0495101 ;  xPos[Index( 154 ,0,2)]= 0.32221 ;
xPos[Index( 154 ,1,0)]= -0.168011 ;  xPos[Index( 154 ,1,1)]= 0.0495101 ;  xPos[Index( 154 ,1,2)]= -0.32221 ;

xPos[Index( 155 ,0,0)]= -0.172261 ;  xPos[Index( 155 ,0,1)]= -0.840007 ;  xPos[Index( 155 ,0,2)]= -0.685784 ;
xPos[Index( 155 ,1,0)]= 0.172261 ;  xPos[Index( 155 ,1,1)]= 0.840007 ;  xPos[Index( 155 ,1,2)]= 0.685784 ;

xPos[Index( 156 ,0,0)]= 1.82717 ;  xPos[Index( 156 ,0,1)]= 0.317694 ;  xPos[Index( 156 ,0,2)]= -2.23187 ;
xPos[Index( 156 ,1,0)]= -1.82717 ;  xPos[Index( 156 ,1,1)]= -0.317694 ;  xPos[Index( 156 ,1,2)]= 2.23187 ;

xPos[Index( 157 ,0,0)]= -0.180038 ;  xPos[Index( 157 ,0,1)]= 0.559712 ;  xPos[Index( 157 ,0,2)]= 0.0161275 ;
xPos[Index( 157 ,1,0)]= 0.180038 ;  xPos[Index( 157 ,1,1)]= -0.559712 ;  xPos[Index( 157 ,1,2)]= -0.0161275 ;

xPos[Index( 158 ,0,0)]= -0.314941 ;  xPos[Index( 158 ,0,1)]= 0.104045 ;  xPos[Index( 158 ,0,2)]= 1.40198 ;
xPos[Index( 158 ,1,0)]= 0.314941 ;  xPos[Index( 158 ,1,1)]= -0.104045 ;  xPos[Index( 158 ,1,2)]= -1.40198 ;

xPos[Index( 159 ,0,0)]= 0.558356 ;  xPos[Index( 159 ,0,1)]= 0.837871 ;  xPos[Index( 159 ,0,2)]= 0.870986 ;
xPos[Index( 159 ,1,0)]= -0.558356 ;  xPos[Index( 159 ,1,1)]= -0.837871 ;  xPos[Index( 159 ,1,2)]= -0.870986 ;

xPos[Index( 160 ,0,0)]= -1.90971 ;  xPos[Index( 160 ,0,1)]= -0.440465 ;  xPos[Index( 160 ,0,2)]= 0.0305316 ;
xPos[Index( 160 ,1,0)]= 1.90971 ;  xPos[Index( 160 ,1,1)]= 0.440465 ;  xPos[Index( 160 ,1,2)]= -0.0305316 ;

xPos[Index( 161 ,0,0)]= -0.205637 ;  xPos[Index( 161 ,0,1)]= -0.31779 ;  xPos[Index( 161 ,0,2)]= -0.623123 ;
xPos[Index( 161 ,1,0)]= 0.205637 ;  xPos[Index( 161 ,1,1)]= 0.31779 ;  xPos[Index( 161 ,1,2)]= 0.623123 ;

xPos[Index( 162 ,0,0)]= 0.160694 ;  xPos[Index( 162 ,0,1)]= -0.337381 ;  xPos[Index( 162 ,0,2)]= -0.430108 ;
xPos[Index( 162 ,1,0)]= -0.160694 ;  xPos[Index( 162 ,1,1)]= 0.337381 ;  xPos[Index( 162 ,1,2)]= 0.430108 ;

xPos[Index( 163 ,0,0)]= -1.4086 ;  xPos[Index( 163 ,0,1)]= 0.174954 ;  xPos[Index( 163 ,0,2)]= -0.866681 ;
xPos[Index( 163 ,1,0)]= 1.4086 ;  xPos[Index( 163 ,1,1)]= -0.174954 ;  xPos[Index( 163 ,1,2)]= 0.866681 ;

xPos[Index( 164 ,0,0)]= 0.239763 ;  xPos[Index( 164 ,0,1)]= -0.592994 ;  xPos[Index( 164 ,0,2)]= -0.836455 ;
xPos[Index( 164 ,1,0)]= -0.239763 ;  xPos[Index( 164 ,1,1)]= 0.592994 ;  xPos[Index( 164 ,1,2)]= 0.836455 ;

xPos[Index( 165 ,0,0)]= 0.657406 ;  xPos[Index( 165 ,0,1)]= 0.767573 ;  xPos[Index( 165 ,0,2)]= -1.10383 ;
xPos[Index( 165 ,1,0)]= -0.657406 ;  xPos[Index( 165 ,1,1)]= -0.767573 ;  xPos[Index( 165 ,1,2)]= 1.10383 ;

xPos[Index( 166 ,0,0)]= -0.606579 ;  xPos[Index( 166 ,0,1)]= -0.445663 ;  xPos[Index( 166 ,0,2)]= 0.577971 ;
xPos[Index( 166 ,1,0)]= 0.606579 ;  xPos[Index( 166 ,1,1)]= 0.445663 ;  xPos[Index( 166 ,1,2)]= -0.577971 ;

xPos[Index( 167 ,0,0)]= -0.373149 ;  xPos[Index( 167 ,0,1)]= 0.45656 ;  xPos[Index( 167 ,0,2)]= 1.08166 ;
xPos[Index( 167 ,1,0)]= 0.373149 ;  xPos[Index( 167 ,1,1)]= -0.45656 ;  xPos[Index( 167 ,1,2)]= -1.08166 ;

xPos[Index( 168 ,0,0)]= 0.382008 ;  xPos[Index( 168 ,0,1)]= -0.314206 ;  xPos[Index( 168 ,0,2)]= -0.771417 ;
xPos[Index( 168 ,1,0)]= -0.382008 ;  xPos[Index( 168 ,1,1)]= 0.314206 ;  xPos[Index( 168 ,1,2)]= 0.771417 ;

xPos[Index( 169 ,0,0)]= 0.581466 ;  xPos[Index( 169 ,0,1)]= 0.129941 ;  xPos[Index( 169 ,0,2)]= 1.13134 ;
xPos[Index( 169 ,1,0)]= -0.581466 ;  xPos[Index( 169 ,1,1)]= -0.129941 ;  xPos[Index( 169 ,1,2)]= -1.13134 ;

xPos[Index( 170 ,0,0)]= 0.88399 ;  xPos[Index( 170 ,0,1)]= 0.314726 ;  xPos[Index( 170 ,0,2)]= 0.218594 ;
xPos[Index( 170 ,1,0)]= -0.88399 ;  xPos[Index( 170 ,1,1)]= -0.314726 ;  xPos[Index( 170 ,1,2)]= -0.218594 ;

xPos[Index( 171 ,0,0)]= -1.07205 ;  xPos[Index( 171 ,0,1)]= 0.473623 ;  xPos[Index( 171 ,0,2)]= -0.710723 ;
xPos[Index( 171 ,1,0)]= 1.07205 ;  xPos[Index( 171 ,1,1)]= -0.473623 ;  xPos[Index( 171 ,1,2)]= 0.710723 ;

xPos[Index( 172 ,0,0)]= -0.987834 ;  xPos[Index( 172 ,0,1)]= 1.4792 ;  xPos[Index( 172 ,0,2)]= 0.547533 ;
xPos[Index( 172 ,1,0)]= 0.987834 ;  xPos[Index( 172 ,1,1)]= -1.4792 ;  xPos[Index( 172 ,1,2)]= -0.547533 ;

xPos[Index( 173 ,0,0)]= 0.252328 ;  xPos[Index( 173 ,0,1)]= 1.38134 ;  xPos[Index( 173 ,0,2)]= -0.0847939 ;
xPos[Index( 173 ,1,0)]= -0.252328 ;  xPos[Index( 173 ,1,1)]= -1.38134 ;  xPos[Index( 173 ,1,2)]= 0.0847939 ;

xPos[Index( 174 ,0,0)]= -0.497142 ;  xPos[Index( 174 ,0,1)]= 0.166167 ;  xPos[Index( 174 ,0,2)]= 0.0344364 ;
xPos[Index( 174 ,1,0)]= 0.497142 ;  xPos[Index( 174 ,1,1)]= -0.166167 ;  xPos[Index( 174 ,1,2)]= -0.0344364 ;

xPos[Index( 175 ,0,0)]= 1.53432 ;  xPos[Index( 175 ,0,1)]= 2.13356 ;  xPos[Index( 175 ,0,2)]= 0.449913 ;
xPos[Index( 175 ,1,0)]= -1.53432 ;  xPos[Index( 175 ,1,1)]= -2.13356 ;  xPos[Index( 175 ,1,2)]= -0.449913 ;

xPos[Index( 176 ,0,0)]= -1.23744 ;  xPos[Index( 176 ,0,1)]= -0.411552 ;  xPos[Index( 176 ,0,2)]= -0.28959 ;
xPos[Index( 176 ,1,0)]= 1.23744 ;  xPos[Index( 176 ,1,1)]= 0.411552 ;  xPos[Index( 176 ,1,2)]= 0.28959 ;

xPos[Index( 177 ,0,0)]= -0.246307 ;  xPos[Index( 177 ,0,1)]= -0.648532 ;  xPos[Index( 177 ,0,2)]= 0.542568 ;
xPos[Index( 177 ,1,0)]= 0.246307 ;  xPos[Index( 177 ,1,1)]= 0.648532 ;  xPos[Index( 177 ,1,2)]= -0.542568 ;

xPos[Index( 178 ,0,0)]= 0.0967455 ;  xPos[Index( 178 ,0,1)]= -1.38019 ;  xPos[Index( 178 ,0,2)]= -0.933828 ;
xPos[Index( 178 ,1,0)]= -0.0967455 ;  xPos[Index( 178 ,1,1)]= 1.38019 ;  xPos[Index( 178 ,1,2)]= 0.933828 ;

xPos[Index( 179 ,0,0)]= -1.35675 ;  xPos[Index( 179 ,0,1)]= 0.14512 ;  xPos[Index( 179 ,0,2)]= 1.44487 ;
xPos[Index( 179 ,1,0)]= 1.35675 ;  xPos[Index( 179 ,1,1)]= -0.14512 ;  xPos[Index( 179 ,1,2)]= -1.44487 ;

xPos[Index( 180 ,0,0)]= -0.340823 ;  xPos[Index( 180 ,0,1)]= 0.36281 ;  xPos[Index( 180 ,0,2)]= 0.457967 ;
xPos[Index( 180 ,1,0)]= 0.340823 ;  xPos[Index( 180 ,1,1)]= -0.36281 ;  xPos[Index( 180 ,1,2)]= -0.457967 ;

xPos[Index( 181 ,0,0)]= -0.638725 ;  xPos[Index( 181 ,0,1)]= 0.0175739 ;  xPos[Index( 181 ,0,2)]= -0.437794 ;
xPos[Index( 181 ,1,0)]= 0.638725 ;  xPos[Index( 181 ,1,1)]= -0.0175739 ;  xPos[Index( 181 ,1,2)]= 0.437794 ;

xPos[Index( 182 ,0,0)]= 0.381549 ;  xPos[Index( 182 ,0,1)]= -1.23059 ;  xPos[Index( 182 ,0,2)]= 1.20198 ;
xPos[Index( 182 ,1,0)]= -0.381549 ;  xPos[Index( 182 ,1,1)]= 1.23059 ;  xPos[Index( 182 ,1,2)]= -1.20198 ;

xPos[Index( 183 ,0,0)]= 0.588321 ;  xPos[Index( 183 ,0,1)]= -0.464379 ;  xPos[Index( 183 ,0,2)]= -0.504514 ;
xPos[Index( 183 ,1,0)]= -0.588321 ;  xPos[Index( 183 ,1,1)]= 0.464379 ;  xPos[Index( 183 ,1,2)]= 0.504514 ;

xPos[Index( 184 ,0,0)]= 1.01003 ;  xPos[Index( 184 ,0,1)]= -0.780854 ;  xPos[Index( 184 ,0,2)]= -0.577152 ;
xPos[Index( 184 ,1,0)]= -1.01003 ;  xPos[Index( 184 ,1,1)]= 0.780854 ;  xPos[Index( 184 ,1,2)]= 0.577152 ;

xPos[Index( 185 ,0,0)]= -0.23854 ;  xPos[Index( 185 ,0,1)]= 0.306711 ;  xPos[Index( 185 ,0,2)]= 1.20668 ;
xPos[Index( 185 ,1,0)]= 0.23854 ;  xPos[Index( 185 ,1,1)]= -0.306711 ;  xPos[Index( 185 ,1,2)]= -1.20668 ;

xPos[Index( 186 ,0,0)]= -0.533227 ;  xPos[Index( 186 ,0,1)]= -0.830859 ;  xPos[Index( 186 ,0,2)]= 0.963152 ;
xPos[Index( 186 ,1,0)]= 0.533227 ;  xPos[Index( 186 ,1,1)]= 0.830859 ;  xPos[Index( 186 ,1,2)]= -0.963152 ;

xPos[Index( 187 ,0,0)]= -0.5221 ;  xPos[Index( 187 ,0,1)]= 0.135643 ;  xPos[Index( 187 ,0,2)]= -0.299168 ;
xPos[Index( 187 ,1,0)]= 0.5221 ;  xPos[Index( 187 ,1,1)]= -0.135643 ;  xPos[Index( 187 ,1,2)]= 0.299168 ;

xPos[Index( 188 ,0,0)]= 0.174419 ;  xPos[Index( 188 ,0,1)]= 0.501992 ;  xPos[Index( 188 ,0,2)]= -0.0525333 ;
xPos[Index( 188 ,1,0)]= -0.174419 ;  xPos[Index( 188 ,1,1)]= -0.501992 ;  xPos[Index( 188 ,1,2)]= 0.0525333 ;

xPos[Index( 189 ,0,0)]= 0.298844 ;  xPos[Index( 189 ,0,1)]= 0.499356 ;  xPos[Index( 189 ,0,2)]= -0.0970517 ;
xPos[Index( 189 ,1,0)]= -0.298844 ;  xPos[Index( 189 ,1,1)]= -0.499356 ;  xPos[Index( 189 ,1,2)]= 0.0970517 ;

xPos[Index( 190 ,0,0)]= -0.597763 ;  xPos[Index( 190 ,0,1)]= -0.0634056 ;  xPos[Index( 190 ,0,2)]= 1.45875 ;
xPos[Index( 190 ,1,0)]= 0.597763 ;  xPos[Index( 190 ,1,1)]= 0.0634056 ;  xPos[Index( 190 ,1,2)]= -1.45875 ;

xPos[Index( 191 ,0,0)]= -0.480874 ;  xPos[Index( 191 ,0,1)]= -0.0471544 ;  xPos[Index( 191 ,0,2)]= 1.07712 ;
xPos[Index( 191 ,1,0)]= 0.480874 ;  xPos[Index( 191 ,1,1)]= 0.0471544 ;  xPos[Index( 191 ,1,2)]= -1.07712 ;

xPos[Index( 192 ,0,0)]= -1.50797 ;  xPos[Index( 192 ,0,1)]= -0.506832 ;  xPos[Index( 192 ,0,2)]= -0.645205 ;
xPos[Index( 192 ,1,0)]= 1.50797 ;  xPos[Index( 192 ,1,1)]= 0.506832 ;  xPos[Index( 192 ,1,2)]= 0.645205 ;

xPos[Index( 193 ,0,0)]= 1.52381 ;  xPos[Index( 193 ,0,1)]= 1.18396 ;  xPos[Index( 193 ,0,2)]= -0.000957711 ;
xPos[Index( 193 ,1,0)]= -1.52381 ;  xPos[Index( 193 ,1,1)]= -1.18396 ;  xPos[Index( 193 ,1,2)]= 0.000957711 ;

xPos[Index( 194 ,0,0)]= 1.3194 ;  xPos[Index( 194 ,0,1)]= 0.828849 ;  xPos[Index( 194 ,0,2)]= 0.152576 ;
xPos[Index( 194 ,1,0)]= -1.3194 ;  xPos[Index( 194 ,1,1)]= -0.828849 ;  xPos[Index( 194 ,1,2)]= -0.152576 ;

xPos[Index( 195 ,0,0)]= -1.22644 ;  xPos[Index( 195 ,0,1)]= -0.019649 ;  xPos[Index( 195 ,0,2)]= -0.765256 ;
xPos[Index( 195 ,1,0)]= 1.22644 ;  xPos[Index( 195 ,1,1)]= 0.019649 ;  xPos[Index( 195 ,1,2)]= 0.765256 ;

xPos[Index( 196 ,0,0)]= 0.246878 ;  xPos[Index( 196 ,0,1)]= -1.00511 ;  xPos[Index( 196 ,0,2)]= -0.417824 ;
xPos[Index( 196 ,1,0)]= -0.246878 ;  xPos[Index( 196 ,1,1)]= 1.00511 ;  xPos[Index( 196 ,1,2)]= 0.417824 ;

xPos[Index( 197 ,0,0)]= -0.644421 ;  xPos[Index( 197 ,0,1)]= -0.74787 ;  xPos[Index( 197 ,0,2)]= -0.151588 ;
xPos[Index( 197 ,1,0)]= 0.644421 ;  xPos[Index( 197 ,1,1)]= 0.74787 ;  xPos[Index( 197 ,1,2)]= 0.151588 ;

xPos[Index( 198 ,0,0)]= 0.31669 ;  xPos[Index( 198 ,0,1)]= 0.638645 ;  xPos[Index( 198 ,0,2)]= -0.955602 ;
xPos[Index( 198 ,1,0)]= -0.31669 ;  xPos[Index( 198 ,1,1)]= -0.638645 ;  xPos[Index( 198 ,1,2)]= 0.955602 ;

xPos[Index( 199 ,0,0)]= 1.02974 ;  xPos[Index( 199 ,0,1)]= -0.18135 ;  xPos[Index( 199 ,0,2)]= 0.619903 ;
xPos[Index( 199 ,1,0)]= -1.02974 ;  xPos[Index( 199 ,1,1)]= 0.18135 ;  xPos[Index( 199 ,1,2)]= -0.619903 ;

xPos[Index( 200 ,0,0)]= 0.951721 ;  xPos[Index( 200 ,0,1)]= -0.34253 ;  xPos[Index( 200 ,0,2)]= 0.656417 ;
xPos[Index( 200 ,1,0)]= -0.951721 ;  xPos[Index( 200 ,1,1)]= 0.34253 ;  xPos[Index( 200 ,1,2)]= -0.656417 ;

xPos[Index( 201 ,0,0)]= 0.475104 ;  xPos[Index( 201 ,0,1)]= 0.295366 ;  xPos[Index( 201 ,0,2)]= -0.527853 ;
xPos[Index( 201 ,1,0)]= -0.475104 ;  xPos[Index( 201 ,1,1)]= -0.295366 ;  xPos[Index( 201 ,1,2)]= 0.527853 ;

xPos[Index( 202 ,0,0)]= -0.133074 ;  xPos[Index( 202 ,0,1)]= 0.278929 ;  xPos[Index( 202 ,0,2)]= -0.313755 ;
xPos[Index( 202 ,1,0)]= 0.133074 ;  xPos[Index( 202 ,1,1)]= -0.278929 ;  xPos[Index( 202 ,1,2)]= 0.313755 ;

xPos[Index( 203 ,0,0)]= -0.832782 ;  xPos[Index( 203 ,0,1)]= -0.326962 ;  xPos[Index( 203 ,0,2)]= -0.375109 ;
xPos[Index( 203 ,1,0)]= 0.832782 ;  xPos[Index( 203 ,1,1)]= 0.326962 ;  xPos[Index( 203 ,1,2)]= 0.375109 ;

xPos[Index( 204 ,0,0)]= 0.358233 ;  xPos[Index( 204 ,0,1)]= -0.0760606 ;  xPos[Index( 204 ,0,2)]= -1.76625 ;
xPos[Index( 204 ,1,0)]= -0.358233 ;  xPos[Index( 204 ,1,1)]= 0.0760606 ;  xPos[Index( 204 ,1,2)]= 1.76625 ;

xPos[Index( 205 ,0,0)]= -0.399516 ;  xPos[Index( 205 ,0,1)]= 0.0163509 ;  xPos[Index( 205 ,0,2)]= -0.346884 ;
xPos[Index( 205 ,1,0)]= 0.399516 ;  xPos[Index( 205 ,1,1)]= -0.0163509 ;  xPos[Index( 205 ,1,2)]= 0.346884 ;

xPos[Index( 206 ,0,0)]= 0.253736 ;  xPos[Index( 206 ,0,1)]= -1.04656 ;  xPos[Index( 206 ,0,2)]= 0.472836 ;
xPos[Index( 206 ,1,0)]= -0.253736 ;  xPos[Index( 206 ,1,1)]= 1.04656 ;  xPos[Index( 206 ,1,2)]= -0.472836 ;

xPos[Index( 207 ,0,0)]= -0.330843 ;  xPos[Index( 207 ,0,1)]= 0.205363 ;  xPos[Index( 207 ,0,2)]= -0.676989 ;
xPos[Index( 207 ,1,0)]= 0.330843 ;  xPos[Index( 207 ,1,1)]= -0.205363 ;  xPos[Index( 207 ,1,2)]= 0.676989 ;

xPos[Index( 208 ,0,0)]= 0.0271521 ;  xPos[Index( 208 ,0,1)]= -0.116627 ;  xPos[Index( 208 ,0,2)]= -0.0246636 ;
xPos[Index( 208 ,1,0)]= -0.0271521 ;  xPos[Index( 208 ,1,1)]= 0.116627 ;  xPos[Index( 208 ,1,2)]= 0.0246636 ;

xPos[Index( 209 ,0,0)]= 0.610848 ;  xPos[Index( 209 ,0,1)]= 0.0314734 ;  xPos[Index( 209 ,0,2)]= 0.59849 ;
xPos[Index( 209 ,1,0)]= -0.610848 ;  xPos[Index( 209 ,1,1)]= -0.0314734 ;  xPos[Index( 209 ,1,2)]= -0.59849 ;

xPos[Index( 210 ,0,0)]= 1.36825 ;  xPos[Index( 210 ,0,1)]= -0.875633 ;  xPos[Index( 210 ,0,2)]= -0.980885 ;
xPos[Index( 210 ,1,0)]= -1.36825 ;  xPos[Index( 210 ,1,1)]= 0.875633 ;  xPos[Index( 210 ,1,2)]= 0.980885 ;

xPos[Index( 211 ,0,0)]= 0.325774 ;  xPos[Index( 211 ,0,1)]= 0.61012 ;  xPos[Index( 211 ,0,2)]= 0.127129 ;
xPos[Index( 211 ,1,0)]= -0.325774 ;  xPos[Index( 211 ,1,1)]= -0.61012 ;  xPos[Index( 211 ,1,2)]= -0.127129 ;

xPos[Index( 212 ,0,0)]= 0.632771 ;  xPos[Index( 212 ,0,1)]= 1.99634 ;  xPos[Index( 212 ,0,2)]= 0.615836 ;
xPos[Index( 212 ,1,0)]= -0.632771 ;  xPos[Index( 212 ,1,1)]= -1.99634 ;  xPos[Index( 212 ,1,2)]= -0.615836 ;

xPos[Index( 213 ,0,0)]= -1.5314 ;  xPos[Index( 213 ,0,1)]= -0.564689 ;  xPos[Index( 213 ,0,2)]= -1.19066 ;
xPos[Index( 213 ,1,0)]= 1.5314 ;  xPos[Index( 213 ,1,1)]= 0.564689 ;  xPos[Index( 213 ,1,2)]= 1.19066 ;

xPos[Index( 214 ,0,0)]= -1.14987 ;  xPos[Index( 214 ,0,1)]= -0.759065 ;  xPos[Index( 214 ,0,2)]= 1.0667 ;
xPos[Index( 214 ,1,0)]= 1.14987 ;  xPos[Index( 214 ,1,1)]= 0.759065 ;  xPos[Index( 214 ,1,2)]= -1.0667 ;

xPos[Index( 215 ,0,0)]= -0.91167 ;  xPos[Index( 215 ,0,1)]= -0.671614 ;  xPos[Index( 215 ,0,2)]= 0.828985 ;
xPos[Index( 215 ,1,0)]= 0.91167 ;  xPos[Index( 215 ,1,1)]= 0.671614 ;  xPos[Index( 215 ,1,2)]= -0.828985 ;

xPos[Index( 216 ,0,0)]= 1.34086 ;  xPos[Index( 216 ,0,1)]= 0.295342 ;  xPos[Index( 216 ,0,2)]= -0.575551 ;
xPos[Index( 216 ,1,0)]= -1.34086 ;  xPos[Index( 216 ,1,1)]= -0.295342 ;  xPos[Index( 216 ,1,2)]= 0.575551 ;

xPos[Index( 217 ,0,0)]= 0.785138 ;  xPos[Index( 217 ,0,1)]= -1.15464 ;  xPos[Index( 217 ,0,2)]= 0.459884 ;
xPos[Index( 217 ,1,0)]= -0.785138 ;  xPos[Index( 217 ,1,1)]= 1.15464 ;  xPos[Index( 217 ,1,2)]= -0.459884 ;

xPos[Index( 218 ,0,0)]= -0.223842 ;  xPos[Index( 218 ,0,1)]= -0.740073 ;  xPos[Index( 218 ,0,2)]= -0.230976 ;
xPos[Index( 218 ,1,0)]= 0.223842 ;  xPos[Index( 218 ,1,1)]= 0.740073 ;  xPos[Index( 218 ,1,2)]= 0.230976 ;

xPos[Index( 219 ,0,0)]= 0.119438 ;  xPos[Index( 219 ,0,1)]= 1.49155 ;  xPos[Index( 219 ,0,2)]= -0.788663 ;
xPos[Index( 219 ,1,0)]= -0.119438 ;  xPos[Index( 219 ,1,1)]= -1.49155 ;  xPos[Index( 219 ,1,2)]= 0.788663 ;

xPos[Index( 220 ,0,0)]= 0.388995 ;  xPos[Index( 220 ,0,1)]= 0.628985 ;  xPos[Index( 220 ,0,2)]= -1.03101 ;
xPos[Index( 220 ,1,0)]= -0.388995 ;  xPos[Index( 220 ,1,1)]= -0.628985 ;  xPos[Index( 220 ,1,2)]= 1.03101 ;

xPos[Index( 221 ,0,0)]= 0.718465 ;  xPos[Index( 221 ,0,1)]= 0.0667586 ;  xPos[Index( 221 ,0,2)]= -0.58562 ;
xPos[Index( 221 ,1,0)]= -0.718465 ;  xPos[Index( 221 ,1,1)]= -0.0667586 ;  xPos[Index( 221 ,1,2)]= 0.58562 ;

xPos[Index( 222 ,0,0)]= 0.207679 ;  xPos[Index( 222 ,0,1)]= 0.172255 ;  xPos[Index( 222 ,0,2)]= 0.768574 ;
xPos[Index( 222 ,1,0)]= -0.207679 ;  xPos[Index( 222 ,1,1)]= -0.172255 ;  xPos[Index( 222 ,1,2)]= -0.768574 ;

xPos[Index( 223 ,0,0)]= 0.0145298 ;  xPos[Index( 223 ,0,1)]= -0.816105 ;  xPos[Index( 223 ,0,2)]= -0.705167 ;
xPos[Index( 223 ,1,0)]= -0.0145298 ;  xPos[Index( 223 ,1,1)]= 0.816105 ;  xPos[Index( 223 ,1,2)]= 0.705167 ;

xPos[Index( 224 ,0,0)]= -1.67793 ;  xPos[Index( 224 ,0,1)]= -0.292148 ;  xPos[Index( 224 ,0,2)]= -0.553029 ;
xPos[Index( 224 ,1,0)]= 1.67793 ;  xPos[Index( 224 ,1,1)]= 0.292148 ;  xPos[Index( 224 ,1,2)]= 0.553029 ;

xPos[Index( 225 ,0,0)]= 0.370195 ;  xPos[Index( 225 ,0,1)]= 1.61751 ;  xPos[Index( 225 ,0,2)]= -0.833511 ;
xPos[Index( 225 ,1,0)]= -0.370195 ;  xPos[Index( 225 ,1,1)]= -1.61751 ;  xPos[Index( 225 ,1,2)]= 0.833511 ;

xPos[Index( 226 ,0,0)]= 0.832443 ;  xPos[Index( 226 ,0,1)]= 0.374054 ;  xPos[Index( 226 ,0,2)]= 0.306807 ;
xPos[Index( 226 ,1,0)]= -0.832443 ;  xPos[Index( 226 ,1,1)]= -0.374054 ;  xPos[Index( 226 ,1,2)]= -0.306807 ;

xPos[Index( 227 ,0,0)]= 1.15733 ;  xPos[Index( 227 ,0,1)]= 1.64813 ;  xPos[Index( 227 ,0,2)]= -1.59621 ;
xPos[Index( 227 ,1,0)]= -1.15733 ;  xPos[Index( 227 ,1,1)]= -1.64813 ;  xPos[Index( 227 ,1,2)]= 1.59621 ;

xPos[Index( 228 ,0,0)]= 0.876629 ;  xPos[Index( 228 ,0,1)]= 1.77066 ;  xPos[Index( 228 ,0,2)]= 0.222946 ;
xPos[Index( 228 ,1,0)]= -0.876629 ;  xPos[Index( 228 ,1,1)]= -1.77066 ;  xPos[Index( 228 ,1,2)]= -0.222946 ;

xPos[Index( 229 ,0,0)]= -0.170604 ;  xPos[Index( 229 ,0,1)]= -0.598126 ;  xPos[Index( 229 ,0,2)]= 0.858045 ;
xPos[Index( 229 ,1,0)]= 0.170604 ;  xPos[Index( 229 ,1,1)]= 0.598126 ;  xPos[Index( 229 ,1,2)]= -0.858045 ;

xPos[Index( 230 ,0,0)]= 0.206848 ;  xPos[Index( 230 ,0,1)]= -0.534066 ;  xPos[Index( 230 ,0,2)]= 1.20336 ;
xPos[Index( 230 ,1,0)]= -0.206848 ;  xPos[Index( 230 ,1,1)]= 0.534066 ;  xPos[Index( 230 ,1,2)]= -1.20336 ;

xPos[Index( 231 ,0,0)]= 0.71849 ;  xPos[Index( 231 ,0,1)]= 1.85836 ;  xPos[Index( 231 ,0,2)]= 0.282124 ;
xPos[Index( 231 ,1,0)]= -0.71849 ;  xPos[Index( 231 ,1,1)]= -1.85836 ;  xPos[Index( 231 ,1,2)]= -0.282124 ;

xPos[Index( 232 ,0,0)]= -0.705981 ;  xPos[Index( 232 ,0,1)]= -1.33981 ;  xPos[Index( 232 ,0,2)]= -0.199673 ;
xPos[Index( 232 ,1,0)]= 0.705981 ;  xPos[Index( 232 ,1,1)]= 1.33981 ;  xPos[Index( 232 ,1,2)]= 0.199673 ;

xPos[Index( 233 ,0,0)]= 0.705137 ;  xPos[Index( 233 ,0,1)]= -0.8363 ;  xPos[Index( 233 ,0,2)]= -1.33285 ;
xPos[Index( 233 ,1,0)]= -0.705137 ;  xPos[Index( 233 ,1,1)]= 0.8363 ;  xPos[Index( 233 ,1,2)]= 1.33285 ;

xPos[Index( 234 ,0,0)]= 1.06549 ;  xPos[Index( 234 ,0,1)]= -1.56453 ;  xPos[Index( 234 ,0,2)]= -0.362061 ;
xPos[Index( 234 ,1,0)]= -1.06549 ;  xPos[Index( 234 ,1,1)]= 1.56453 ;  xPos[Index( 234 ,1,2)]= 0.362061 ;

xPos[Index( 235 ,0,0)]= -0.106038 ;  xPos[Index( 235 ,0,1)]= -0.349049 ;  xPos[Index( 235 ,0,2)]= -0.454154 ;
xPos[Index( 235 ,1,0)]= 0.106038 ;  xPos[Index( 235 ,1,1)]= 0.349049 ;  xPos[Index( 235 ,1,2)]= 0.454154 ;

xPos[Index( 236 ,0,0)]= -1.9951 ;  xPos[Index( 236 ,0,1)]= -0.386183 ;  xPos[Index( 236 ,0,2)]= -0.576074 ;
xPos[Index( 236 ,1,0)]= 1.9951 ;  xPos[Index( 236 ,1,1)]= 0.386183 ;  xPos[Index( 236 ,1,2)]= 0.576074 ;

xPos[Index( 237 ,0,0)]= -1.90986 ;  xPos[Index( 237 ,0,1)]= 0.0365954 ;  xPos[Index( 237 ,0,2)]= 0.360527 ;
xPos[Index( 237 ,1,0)]= 1.90986 ;  xPos[Index( 237 ,1,1)]= -0.0365954 ;  xPos[Index( 237 ,1,2)]= -0.360527 ;

xPos[Index( 238 ,0,0)]= 0.371409 ;  xPos[Index( 238 ,0,1)]= -0.08474 ;  xPos[Index( 238 ,0,2)]= -0.199181 ;
xPos[Index( 238 ,1,0)]= -0.371409 ;  xPos[Index( 238 ,1,1)]= 0.08474 ;  xPos[Index( 238 ,1,2)]= 0.199181 ;

xPos[Index( 239 ,0,0)]= 0.55845 ;  xPos[Index( 239 ,0,1)]= -0.561699 ;  xPos[Index( 239 ,0,2)]= -0.469898 ;
xPos[Index( 239 ,1,0)]= -0.55845 ;  xPos[Index( 239 ,1,1)]= 0.561699 ;  xPos[Index( 239 ,1,2)]= 0.469898 ;

xPos[Index( 240 ,0,0)]= 0.251528 ;  xPos[Index( 240 ,0,1)]= -0.706647 ;  xPos[Index( 240 ,0,2)]= 0.0968944 ;
xPos[Index( 240 ,1,0)]= -0.251528 ;  xPos[Index( 240 ,1,1)]= 0.706647 ;  xPos[Index( 240 ,1,2)]= -0.0968944 ;

xPos[Index( 241 ,0,0)]= -0.465533 ;  xPos[Index( 241 ,0,1)]= 0.230659 ;  xPos[Index( 241 ,0,2)]= -1.38663 ;
xPos[Index( 241 ,1,0)]= 0.465533 ;  xPos[Index( 241 ,1,1)]= -0.230659 ;  xPos[Index( 241 ,1,2)]= 1.38663 ;

xPos[Index( 242 ,0,0)]= 0.282631 ;  xPos[Index( 242 ,0,1)]= 0.0520884 ;  xPos[Index( 242 ,0,2)]= 1.27916 ;
xPos[Index( 242 ,1,0)]= -0.282631 ;  xPos[Index( 242 ,1,1)]= -0.0520884 ;  xPos[Index( 242 ,1,2)]= -1.27916 ;

xPos[Index( 243 ,0,0)]= -0.192481 ;  xPos[Index( 243 ,0,1)]= -0.0417057 ;  xPos[Index( 243 ,0,2)]= -0.791409 ;
xPos[Index( 243 ,1,0)]= 0.192481 ;  xPos[Index( 243 ,1,1)]= 0.0417057 ;  xPos[Index( 243 ,1,2)]= 0.791409 ;

xPos[Index( 244 ,0,0)]= 0.0692963 ;  xPos[Index( 244 ,0,1)]= -0.327577 ;  xPos[Index( 244 ,0,2)]= 0.203129 ;
xPos[Index( 244 ,1,0)]= -0.0692963 ;  xPos[Index( 244 ,1,1)]= 0.327577 ;  xPos[Index( 244 ,1,2)]= -0.203129 ;

xPos[Index( 245 ,0,0)]= -0.534298 ;  xPos[Index( 245 ,0,1)]= 0.133085 ;  xPos[Index( 245 ,0,2)]= 0.789352 ;
xPos[Index( 245 ,1,0)]= 0.534298 ;  xPos[Index( 245 ,1,1)]= -0.133085 ;  xPos[Index( 245 ,1,2)]= -0.789352 ;

xPos[Index( 246 ,0,0)]= -0.189469 ;  xPos[Index( 246 ,0,1)]= -0.293385 ;  xPos[Index( 246 ,0,2)]= -0.0316542 ;
xPos[Index( 246 ,1,0)]= 0.189469 ;  xPos[Index( 246 ,1,1)]= 0.293385 ;  xPos[Index( 246 ,1,2)]= 0.0316542 ;

xPos[Index( 247 ,0,0)]= 0.316422 ;  xPos[Index( 247 ,0,1)]= -0.418745 ;  xPos[Index( 247 ,0,2)]= 0.732962 ;
xPos[Index( 247 ,1,0)]= -0.316422 ;  xPos[Index( 247 ,1,1)]= 0.418745 ;  xPos[Index( 247 ,1,2)]= -0.732962 ;

xPos[Index( 248 ,0,0)]= 0.0187916 ;  xPos[Index( 248 ,0,1)]= 0.436296 ;  xPos[Index( 248 ,0,2)]= -0.371419 ;
xPos[Index( 248 ,1,0)]= -0.0187916 ;  xPos[Index( 248 ,1,1)]= -0.436296 ;  xPos[Index( 248 ,1,2)]= 0.371419 ;

xPos[Index( 249 ,0,0)]= -1.18374 ;  xPos[Index( 249 ,0,1)]= 0.749769 ;  xPos[Index( 249 ,0,2)]= 0.707003 ;
xPos[Index( 249 ,1,0)]= 1.18374 ;  xPos[Index( 249 ,1,1)]= -0.749769 ;  xPos[Index( 249 ,1,2)]= -0.707003 ;

xPos[Index( 250 ,0,0)]= 0.993151 ;  xPos[Index( 250 ,0,1)]= -0.579624 ;  xPos[Index( 250 ,0,2)]= 0.624269 ;
xPos[Index( 250 ,1,0)]= -0.993151 ;  xPos[Index( 250 ,1,1)]= 0.579624 ;  xPos[Index( 250 ,1,2)]= -0.624269 ;

xPos[Index( 251 ,0,0)]= 0.340213 ;  xPos[Index( 251 ,0,1)]= -0.237026 ;  xPos[Index( 251 ,0,2)]= -0.0543208 ;
xPos[Index( 251 ,1,0)]= -0.340213 ;  xPos[Index( 251 ,1,1)]= 0.237026 ;  xPos[Index( 251 ,1,2)]= 0.0543208 ;

xPos[Index( 252 ,0,0)]= 0.500547 ;  xPos[Index( 252 ,0,1)]= -0.0913322 ;  xPos[Index( 252 ,0,2)]= -0.576757 ;
xPos[Index( 252 ,1,0)]= -0.500547 ;  xPos[Index( 252 ,1,1)]= 0.0913322 ;  xPos[Index( 252 ,1,2)]= 0.576757 ;

xPos[Index( 253 ,0,0)]= 1.05972 ;  xPos[Index( 253 ,0,1)]= 0.0130157 ;  xPos[Index( 253 ,0,2)]= 0.709781 ;
xPos[Index( 253 ,1,0)]= -1.05972 ;  xPos[Index( 253 ,1,1)]= -0.0130157 ;  xPos[Index( 253 ,1,2)]= -0.709781 ;

xPos[Index( 254 ,0,0)]= 1.074 ;  xPos[Index( 254 ,0,1)]= -0.715646 ;  xPos[Index( 254 ,0,2)]= 0.384867 ;
xPos[Index( 254 ,1,0)]= -1.074 ;  xPos[Index( 254 ,1,1)]= 0.715646 ;  xPos[Index( 254 ,1,2)]= -0.384867 ;

xPos[Index( 255 ,0,0)]= 0.357743 ;  xPos[Index( 255 ,0,1)]= 0.581659 ;  xPos[Index( 255 ,0,2)]= 0.0136315 ;
xPos[Index( 255 ,1,0)]= -0.357743 ;  xPos[Index( 255 ,1,1)]= -0.581659 ;  xPos[Index( 255 ,1,2)]= -0.0136315 ;

xPos[Index( 256 ,0,0)]= 0.699119 ;  xPos[Index( 256 ,0,1)]= 1.04708 ;  xPos[Index( 256 ,0,2)]= -0.420968 ;
xPos[Index( 256 ,1,0)]= -0.699119 ;  xPos[Index( 256 ,1,1)]= -1.04708 ;  xPos[Index( 256 ,1,2)]= 0.420968 ;

xPos[Index( 257 ,0,0)]= -2.86222 ;  xPos[Index( 257 ,0,1)]= -0.286125 ;  xPos[Index( 257 ,0,2)]= 0.412544 ;
xPos[Index( 257 ,1,0)]= 2.86222 ;  xPos[Index( 257 ,1,1)]= 0.286125 ;  xPos[Index( 257 ,1,2)]= -0.412544 ;

xPos[Index( 258 ,0,0)]= 0.400048 ;  xPos[Index( 258 ,0,1)]= 1.02611 ;  xPos[Index( 258 ,0,2)]= 1.75852 ;
xPos[Index( 258 ,1,0)]= -0.400048 ;  xPos[Index( 258 ,1,1)]= -1.02611 ;  xPos[Index( 258 ,1,2)]= -1.75852 ;

xPos[Index( 259 ,0,0)]= -0.199191 ;  xPos[Index( 259 ,0,1)]= 1.04292 ;  xPos[Index( 259 ,0,2)]= 0.0129566 ;
xPos[Index( 259 ,1,0)]= 0.199191 ;  xPos[Index( 259 ,1,1)]= -1.04292 ;  xPos[Index( 259 ,1,2)]= -0.0129566 ;

xPos[Index( 260 ,0,0)]= -0.747683 ;  xPos[Index( 260 ,0,1)]= 2.45087 ;  xPos[Index( 260 ,0,2)]= -0.326639 ;
xPos[Index( 260 ,1,0)]= 0.747683 ;  xPos[Index( 260 ,1,1)]= -2.45087 ;  xPos[Index( 260 ,1,2)]= 0.326639 ;

xPos[Index( 261 ,0,0)]= 0.00724341 ;  xPos[Index( 261 ,0,1)]= -0.08995 ;  xPos[Index( 261 ,0,2)]= -0.842643 ;
xPos[Index( 261 ,1,0)]= -0.00724341 ;  xPos[Index( 261 ,1,1)]= 0.08995 ;  xPos[Index( 261 ,1,2)]= 0.842643 ;

xPos[Index( 262 ,0,0)]= 0.363899 ;  xPos[Index( 262 ,0,1)]= 1.41888 ;  xPos[Index( 262 ,0,2)]= 0.129183 ;
xPos[Index( 262 ,1,0)]= -0.363899 ;  xPos[Index( 262 ,1,1)]= -1.41888 ;  xPos[Index( 262 ,1,2)]= -0.129183 ;

xPos[Index( 263 ,0,0)]= 1.40371 ;  xPos[Index( 263 ,0,1)]= -0.117491 ;  xPos[Index( 263 ,0,2)]= 0.809557 ;
xPos[Index( 263 ,1,0)]= -1.40371 ;  xPos[Index( 263 ,1,1)]= 0.117491 ;  xPos[Index( 263 ,1,2)]= -0.809557 ;

xPos[Index( 264 ,0,0)]= 0.30487 ;  xPos[Index( 264 ,0,1)]= 0.600831 ;  xPos[Index( 264 ,0,2)]= -0.226584 ;
xPos[Index( 264 ,1,0)]= -0.30487 ;  xPos[Index( 264 ,1,1)]= -0.600831 ;  xPos[Index( 264 ,1,2)]= 0.226584 ;

xPos[Index( 265 ,0,0)]= 0.371481 ;  xPos[Index( 265 ,0,1)]= -1.11034 ;  xPos[Index( 265 ,0,2)]= -0.208249 ;
xPos[Index( 265 ,1,0)]= -0.371481 ;  xPos[Index( 265 ,1,1)]= 1.11034 ;  xPos[Index( 265 ,1,2)]= 0.208249 ;

xPos[Index( 266 ,0,0)]= -0.0847622 ;  xPos[Index( 266 ,0,1)]= 1.7548 ;  xPos[Index( 266 ,0,2)]= 0.450134 ;
xPos[Index( 266 ,1,0)]= 0.0847622 ;  xPos[Index( 266 ,1,1)]= -1.7548 ;  xPos[Index( 266 ,1,2)]= -0.450134 ;

xPos[Index( 267 ,0,0)]= 0.02427 ;  xPos[Index( 267 ,0,1)]= -0.709892 ;  xPos[Index( 267 ,0,2)]= -0.725248 ;
xPos[Index( 267 ,1,0)]= -0.02427 ;  xPos[Index( 267 ,1,1)]= 0.709892 ;  xPos[Index( 267 ,1,2)]= 0.725248 ;

xPos[Index( 268 ,0,0)]= -0.392126 ;  xPos[Index( 268 ,0,1)]= -0.282515 ;  xPos[Index( 268 ,0,2)]= -1.00675 ;
xPos[Index( 268 ,1,0)]= 0.392126 ;  xPos[Index( 268 ,1,1)]= 0.282515 ;  xPos[Index( 268 ,1,2)]= 1.00675 ;

xPos[Index( 269 ,0,0)]= 1.13878 ;  xPos[Index( 269 ,0,1)]= -0.786207 ;  xPos[Index( 269 ,0,2)]= 1.10289 ;
xPos[Index( 269 ,1,0)]= -1.13878 ;  xPos[Index( 269 ,1,1)]= 0.786207 ;  xPos[Index( 269 ,1,2)]= -1.10289 ;

xPos[Index( 270 ,0,0)]= -0.0157036 ;  xPos[Index( 270 ,0,1)]= 0.492728 ;  xPos[Index( 270 ,0,2)]= 0.599167 ;
xPos[Index( 270 ,1,0)]= 0.0157036 ;  xPos[Index( 270 ,1,1)]= -0.492728 ;  xPos[Index( 270 ,1,2)]= -0.599167 ;

xPos[Index( 271 ,0,0)]= 0.64112 ;  xPos[Index( 271 ,0,1)]= -0.628808 ;  xPos[Index( 271 ,0,2)]= -1.83153 ;
xPos[Index( 271 ,1,0)]= -0.64112 ;  xPos[Index( 271 ,1,1)]= 0.628808 ;  xPos[Index( 271 ,1,2)]= 1.83153 ;

xPos[Index( 272 ,0,0)]= -1.45356 ;  xPos[Index( 272 ,0,1)]= -0.609152 ;  xPos[Index( 272 ,0,2)]= -1.08851 ;
xPos[Index( 272 ,1,0)]= 1.45356 ;  xPos[Index( 272 ,1,1)]= 0.609152 ;  xPos[Index( 272 ,1,2)]= 1.08851 ;

xPos[Index( 273 ,0,0)]= 0.718545 ;  xPos[Index( 273 ,0,1)]= 0.88381 ;  xPos[Index( 273 ,0,2)]= 0.776355 ;
xPos[Index( 273 ,1,0)]= -0.718545 ;  xPos[Index( 273 ,1,1)]= -0.88381 ;  xPos[Index( 273 ,1,2)]= -0.776355 ;

xPos[Index( 274 ,0,0)]= -1.26801 ;  xPos[Index( 274 ,0,1)]= 1.76036 ;  xPos[Index( 274 ,0,2)]= -0.701731 ;
xPos[Index( 274 ,1,0)]= 1.26801 ;  xPos[Index( 274 ,1,1)]= -1.76036 ;  xPos[Index( 274 ,1,2)]= 0.701731 ;

xPos[Index( 275 ,0,0)]= -0.935127 ;  xPos[Index( 275 ,0,1)]= 0.405728 ;  xPos[Index( 275 ,0,2)]= -0.534184 ;
xPos[Index( 275 ,1,0)]= 0.935127 ;  xPos[Index( 275 ,1,1)]= -0.405728 ;  xPos[Index( 275 ,1,2)]= 0.534184 ;

xPos[Index( 276 ,0,0)]= -0.232139 ;  xPos[Index( 276 ,0,1)]= -0.233082 ;  xPos[Index( 276 ,0,2)]= -0.343193 ;
xPos[Index( 276 ,1,0)]= 0.232139 ;  xPos[Index( 276 ,1,1)]= 0.233082 ;  xPos[Index( 276 ,1,2)]= 0.343193 ;

xPos[Index( 277 ,0,0)]= -1.9443 ;  xPos[Index( 277 ,0,1)]= 0.939606 ;  xPos[Index( 277 ,0,2)]= 0.165486 ;
xPos[Index( 277 ,1,0)]= 1.9443 ;  xPos[Index( 277 ,1,1)]= -0.939606 ;  xPos[Index( 277 ,1,2)]= -0.165486 ;

xPos[Index( 278 ,0,0)]= -0.230191 ;  xPos[Index( 278 ,0,1)]= 1.00457 ;  xPos[Index( 278 ,0,2)]= 0.910822 ;
xPos[Index( 278 ,1,0)]= 0.230191 ;  xPos[Index( 278 ,1,1)]= -1.00457 ;  xPos[Index( 278 ,1,2)]= -0.910822 ;

xPos[Index( 279 ,0,0)]= -0.209377 ;  xPos[Index( 279 ,0,1)]= -0.835872 ;  xPos[Index( 279 ,0,2)]= -0.184592 ;
xPos[Index( 279 ,1,0)]= 0.209377 ;  xPos[Index( 279 ,1,1)]= 0.835872 ;  xPos[Index( 279 ,1,2)]= 0.184592 ;

xPos[Index( 280 ,0,0)]= -0.169375 ;  xPos[Index( 280 ,0,1)]= -0.349862 ;  xPos[Index( 280 ,0,2)]= -0.167825 ;
xPos[Index( 280 ,1,0)]= 0.169375 ;  xPos[Index( 280 ,1,1)]= 0.349862 ;  xPos[Index( 280 ,1,2)]= 0.167825 ;

xPos[Index( 281 ,0,0)]= 0.285599 ;  xPos[Index( 281 ,0,1)]= 0.41747 ;  xPos[Index( 281 ,0,2)]= 0.617286 ;
xPos[Index( 281 ,1,0)]= -0.285599 ;  xPos[Index( 281 ,1,1)]= -0.41747 ;  xPos[Index( 281 ,1,2)]= -0.617286 ;

xPos[Index( 282 ,0,0)]= 0.887645 ;  xPos[Index( 282 ,0,1)]= -0.484872 ;  xPos[Index( 282 ,0,2)]= -3.50041 ;
xPos[Index( 282 ,1,0)]= -0.887645 ;  xPos[Index( 282 ,1,1)]= 0.484872 ;  xPos[Index( 282 ,1,2)]= 3.50041 ;

xPos[Index( 283 ,0,0)]= 1.09103 ;  xPos[Index( 283 ,0,1)]= 0.197623 ;  xPos[Index( 283 ,0,2)]= 0.676149 ;
xPos[Index( 283 ,1,0)]= -1.09103 ;  xPos[Index( 283 ,1,1)]= -0.197623 ;  xPos[Index( 283 ,1,2)]= -0.676149 ;

xPos[Index( 284 ,0,0)]= 0.134721 ;  xPos[Index( 284 ,0,1)]= 1.21679 ;  xPos[Index( 284 ,0,2)]= 0.98895 ;
xPos[Index( 284 ,1,0)]= -0.134721 ;  xPos[Index( 284 ,1,1)]= -1.21679 ;  xPos[Index( 284 ,1,2)]= -0.98895 ;

xPos[Index( 285 ,0,0)]= -0.107129 ;  xPos[Index( 285 ,0,1)]= -1.44548 ;  xPos[Index( 285 ,0,2)]= 1.61537 ;
xPos[Index( 285 ,1,0)]= 0.107129 ;  xPos[Index( 285 ,1,1)]= 1.44548 ;  xPos[Index( 285 ,1,2)]= -1.61537 ;

xPos[Index( 286 ,0,0)]= 0.183634 ;  xPos[Index( 286 ,0,1)]= -0.820375 ;  xPos[Index( 286 ,0,2)]= -0.973422 ;
xPos[Index( 286 ,1,0)]= -0.183634 ;  xPos[Index( 286 ,1,1)]= 0.820375 ;  xPos[Index( 286 ,1,2)]= 0.973422 ;

xPos[Index( 287 ,0,0)]= -0.490311 ;  xPos[Index( 287 ,0,1)]= 0.589601 ;  xPos[Index( 287 ,0,2)]= 0.00208554 ;
xPos[Index( 287 ,1,0)]= 0.490311 ;  xPos[Index( 287 ,1,1)]= -0.589601 ;  xPos[Index( 287 ,1,2)]= -0.00208554 ;

xPos[Index( 288 ,0,0)]= 0.63639 ;  xPos[Index( 288 ,0,1)]= 1.77633 ;  xPos[Index( 288 ,0,2)]= -0.0336308 ;
xPos[Index( 288 ,1,0)]= -0.63639 ;  xPos[Index( 288 ,1,1)]= -1.77633 ;  xPos[Index( 288 ,1,2)]= 0.0336308 ;

xPos[Index( 289 ,0,0)]= -0.466598 ;  xPos[Index( 289 ,0,1)]= 0.279133 ;  xPos[Index( 289 ,0,2)]= 1.08751 ;
xPos[Index( 289 ,1,0)]= 0.466598 ;  xPos[Index( 289 ,1,1)]= -0.279133 ;  xPos[Index( 289 ,1,2)]= -1.08751 ;

xPos[Index( 290 ,0,0)]= -0.0956189 ;  xPos[Index( 290 ,0,1)]= -2.18041 ;  xPos[Index( 290 ,0,2)]= 0.840182 ;
xPos[Index( 290 ,1,0)]= 0.0956189 ;  xPos[Index( 290 ,1,1)]= 2.18041 ;  xPos[Index( 290 ,1,2)]= -0.840182 ;

xPos[Index( 291 ,0,0)]= 1.137 ;  xPos[Index( 291 ,0,1)]= 0.346984 ;  xPos[Index( 291 ,0,2)]= 0.631821 ;
xPos[Index( 291 ,1,0)]= -1.137 ;  xPos[Index( 291 ,1,1)]= -0.346984 ;  xPos[Index( 291 ,1,2)]= -0.631821 ;

xPos[Index( 292 ,0,0)]= 0.805067 ;  xPos[Index( 292 ,0,1)]= -0.472507 ;  xPos[Index( 292 ,0,2)]= 0.221027 ;
xPos[Index( 292 ,1,0)]= -0.805067 ;  xPos[Index( 292 ,1,1)]= 0.472507 ;  xPos[Index( 292 ,1,2)]= -0.221027 ;

xPos[Index( 293 ,0,0)]= -1.55257 ;  xPos[Index( 293 ,0,1)]= -0.566448 ;  xPos[Index( 293 ,0,2)]= -0.357319 ;
xPos[Index( 293 ,1,0)]= 1.55257 ;  xPos[Index( 293 ,1,1)]= 0.566448 ;  xPos[Index( 293 ,1,2)]= 0.357319 ;

xPos[Index( 294 ,0,0)]= -0.114832 ;  xPos[Index( 294 ,0,1)]= 0.0526071 ;  xPos[Index( 294 ,0,2)]= -0.25937 ;
xPos[Index( 294 ,1,0)]= 0.114832 ;  xPos[Index( 294 ,1,1)]= -0.0526071 ;  xPos[Index( 294 ,1,2)]= 0.25937 ;

xPos[Index( 295 ,0,0)]= -0.0728052 ;  xPos[Index( 295 ,0,1)]= 0.058356 ;  xPos[Index( 295 ,0,2)]= 0.374395 ;
xPos[Index( 295 ,1,0)]= 0.0728052 ;  xPos[Index( 295 ,1,1)]= -0.058356 ;  xPos[Index( 295 ,1,2)]= -0.374395 ;

xPos[Index( 296 ,0,0)]= -0.653151 ;  xPos[Index( 296 ,0,1)]= -0.157136 ;  xPos[Index( 296 ,0,2)]= -1.41095 ;
xPos[Index( 296 ,1,0)]= 0.653151 ;  xPos[Index( 296 ,1,1)]= 0.157136 ;  xPos[Index( 296 ,1,2)]= 1.41095 ;

xPos[Index( 297 ,0,0)]= 1.32611 ;  xPos[Index( 297 ,0,1)]= -0.106171 ;  xPos[Index( 297 ,0,2)]= 0.446031 ;
xPos[Index( 297 ,1,0)]= -1.32611 ;  xPos[Index( 297 ,1,1)]= 0.106171 ;  xPos[Index( 297 ,1,2)]= -0.446031 ;

xPos[Index( 298 ,0,0)]= -0.82467 ;  xPos[Index( 298 ,0,1)]= 0.125206 ;  xPos[Index( 298 ,0,2)]= -0.373355 ;
xPos[Index( 298 ,1,0)]= 0.82467 ;  xPos[Index( 298 ,1,1)]= -0.125206 ;  xPos[Index( 298 ,1,2)]= 0.373355 ;

xPos[Index( 299 ,0,0)]= -0.446996 ;  xPos[Index( 299 ,0,1)]= -0.377854 ;  xPos[Index( 299 ,0,2)]= -0.599044 ;
xPos[Index( 299 ,1,0)]= 0.446996 ;  xPos[Index( 299 ,1,1)]= 0.377854 ;  xPos[Index( 299 ,1,2)]= 0.599044 ;

}

double uni_H2_rn(){return drand48();}

void GetNucleonPositions(double *x1,double *x2){
  int Sample=int(NumberOfSamples*uni_H2_rn());
  x1[0]=xPos[Index(Sample,0,0)];  x1[1]=xPos[Index(Sample,0,1)];  x1[2]=xPos[Index(Sample,0,2)];
  x2[0]=xPos[Index(Sample,1,0)];  x2[1]=xPos[Index(Sample,1,1)];  x2[2]=xPos[Index(Sample,1,2)];
}

}
