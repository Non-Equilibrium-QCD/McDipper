/* Copyright (c) 2023 Oscar Garcia-Montero
 * All rights reserved. */

#ifndef TABLE_IO_H
#define TABLE_IO_H
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

/* Shape checks for the tabulated inputs.
 *
 * Tables are written by one run and read back by a later one, often after
 * being merged by hand from per-rapidity-slice jobs. Three ways that goes
 * wrong used to be unreported: a missing file made fscanf dereference a null
 * FILE* and segfault, a truncated file gave a message naming neither the file
 * nor the expected size, and a file holding more rows than its grid declares
 * was read only as far as the grid and the rest silently dropped.
 *
 * The last one is the reason these exist. A merge that duplicates or omits a
 * slice leaves a table whose row count disagrees with the config.yaml sitting
 * next to it, and without a count the run just uses the wrong rows. */

namespace TableIO{

  inline FILE * open_or_exit(const std::string & path){
    FILE * f = fopen(path.c_str(),"r");
    if(f==NULL){
      std::cerr<<"[ Error::Table ]: Cannot open "<<path<<std::endl;
      std::cerr<<"                  Delete the set directory to force a rebuild."<<std::endl;
      exit(EXIT_FAILURE);
    }
    return f;
  }

  inline void short_table(const std::string & path, long row, long expected){
    std::cerr<<"[ Error::Table ]: "<<path<<" ran out after "<<row<<" of the "
             <<expected<<" rows its grid declares."<<std::endl;
    std::cerr<<"                  The table is truncated - most often the job writing it was killed."<<std::endl;
    std::cerr<<"                  Delete the set directory to force a rebuild."<<std::endl;
    exit(EXIT_FAILURE);
  }

  inline void assert_at_end(FILE * f, const std::string & path, long expected){
    double junk;
    if(fscanf(f,"%lf",&junk)==1){
      std::cerr<<"[ Error::Table ]: "<<path<<" holds more than the "<<expected
               <<" rows its grid declares."<<std::endl;
      std::cerr<<"                  If the table was merged from per-slice jobs, a slice is duplicated"<<std::endl;
      std::cerr<<"                  or the grid recorded in config.yaml does not match the data."<<std::endl;
      exit(EXIT_FAILURE);
    }
  }

}

#endif /* table_io.h */
