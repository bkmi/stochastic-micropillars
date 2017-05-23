#ifndef KATANA_GET_SPECS__
#define KATANA_GET_SPECS__


// hier muss wohl kein include rein =)
//#include <cmath>
//#include <iostream>
//#include <vector>
//#include <complex.h>
//#include "fftw3.h"
//#include <fstream>
#include <string>


namespace katana
{
using Time=double;
void get_power_spec(std::vector<double> &timeseries, std::vector<double> &powerspec);
void get_optical_spec(std::vector<std::complex<double>> &timeseries, std::vector<double> &spec);
void get_fourier_spec(std::vector<std::complex<double>> &timeseries, std::vector<std::complex<double>> &spec);

//Dump shit

void dump_power_spec(std::vector<double> &powerspec, Time dt, int maxN,  std::string fourieroutname); 
void dump_optical_spec(std::vector<double> &spec, Time dt, int maxN, std::string fourieroutname); 
void dump_fourer_spec(std::vector<std::complex<double>> &spec, Time dt, int maxN, std::string fourieroutname); 


}

#endif
