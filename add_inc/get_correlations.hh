#ifndef KATANA_GET_CORRELATIONS__
#define KATANA_GET_CORRELATIONS__


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


/*############EXPLANATION TO USAGE!##############
################################################# 
Die Funktionen get_autocorr sind eigentlich hidden und werden nicht benötigt.
Es kann eigentlich immer get_g2, oder get_g1 benutzt werden.
Die Funktion "normalize correlations" nützt dann, wenn die Funktionen nicht im
Sinne von Lasern benutzt werden.
get_auto(cross)corr benutzt werden, entspricht der Rückgabewert dem mean_value der (complexen) Zeitserie, (evtl auch 0!)
#################################################
################################################# 

*/


// Get the correlation functions and normalize
double get_autocorr(std::vector<std::complex<double>> &timeseries, std::vector<std::complex<double>> &autocorrelation);
double get_autocorr_real(std::vector<double> &timeseries, std::vector<double> &autocorrelation);
double get_autocorr(std::vector<double> &timeseries, std::vector<double> &autocorrelation); //Overload of get_autocorr with exact definition of get_autocorr_real 

std::vector<double> get_crosscorr(std::vector<std::complex<double>> &timeseries1, std::vector<std::complex<double>> &timeseries2, std::vector<std::complex<double>> &crosscorrelation);
std::vector<double> get_crosscorr_real(std::vector<double> &timeseries1, std::vector<double> &timeseries2, std::vector<double> &crosscorrelation);
std::vector<double> get_crosscorr(std::vector<double> &timeseries1, std::vector<double> &timeseries2, std::vector<double> &crosscorrelation); //Same as _real function. Overload for simplicity.

void normalize_correlation(std::vector<double> &corrleation, double normalizer);

// Get the correlation functions unnormalized but return the vector
std::vector<std::complex<double>> get_autocorr(std::vector<std::complex<double>> &timeseries);
std::vector<double> get_autocorr_real(std::vector<double> &timeseries);
std::vector<double> get_autocorr(std::vector<double> &timeseries); //Overload of get_autocorr with exact definition of get_autocorr_real 
/*
std::vector<double> get_crosscorr(std::vector<std::complex<double>> &timeseries1, std::vector<std::complex<double>> &timeseries2, std::vector<std::complex<double>> &crosscorrelation);
std::vector<double> get_crosscorr_real(std::vector<double> &timeseries1, std::vector<double> &timeseries2, std::vector<double> &crosscorrelation);
std::vector<double> get_crosscorr(std::vector<double> &timeseries1, std::vector<double> &timeseries2, std::vector<double> &crosscorrelation); //Same as _real function. Overload for simplicity.
*/


//get g2 auto and cross_correlations 

void get_g2(std::vector<double> &timeseries, std::vector<double> &autocorrelation); 								//Auto g2 correlation
void get_g2(std::vector<double> &timeseries1, std::vector<double> &timeseries2, std::vector<double> &correlation); 		//Cross g2 correlation

//get g2 auto and cross correlation and return the corresponding vector instead of using a void
std::vector<double> get_g2(std::vector<double> &timeseries); 								//Auto g2 correlation
std::vector<double> get_g2_cc(std::vector<double> &timeseries1, std::vector<double> &timeseries2); 		//Cross g2 correlation

// Print functions to file
void print_correlation(std::vector<double> &correlation, Time dt, std::string traceoutname, std::size_t maxN, std::size_t every, std::string unit);
void print_correlation(std::vector<std::complex<double>> &correlation, Time dt, std::string traceoutname, std::size_t maxN, std::size_t every, std::string unit);
void print_correlations(std::vector<std::vector<double>> &correlation, Time dt, std::string correlationoutname, std::string correlationheader, std::size_t maxN, int every, unsigned int precision=9);


}

#endif
