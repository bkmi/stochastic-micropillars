#ifndef KATANA_ARMA_WRAP__
#define KATANA_ARMA_WRAP__


// hier muss wohl kein include rein =)
//#include <cmath>
//#include <iostream>
//#include <vector>
//#include "fftw3.h"
//#include <fstream>
#include <string>
#include <complex.h>
#include <armadillo>

namespace katana
{

void add_vec_to_matrix(arma::mat &A, std::vector<double> vector);
std::vector<std::vector<double>> mat_to_stdvecvec(arma::mat &A); 
arma::mat stdvecvec_to_mat(std::vector<std::vector<double>> &vecvec);
}

#endif
