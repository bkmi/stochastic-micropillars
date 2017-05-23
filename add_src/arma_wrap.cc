#include <vector>
#include <armadillo>
#include "arma_wrap.hh"

namespace katana
{

void add_vec_to_matrix(arma::mat &A, std::vector<double> vector)
   {	if(size(A)(0)==0)
		{std::cout << "Matrix leer, initialisiere mit Vektor der Länge " << vector.size() << std::endl;
		 A.set_size(vector.size(),1);
		 arma::colvec vectorcopy=arma::conv_to<arma::colvec>::from(vector);
		 A.insert_cols(1,vectorcopy);
 		 A.shed_col(0);
	 	 return;
		}
	if(vector.size()!=size(A)(0)) 
		{std::cout << "Vector und Matrix passen nicht zusammen!" << std::endl;
		 std::cout << "Matrix hat Zeilenlänge " << size(A)(0) << ", aber Vektor ist " << vector.size() << " Zeilen lang! ABORT!" << std::endl;
		 return;
		};
	int cols=size(A)(1);	
        arma::colvec vectorcopy=arma::conv_to<arma::colvec>::from(vector);
	A.insert_cols(cols, vectorcopy);
   }

std::vector<std::vector<double>> mat_to_stdvecvec(arma::mat &A)
   {
    std::vector<std::vector<double>> V(A.n_rows);
    for (size_t i = 0; i < A.n_rows; ++i)
	 { V[i] = arma::conv_to<std::vector<double>>::from(A.row(i));
	 };
    return V;
   }

arma::mat stdvecvec_to_mat(std::vector<std::vector<double>> &vecvec)
	{ int rows=vecvec.size(), cols=vecvec[0].size();
	  arma::mat MAT(1, cols);
          for (size_t i=0; i<rows; ++i)
	 	{ arma::rowvec rowvector=arma::conv_to<arma::rowvec>::from(vecvec[i]);
		  MAT.insert_rows(i, rowvector);		
		}
	  MAT.shed_row(1);
 	  return MAT;
	}

}
