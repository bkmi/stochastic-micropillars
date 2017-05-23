#ifndef KATANA_AVG_CONTAINER__
#define KATANA_AVG_CONTAINER__

#include <cstdlib>
#include <vector>
#include <iostream>
#include <armadillo>

namespace katana
{

class container_1d
{ 
	int counter=0;	
	std::vector<double> sum_state;
	public:
		container_1d();
		container_1d(std::vector<double> init_state);
		void add(std::vector<double> toadd_state);
		void get_avg(std::vector<double> &outvec);
		std::vector<double> get_avg();
		void clear();
		int get_counter();
};

class container_2d
{ 
	int counter=0;	
	arma::mat sum_state;
	arma::mat stddev_state;
	public:
		container_2d();
		container_2d(std::vector<std::vector<double>> init_state);
		void add(std::vector<std::vector<double>> toadd_state);
		void get_avg(std::vector<double> &outvec);
		std::vector<std::vector<double>> get_avg();
		void clear();
		int get_counter();
};

}



#endif
