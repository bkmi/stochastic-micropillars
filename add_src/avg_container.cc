#include "avg_container.hh"
#include "arma_wrap.hh"
#include <armadillo>
#include <cstdlib>
#include <vector>

namespace katana
{

container_1d::container_1d(){counter=0;};
container_1d::container_1d(std::vector<double>init_state)	
	{
	sum_state.resize(init_state.size());
	for (std::size_t i = 0; i < init_state.size(); i++)  sum_state[i] = init_state[i];
	counter=1;
	}	

void container_1d::add(std::vector<double> add_state)
	{
	if(counter==0) {
		 sum_state.resize(add_state.size());
		 for (std::size_t i = 0; i < add_state.size(); i++)  sum_state[i] = add_state[i];
		}
	else{	for (std::size_t i = 0; i < add_state.size(); i++)  sum_state[i] += add_state[i];};
	counter++;
	}

void container_1d::get_avg(std::vector<double> &outvec)
	{	
	  outvec.resize(sum_state.size());
	  for (std::size_t i = 0; i < sum_state.size(); i++) outvec[i] = sum_state[i]/counter;
	}
std::vector<double> container_1d::get_avg()
	{ std::vector<double> outvec(sum_state.size());
	  for (std::size_t i = 0; i < sum_state.size(); i++) outvec[i] = sum_state[i]/counter;
	  return outvec;
	}


void container_1d::clear()
	{
	 sum_state.clear();
	 counter=0;
	}

container_2d::container_2d(){counter=0;};
container_2d::container_2d(std::vector<std::vector<double>>init_state)	
	{
//	sum_state.resize(init_state.size(), init_state[0].size());
	//for (std::size_t i = 0; i < init_state.size(); i++) sum_state(i) = arma::conv_to 
//	counter=1;
	}	

void container_2d::add(std::vector<std::vector<double>> add_state)
	{
	//if(counter==0) {
//		 sum_state.resize(add_state.size());
//		 for (std::size_t i = 0; i < add_state.size(); i++)  sum_state[i] = add_state[i];
//		}
//	else{	for (std::size_t i = 0; i < add_state.size(); i++)	{
//			  for (std::size_t j=0; j<add_state[i].size(); j++) sum_state[i][j] = sum_state[i][j]+ add_state[i][j];};
//	counter++;
	}

void container_2d::get_avg(std::vector<double> &outvec)
	{	
	}

std::vector<std::vector<double>> container_2d::get_avg()
	{std::vector<std::vector<double>> returnvec;

	return returnvec;	
	}

void container_2d::clear()
	{
	 sum_state.clear();
	 counter=0;
	}



}
