#ifndef ROEHM_ASSIGN_STATE__
#define ROEHM_ASSIGN_STATE__  	//include guard for the precompiler


#include <array>
#include <type_traits>
#include <vector>

DECREPIT

namespace roehm 
{	//initialize a new object of the same type as the argument. We need this, because we do not always know wether we are dealing with vectors, complex numbers or something else...

	//for arithmetic types (int, double, ...)

	template <typename T> 
	//typename std::enable_if<std::is_arithmetic<T>::value, T>::type 	//this line defines the return type as the same type as T, but only of T is arithmetic (thanks to the enable_if_c template function)
	void assign_state(const T &x0, T &x1)
	{
	return x0;	//value for a newly created arithmetic variable
	}
	
	//overloading for arrays
	
	template <typename T, size_t N>
	void assign_state(const std::array<T, N> &x0, std::array<T, N> T &x1) -> std::array<T, N> 
	{
	return x0; //returns a new array of Size N with arguments of type T, all of value 0
	}
	
	//overloading for vectors 
	
	template <typename T>
	void assign_state(const std::vector<T> &x0, std::vector<T> &x1) -> std::vector<T> 
	{
	return x0;
	}
	
	
	} // end namespace roehm
	
#endif /* end of include guard: ROEHM_ASSIGN_STATE__ */
