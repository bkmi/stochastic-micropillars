#ifndef ROEHM_INITIALIZE_STATE__
#define ROEHM_INITIALIZE_STATE__  	//include guard for the precompiler


#include <array>
#include <type_traits>
#include <vector>

namespace roehm 
{	//initialize a new object of the same type as the argument. We need this, because we do not always know wether we are dealing with vectors, complex numbers or something else...

	//for arithmetic types (int, double, ...)

	template <typename T> 
	typename std::enable_if<std::is_arithmetic<T>::value, T>::type 	//this line defines the return type as the same type as T, but only of T is arithmetic (thanks to the enable_if_c template function)
		initialize_state(const T &x0)
	{
	return 0.0;	//value for a newly created arithmetic variable
	}
	
	//overloading for arrays
	
	template <typename T, size_t N>
	auto initialize_state(const std::array<T, N> &x0) -> std::array<T, N> 
	{
	return std::array<T, N>(); //returns a new array of Size N with arguments of type T, all of value 0
	}
	
	//overloading for vectors 
	
	template <typename T>
	auto initialize_state(const std::vector<T> &x0) -> std::vector<T> 
	{
	return std::vector<T>(x0.size());
	}
	
	
	} // end namespace roehm
	
#endif /* end of include guard: ROEHM_INITIALIZE_STATE__ */
