#ifndef ROEHM_SOLVER_ADD_OPERATOR__
#define ROEHM_SOLVER_ADD_OPERATOR__  	//include guard for the precompiler

#include <vector>
//obsolete header for now


namespace roehm
{
	
	//define a general addition for two variables (think: integer, double, complex...) where the first entry gets rewritten

	template<typename T>  
	void add_operator(T* a, T* b)
	{
		*a = *a + *b;		//Take value at the adress of a and add to it the vale at the adress of b
	}

	template<typename T>  
	void add_operator(T* a, const T* b)			//should the second pointer given be considered a constant, we will just use the same code again (probably a smarter way?)
	{
		*a = *a + *b;		
	}
		
	//define a specialized rule if vectors are involved, possibly slow for large vectors?
	//vector pointers

	template<typename T>
	void add_operator(std::vector<T>* a, std::vector<T>* b)
	{		//for performance reason I do not check if the vectors have the same size. Lets just hope the compiler will throw an error in all relevant cases.

			for (std::size_t i = 0; i < a->size(); i++) a->at(i) = a->at(i) + b->at(i);  //for-loop for adding, std::size_t is the unsigned integer type that size() returns
	}
		
	template<typename T>
	void add_operator(std::vector<T>* a, const std::vector<T>* b)
	{		//for performance reason I do not check if the vectors have the same size. Lets just hope the compiler will throw an error in all relevant cases.

			for (std::size_t i = 0; i < a->size(); i++) a->at(i) = a->at(i) + b->at(i);  //for loop for adding
	}
		
	
	
	
		
	
	
	
}

#endif  //End of include guard for ROEHM_SOLVER_ADD_OPERATOR__
