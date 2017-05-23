#ifndef ROEHM_SOLVER_TRANSFORM__
#define ROEHM_SOLVER_TRANSFORM__  	//include guard for the precompiler

#include <algorithm>
#include <type_traits>

namespace roehm 
{	
	
	//define a general mechanism to add/multiply + add several inputs - independent of wether they are vectors or arithmetic types (int, double...)
	//however, the compiler must know at compile time, so it will only choose the correct interpretation of transform


	//template for two arithemtic inputs which are transformed into an output via Operation
	template <typename InputScalar, typename OutputScalar, typename Operation, typename std::enable_if<std::is_arithmetic<InputScalar>::value>::type* = nullptr>
	void transform(const InputScalar &src1, const InputScalar &src2, OutputScalar &dest, Operation op)
	{
	dest = op(src1, src2);
	}
	
	
	//template for vector input/output
	template <typename InputRARange, typename OutputRARange, typename Operation, typename std::enable_if<!std::is_arithmetic<InputRARange>::value>::type* = nullptr>
	void transform(const InputRARange &src1, const InputRARange &src2, OutputRARange &dest, Operation op) 
	{
	auto n = std::min(src1.size(), src2.size());
		for (size_t i = 0; i < n; ++i) 				//do the same for alle elements of the vector
		{
			dest[i] = op(src1[i], src2[i]);
		}
	}
	
	//template <
		//typename InputScalar, typename OutputScalar, typename Operation,
		//enable_if<std::is_arithmetic<InputScalar>::value> = detail::dummy_enabler>
	//auto transform(const InputScalar &src1, const InputScalar &src2,
				//const InputScalar &src3, const InputScalar &src4,
				//OutputScalar &dest, Operation op) -> void {
	//dest = op(src1, src2, src3, src4);
	//}
	
	//template <
		//typename InputRARange, typename OutputRARange, typename Operation,
		//enable_if<!std::is_arithmetic<InputRARange>::value> = detail::dummy_enabler>
	//auto transform(const InputRARange &src1, const InputRARange &src2,
				//const InputRARange &src3, const InputRARange &src4,
				//OutputRARange &dest, Operation op) -> void {
	//auto n = std::min(src1.size(), src2.size());
	//for (size_t i = 0; i < n; ++i) {
		//dest[i] = op(src1[i], src2[i], src3[i], src4[i]);
	//}
	//}

} 

#endif // ifndef ROEHM_SOLVER_TRANSFORM__

