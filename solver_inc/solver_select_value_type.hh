#ifndef ROEHM_SELECT_VALUE_TYPE__
#define ROEHM_SELECT_VALUE_TYPE__

	//adopted from oliver esser = lego

#include <type_traits>

namespace roehm
{
	namespace detail 
	{
	
		template <bool IsFloatingPoint, class State> 
		struct select_value_type_impl 
		{
			using type = State;
		};
		
		template <class State>
		struct select_value_type_impl<false, State> 
		{
			using type = typename State::value_type;
		};
	
	} // end namespace roehm::detail

	template <class State> 
	struct select_value_type 
	{
		using type = typename detail::select_value_type_impl<std::is_floating_point<State>::value, State>::type;
	};

}

#endif // ifndef ROEHM_SELECT_VALUE_TYPE__
