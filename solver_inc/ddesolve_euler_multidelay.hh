//DDE Euler based Multi Delay Solver 2016/01/13

//Ver 0.1

#ifndef ROEHM_DDESOLVE_EULER_MULTIDELAY__
#define ROEHM_DDESOLVE_EULER_MULTIDELAY__ 	//include guard for the precompiler


#include "solver_errorchecks.hh"
#include "solver_transform.hh"
#include "solver_initialize_state.hh"
#include "solver_select_value_type.hh"
#include <iostream>
#include <vector>
#include <cassert>

namespace roehm
{

/*ddesolve_euler_multidelay(1, 2, 3, 4, 5, 6, 7, 8) needs: 
*1. a (void) function, here called system, which has as arguments [State X, Statevector X_tau, Time t and State dxdt] (in that order) and then writes into dxdt after solving a set of time-discretized differential equations dependent on State X, Statevector X_tau and Time t 
*2. starting State (technically a pointer): &X
*3. pointer to the delay state vector: &X_tau_vector  which should be FILLED, otherwise it is ASSUMED TO BE EQUAL TO X_INIT IN EVERY ENTRY!!! and its length must be tau in timesteps! [Yes! not +1!]
*4. array of delay times: &tau_list_vector; first entry will correspond to the first entry of Statevector X_tau
*5. starting Time of the integration: t_init 
*6. ending Time of the integration: t_end 
*7. Time step: dt
*8. a function to be executed after every timestep dt called afterstepfunction, dependent on States &X, &X_tau, Time T and dxdt 
*/

template <typename System, typename State, typename Statevector, typename Time, typename AfterStepFunction>	
void ddesolve_euler_multidelay(System system, State &x_init, Statevector &x_tau_vector, std::vector<size_t> &tau_list_vector, const Time t_init, const Time t_end, const Time dt, AfterStepFunction afterstepfunction)
	{	
		using Value = typename select_value_type<State>::type;
		size_t Hist_len = x_tau_vector.size();
		size_t Number_of_delays = tau_list_vector.size();
				
		/*error checks should enter here */
		
		{	//Check correct delay array size
			size_t maxdelay = 0;
			for(size_t delaynumber = 0; delaynumber < Number_of_delays; delaynumber++)
			{
				if ( tau_list_vector[delaynumber] > maxdelay)
				{
					maxdelay = tau_list_vector[delaynumber];
				}
			}		
			assert(maxdelay ==  Hist_len);
		}
		
		
		//Initialize State and Time variables used for calculations during time integration
		State roehm_x_dynamic = roehm::initialize_state(x_init);	//adopt the values for starting		
		std::vector<State> roehm_x_tau(Number_of_delays, x_init);	//create our vector of delayed states, that needs to go into the system function; use x_init as initiliazer, because otherwise our constructor does not seem to know what type of vector to use
		Time roehm_t_dynamic = roehm::initialize_state(t_init);
		
		
		roehm::transform(roehm_x_dynamic, x_init, roehm_x_dynamic, 		//from  "solver_transform.hh"
		[&](Value x0_i, Value x1_i) { return x1_i;});		//add values of x_init to roehm_x_dynamic so that it carries the initial conditions, all in the name of not using the global variable x_init...
		roehm_t_dynamic = t_init;  // Now adopt the values for starting		
		
		//Initiliaze derivative
		State roehm_dxdt_dynamic = roehm::initialize_state(x_init);
		
		//delay-counter
		size_t k = 0; 
				
		for (roehm_t_dynamic = t_init; roehm_t_dynamic < t_end;)
			{	
				
				
				for(size_t delaynumber = 0; delaynumber < Number_of_delays; delaynumber++)
				{	
					size_t Correct_spot = (k - tau_list_vector[delaynumber] + Hist_len)%Hist_len;  //can modulus use negative numbers correctly in c++?
					
					roehm_x_tau[delaynumber] =  x_tau_vector[Correct_spot];
				}
				//read in current tau-values
												
				system(roehm_x_dynamic, roehm_x_tau, roehm_t_dynamic, roehm_dxdt_dynamic);  //evaluate differential equations 
												
				x_tau_vector[k] = roehm_x_dynamic;						//record the old state
	
				roehm::transform(roehm_x_dynamic, roehm_dxdt_dynamic, roehm_x_dynamic, 		//from  "solver_transform.hh"
			[&](Value x0_i, Value dxdt1_i) {   return x0_i + dxdt1_i * dt;});		//add differential to state variable, we inline-define a lambda-function to save space

				roehm_t_dynamic += dt; // better spot for the increment than the increment slot in the for-loop, this way the output happens at the (in my opinion) most logical time  //update 07/2015: really? I'm no longer sure...	
												
				afterstepfunction(roehm_x_dynamic, roehm_x_tau, roehm_t_dynamic, roehm_dxdt_dynamic);	//evaluate after-step-function
				
				k += 1; //increase counter for delay array by 1
		
				k = k%Hist_len;		//if outside of delayarraylimits, decrease by Hist_len
					
			}
						
		//to ensure a smooth transition to the next call of the dde solver, we will set all external variables in a state compatible with a new call	
		x_init = roehm_x_dynamic;	//the input state variables will be set to their last value, making all our efforts to not use it void. well, not ALL of them. We are outside of our Ouput routine, at least.
		
		if (k != 0 )	//if the counter for the delay array is not 0, the delay array is not ordered properly!
		{
			Statevector roehm_sort_delay_vector = roehm::initialize_state(x_tau_vector);
			
			for(int u = 0; u < Hist_len; u++)
			{
				roehm_sort_delay_vector[u] = x_tau_vector[k];
				
				k += 1; //increase counter for delay array by 1
		
				k = k%Hist_len;		//if outside of delayarraylimits, decrease by Hist_len				
			}
			x_tau_vector = roehm_sort_delay_vector;				
		}	
		
		
		
	}

}


#endif //end of include guard ROEHM_DDESOLVE_EULER_MULTIDELAY__
