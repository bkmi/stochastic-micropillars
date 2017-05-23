//DDE Euler based Single Delay Solver 2015/05/03

//Ver 0.6

#ifndef ROEHM_DDESOLVE_EULER_SINGLEDELAY__
#define ROEHM_DDESOLVE_EULER_SINGLEDELAY__ 	//include guard for the precompiler


#include "solver_errorchecks.hh"
#include "solver_transform.hh"
#include "solver_initialize_state.hh"
#include "solver_select_value_type.hh"
#include <iostream>
#include <vector>

namespace roehm
{

/*ddesolve_euler_singledelay(1, 2, 3, 4, 5, 6, 7) needs: 
*1. a (void) function, here called system, which has as arguments [State X, State X_tau, Time t and State dxdt] (in that order) and then writes into dxdt after solving a set of time-discretized differential equations dependent on State X, X_tau and Time t 
*2. starting State (technically a pointer): &X
*3. delay state vector pointer: &X_tau_vector
*4. starting Time of the integration: t_init 
*5. ending Time of the integration: t_end 
*6. Time step: dt
*7. a function to be executed after every timestep dt called afterstepfunction, dependent on States &X, &X_tau, Time T and dxdt 
*/

template <typename System, typename State, typename Statevector, typename Time, typename AfterStepFunction>	
void ddesolve_euler_singledelay(System system, State &x_init, Statevector &x_tau_vector, const Time t_init, const Time t_end, const Time dt, AfterStepFunction afterstepfunction)
	{	
		using Value = typename select_value_type<State>::type;
		
		/*error checks should enter here */
	
		//Initialize State and Time variables used for calculations during time integration
		State roehm_x_dynamic = roehm::initialize_state(x_init);	//adopt the values for starting		
		State roehm_x_tau = roehm::initialize_state(x_init);		
		Time roehm_t_dynamic = roehm::initialize_state(t_init);
		
		
		roehm::transform(roehm_x_dynamic, x_init, roehm_x_dynamic, 		//from  "solver_transform.hh"
		[&](Value x0_i, Value x1_i) { return x1_i;});		//add values of x_init to roehm_x_dynamic so that it carries the initial conditions, all in the name of not using the global variable x_init...
		roehm_t_dynamic = t_init;  // Now adopt the values for starting		
		
		//Initiliaze derivative
		State roehm_dxdt_dynamic = roehm::initialize_state(x_init);
		
		//delay-counter
		size_t k = 0; 
		size_t hist_len = x_tau_vector.size();
		
		//for-loop 
		
		for (roehm_t_dynamic = t_init; roehm_t_dynamic < t_end;)
			{			
				roehm_x_tau = x_tau_vector[k];	 //read in current tau-value
								
				system(roehm_x_dynamic, roehm_x_tau, roehm_t_dynamic, roehm_dxdt_dynamic);  //evaluate differential equations 
								
				x_tau_vector[k] = roehm_x_dynamic;						//record the old state
	
				roehm::transform(roehm_x_dynamic, roehm_dxdt_dynamic, roehm_x_dynamic, 		//from  "solver_transform.hh"
			[&](Value x0_i, Value dxdt1_i) {   return x0_i + dxdt1_i * dt;});		//add differential to state variable, we inline-define a lambda-function to save space

				roehm_t_dynamic += dt; // better spot for the increment than the increment slot in the for-loop, this way the output happens at the (in my opinion) most logical time  //update 07/2015: really? I'm no longer sure...	
								
				afterstepfunction(roehm_x_dynamic, roehm_x_tau, roehm_t_dynamic, roehm_dxdt_dynamic);	//evaluate after-step-function
		
				k += 1; //increase counter for delay array by 1
		
				k = k%hist_len;		//if outside of delayarraylimits, decrease by hist_len
					
			}
			
			
			
		//to ensure a smooth transition to the next call of the dde solver, we will set all external variables in a state compatible with a new call	
		x_init = roehm_x_dynamic;	//the input state variables will be set to their last value
		
		if (k != 0 )	//if the counter for the delay array is not 0, the delay array is not ordered properly!
		{
			Statevector roehm_sort_delay_vector = roehm::initialize_state(x_tau_vector);
			
			for(int u = 0; u < hist_len; u++)
			{
				roehm_sort_delay_vector[u] = x_tau_vector[k];
				
				k += 1; //increase counter for delay array by 1
		
				k = k%hist_len;		//if outside of delayarraylimits, decrease by hist_len				
			}
			x_tau_vector = roehm_sort_delay_vector;				
		}	
		
		
		
	}

}


#endif //end of include guard ROEHM_DDESOLVE_EULER_SINGLEDELAY
