#include <cstdlib>
#include <sys/time.h>
#include <climits>
#include <sstream>
#include <iostream>
#include "set_seed.hh"

namespace katana{

Seed::Seed(){}
void Seed::set_seed()
		{
			struct timeval t;
			gettimeofday(&t, NULL);
		
			unsigned int time_sec = t.tv_sec;     // Current seconds since epoc (midnight, jan 1, 1970)
			unsigned int time_musec = t.tv_usec; // Number of microseconds this second
			unsigned long long int rawseed;			// long long is sufficient in every case. for AMD64 same as u long int
			unsigned int seed;
			
		
			std::stringstream seedstream;
			std::string seedstring;
									
			seedstream << time_sec << time_musec; 
			seedstream >> rawseed; // time in sec+musec into seedstring
			
			seed= (unsigned int) rawseed % UINT_MAX;
			std::srand(seed);				
							
		}
void Seed::set_seed(int seed)
		{ 
			srand(seed);
		}


}
