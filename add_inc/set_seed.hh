#ifndef KATANA_SET_SEED__
#define KATANA_SET_SEED__


#include <cstdlib>
#include <sys/time.h>
#include <climits>
#include <sstream>
#include <iostream>

namespace katana
{

class Seed
{ 	public:
		Seed();
		void set_seed();
		void set_seed(int );
};
}

#endif
