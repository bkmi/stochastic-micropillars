#include <complex>
#include "gwNoise.hh"


// BOX MUELLER VERFAHREN

namespace katana
{

std::complex<double> gwNoise()
{	double spont_phi= (rand()%RAND_MAX)/(double) RAND_MAX *2.0*3.14159265359;
	double spont_amp= sqrt(-2* log((rand()%RAND_MAX+1)/(double) RAND_MAX)); //gibt nur die Zufallszahl nach normierter Gaussverteilungi
	return  {spont_amp*cos(spont_phi), spont_amp*sin(spont_phi)};
}

}
