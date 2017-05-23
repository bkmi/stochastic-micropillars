#ifndef KATANA_PRINT_TRACE__
#define KATANA_PRINT TRACE__


// hier muss wohl kein include rein =)
//#include <cmath>
//#include <iostream>
//#include <vector>
//#include "fftw3.h"
//#include <fstream>
#include <string>
#include <complex.h>

namespace katana
{
using Time=double;
void dump_trace(std::vector<std::complex<double>> &timetrace, Time dt, std::string traceoutname, std::string traceheader, int every, unsigned int precision=9);
void dump_trace(std::vector<std::vector<double>> &timetrace, Time dt, std::string traceoutname, std::string traceheader, int every, unsigned int precision=9);
void dump_trace(std::vector<double> &timetrace, Time dt, std::string traceoutname, std::string traceheader, int every, unsigned int precision=9);
}

#endif 
