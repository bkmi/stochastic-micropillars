#include <cmath>
#include <iostream>
#include <vector>
#include <complex.h>
#include <fstream>
#include <string>
#include "print_trace.hh"

namespace katana 
{
using Time=double;

void dump_trace(std::vector<std::complex<double>> &timetrace, Time dt, std::string traceoutname, std::string traceheader, int every, unsigned int precision)
	{	std::ofstream traceout (traceoutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
		traceout.precision(precision);
		traceout << traceheader << std::endl; 
		for(int iter= 0; iter<timetrace.size(); iter+=every )  
			{traceout << iter*dt << "\t" << timetrace[iter].real() << "\t" << timetrace[iter].imag() << std::endl;
			}				
		traceout.close();	
	}
void dump_trace(std::vector<std::vector<double>> &timetrace, Time dt, std::string traceoutname, std::string traceheader, int every, unsigned int precision)
	{	std::ofstream traceout (traceoutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
		traceout.precision(precision);
		traceout << traceheader << std::endl; 
		for(int iter= 0; iter<timetrace.size(); iter+=every )  
			{traceout << iter*dt;
				for(int elem=0; elem<timetrace[0].size(); elem++)
					{traceout << "\t" << timetrace[iter][elem];}
				traceout << std::endl;
			}					
		traceout.close();
	
	}


void dump_trace(std::vector<double> &timetrace, Time dt, std::string traceoutname, std::string traceheader, int every, unsigned int precision)
	{	std::ofstream traceout (traceoutname, std::ios::out | std::ios::trunc);      //Output File (ofstream) overwrite
		traceout.precision(precision);
		traceout << traceheader << std::endl; 
		for(int iter= 0; iter<timetrace.size(); iter+=every )  
			{traceout << iter*dt;
			 traceout << "\t" << timetrace[iter] << std::endl;
			}					
		traceout.close();
	
	}



}
