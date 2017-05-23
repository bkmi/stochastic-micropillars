#include "get_params.hh"
#include <algorithm>
#include <string>
#include <cstdlib>
#include <iostream>
namespace  katana
{

double getCmdOption(char ** begin, char ** end, const std::string & option, double Default)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
	char* b;
	double d;
	d = strtod(*itr, &b);
	if (0 == d && *itr == b) 
		{	std::cout << "Input of Option ''" << option << "'' was wrong. Setting to default: " << option << " " << Default << std::endl;      // error handling.
   			return Default;
		}
	std::cout << "Set Option: "<< option << " " << d << std::endl;
        return d;

    }
    return Default;
}

std::string getCmdOption(char ** begin, char ** end, const std::string & option, std::string Default)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
	std::cout << "Set Option: "<< option << " = " << *itr << std::endl;
        return *itr;

    }
    return Default;
}


bool getCmdOption_bool(char ** begin, char ** end, const std::string & option, bool Default)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end)
    {	bool Val=!Default;
	std::cout << "Set Option: "<< option << " = " << Val << std::endl;
        return Val;

    }
    return Default;
}







bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}


//char* getCmdOption(char ** begin, char ** end, const std::string & option)
//{
//    char ** itr = std::find(begin, end, option);
//    if (itr != end && ++itr != end)
//    {
//        return *itr;
//    }
//    return 0;
//}
//
//bool cmdOptionExists(char** begin, char** end, const std::string& option)
//{
//    return std::find(begin, end, option) != end;
//}
}
