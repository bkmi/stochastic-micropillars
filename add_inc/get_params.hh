#ifndef KATANA_GET_PARAMS__
#define KATANA_GET_PARAMS__


#include <string>


namespace katana
{

double getCmdOption(char ** begin, char ** end, const std::string & option, double Default);
bool getCmdOption_bool(char ** begin, char ** end, const std::string & option, bool Default);  // Switches to non Default
std::string getCmdOption(char ** begin, char ** end, const std::string & option, std::string Default);


bool cmdOptionExists(char** begin, char** end, const std::string& option);

}
#endif
