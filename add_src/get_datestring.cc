#include <string>
#include <sys/time.h>
#include <iostream>
namespace katana

{
std::string getDateString()
	{
		time_t t=time(0);
		struct tm * now = localtime( & t );
		   std::string date= std::to_string(now->tm_year + 1900) + '-' + std::to_string(now->tm_mon + 1) + '-' + std::to_string(now->tm_mday) + '/';  
		return date;
	};

}

