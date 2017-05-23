#include <cstdio>
#include <iostream>
#include "remove.hh"

namespace katana
{

void rm(std::string filename)
{
 if(std::remove(filename.c_str())!=0)
	{ std::cout << "ERROR deleting file, maybe it doesn't exist or wrong path"<< std::endl;}
}

}


