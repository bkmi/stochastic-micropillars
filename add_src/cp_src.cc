#include <fstream>
#include <string>
#include "cp_src.hh"

//----- CP Source Code Main to data/date Folder --- //

namespace katana{

void cp_src(std::string sourcecode, std::string path)
{
	std::ifstream source(sourcecode.c_str(), std::ios::binary);
	std::ofstream dest((path+sourcecode).c_str(), std::ios::binary);
	dest << source.rdbuf();

	source.close();
	dest.close();
}

}	
//----- END CP ------------------------------------ //


