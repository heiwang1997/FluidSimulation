#include "stdafx.h"
#include "config.h"
#include <iostream>
#include <fstream>

#ifdef _M_X64
#pragma comment (lib, "libconfig-1.3.2/lib64/libconfig++.lib")
#else
#pragma comment (lib, "libconfig-1.3.2/lib/libconfig++.lib")
#endif



Config::Config(std::string filename) {
	this->filename = filename;
	try {
		cfg_.readFile(filename.c_str());
	} catch (libconfig::FileIOException e) {
		std::cout <<filename<<" not found."<<std::endl;
	} catch (libconfig::ParseException e) {
		//LOG4CXX_INFO(logger, e.what());
		std::cout <<filename<<" parse error: "<<e.getError()<<" at line "<<e.getLine()<<std::endl;
	}
}

Config::~Config() {
	return;
}
