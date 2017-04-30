#pragma once

#include <string>

class Field;
class TimeStepController;

class io
{
public:
	static void dumpSolverToFile(const std::string& filename, Field* rhoFiled, Field* vxField,
		Field* vyField, Field* vzField, Field* thetaField, TimeStepController* step);
	static void loadSolverFromFile(const std::string& filename, Field* rhoFiled, Field* vxField,
		Field* vyField, Field* vzField, Field* thetaField, TimeStepController* step);
};

