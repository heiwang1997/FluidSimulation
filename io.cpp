#include "stdafx.h"
#include "io.h"
#include "Field.h"
#include "TimeStepController.h"
#include <fstream>

void io::dumpSolverToFile(const std::string & filename, Field * rhoField, Field * vxField,
	Field * vyField, Field * vzField, Field* thetaField, TimeStepController * step)
{
	std::ofstream fout(filename, std::ios::binary);
	rhoField->dumpFieldToFile(fout);
	vxField->dumpFieldToFile(fout);
	vyField->dumpFieldToFile(fout);
	vzField->dumpFieldToFile(fout);
	thetaField->dumpFieldToFile(fout);
	step->dumpToFile(fout);
	fout.close();
}

void io::loadSolverFromFile(const std::string & filename, Field * rhoField, Field * vxField,
	Field * vyField, Field * vzField, Field* thetaField, TimeStepController * step)
{
	LOG(INFO) << "Loading model from " << filename;
	std::ifstream fin(filename, std::ios::binary);
	rhoField->loadFiledFromFile(fin);
	vxField->loadFiledFromFile(fin);
	vyField->loadFiledFromFile(fin);
	vzField->loadFiledFromFile(fin);
	thetaField->loadFiledFromFile(fin);
	step->loadFromFile(fin);
	// Test Model Compatibility
	CHECK_EQ(fin.eof(), false) << "Missing values in model!";
	unsigned char probe; fin >> probe;
	CHECK_EQ(fin.eof(), true) << "Superfluous values in model!";
	LOG(INFO) << "Model integrity check passed!";
	fin.close();
}
