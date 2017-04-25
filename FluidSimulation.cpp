#include "stdafx.h"
#include "io.h"
#include "meta.h"
#include "Field.h"
#include "config.h"
#include "SimpleSolver.h"
#include "ThermalSolver.h"
#include "TimeStepController.h"
#include <direct.h>

const char* defaultConfigFilename = "thermal.cfg";

Field* getInitRhoField(int resX, int resY, int resZ, real dx, real ld, real vd);

Field* getSemiLiquidField(int resX, int resY, int resZ, real ld, real vd, real lperc);

Field* getSingleBubbleRhoField(int resX, int resY, int resZ, real dx, real ld, real vd);

Field* getInitThetaField(int resX, int resY, int resZ, real startTheta);

int main(int argc, char** argv)
{
	FLAGS_alsologtostderr = 1;
	FLAGS_colorlogtostderr = 1;
	google::InitGoogleLogging((const char*)argv[0]);
	Config* config;
	Field* initRhoField = 0, *initVxField = 0, *initVyField = 0, *initVzField = 0,
		*initThetaField = 0;
	TimeStepController* timeStep = 0;
	if (argc <= 2) {
		if (argc == 1) {
			LOG(WARNING) << "Configuration file not specified in argv[1], default to " <<
				defaultConfigFilename;
			config = new Config(defaultConfigFilename);
		}
		else {
			LOG(INFO) << "Starting from default rho field and all zero velocity fields";
			LOG(INFO) << "Now loading " << argv[1];
			config = new Config(argv[1]);
		}

		//initRhoField = getInitRhoField(config->resX(), config->resY(), config->resZ(),
		//	config->h(), config->vdwLiquidRho(), config->vdwVaporRho());

		//initRhoField = getSingleBubbleRhoField(config->resX(), config->resY(), config->resZ(),
		//	config->h(), config->vdwLiquidRho(), config->vdwVaporRho());

		initRhoField = getSemiLiquidField(config->resX(), config->resY(), config->resZ(),
			config->vdwLiquidRho(), config->vdwVaporRho(), 1.0f);

		initThetaField = getInitThetaField(config->resX(), config->resY(), config->resZ(),
			config->startTheta());

		timeStep = new TimeStepController(
			config->totalFrame(), config->gridFPS(), config->gridDt());
	}
	else if (argc == 3) {
		LOG(FATAL) << "Not supported.";
		return -1;

		config = new Config(argv[1]);
		initRhoField = new Field(config->resX(), config->resY(), config->resZ());
		initVxField = new Field(config->resX() + 1, config->resY(), config->resZ());
		initVyField = new Field(config->resX(), config->resY() + 1, config->resZ());
		initVzField = new Field(config->resX(), config->resY(), config->resZ() + 1);
		timeStep = new TimeStepController(1.0f, 1, 1.0f);
		io::loadSolverFromFile(argv[2], initRhoField, initVxField, initVyField, initVzField, timeStep);
	}
	else {
		std::cout << "Usage: FluidSimulation.exe <config file> <model file>" << std::endl;
		LOG(FATAL) << "Startup parameter not supported.";
		return 1;
	}
	
	// Copy configuration file into Snapshot folder.
	{
		std::string baseFolder = config->snapshotOutputDir() + "/" + config->runName() + "/";
		_mkdir(baseFolder.c_str());
		std::string configFilename = baseFolder + "config.cfg";
		std::ifstream src(config->filename, std::ios::binary);
		std::ofstream dst(configFilename, std::ios::binary);
		dst << src.rdbuf();
		src.close();
		dst.close();
		LOG(INFO) << "Configuration file copied to " << configFilename;
	}

	//SimpleSolver* simpleSolver = new SimpleSolver(config, initRhoField, initVxField, initVyField, initVzField);
	//simpleSolver->run(timeStep);

	//delete simpleSolver;
	ThermalSolver* thermalSolver = new ThermalSolver(config, initRhoField, 
		initVxField, initVyField, initVzField, initThetaField);
	thermalSolver->run(timeStep);

	delete thermalSolver;
	delete initRhoField;
	if (initThetaField) delete initThetaField;
	delete config;
	delete timeStep;
    return 0;
}

Field* getInitRhoField(int resX, int resY, int resZ, real dx, real ld, real vd) {
	Field* result = new Field(resX, resY, resZ);
	real* rho = result->content;

	real xTotal = dx * resX;
	real yTotal = dx * resY;
	real zTotal = dx * resZ;

	Vec3f bubble1 = Vec3f(0.335f /*0.41f*/ * xTotal, 0.50f * yTotal, 0.50f * zTotal);// *_xRes;
	Vec3f bubble2 = Vec3f(0.665f /*0.67f*/ * xTotal, 0.50f * yTotal, 0.50f * zTotal);// *_yRes;
	float Rb1 = 0.16f * xTotal;// 0.16;
	float Rb2 = 0.16f * xTotal;// 0.08;

	double totalmass = 0;

	real liquidDensity = ld, vaporDensity = vd;
	real interfaceScalingFactor = 200.0f;

	bool regularizedInterface = true;

	for (int z = 0; z < resZ; z++)
		for (int y = 0; y < resY; y++)
			for (int x = 0; x < resX; x++)
			{
				int index = result->getIndex(x, y, z);
				Vec3f gc = Vec3f(x + 0.5f, y + 0.5f, z + 0.5f) * dx;

				if (regularizedInterface) {
					real bubble1_dis = mag(gc - bubble1) - Rb1;
					real bubble2_dis = mag(gc - bubble2) - Rb2;

					rho[index] = vaporDensity + (liquidDensity - vaporDensity) / 2 * (
						tanh(interfaceScalingFactor * bubble1_dis) +
						tanh(interfaceScalingFactor * bubble2_dis)
						);
					/*
					rho[index] = vaporDensity + (liquidDensity - vaporDensity) / 2 *
					(1 + tanh(interfaceScalingFactor * bubble1_dis));
					*/
				}
				totalmass += rho[index];
			}
	printf("initial total mass: %f\n", totalmass);
	return result;
}

Field * getSemiLiquidField(int resX, int resY, int resZ, real ld, real vd, real lperc)
{
	Field* result = new Field(resX, resY, resZ);
	real* rho = result->content;

	real interfaceScalingFactor = 200.0f;

	float divPlaneY = lperc * resY;

	for (int z = 0; z < resZ; z++)
		for (int y = 0; y < resY; y++)
			for (int x = 0; x < resX; x++)
			{
				int index = result->getIndex(x, y, z);
				rho[index] = vd + (ld - vd) / 2 * (1 + tanh(interfaceScalingFactor * (divPlaneY - y)));
			}
	LOG(INFO) << "Initialized with SEMI-LIQUID " << lperc;
	return result;
}

Field* getSingleBubbleRhoField(int resX, int resY, int resZ, real dx, real ld, real vd) {
	Field* result = new Field(resX, resY, resZ);
	real* rho = result->content;

	real xTotal = dx * resX;
	real yTotal = dx * resY;
	real zTotal = dx * resZ;

	Vec3f bubble1 = Vec3f(0.5f /*0.41f*/ * xTotal, 0.50f * yTotal, 0.50f * zTotal);// *_xRes;
	float Rb1 = 0.16f * xTotal;// 0.16;

	double totalmass = 0;

	real liquidDensity = ld, vaporDensity = vd;
	real interfaceScalingFactor = 200.0f;

	bool regularizedInterface = true;

	for (int z = 0; z < resZ; z++)
		for (int y = 0; y < resY; y++)
			for (int x = 0; x < resX; x++)
			{
				int index = result->getIndex(x, y, z);
				Vec3f gc = Vec3f(x + 0.5f, y + 0.5f, z + 0.5f) * dx;

				if (regularizedInterface) {
					real bubble1_dis = mag(gc - bubble1) - Rb1;
					rho[index] = vaporDensity + (liquidDensity - vaporDensity) / 2 *
					(1 + tanh(interfaceScalingFactor * bubble1_dis));
				}
				totalmass += rho[index];
			}
	printf("initial total mass: %f\n", totalmass);
	return result;
}

Field * getInitThetaField(int resX, int resY, int resZ, real startTheta)
{
	Field* result = new Field(resX, resY, resZ);
	real* theta = result->content;

	int totalCells = resX * resY * resZ;
	for (int i = 0; i < totalCells; ++i) {
		theta[i] = startTheta;
	}

	return result;
}
