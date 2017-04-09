#include "stdafx.h"
#include "io.h"
#include "meta.h"
#include "Field.h"
#include "config.h"
#include "SimpleSolver.h"
#include "TimeStepController.h"

const char* defaultConfigFilename = "isothermal.cfg";

Field* getInitRhoField(int resX, int resY, int resZ, real dx, real ld, real vd);

int main(int argc, char** argv)
{
	FLAGS_alsologtostderr = 1;
	FLAGS_colorlogtostderr = 1;
	google::InitGoogleLogging((const char*)argv[0]);
	Config* config;
	Field* initRhoField = 0, *initVxField = 0, *initVyField = 0, *initVzField = 0;
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
		initRhoField = getInitRhoField(config->resX(), config->resY(), config->resZ(),
			config->h(), config->vdwLiquidRho(), config->vdwVaporRho());
		timeStep = new TimeStepController(
			config->totalFrame(), config->gridFPS(), config->gridDt());
	}
	else if (argc == 3) {
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
	
	SimpleSolver* simpleSolver = new SimpleSolver(config, initRhoField, initVxField, initVyField, initVzField);
	simpleSolver->run(timeStep);

	delete simpleSolver;
	delete initRhoField;
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

	Vec3f bubble1 = Vec3f(0.457f /*0.41f*/ * xTotal, 0.50f * yTotal, 0.50f * zTotal);// *_xRes;
	Vec3f bubble2 = Vec3f(0.543f /*0.67f*/ * xTotal, 0.50f * yTotal, 0.50f * zTotal);// *_yRes;
	float Rb1 = 0.04f * xTotal;// 0.16;
	float Rb2 = 0.04f * xTotal;// 0.08;
							  /*
							  // Rotate 30 test
							  real rotate_angle = -0.5236f;
							  Vec3f rotate_center = Vec3f(0.50, 0.50, 0.50) * xTotal;
							  Vec3f bubble1_rotate_vec = bubble1 - rotate_center;
							  bubble1 = rotate_center + Vec3f(cos(rotate_angle) * bubble1_rotate_vec[0] - sin(rotate_angle) * bubble1_rotate_vec[1],
							  sin(rotate_angle) * bubble1_rotate_vec[0] + cos(rotate_angle) * bubble1_rotate_vec[1], 0.0);
							  Vec3f bubble2_rotate_vec = bubble2 - rotate_center;
							  bubble2 = rotate_center + Vec3f(cos(rotate_angle) * bubble2_rotate_vec[0] - sin(rotate_angle) * bubble2_rotate_vec[1],
							  sin(rotate_angle) * bubble2_rotate_vec[0] + cos(rotate_angle) * bubble2_rotate_vec[1], 0.0);
							  // End rotate
							  */
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
