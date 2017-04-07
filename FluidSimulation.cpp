#include "stdafx.h"
#include "meta.h"
#include "Field.h"
#include "config.h"
#include "SimpleSolver.h"

const char* defaultConfigFilename = "isothermal.cfg";

Field* getInitRhoField(int resX, int resY, int resZ, real dx, real ld, real vd);

int main(int argc, char** argv)
{
	FLAGS_alsologtostderr = 1;
	FLAGS_colorlogtostderr = 1;
	google::InitGoogleLogging((const char*)argv[0]);
	Config* config;
	if (argc == 1) {
		LOG(INFO) << "Configuration file not specified in argv[1], default to " <<
			defaultConfigFilename;
		config = new Config(defaultConfigFilename);
	}
	else {
		config = new Config(argv[1]);
	}
	Field* initRhoField = getInitRhoField(config->resX(), config->resY(), config->resZ(),
		config->h(), config->vdwLiquidRho(), config->vdwVaporRho());
	SimpleSolver* simpleSolver = new SimpleSolver(config, initRhoField, 0, 0, 0);
	simpleSolver->run();
	delete simpleSolver;
	delete initRhoField;
	delete config;
    return 0;
}

Field* getInitRhoField(int resX, int resY, int resZ, real dx, real ld, real vd) {
	Field* result = new Field(resX, resY, resZ);
	real* rho = result->content;

	real xTotal = dx * resX;
	real yTotal = dx * resY;
	real zTotal = dx * resZ;

	Vec3f bubble1 = Vec3f(0.41f, 0.50f, 0.50f) * xTotal;// *_xRes;
	Vec3f bubble2 = Vec3f(0.66f /*0.67f*/, 0.50f, 0.50f) * yTotal;// *_yRes;
	float Rb1 = 0.16f * xTotal;// *_xRes;
	float Rb2 = 0.08f * yTotal;// *_xRes;
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
