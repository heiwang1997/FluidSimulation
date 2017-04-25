#pragma once

#include <string>
#include <map>
#include <stdio.h>
#include "libconfig-1.3.2\include\libconfig.hpp"

class Config {
public:
	Config(std::string filename);
	~Config();

	std::string filename;
	// Getters begin
	std::string runName() { return (const char*)cfg_.lookup("runName"); }

	int resX() {return cfg_.lookup("grid.sizeX");}
	int resY() {return cfg_.lookup("grid.sizeY");}
	int resZ() {return cfg_.lookup("grid.sizeZ");}
	float h() { return cfg_.lookup("grid.h"); }

	int totalFrame() {return cfg_.lookup("timing.totalFrame");}
	int gridFPS() {return cfg_.lookup("timing.FPS");}
	float gridDt() {return cfg_.lookup("timing.dt");}
	
	float vdwPA() { return cfg_.lookup("vdwFluid.a"); }
	float vdwPB() { return cfg_.lookup("vdwFluid.b"); }
	float vdwTheta() { return cfg_.lookup("vdwFluid.theta"); }
	float vdwLiquidRho() { return cfg_.lookup("vdwFluid.liquidRho"); }
	float vdwVaporRho() { return cfg_.lookup("vdwFluid.vaporRho"); }
	float vdwPWE() { return cfg_.lookup("vdwFluid.we"); }
	float vdwCv() { return cfg_.lookup("vdwFluid.cv"); }

	float velConvergeTol() { return cfg_.lookup("simpleAlgorithm.velocityConvergenceTolerance"); }
	float rhoConvergeTol() { return cfg_.lookup("simpleAlgorithm.rhoConvergenceTolerance"); }
	float rhoRelaxCoef() { return cfg_.lookup("simpleAlgorithm.rhoRelaxCoefficient"); }

	float gravity() { return cfg_.lookup("envConstants.gravity"); }
	float heatDiffuse() { return cfg_.lookup("envConstants.heatDiffuseSpeed"); }

	float startTheta() { return cfg_.lookup("boilingProcess.startTheta"); }
	float targetTheta() { return cfg_.lookup("boilingProcess.targetTheta"); }
	float heatSpeed() { return cfg_.lookup("boilingProcess.heatSpeed"); }
	int heaterIntervalX() { return cfg_.lookup("boilingProcess.heaterPositions.intervalX"); }
	int heaterIntervalY() { return cfg_.lookup("boilingProcess.heaterPositions.intervalY"); }
	int heaterIntervalZ() { return cfg_.lookup("boilingProcess.heaterPositions.intervalZ"); }

	std::string fieldOutputDir() { return (const char*)cfg_.lookup("debug.outputDir"); }
	std::string snapshotOutputDir() { return (const char*)cfg_.lookup("debug.snapshotDir"); }
	int snapshotInterval() { return cfg_.lookup("debug.snapshotInterval"); }
private:
	libconfig::Config cfg_;
};