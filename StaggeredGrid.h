#pragma once

#include "TimeStepController.h"
#include "vec.h"

class Config;

enum CELL_TYPE {
	FLUID = 0,
	SOLID = 1,
	AIR = 2,
};

enum LEVEL_SET_INFO {
	UNKNOWN = 0,
	IMMEDIATE_NEIGHBOR = 1,
	KNOWN = 2,
};

enum SIMULATION_MODE {
	SMOKE = 0,
	WATER = 1,
};

class StaggeredGrid {
public:
	StaggeredGrid(Config *config);
	~StaggeredGrid();

	StaggeredGrid(bool test);

	void step(real _dt);
	void stepLevelSet(real _dt);
	void setLeftBoundaryVelocity(int v);
	void addSingleExplosion();
	void dumpPreview(real *field, int No);
	void dumpSlicePreview(int No, real * field, int z, float scale);
	void dumpSlicePreview(int No, CELL_TYPE * field, int z);
	void addWater(); // fill the [left, rightmost] box region with water.
	void test();
	void run();
	void runWater();

private:
	SIMULATION_MODE	mode = SMOKE;

	// grid dimension
	int resX, resY, resZ;
	int totalCells, slabSize;
	int totalVX, totalVY, totalVZ;
	
	// simulation control
	real totalTime;
	int totalSteps;
	real dx, dt, buoyancy;
	real threshold;
	int iterations;

	TimeStepController *timeStep;

	// fields
	real *velocityX, *velocityY, *velocityZ;
	real *velocityXOld, *velocityYOld, *velocityZOld;
	real *concentration, *concentrationOld;
	real *heat, *heatOld;
	real *forceX, *forceY, *forceZ;
	real *signedDistanceField, *signedDistanceFieldOld;

	real gravity;
	
	CELL_TYPE *cellType;
	LEVEL_SET_INFO *levelSetInfo;
	Vec3f **closestPoint;
	int *closestPointIndex;

	int leftBoundaryVelocity = 0;
	real leftBoundaryPosition = 0;

	Config *config;

private:
	//methods
	void applyForce(real dt, real *vx, real *vy, real *vz);
	void addGravity();
	void addBuoyancy(real *heat);
	void advectMacCormack(real dt);

	void advectFieldMacCormack(const real dt, real *vx, real *vy, real *vz, real *fieldOld, real*fieldNew);
	void advectFieldSemiLagrange(const real dt, real *vx, real *vy, real*vz, real *fieldOld, real *fieldNew, bool clampExtrema);
	void clampOutsideRays(const real dt, real *vx, real *vy, real *vz, real* field, const real *interimField);

	void advectVelocitySemiLagrange(real dt, real *vx, real *vy, real *vz, real *vxNew, real *vyNew, real *vzNew);
	void advectVelocitySemiLagrange(const real dt, real * vxBackground, real * vyBackground, real * vzBackground, real * vxInterim, real * vyInterim, real * vzInterim, real * vxNew, real * vyNew, real * vzNew);
	void project();

	void genA(real *Adiag, real *Aplusi, real *Aplusj, real *Aplusk, real scale); // tested
	void genPrecon(real *Adiag, real *Aplusi, real *Aplusj, real *Aplusk, real *precon); //tested
	void applyPrecon(real *Aplusi, real *Aplusj, real *Aplusk, real *precon, real *r, real *pcg_z);
	void applyA(real *Adiag, real *Aplusi, real *Aplusj, real *Aplusk, real *s, real *t); // tested
	void solvePressure(real *Adiag, real *Aplusi, real *Aplusj, real *Aplusk, real *precon, real *pressure, real *divergence);
	
	
	void updateCellType(real *sdf);

	// from a distorted signed distance, compute a new signed distance.
	// input: SDF old; output: SDF new;
	void computeSignedDistanceField(real * sdf);

	// input: a sdf, output: mark cells near interface and give their closet point on the surface.
	// called by computeSDF();
	void findSurfaceCells(real *sdf, real *usdf);

	int isValidVelocity(int i, int j, int k, int direction);

	/* read the vXYZ field */
	void getAveragedVelocityAt(real x, real y, real z, real &vx, real &vy, real &vz);
	
	/* implied input: velocityXYZ, their values are modified*/
	void extrapolateVelocity();

	// generate A from a signed distance field
	void genA(real *signedDistanceField);
};

