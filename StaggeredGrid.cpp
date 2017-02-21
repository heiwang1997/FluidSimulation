#include "stdafx.h"
#include "StaggeredGrid.h"
#include "FieldManipulation.h"
#include "config.h"
#include "vec.h"
#include "mkl_cblas.h"
#include <string>


using std::cout;
using std::endl;

StaggeredGrid::StaggeredGrid(Config *_config) {
	config = _config;
	resX = config->resX();
	resY = config->resY();
	resZ = config->resZ();
	// set simulation constants
	buoyancy = config->buoyancy();
	gravity = -1;
	totalTime = 0.0f;
	totalSteps = 0;
	threshold = 1e-4f;
	iterations = config->iteration();

	dx = config->h();

	totalCells = resX * resY * resZ;
	slabSize = resX * resY;

	// scalar field stored in the center of each cell
	concentration = (real*)calloc(totalCells, sizeof(real));
	concentrationOld = (real*)calloc(totalCells, sizeof(real));
	heat = (real*)calloc(totalCells, sizeof(real));
	heatOld = (real*)calloc(totalCells, sizeof(real));
	rho = (real*)calloc(totalCells, sizeof(real));
	forceX = (real*)calloc(totalCells, sizeof(real));
	forceY = (real*)calloc(totalCells, sizeof(real));
	forceZ = (real*)calloc(totalCells, sizeof(real));

	// u stored in staggered positions
	totalVX = (resX + 1) * resY * resZ;
	totalVY = resX * (resY + 1) * resZ;
	totalVZ = resX * resY * (resZ + 1);
	velocityX = (real*)calloc(totalVX, sizeof(real));
	velocityXOld = (real*)calloc(totalVX, sizeof(real));
	velocityY = (real*)calloc(totalVY, sizeof(real));
	velocityYOld = (real*)calloc(totalVY, sizeof(real));
	velocityZ = (real*)calloc(totalVZ, sizeof(real));
	velocityZOld = (real*)calloc(totalVZ, sizeof(real));

	cellType = (CELL_TYPE*)calloc(totalCells, sizeof(CELL_TYPE));
	levelSetInfo = (LEVEL_SET_INFO*)calloc(totalCells, sizeof(LEVEL_SET_INFO));
	closestPoint = (Vec3f**)calloc(totalCells, sizeof(Vec3f*));
	for (int i = 0; i < totalCells; i++) {
		closestPoint[i] = new Vec3f(-1, -2, -3); // value set for debugging only.
	}
	signedDistanceField = (real*)calloc(totalCells, sizeof(real));
	signedDistanceFieldOld = (real*)calloc(totalCells, sizeof(real));
	closestPointIndex = (int*)calloc(totalCells, sizeof(int));
}

StaggeredGrid::~StaggeredGrid() {
	delete concentration;
	delete concentrationOld;
	delete heat;
	delete heatOld;
	delete forceX;
	delete forceY;
	delete forceZ;

	delete velocityX;
	delete velocityXOld;
	delete velocityY;
	delete velocityYOld;
	delete velocityZ;
	delete velocityZOld;

	delete cellType;
}

void StaggeredGrid::addSingleExplosion() {
	float xTotal = dx * resX;
	float yTotal = dx * resY;
	float zTotal = dx * resZ;

	float heighMin = 0.05f;
	float heighMax = 0.10f;
	float sourceRadius = 0.05f;

	for (int z = 0; z < resZ; z++) {
		for (int y = (int)(heighMin / dx); y <= (int)(heighMax / dx); y++) {
			for (int x = 0; x < resX; x++) {
				float xLength = (x + 0.5) * dx - xTotal * 0.5f;
				float zLength = (z + 0.5) * dx - zTotal * 0.5f;
				float radius = sqrtf(xLength * xLength + zLength * zLength);
				if (radius < sourceRadius) {
					int index = x + y * resX + z * slabSize;
					concentration[index] = 1.0f;
					heat[index] = 1.0f;
				}
			}
		}
	}
}

void StaggeredGrid::dumpPreview(real *field, int No) {
	const int nitems = resX * resY;
	const int otherDir = 2;
	float *buf = new float[nitems];
	float div = 1. / (float)resX; // normalize for shorter sides, old: res[otherDir];
	div *= 4.; //slightly increase contrast
	for (int i = 0; i<nitems; i++)
		buf[i] = 0.0f;
	// float scale = 5.0f;
	float scale = 100.0f;
	for (int x = 0; x < resX; x++) {
		for (int y = 0; y < resY; y++) {
			for (int z = 0; z < resZ; z++) {
				const int index = getIndex(x, y, z, resX, resY, resZ);
				const int bufindex = x + y * resX;
				buf[bufindex] += field[index] * div * scale;
			}
		}
	}
	dumpNumberedPNG(No, "result/front", buf, resX, resY);
	delete[] buf;
}

void StaggeredGrid::dumpSlicePreview(int No, real * field, int z, float scale) {
	const int nitems = resX * resY;
	const int otherDir = 2;
	float *buf = new float[nitems];
	for (int i = 0; i<nitems; i++)
		buf[i] = 0.0f;
	for (int x = 0; x < resX; x++) {
		for (int y = 0; y < resY; y++) {
			const int index = getIndex(x, y, z, resX, resY, resZ);
			const int bufindex = x + y * resX;
			buf[bufindex] = field[index] * scale;
		}
	}
	dumpNumberedPNG(No, "result/slice", buf, resX, resY);
	delete[] buf;
}

void StaggeredGrid::dumpSlicePreview(int No, CELL_TYPE * field, int z) {
	const int nitems = resX * resY;
	const int otherDir = 2;
	float *buf = new float[nitems];
	for (int i = 0; i<nitems; i++)
		buf[i] = 0.0f;
	for (int x = 0; x < resX; x++) {
		for (int y = 0; y < resY; y++) {
			const int index = getIndex(x, y, z, resX, resY, resZ);
			const int bufindex = x + y * resX;
			if (field[index] == SOLID) {
				buf[bufindex] = 1.0f;
			}
			else if (field[index] == FLUID) {
				buf[bufindex] = -1.0f;
			}
			else {
				buf[bufindex] = 0;
			}
		}
	}
	dumpNumberedPNG(No, "result/slice", buf, resX, resY);
	delete[] buf;
}

void StaggeredGrid::addWater() {
	mode = WATER;
	setBorderType(cellType, resX, resY, resZ, SOLID);
	for (int i = 1; i < resX - 1; i++) {
		for (int j = 1; j < resY - 1; j++) {
			for (int k = 1; k < resZ - 1; k++) {
				int index = getIndex(i, j, k, resX, resY, resZ);
				signedDistanceField[index] = j - i;
			}
		}
	}
	computeSignedDistanceField(signedDistanceField);
}

void StaggeredGrid::run() {
	timeStep = new TimeStepController(config->totalFrame(), config->gridFPS(), config->gridDt());

	addSingleExplosion();

	while (!timeStep->isFinished()) {
		step(timeStep->getStepDt());
		int fmCnt;
		if (timeStep->isFrameTime(fmCnt)) {
			dumpPreview(concentration, totalSteps);
			addSingleExplosion();
		}
		std::cout << "This step elapsed *** seconds." << std::endl;
	}
	delete timeStep;
	std::cout << "Total elapsed *** seconds." << std::endl;
}

void StaggeredGrid::runWater() {
	timeStep = new TimeStepController(config->totalFrame(), config->gridFPS(), config->gridDt());

	cout << "Simulation started." << endl;
	addWater();

	while (!timeStep->isFinished()) {
		stepLevelSet(timeStep->getStepDt());
		int fmCnt;
		if (timeStep->isFrameTime(fmCnt)) {
#ifndef _DEBUG
			dumpSlicePreview(fmCnt, cellType, resZ / 2);
#endif
		}
		std::cout << "Step finished" << std::endl;
	}
	delete timeStep;
}

void StaggeredGrid::step(real _dt) {
	// wipe forces
	for (int i = 0; i < totalCells; i++)
		forceX[i] = forceY[i] = forceZ[i] = 0.0f;

	// wipe boundaries
	setBorderValue(concentration, resX, resY, resZ, 0);
	setBorderValue(heat, resX, resY, resZ, 0);

	// advect everything
	advectMacCormack(_dt);


	advectVelocitySemiLagrange(_dt / dx, velocityX, velocityY, velocityZ, velocityXOld, velocityYOld, velocityZOld);
	real *temp;
	temp = velocityX; velocityX = velocityXOld; velocityXOld = temp;
	temp = velocityY; velocityY = velocityYOld; velocityYOld = temp;
	temp = velocityZ; velocityZ = velocityZOld; velocityZOld = temp;

	std::cout << "advection used *** seconds." << std::endl;

	// run the solvers
	addBuoyancy(heat);

	applyForce(_dt, velocityX, velocityY, velocityZ);

	project();

	totalTime += _dt;
	totalSteps++;
}

void StaggeredGrid::stepLevelSet(real _dt) {
	// wipe forces
	for (int i = 0; i < totalCells; i++)
		forceX[i] = forceY[i] = forceZ[i] = 0.0f;
	cout << "initalize SDF" << endl;
	computeSignedDistanceField(signedDistanceField);

	cout << "extrapolate velocity field" << endl;
	extrapolateVelocity();

	real *vxInterim = new real[totalVX];
	real *vyInterim = new real[totalVY];
	real *vzInterim = new real[totalVZ];
	memcpy(velocityXOld, velocityX, totalVX * sizeof(real));
	memcpy(velocityYOld, velocityY, totalVY * sizeof(real));
	memcpy(velocityZOld, velocityZ, totalVZ * sizeof(real));
	memcpy(vxInterim, velocityX, totalVX * sizeof(real));
	memcpy(vyInterim, velocityY, totalVY * sizeof(real));
	memcpy(vzInterim, velocityZ, totalVZ * sizeof(real));

	addGravity();
	applyForce(_dt, vxInterim, vyInterim, vzInterim);
	cout << "advect SDF" << endl;
	// advectSDF in v_old
	real *temp;
	temp = signedDistanceField; signedDistanceField = signedDistanceFieldOld; signedDistanceFieldOld = temp;
	advectFieldMacCormack(_dt / dx, velocityXOld, velocityYOld, velocityZOld, signedDistanceFieldOld, signedDistanceField);
	computeSignedDistanceField(signedDistanceField);

	cout << "advect velocity" << endl;
	// advectVelocitySemiLagrange(_dt / dx, velocityXOld, velocityYOld, velocityZOld, velocityXOld, velocityYOld, velocityZOld);
	advectVelocitySemiLagrange(_dt / dx, velocityXOld, velocityYOld, velocityZOld,
		vxInterim, vyInterim, vzInterim,
		velocityX, velocityY, velocityZ);

	cout << "project" << endl;
	// run the solvers
	project();

	totalTime += _dt;
	totalSteps++;
	delete[] vxInterim;
	delete[] vyInterim;
	delete[] vzInterim;
}


void StaggeredGrid::applyForce(real dt, real *vx, real *vy, real *vz) {
	for (int z = 0; z < resZ; z++) {
		for (int y = 0; y < resY; y++) {
			for (int x = 0; x < resX; x++) {
				int index = getIndex(x, y, z, resX, resY, resZ);
				if (cellType[index] != FLUID) {
					continue;
				}
				int indexLeft = getIndex(x, y, z, resX, resY, resZ);
				int indexRight = getIndex(x + 1, y, z, resX, resY, resZ);
				int indexBelow = getIndex(x, y, z, resX, resY, resZ);
				int indexAbove = getIndex(x, y + 1, z, resX, resY, resZ);
				int indexFront = getIndex(x, y, z, resX, resY, resZ);
				int indexBack = getIndex(x, y, z + 1, resX, resY, resZ);

				if (cellType[indexLeft] != SOLID) {
					vx[getIndex(x, y, z, resX + 1, resY, resZ)] += dt * forceX[index] * 0.5f;
				}
				if (cellType[indexRight] != SOLID) {
					vx[getIndex(x + 1, y, z, resX + 1, resY, resZ)] += dt * forceX[index] * 0.5f;
				}

				if (cellType[indexBelow] != SOLID) {
					vy[getIndex(x, y, z, resX, resY + 1, resZ)] += dt * forceY[index] * 0.5f;
				}
				if (cellType[indexAbove] != SOLID) {
					vy[getIndex(x, y + 1, z, resX, resY + 1, resZ)] += dt * forceY[index] * 0.5f;
				}

				if (cellType[indexFront] != SOLID) {
					vz[getIndex(x, y, z, resX, resY, resZ + 1)] += dt * forceZ[index] * 0.5f;
				}
				if (cellType[indexBack] != SOLID) {
					vz[getIndex(x, y, z + 1, resX, resY, resZ + 1)] += dt * forceZ[index] * 0.5f;
				}
			}
		}
	}
}

void StaggeredGrid::addGravity() {
	for (int i = 0; i < totalCells; i++) {
		forceY[i] += gravity;
	}
}

void StaggeredGrid::addBuoyancy(real *field)
{
	int index = 0;

	real beta = buoyancy;
	if (beta == 0.)
		return;
	for (int z = 0; z < resZ; z++)
		for (int y = 0; y < resY; y++)
			for (int x = 0; x < resX; x++, index++)
				//for heat buoyancy
				//0.5f = ambient temperature , alpha = 0.1
				forceY[index] += beta * (field[index] - 0.5f) - 0.1 * concentration[index];
}

void StaggeredGrid::advectMacCormack(real _dt) {
	real *temp;
	temp = concentration; concentration = concentrationOld; concentrationOld = temp;
	temp = heat; heat = heatOld; heatOld = temp;

	/* space change */
	const float dt0 = _dt / dx;
	advectFieldMacCormack(dt0, velocityX, velocityY, velocityZ, concentrationOld, concentration);
	advectFieldMacCormack(dt0, velocityX, velocityY, velocityZ, heatOld, heat);

	setBorderValue(concentration, resX, resY, resZ, 0.0f);
	setBorderValue(heat, resX, resY, resZ, 0.0f);
}

void StaggeredGrid::advectFieldMacCormack(const real dt, real *vx, real *vy, real *vz, real *fieldOld, real *fieldNew) {

	float* phiHatN = (real*)calloc(totalCells, sizeof(real));
	float* phiHatN1 = (real*)calloc(totalCells, sizeof(real));

	for (int i = 0; i < totalCells; i++)
		phiHatN[i] = phiHatN1[i] = fieldOld[i];

	float *phiN = fieldOld;
	float *phiN1 = fieldNew;

	// phiHatN1 = A(phiN)
	advectFieldSemiLagrange(dt, velocityX, velocityY, velocityZ, phiN, phiHatN1, false);
	// phiHatN = A^R(phiHatN1)
	advectFieldSemiLagrange(-1.0*dt, velocityX, velocityY, velocityZ, phiHatN1, phiHatN, false);

	// phiN1 = phiHatN1 + (phiN - phiHatN) / 2
	const int border = 0;
	for (int z = border; z < resX - border; z++)
		for (int y = border; y < resY - border; y++)
			for (int x = border; x < resZ - border; x++) {
				int index = getIndex(x, y, z, resX, resY, resZ);
				phiN1[index] = phiHatN1[index] + (phiN[index] - phiHatN[index]) * 0.50f;
			}

	// clamp any newly created extrema
	bool clampExtrema = true;
	advectFieldSemiLagrange(dt, velocityX, velocityY, velocityZ, fieldOld, fieldNew, clampExtrema);

	// if the error estimate was bad, revert to first order
	clampOutsideRays(dt, velocityX, velocityY, velocityZ, fieldNew, phiHatN1);

	delete phiHatN;
	delete phiHatN1;
}

// tested for advection
void StaggeredGrid::advectFieldSemiLagrange(const real dt, real *vx, real *vy, real *vz,
	real *fieldOld, real *fieldNew, bool clampExtrema = false) {
	// scale dt up to grid resolution
	for (int z = 1; z < resZ - 1; z++) {
		for (int y = 1; y < resY - 1; y++) {
			for (int x = 1; x < resX - 1; x++) {
				int index = getIndex(x, y, z, resX, resY, resZ);

				int vxLeftIndex = getIndex(x, y, z, resX + 1, resY, resZ);
				int vxRightIndex = getIndex(x + 1, y, z, resX + 1, resY, resZ);
				real xTrace = x - dt * (vx[vxLeftIndex] + vx[vxRightIndex]) * 0.5f;

				int vyBelowIndex = getIndex(x, y, z, resX, resY + 1, resZ);
				int vyAboveIndex = getIndex(x, y + 1, z, resX, resY + 1, resZ);
				real yTrace = y - dt * (vy[vyBelowIndex] + vy[vyAboveIndex]) * 0.5f;

				int vzFrontIndex = getIndex(x, y, z, resX, resY, resZ + 1);
				int vzBackIndex = getIndex(x, y, z + 1, resX, resY, resZ + 1);
				real zTrace = z - dt * (vz[vzFrontIndex] + vz[vzBackIndex]) * 0.5f;
				// backtrace

				// clamp backtrace to grid boundaries
				if (!loopBoundary) {
					if (xTrace < 0.5) xTrace = 0.5;
					if (xTrace > resX - 1.5) xTrace = resX - 1.5;
					if (yTrace < 0.5) yTrace = 0.5;
					if (yTrace > resY - 1.5) yTrace = resY - 1.5;
					if (zTrace < 0.5) zTrace = 0.5;
					if (zTrace > resZ - 1.5) zTrace = resZ - 1.5;
				}

#ifndef IFLOOR
				// locate neighbors to interpolate
				const int x0 = (int)xTrace;
				const int x1 = x0 + 1;
				const int y0 = (int)yTrace;
				const int y1 = y0 + 1;
				const int z0 = (int)zTrace;
				const int z1 = z0 + 1;

				// get interpolation weights
				const float s1 = xTrace - (int)xTrace;
				const float s0 = 1.0f - s1;
				const float t1 = yTrace - (int)yTrace;
				const float t0 = 1.0f - t1;
				const float u1 = zTrace - (int)zTrace;
				const float u0 = 1.0f - u1;
#else

				// locate neighbors to interpolate
				// due to unexpectedly large value of v*, the trace back position may be far away from [0, resX], 
				const int x0 = (ifloor(xTrace) + resX) % resX;
				const int x1 = (x0 + 1) % resX;
				const int y0 = (ifloor(yTrace) + resY) % resY;
				const int y1 = (y0 + 1) % resY;
				const int z0 = (ifloor(zTrace) + resZ) % resZ;
				const int z1 = (z0 + 1) % resZ;

				// get interpolation weights
				const float s1 = xTrace - floor(xTrace);
				const float s0 = 1.0f - s1;
				const float t1 = yTrace - floor(yTrace);
				const float t0 = 1.0f - t1;
				const float u1 = zTrace - floor(zTrace);
				const float u0 = 1.0f - u1;
#endif 
				const int i000 = getIndex(x0, y0, z0, resX, resY, resZ);
				const int i010 = getIndex(x0, y1, z0, resX, resY, resZ);
				const int i100 = getIndex(x1, y0, z0, resX, resY, resZ);
				const int i110 = getIndex(x1, y1, z0, resX, resY, resZ);
				const int i001 = getIndex(x0, y0, z1, resX, resY, resZ);
				const int i011 = getIndex(x0, y1, z1, resX, resY, resZ);
				const int i101 = getIndex(x1, y0, z1, resX, resY, resZ);
				const int i111 = getIndex(x1, y1, z1, resX, resY, resZ);

				// interpolate
				if (!clampExtrema) {
					fieldNew[index] = u0 * (s0 * (t0 * fieldOld[i000] +
						t1 * fieldOld[i010]) +
						s1 * (t0 * fieldOld[i100] +
							t1 * fieldOld[i110])) +
						u1 * (s0 * (t0 * fieldOld[i001] +
							t1 * fieldOld[i011]) +
							s1 * (t0 * fieldOld[i101] +
								t1 * fieldOld[i111]));
				}
				else {
					real minValue = min(fieldOld[i000], fieldOld[i001], fieldOld[i010], fieldOld[i011],
						fieldOld[i100], fieldOld[i101], fieldOld[i110], fieldOld[i111]);
					real maxValue = max(fieldOld[i000], fieldOld[i001], fieldOld[i010], fieldOld[i011],
						fieldOld[i100], fieldOld[i101], fieldOld[i110], fieldOld[i111]);
					fieldNew[index] = (fieldNew[index] > maxValue) ? maxValue : fieldNew[index];
					fieldNew[index] = (fieldNew[index] < minValue) ? minValue : fieldNew[index];
				}
			}
		}
	}
}

void StaggeredGrid::clampOutsideRays(const real dt, real *vx, real *vy, real *vz, real *fieldNew, const real *interimField) {
	for (int z = 1; z < resZ - 1; z++) {
		for (int y = 1; y < resY - 1; y++) {
			for (int x = 1; x < resX - 1; x++) {
				const int index = getIndex(x, y, z, resX, resY, resZ);
				// backtrace

				int vxLeftIndex = getIndex(x, y, z, resX + 1, resY, resZ);
				int vxRightIndex = getIndex(x + 1, y, z, resX + 1, resY, resZ);
				real vxdt = dt * (vx[vxLeftIndex] + vx[vxRightIndex]) * 0.5f;

				int vyBelowIndex = getIndex(x, y, z, resX, resY + 1, resZ);
				int vyAboveIndex = getIndex(x, y + 1, z, resX, resY + 1, resZ);
				real vydt = dt * (vy[vyBelowIndex] + vy[vyAboveIndex]) * 0.5f;

				int vzFrontIndex = getIndex(x, y, z, resX, resY, resZ + 1);
				int vzBackIndex = getIndex(x, y, z + 1, resX, resY, resZ + 1);
				real vzdt = dt * (vz[vzFrontIndex] + vz[vzBackIndex]) * 0.5f;

				float xBackward = x + vxdt;
				float yBackward = y + vydt;
				float zBackward = z + vzdt;
				float xTrace = x - vxdt;
				float yTrace = y - vydt;
				float zTrace = z - vzdt;

				// see if it goes outside the boundaries
				// todo
				bool hasObstacle =
					(zTrace < 1.0f) || (zTrace > resZ - 2.0f) ||
					(yTrace < 1.0f) || (yTrace > resY - 2.0f) ||
					(xTrace < 1.0f) || (xTrace > resX - 2.0f) ||
					(zBackward < 1.0f) || (zBackward > resZ - 2.0f) ||
					(yBackward < 1.0f) || (yBackward > resY - 2.0f) ||
					(xBackward < 1.0f) || (xBackward > resX - 2.0f);
				// end todo
				// reuse old advection instead of doing another one...
				if (hasObstacle) {
					fieldNew[index] = interimField[index];
				}
			}
		}
	}
}



/* advect vx, vy, vz inside itself, new values are stored in vNew */
void StaggeredGrid::advectVelocitySemiLagrange(const real dt, real * vx, real * vy, real * vz, real * vxNew, real * vyNew, real * vzNew) {
	advectVelocitySemiLagrange(dt, vx, vy, vz, vx, vy, vz, vxNew, vyNew, vzNew);
}

/* advect vx, vy, vz inside background value, new values are stored in vNew */
void StaggeredGrid::advectVelocitySemiLagrange(const real dt,
	real * vxBackground, real * vyBackground, real * vzBackground,
	real * vxInterim, real *vyInterim, real *vzInterim,
	real * vxNew, real * vyNew, real * vzNew) {
	// velocity X
	for (int z = 1; z < resZ - 1; z++) {
		for (int y = 1; y < resY - 1; y++) {
			for (int x = 2; x < resX - 1; x++) {
				int index = getIndex(x, y, z, resX + 1, resY, resZ);
				float velx = vxBackground[index];
				float vely = (
					vyBackground[getIndex(x - 1, y, z, resX, resY + 1, resZ)] +
					vyBackground[getIndex(x - 1, y + 1, z, resX, resY + 1, resZ)] +
					vyBackground[getIndex(x, y, z, resX, resY + 1, resZ)] +
					vyBackground[getIndex(x, y + 1, z, resX, resY + 1, resZ)])*0.25f;
				float velz = (
					vzBackground[getIndex(x, y, z, resX, resY, resZ + 1)] +
					vzBackground[getIndex(x - 1, y, z, resX, resY, resZ + 1)] +
					vzBackground[getIndex(x, y, z + 1, resX, resY, resZ + 1)] +
					vzBackground[getIndex(x - 1, y, z + 1, resX, resY, resZ + 1)])*0.25f;

				// backtrace
				float xTrace = x - dt * velx;
				float yTrace = y - dt * vely;
				float zTrace = z - dt * velz;

				// clamp backtrace to grid boundaries
				if (!loopBoundary) {
					if (xTrace < 1.0) xTrace = 1.0f;
					if (xTrace > resX - 1.0) xTrace = resX - 1.0f;
					if (yTrace < 0.5f) yTrace = 0.5f;
					if (yTrace > resY - 1.5) yTrace = resY - 1.5f;
					if (zTrace < 0.5f) zTrace = 0.5f;
					if (zTrace > resZ - 1.5) zTrace = resZ - 1.5f;
				}
				else {
					assert(false);
				}
#ifndef IFLOOR
				// locate neighbors to interpolate
				const int x0 = (int)xTrace;
				const int x1 = x0 + 1;
				const int y0 = (int)yTrace;
				const int y1 = y0 + 1;
				const int z0 = (int)zTrace;
				const int z1 = z0 + 1;

				// get interpolation weights
				const float s1 = xTrace - (int)xTrace;
				const float s0 = 1.0f - s1;
				const float t1 = yTrace - (int)yTrace;
				const float t0 = 1.0f - t1;
				const float u1 = zTrace - (int)zTrace;
				const float u0 = 1.0f - u1;

#else
				// locate neighbors to interpolate
				const int x0 = (ifloor(xTrace) + resX) % resX;
				const int x1 = (x0 + 1) % resX;
				const int y0 = (ifloor(yTrace) + resY) % resY;
				const int y1 = (y0 + 1) % resY;
				const int z0 = (ifloor(zTrace) + resZ) % resZ;
				const int z1 = (z0 + 1) % resZ;

				// get interpolation weights
				const float s1 = xTrace - floor(xTrace);
				const float s0 = 1.0f - s1;
				const float t1 = yTrace - floor(yTrace);
				const float t0 = 1.0f - t1;
				const float u1 = zTrace - floor(zTrace);
				const float u0 = 1.0f - u1;

#endif
				const int i000 = getIndex(x0, y0, z0, resX + 1, resY, resZ);
				const int i010 = getIndex(x0, y1, z0, resX + 1, resY, resZ);
				const int i100 = getIndex(x1, y0, z0, resX + 1, resY, resZ);
				const int i110 = getIndex(x1, y1, z0, resX + 1, resY, resZ);
				const int i001 = getIndex(x0, y0, z1, resX + 1, resY, resZ);
				const int i011 = getIndex(x0, y1, z1, resX + 1, resY, resZ);
				const int i101 = getIndex(x1, y0, z1, resX + 1, resY, resZ);
				const int i111 = getIndex(x1, y1, z1, resX + 1, resY, resZ);

				// interpolate

				vxNew[index] = u0 * (s0 * (t0 * vxInterim[i000] +
					t1 * vxInterim[i010]) +
					s1 * (t0 * vxInterim[i100] +
						t1 * vxInterim[i110])) +
					u1 * (s0 * (t0 * vxInterim[i001] +
						t1 * vxInterim[i011]) +
						s1 * (t0 * vxInterim[i101] +
							t1 * vxInterim[i111]));

			}
		}
	}
	// velocity Y
	for (int z = 1; z < resZ - 1; z++) {
		for (int y = 2; y < resY - 1; y++) {
			for (int x = 1; x < resX - 1; x++) {
				int index = getIndex(x, y, z, resX, resY + 1, resZ);
				float velx = (
					vxBackground[getIndex(x, y - 1, z, resX + 1, resY, resZ)] +
					vxBackground[getIndex(x, y, z, resX + 1, resY, resZ)] +
					vxBackground[getIndex(x + 1, y - 1, z, resX + 1, resY, resZ)] +
					vxBackground[getIndex(x + 1, y, z, resX + 1, resY, resZ)])*0.25f;
				float vely = vyBackground[index];
				float velz = (
					vzBackground[getIndex(x, y, z, resX, resY, resZ + 1)] +
					vzBackground[getIndex(x, y, z + 1, resX, resY, resZ + 1)] +
					vzBackground[getIndex(x, y - 1, z, resX, resY, resZ + 1)] +
					vzBackground[getIndex(x, y - 1, z + 1, resX, resY, resZ + 1)]) * 0.25f;

				// backtrace
				float xTrace = x - dt * velx;
				float yTrace = y - dt * vely;
				float zTrace = z - dt * velz;

				//// clamp backtrace to grid boundaries
				if (!loopBoundary) {
					if (xTrace < 0.5f) xTrace = 0.5f;
					if (xTrace > resX - 1.5) xTrace = resX - 1.5f;
					if (yTrace < 1.0f) yTrace = 1.0f;
					if (yTrace > resY - 1.0f) yTrace = resY - 1.0f;
					if (zTrace < 0.5f) zTrace = 0.5f;
					if (zTrace > resZ - 1.5) zTrace = resZ - 1.5f;
				}
#ifndef IFLOOR
				//// locate neighbors to interpolate
				const int x0 = (int)xTrace;
				const int x1 = x0 + 1;
				const int y0 = (int)yTrace;
				const int y1 = y0 + 1;
				const int z0 = (int)zTrace;
				const int z1 = z0 + 1;

				// get interpolation weights
				const float s1 = xTrace - (int)xTrace;
				const float s0 = 1.0f - s1;
				const float t1 = yTrace - (int)yTrace;
				const float t0 = 1.0f - t1;
				const float u1 = zTrace - (int)zTrace;
				const float u0 = 1.0f - u1;
#else
				const int x0 = (ifloor(xTrace) + resX) % resX;
				const int x1 = (x0 + 1) % resX;
				const int y0 = (ifloor(yTrace) + resY) % resY;
				const int y1 = (y0 + 1) % resY;
				const int z0 = (ifloor(zTrace) + resZ) % resZ;
				const int z1 = (z0 + 1) % resZ;

				// get interpolation weights
				const float s1 = xTrace - floor(xTrace);
				const float s0 = 1.0f - s1;
				const float t1 = yTrace - floor(yTrace);
				const float t0 = 1.0f - t1;
				const float u1 = zTrace - floor(zTrace);
				const float u0 = 1.0f - u1;
#endif
				const int i000 = getIndex(x0, y0, z0, resX, resY + 1, resZ);
				const int i010 = getIndex(x0, y1, z0, resX, resY + 1, resZ);
				const int i100 = getIndex(x1, y0, z0, resX, resY + 1, resZ);
				const int i110 = getIndex(x1, y1, z0, resX, resY + 1, resZ);
				const int i001 = getIndex(x0, y0, z1, resX, resY + 1, resZ);
				const int i011 = getIndex(x0, y1, z1, resX, resY + 1, resZ);
				const int i101 = getIndex(x1, y0, z1, resX, resY + 1, resZ);
				const int i111 = getIndex(x1, y1, z1, resX, resY + 1, resZ);

				// interpolate
				vyNew[index] = u0 * (s0 * (t0 * vyInterim[i000] +
					t1 * vyInterim[i010]) +
					s1 * (t0 * vyInterim[i100] +
						t1 * vyInterim[i110])) +
					u1 * (s0 * (t0 * vyInterim[i001] +
						t1 * vyInterim[i011]) +
						s1 * (t0 * vyInterim[i101] +
							t1 * vyInterim[i111]));
			}
		}
	}
	// velocity Z
	for (int z = 2; z < resZ - 1; z++) {
		for (int y = 1; y < resY - 1; y++) {
			for (int x = 1; x < resX - 1; x++) {
				int index = getIndex(x, y, z, resX, resY, resZ + 1);
				real velx = (
					vxBackground[getIndex(x, y, z - 1, resX + 1, resY, resZ)] +
					vxBackground[getIndex(x, y, z, resX + 1, resY, resZ)] +
					vxBackground[getIndex(x + 1, y, z - 1, resX + 1, resY, resZ)] +
					vxBackground[getIndex(x + 1, y, z, resX + 1, resY, resZ)]) *0.25f;
				real vely = (
					vyBackground[getIndex(x, y, z - 1, resX, resY + 1, resZ)] +
					vyBackground[getIndex(x, y + 1, z - 1, resX, resY + 1, resZ)] +
					vyBackground[getIndex(x, y, z, resX, resY + 1, resZ)] +
					vyBackground[getIndex(x, y + 1, z, resX, resY + 1, resZ)]) * 0.25f;
				real velz = vzBackground[index];

				// backtrace
				float xTrace = x - dt * velx;
				float yTrace = y - dt * vely;
				float zTrace = z - dt * velz;

				// clamp backtrace to grid boundaries
				if (!loopBoundary) {
					if (xTrace < 0.5f) xTrace = 0.5f;
					if (xTrace > resX - 1.5) xTrace = resX - 1.5f;
					if (yTrace < 0.5f) yTrace = 0.5f;
					if (yTrace > resY - 1.5) yTrace = resY - 1.5f;
					if (zTrace < 1.0f) zTrace = 1.0f;
					if (zTrace > resZ - 1.0) zTrace = resZ - 1.0f;
				}
#ifndef IFLOOR
				// locate neighbors to interpolate
				const int x0 = (int)xTrace;
				const int x1 = x0 + 1;
				const int y0 = (int)yTrace;
				const int y1 = y0 + 1;
				const int z0 = (int)zTrace;
				const int z1 = z0 + 1;

				// get interpolation weights
				const float s1 = xTrace - (int)xTrace;
				const float s0 = 1.0f - s1;
				const float t1 = yTrace - (int)yTrace;
				const float t0 = 1.0f - t1;
				const float u1 = zTrace - (int)zTrace;
				const float u0 = 1.0f - u1;

#else
				const int x0 = (ifloor(xTrace) + resX) % resX;
				const int x1 = (x0 + 1) % resX;
				const int y0 = (ifloor(yTrace) + resY) % resY;
				const int y1 = (y0 + 1) % resY;
				const int z0 = (ifloor(zTrace) + resZ) % resZ;
				const int z1 = (z0 + 1) % resZ;

				// get interpolation weights
				const float s1 = xTrace - floor(xTrace);
				const float s0 = 1.0f - s1;
				const float t1 = yTrace - floor(yTrace);
				const float t0 = 1.0f - t1;
				const float u1 = zTrace - floor(zTrace);
				const float u0 = 1.0f - u1;
#endif
				const int i000 = getIndex(x0, y0, z0, resX, resY, resZ + 1);
				const int i010 = getIndex(x0, y1, z0, resX, resY, resZ + 1);
				const int i100 = getIndex(x1, y0, z0, resX, resY, resZ + 1);
				const int i110 = getIndex(x1, y1, z0, resX, resY, resZ + 1);
				const int i001 = getIndex(x0, y0, z1, resX, resY, resZ + 1);
				const int i011 = getIndex(x0, y1, z1, resX, resY, resZ + 1);
				const int i101 = getIndex(x1, y0, z1, resX, resY, resZ + 1);
				const int i111 = getIndex(x1, y1, z1, resX, resY, resZ + 1);

				// interpolate
				// int index = (*this.*idx)(x, y, z);

				vzNew[index] = u0 * (s0 * (t0 * vzInterim[i000] +
					t1 * vzInterim[i010]) +
					s1 * (t0 * vzInterim[i100] +
						t1 * vzInterim[i110])) +
					u1 * (s0 * (t0 * vzInterim[i001] +
						t1 * vzInterim[i011]) +
						s1 * (t0 * vzInterim[i101] +
							t1 * vzInterim[i111]));
			}
		}
	}
}


void StaggeredGrid::project() {

	unsigned char *cellFlag = new unsigned char[totalCells];
	real *precon = new float[totalCells];
	real *pressure = new float[totalCells];

	memset(cellFlag, 0, sizeof(unsigned char)*totalCells);
	memset(precon, 0, sizeof(float)*totalCells);
	memset(pressure, 0, sizeof(float)*totalCells);

	setBorderType(cellType, resX, resY, resZ, SOLID);

	for (int x = 0; x < resX; x++) {
		for (int y = 0; y < resY; y++) {
			velocityZ[getIndex(x, y, 0, resX, resY, resZ + 1)] = 0;
			velocityZ[getIndex(x, y, 1, resX, resY, resZ + 1)] = 0;
			velocityZ[getIndex(x, y, resZ, resX, resY, resZ + 1)] = 0;
			velocityZ[getIndex(x, y, resZ - 1, resX, resY, resZ + 1)] = 0;

			velocityX[getIndex(x, y, 0, resX + 1, resY, resZ)] = velocityX[getIndex(x, y, 1, resX + 1, resY, resZ)];
			velocityY[getIndex(x, y, 0, resX, resY + 1, resZ)] = velocityY[getIndex(x, y, 1, resX, resY + 1, resZ)];
			velocityX[getIndex(x, y, resZ - 1, resX + 1, resY, resZ)] = velocityX[getIndex(x, y, resZ - 2, resX + 1, resY, resZ)];
			velocityY[getIndex(x, y, resZ - 1, resX, resY + 1, resZ)] = velocityY[getIndex(x, y, resZ - 2, resX, resY + 1, resZ)];
		}
	}
	for (int x = 0; x < resX; x++) {
		for (int z = 0; z < resZ; z++) {
			velocityY[getIndex(x, 0, z, resX, resY + 1, resZ)] = 0;
			velocityY[getIndex(x, 1, z, resX, resY + 1, resZ)] = 0;
			velocityY[getIndex(x, resY, z, resX, resY + 1, resZ)] = 0;
			velocityY[getIndex(x, resY - 1, z, resX, resY + 1, resZ)] = 0;

			velocityX[getIndex(x, 0, z, resX + 1, resY, resZ)] = velocityX[getIndex(x, 1, z, resX + 1, resY, resZ)];
			velocityZ[getIndex(x, 0, z, resX, resY, resZ + 1)] = velocityZ[getIndex(x, 1, z, resX, resY, resZ + 1)];
			velocityX[getIndex(x, resY - 1, z, resX + 1, resY, resZ)] = velocityX[getIndex(x, resY - 2, z, resX + 1, resY, resZ)];
			velocityZ[getIndex(x, resY - 1, z, resX, resY, resZ + 1)] = velocityZ[getIndex(x, resY - 2, z, resX, resY, resZ + 1)];
		}
	}
	for (int y = 0; y < resY; y++) {
		for (int z = 0; z < resZ; z++) {
			velocityX[getIndex(0, y, z, resX + 1, resY, resZ)] = 0;
			velocityX[getIndex(1, y, z, resX + 1, resY, resZ)] = 0;
			velocityX[getIndex(resX, y, z, resX + 1, resY, resZ)] = 0;
			velocityX[getIndex(resX - 1, y, z, resX + 1, resY, resZ)] = 0;

			velocityY[getIndex(0, y, z, resX, resY + 1, resZ)] = velocityY[getIndex(1, y, z, resX, resY + 1, resZ)];
			velocityZ[getIndex(0, y, z, resX, resY, resZ + 1)] = velocityZ[getIndex(1, y, z, resX, resY, resZ + 1)];
			velocityY[getIndex(resX - 1, y, z, resX, resY + 1, resZ)] = velocityY[getIndex(resX - 2, y, z, resX, resY + 1, resZ)];
			velocityZ[getIndex(resX - 1, y, z, resX, resY, resZ + 1)] = velocityZ[getIndex(resX - 2, y, z, resX, resY, resZ + 1)];
		}
	}
	// end set boundary

	real *Adiag = new real[totalCells]; memset(Adiag, 0, totalCells * sizeof(real));
	real *Aplusi = new real[totalCells]; memset(Aplusi, 0, totalCells * sizeof(real));
	real *Aplusj = new real[totalCells]; memset(Aplusj, 0, totalCells * sizeof(real));
	real *Aplusk = new real[totalCells]; memset(Aplusk, 0, totalCells * sizeof(real));
	cout << Adiag[totalCells - 1] << endl;
	genA(Adiag, Aplusi, Aplusj, Aplusk, 1.0f);
	// Boundaries are set to SOLID. Different from lodsmoke.
	real *preconA = new real[totalCells]; memset(preconA, 0, totalCells);
	genPrecon(Adiag, Aplusi, Aplusj, Aplusk, precon);
	//// get pre-conditioner
	real *divergence = new real[totalCells];
	memset(divergence, 0, sizeof(real)*totalCells);
	// calculate divergence
	for (int z = 1; z < resZ - 1; z++) {
		for (int y = 1; y < resY - 1; y++) {
			for (int x = 1; x < resX - 1; x++) {
				int index = getIndex(x, y, z, resX, resY, resZ);
				int indexRight = getIndex(x + 1, y, z, resX + 1, resY, resZ);
				int indexLeft = getIndex(x, y, z, resX + 1, resY, resZ);
				int indexAbove = getIndex(x, y + 1, z, resX, resY + 1, resZ);
				int indexBelow = getIndex(x, y, z, resX, resY + 1, resZ);
				int indexFront = getIndex(x, y, z, resX, resY, resZ + 1);
				int indexBack = getIndex(x, y, z + 1, resX, resY, resZ + 1);
				if (cellType[index] != FLUID) {
					divergence[index] = 0.0f;
				}
				else {
					divergence[index] = -(
						velocityX[indexRight]
						- velocityX[indexLeft]
						+ velocityY[indexAbove]
						- velocityY[indexBelow]
						+ velocityZ[indexBack]
						- velocityZ[indexFront]
						);
				}
			}
		}
	}

	double tt = cblas_dsdot(slabSize*resZ, divergence, 1, divergence, 1);
	std::cout << "2-norm of div :" << tt << std::endl;

	// solve Poisson equation
	// laplace(p*dt/rho) = div(u)
	// pressure = [0,...,0] now
	solvePressure(Adiag, Aplusi, Aplusj, Aplusk, precon, pressure, divergence);

	int index = slabSize + resX + 1;
	for (int z = 1; z < resZ - 1; z++, index += 2 * resX) {
		for (int y = 1; y < resY - 1; y++, index += 2) {
			for (int x = 1; x < resX - 1; x++, index++) {
				if (cellType[index] == FLUID) {
					velocityX[getIndex(x, y, z, resX + 1, resY, resZ)] -= pressure[index];
					velocityX[getIndex(x + 1, y, z, resX + 1, resY, resZ)] += pressure[index];
					velocityY[getIndex(x, y, z, resX, resY + 1, resZ)] -= pressure[index];
					velocityY[getIndex(x, y + 1, z, resX, resY + 1, resZ)] += pressure[index];
					velocityZ[getIndex(x, y, z, resX, resY, resZ)] -= pressure[index];
					velocityZ[getIndex(x, y, z + 1, resX, resY, resZ + 1)] += pressure[index];
				}
			}
		}
	}
	index = 0;
	for (int z = 0; z < resZ; z++) {
		for (int y = 0; y < resY; y++) {
			for (int x = 0; x < resX; x++, index++) {
				if (cellType[index] == SOLID) {
					velocityX[getIndex(x, y, z, resX + 1, resY, resZ)] = 0;
					velocityX[getIndex(x + 1, y, z, resX + 1, resY, resZ)] = 0;
					velocityY[getIndex(x, y, z, resX, resY + 1, resZ)] = 0;
					velocityY[getIndex(x, y + 1, z, resX, resY + 1, resZ)] = 0;
					velocityZ[getIndex(x, y, z, resX, resY, resZ + 1)] = 0;
					velocityZ[getIndex(x, y, z + 1, resX, resY, resZ + 1)] = 0;
				}
			}
		}
	}



	delete[] Adiag;
	delete[] Aplusi;
	delete[] Aplusj;
	delete[] Aplusk;
	delete[] preconA;
	delete[] divergence;
	delete[] precon;
	delete[] pressure;
}

void StaggeredGrid::genA(real * Adiag, real * Aplusi, real * Aplusj, real * Aplusk, real scale) {
	for (int i = 0; i < resX; i++) {
		for (int j = 0; j < resY; j++) {
			for (int k = 0; k < resZ; k++) {
				int index = getIndex(i, j, k, resX, resY, resZ);
				int rightIndex = index + 1;
				int upIndex = index + resX;
				int backIndex = index + slabSize;
				if (cellType[index] == FLUID && cellType[rightIndex] == FLUID) {
					Adiag[index] += scale;
					Adiag[rightIndex] += scale;
					Aplusi[index] = -scale;
				}
				else if (cellType[index] == FLUID && cellType[rightIndex] == AIR) {
					if (mode == WATER) {
						real mag = (signedDistanceField[index + 1] - signedDistanceField[index]) / signedDistanceField[index];
						assert(mag <= 0);
						if (mag < -1000) mag = -1000;
						Adiag[index] += -mag * scale;
					}
					else {
						Adiag[index] += scale;
					}
				}

				if (cellType[index] == FLUID && cellType[upIndex] == FLUID) {
					Adiag[index] += scale;
					Adiag[upIndex] += scale;
					Aplusj[index] = -scale;
				}
				else if (cellType[index] == FLUID && cellType[upIndex] == AIR) {
					if (mode = WATER) {
						real mag = (signedDistanceField[upIndex] - signedDistanceField[index]) / signedDistanceField[index];
						assert(mag <= 0);
						if (mag < -1000) mag = -1000;
						Adiag[index] += -mag * scale;
					}
					else {
						Adiag[index] += scale;
					}
				}

				if (cellType[index] == FLUID && cellType[backIndex] == FLUID) {
					Adiag[index] += scale;
					Adiag[backIndex] += scale;
					Aplusk[index] = -scale;
				}
				else if (cellType[index] == FLUID && cellType[backIndex] == AIR) {
					if (mode = WATER) {
						real mag = (signedDistanceField[backIndex] - signedDistanceField[index]) / signedDistanceField[index];
						assert(mag <= 0);
						if (mag < -1000) mag = -1000;
						Adiag[index] += -mag * scale;
					}
					else {
						Adiag[index] += scale;
					}
				}
			}
		}
	}
}

real sq(real t) {
	return t * t;
}

void StaggeredGrid::genPrecon(real * Adiag, real * Aplusi, real * Aplusj, real * Aplusk, real * precon) {
	double tau = 0.97, sigma = 0.25;
	for (int i = 0; i < resX; i++) {
		for (int j = 0; j < resY; j++) {
			for (int k = 0; k < resZ; k++) {
				int index = getIndex(i, j, k, resX, resY, resZ);
				int iMinus = index - 1;
				int jMinus = index - resX;
				int kMinus = index - slabSize;
				if (cellType[index] == FLUID) {
					double e = Adiag[index]
						- sq(Aplusi[iMinus] * precon[iMinus])
						- sq(Aplusj[jMinus] * precon[jMinus])
						- sq(Aplusk[kMinus] * precon[kMinus])
						- tau * (
							Aplusi[iMinus] * (Aplusj[iMinus] + Aplusk[iMinus]) * sq(precon[iMinus])
							+ Aplusj[jMinus] * (Aplusi[jMinus] + Aplusk[jMinus]) * sq(precon[jMinus])
							+ Aplusk[kMinus] * (Aplusi[kMinus] + Aplusj[kMinus]) * sq(precon[kMinus])
							);
					if (e < sigma) {
						e = Adiag[index];
					}
					precon[index] = 1 / sqrt(e);
				}
				else {
					precon[index] = 0.0f;
				}
			}
		}
	}
}

void StaggeredGrid::applyPrecon(real * Aplusi, real * Aplusj, real * Aplusk, real * precon, real * r, real * pcg_z) {
	int index, x, y, z;
	index = slabSize + resX + 1;
	float *q = (float *)calloc(totalCells, sizeof(float));
	for (z = 1; z < resZ - 1; z++, index += 2 * resX)
		for (y = 1; y < resY - 1; y++, index += 2)
			for (x = 1; x < resX - 1; x++, index++) {
				if (cellType[index] != FLUID)
					continue;
				float t = r[index] -
					Aplusi[index - 1] * precon[index - 1] * q[index - 1] -
					Aplusj[index - resX] * precon[index - resX] * q[index - resX] -
					Aplusk[index - slabSize] * precon[index - slabSize] * q[index - slabSize];
				q[index] = t * precon[index];
			}
	index = resX - 2 + (resY - 2)*resX + (resZ - 2)*slabSize;
	for (z = resZ - 2; z > 0; z--, index -= 2 * resX) {
		for (y = resY - 2; y > 0; y--, index -= 2) {
			for (x = resX - 2; x > 0; x--, index--) {
				if (cellType[index] != FLUID)
					continue;
				float t = q[index] -
					Aplusi[index] * precon[index] * pcg_z[index + 1] -
					Aplusj[index] * precon[index] * pcg_z[index + resX] -
					Aplusk[index] * precon[index] * pcg_z[index + slabSize];
				pcg_z[index] = t * precon[index];
			}
		}
	}
	free(q);
}

void StaggeredGrid::applyA(real * Adiag, real * Aplusi, real * Aplusj, real * Aplusk, real * s, real * t) {
	for (int i = 0; i < resX; i++) {
		for (int j = 0; j < resY; j++) {
			for (int k = 0; k < resZ; k++) {
				int index = getIndex(i, j, k, resX, resY, resZ);
				if (cellType[index] != FLUID) {
					t[index] = 0.0f;
					continue;
				}
				t[index] = Adiag[index] * s[index] +
					Aplusi[index] * s[index + 1] +
					Aplusj[index] * s[index + resX] +
					Aplusk[index] * s[index + slabSize] +
					Aplusi[index - 1] * s[index - 1] +
					Aplusj[index - resX] * s[index - resX] +
					Aplusk[index - slabSize] * s[index - slabSize];
			}
		}
	}
}

void StaggeredGrid::solvePressure(real * Adiag, real * Aplusi, real * Aplusj, real * Aplusk, real * precon, real * pressure, real * divergence) {
	float *x = pressure, *b = divergence;
	int i = 0;
	float *direction = (float *)calloc(totalCells, sizeof(float));
	float *z = (float *)calloc(totalCells, sizeof(float));
	float *residual = (float *)calloc(totalCells, sizeof(float));
	//1.  r = b - Ax
	applyA(Adiag, Aplusi, Aplusj, Aplusk, x, residual);
	cblas_saxpy(totalCells, -1.0f, b, 1, residual, 1);
	cblas_sscal(totalCells, -1.0f, residual, 1);
	//2. z = Mr, s = z
	// applyPrecon(cellFlag, precon, _residual, _z);
	applyPrecon(Aplusi, Aplusj, Aplusk, precon, residual, z);
	cblas_scopy(totalCells, z, 1, direction, 1);
	//3. sigma = z'r
	double sigma = cblas_dsdot(totalCells, z, 1, residual, 1);
	//3.  error = max|r|
	float err = fabs(residual[cblas_isamax(totalCells, residual, 1)]);
	static int cnt = 0;
	while (err>threshold && i<iterations) {
		// 4. z = Ad
		//tic.tictac();
		// applyA(cellFlag, _direction_, _z);
		applyA(Adiag, Aplusi, Aplusj, Aplusk, direction, z);
		// 5. alpha = sigma / d'z
		double dz = cblas_dsdot(totalCells, direction, 1, z, 1);
		double alpha = sigma / dz;
		// 6. x = x + alpha * d
		cblas_saxpy(totalCells, (float)alpha, direction, 1, x, 1);
		// 7. r = r - alpha * q
		cblas_saxpy(totalCells, (float)-alpha, z, 1, residual, 1);

		err = cblas_snrm2(totalCells, residual, 1);
		err = min(err, fabs(residual[cblas_isamax(totalCells, residual, 1)]));
		if (err <= threshold) break;
		// z = Mr
		//tic.tictac();
		// applyPrecon(cellFlag, precon, residual, z);
		applyPrecon(Aplusi, Aplusj, Aplusk, precon, residual, z);
		//cout <<(cnt)<<" applyPrecon: "<<tic.tictac()<<endl;
		// 10. sigma_new = z'r
		double sigma_new = cblas_dsdot(totalCells, z, 1, residual, 1);
		// beta = sigma_new / sigma
		double beta = sigma_new / sigma;
		// 11. d = z + beta * d

		cblas_sscal(totalCells, (float)beta, direction, 1);
		cblas_saxpy(totalCells, 1.0f, z, 1, direction, 1);
		// sigma = sigma_new
		sigma = sigma_new;
		i++;
	}

	std::cout << i << " iterations converged to " << err << " time used *** " << std::endl;

	free(direction);
	free(z);
	free(residual);
}

// update cell type according to sdf.
void StaggeredGrid::updateCellType(real * sdf) {
	for (int i = 0; i < totalCells; i++) {
		if (cellType[i] != SOLID) {
			if (sdf[i] > 0) {
				cellType[i] = AIR;
			}
			else {
				cellType[i] = FLUID;
			}
		}
	}
}


#include<set>
#include<vector>

struct IndexDistanceWrapper {
	IndexDistanceWrapper(int _i, real _d) : index(_i), distance(_d) { assert(_d >= 0); }
	int index = -1;
	real distance = -1.;
};

bool operator<(const IndexDistanceWrapper &lhs, const IndexDistanceWrapper &rhs) {
	if (lhs.distance < rhs.distance)
		return true;
	if (lhs.distance > rhs.distance)
		return false;
	if (lhs.index < rhs.index)
		return true;
	if (lhs.index > rhs.index)
		return false;
	// indentical
	return false;
	// cout << "operator " << lhs.index << ", " << lhs.distance << " : " << rhs.index << ", " << rhs.distance << endl;
	// assert(false);
}

real norm(Vec3f &vec) {
	real ans = sqrt(vec.v[0] * vec.v[0] + vec.v[1] * vec.v[1] + vec.v[2] * vec.v[2]);
	assert(ans > 0);
	return ans;
}

void vecCopy(Vec3f *s, Vec3f *t) {
	t->v[0] = s->v[0];
	t->v[1] = s->v[1];
	t->v[2] = s->v[2];
}

/* re-assign sdf according to the interface it implies */
void StaggeredGrid::computeSignedDistanceField(real *sdf) {
	for (int i = 0; i < totalCells; i++) {
		levelSetInfo[i] = UNKNOWN;
	}
	real *usdf = new real[totalCells];
	for (int i = 0; i < totalCells; i++) {
		usdf[i] = 10000.0f;
	}
	findSurfaceCells(sdf, usdf);
	updateCellType(sdf);

	std::set<IndexDistanceWrapper> priorityQueue;
	// construct a priority queue
	for (int i = 0; i < resX; i++) {
		for (int j = 0; j < resY; j++) {
			for (int k = 0; k < resZ; k++) {
				int index = getIndex(i, j, k, resX, resY, resZ);
				if (levelSetInfo[index] != IMMEDIATE_NEIGHBOR) {
					continue;
				}
				// real shortestDistance = 10000.0f;
				Vec3f pos(i, j, k);
				if (levelSetInfo[index - 1] == KNOWN) {
					real distance = norm(*closestPoint[index - 1] - pos);
					if (distance < usdf[index]) {
						usdf[index] = distance;
						closestPointIndex[index] = index - 1;
						assert(cellType[closestPointIndex[index]] != SOLID);
						vecCopy(closestPoint[index - 1], closestPoint[index]);
					}
				}
				if (levelSetInfo[index + 1] == KNOWN) {
					real distance = norm(*closestPoint[index + 1] - pos);
					if (distance < usdf[index]) {
						usdf[index] = distance;
						closestPointIndex[index] = index + 1;
						vecCopy(closestPoint[index + 1], closestPoint[index]);
						assert(cellType[closestPointIndex[index]] != SOLID);
					}
				}
				if (levelSetInfo[index - resX] == KNOWN) {
					real distance = norm(*closestPoint[index - resX] - pos);
					if (distance < usdf[index]) {
						usdf[index] = distance;
						closestPointIndex[index] = index - resX;
						vecCopy(closestPoint[index - resX], closestPoint[index]);
						assert(cellType[closestPointIndex[index]] != SOLID);
					}
				}
				if (levelSetInfo[index + resX] == KNOWN) {
					real distance = norm(*closestPoint[index + resX] - pos);
					if (distance < usdf[index]) {
						usdf[index] = distance;
						closestPointIndex[index] = index + resX;
						vecCopy(closestPoint[index + resX], closestPoint[index]);
						assert(cellType[closestPointIndex[index]] != SOLID);
					}
				}
				if (levelSetInfo[index - slabSize] == KNOWN) {
					real distance = norm(*closestPoint[index - slabSize] - pos);
					if (distance < usdf[index]) {
						usdf[index] = distance;
						closestPointIndex[index] = index - slabSize;
						vecCopy(closestPoint[index - slabSize], closestPoint[index]);
						assert(cellType[closestPointIndex[index]] != SOLID);
					}
				}
				if (levelSetInfo[index + slabSize] == KNOWN) {
					real distance = norm(*closestPoint[index + slabSize] - pos);
					if (distance < usdf[index]) {
						usdf[index] = distance;
						closestPointIndex[index] = index + slabSize;
						vecCopy(closestPoint[index + slabSize], closestPoint[index]);
						assert(cellType[closestPointIndex[index]] != SOLID);
					}
				}
				priorityQueue.insert(IndexDistanceWrapper(index, usdf[index]));
			}
		}
	}

	if (priorityQueue.empty()) {
		cout << "No interface" << endl;
		assert(false);
	}

	// at this point, every immediate neighbor has a well defined distance, a closest index and a closest point
	while (!priorityQueue.empty()) {
		IndexDistanceWrapper dequeueItem = *priorityQueue.begin();
		priorityQueue.erase(priorityQueue.begin());
		int index = dequeueItem.index;
		assert(usdf[index] == dequeueItem.distance);
		usdf[index] = dequeueItem.distance;
		levelSetInfo[index] = KNOWN;
		assert(cellType[index] != SOLID);

		// for each of its neighbors
		int neighborIndices[] = { index - 1, index + 1, index - resX, index + resX, index - slabSize, index + slabSize };
		// std::vector<int> neighborIndices
		for (int n = 0; n < 6; n++) {
			int neighborIndex = neighborIndices[n];
			if (cellType[neighborIndex] != SOLID && levelSetInfo[neighborIndex] != KNOWN) {
				int i, j, k;
				getPos(neighborIndex, resX, resY, resZ, i, j, k);
				if (levelSetInfo[neighborIndex] == IMMEDIATE_NEIGHBOR) {
					// update its  shortest distance and point and (maybe) re-insert it in the priority queue
					if (closestPointIndex[neighborIndex] == index) {
						assert(false);
					}
					else {
						assert(usdf[neighborIndex] < 9999.0f);
						real newDistance = norm(Vec3f(i, j, k) - *closestPoint[index]);
						if (newDistance < usdf[neighborIndex]) {
							auto iter = priorityQueue.find(IndexDistanceWrapper(neighborIndex, usdf[neighborIndex]));
							assert(iter != priorityQueue.end());
							priorityQueue.erase(iter);

							usdf[neighborIndex] = newDistance;
							closestPointIndex[neighborIndex] = closestPointIndex[index];
							vecCopy(closestPoint[index], closestPoint[neighborIndex]);
							assert(cellType[closestPointIndex[neighborIndex]] != SOLID);

							priorityQueue.insert(IndexDistanceWrapper(neighborIndex, usdf[neighborIndex]));
						}
						else {
							;
							// leave as it is;
						}
					}
				}
				else if (levelSetInfo[neighborIndex] == UNKNOWN) {
					// add its information and insert it to priority queue
					assert(usdf[neighborIndex] > 9999.0f);
					real newDistance = norm(Vec3f(i, j, k) - *closestPoint[index]);
					usdf[neighborIndex] = newDistance;
					closestPointIndex[neighborIndex] = closestPointIndex[index];
					vecCopy(closestPoint[index], closestPoint[neighborIndex]);
					assert(cellType[closestPointIndex[index]] != SOLID);

					priorityQueue.insert(IndexDistanceWrapper(neighborIndex, usdf[neighborIndex]));
					levelSetInfo[neighborIndex] = IMMEDIATE_NEIGHBOR;
				}
				else {
					assert(false);
				}
			}
		}

		real narrowBand = 555.0f;
		if (dequeueItem.distance > narrowBand) {
			// break;
			;
		}
	}
	cout << signedDistanceField[getIndex(0, 0, 3, resX, resY, resZ)] << endl;
	for (int i = 0; i < totalCells; i++) {
		assert(usdf[i] >= 0);
		if (cellType[i] == AIR) {
			sdf[i] = usdf[i];
		}
		else if (cellType[i] == FLUID) {
			sdf[i] = -usdf[i];
		}
		else {
			// solid
			;
		}
	}
	// dumpSlicePreview(11, usdf, 3, 0.1);
	assert(sdf == signedDistanceField);
	// dumpSlicePreview(12, sdf, 3, 1);
}

/* surface cells are marked as KNOWN, their unsigned shortest distance to surface is stored in usdf,
and corresponding surface point is stored in cloestPoint[i] */
void StaggeredGrid::findSurfaceCells(real *sdf, real *usdf) {
	bool interfaceFound = false;
	for (int i = 0; i < resX; i++) {
		for (int j = 0; j < resY; j++) {
			for (int k = 0; k < resZ; k++) {
				int index = getIndex(i, j, k, resX, resY, resZ);
				if (cellType[index] == SOLID) {
					continue;
				}
				bool isSurfacePoint = false;
				real shortestDistance = 10000.0;
				Vec3f *cp = closestPoint[index];
				// i, j, k fluid cell, i-1, j, k aircell
				if (cellType[index - 1] != SOLID && sdf[index - 1] * sdf[index] <= 0) {
					real distance = abs(sdf[index] / (sdf[index] - sdf[index - 1]));
					if (distance < shortestDistance) {
						shortestDistance = distance;
						isSurfacePoint = true;
						cp->v[0] = i - distance;
						cp->v[1] = j;
						cp->v[2] = k;
					}
				}
				if (cellType[index + 1] != SOLID && sdf[index + 1] * sdf[index] <= 0) {
					real distance = abs(sdf[index] / (sdf[index] - sdf[index + 1]));
					if (distance < shortestDistance) {
						shortestDistance = distance;
						isSurfacePoint = true;
						cp->v[0] = i + distance;
						cp->v[1] = j;
						cp->v[2] = k;
					}
				}
				if (cellType[index - resX] != SOLID && sdf[index - resX] * sdf[index] <= 0) {
					real distance = abs(sdf[index] / (sdf[index] - sdf[index - resX]));
					if (distance < shortestDistance) {
						shortestDistance = distance;
						isSurfacePoint = true;
						cp->v[0] = i;
						cp->v[1] = j - distance;
						cp->v[2] = k;
					}
				}
				if (cellType[index + resX] != SOLID && sdf[index + resX] * sdf[index] <= 0) {
					real distance = abs(sdf[index] / (sdf[index] - sdf[index + resX]));
					if (distance < shortestDistance) {
						shortestDistance = distance;
						isSurfacePoint = true;
						cp->v[0] = i;
						cp->v[1] = j + distance;
						cp->v[2] = k;
					}
				}
				if (cellType[index - slabSize] != SOLID && sdf[index - slabSize] * sdf[index] <= 0) {
					real distance = abs(sdf[index] / (sdf[index] - sdf[index - slabSize]));
					if (distance < shortestDistance) {
						shortestDistance = distance;
						isSurfacePoint = true;
						cp->v[0] = i;
						cp->v[1] = j;
						cp->v[2] = k - distance;
					}
				}
				if (cellType[index + slabSize] != SOLID && sdf[index + slabSize] * sdf[index] <= 0) {
					real distance = abs(sdf[index] / (sdf[index] - sdf[index + slabSize]));
					if (distance < shortestDistance) {
						shortestDistance = distance;
						isSurfacePoint = true;
						cp->v[0] = i;
						cp->v[1] = j;
						cp->v[2] = k + distance;
					}
				}
				if (isSurfacePoint) {
					interfaceFound = true;
					levelSetInfo[index] = KNOWN;
					usdf[index] = shortestDistance;
					closestPointIndex[index] = -1; // value for debugging purpose.
					if (cellType[index - 1] != SOLID && levelSetInfo[index - 1] != KNOWN)
						levelSetInfo[index - 1] = IMMEDIATE_NEIGHBOR;
					if (cellType[index + 1] != SOLID && levelSetInfo[index + 1] != KNOWN)
						levelSetInfo[index + 1] = IMMEDIATE_NEIGHBOR;
					if (cellType[index - resX] != SOLID && levelSetInfo[index - resX] != KNOWN)
						levelSetInfo[index - resX] = IMMEDIATE_NEIGHBOR;
					if (cellType[index + resX] != SOLID && levelSetInfo[index + resX] != KNOWN)
						levelSetInfo[index + resX] = IMMEDIATE_NEIGHBOR;
					if (cellType[index - slabSize] != SOLID && levelSetInfo[index - slabSize] != KNOWN)
						levelSetInfo[index - slabSize] = IMMEDIATE_NEIGHBOR;
					if (cellType[index + slabSize] != SOLID && levelSetInfo[index + slabSize] != KNOWN)
						levelSetInfo[index + slabSize] = IMMEDIATE_NEIGHBOR;
				}
			}
		}
	}
	assert(interfaceFound);
}

int StaggeredGrid::isValidVelocity(int i, int j, int k, int direction) {
	switch (direction) {
	case 0:
		// vx
		return (cellType[getIndex(i - 1, j, k, resX, resY, resZ)] == FLUID || cellType[getIndex(i, j, k, resX, resY, resZ)] == FLUID)
			? 1 : 0;
	case 1:
		return (cellType[getIndex(i, j - 1, k, resX, resY, resZ)] == FLUID || cellType[getIndex(i, j, k, resX, resY, resZ)] == FLUID)
			? 1 : 0;
	case 2:
		return (cellType[getIndex(i, j, k - 1, resX, resY, resZ)] == FLUID || cellType[getIndex(i, j, k, resX, resY, resZ)] == FLUID)
			? 1 : 0;
	default:
		assert(false);
		break;
	}
}

void StaggeredGrid::getAveragedVelocityAt(real x, real y, real z, real &vx, real &vy, real &vz) {
	real xInX = x + 0.5;
	real yInY = y + 0.5;
	real zInZ = z + 0.5;
	{
		// compute vx;
		// indices
		int xLeft = (int)xInX;
		int xRight = xLeft + 1;
		int xBelow = (int)y;
		int xAbove = xBelow + 1;
		int xFront = (int)z;
		int xBack = xFront + 1;
		// weights
		real r1 = xInX - xLeft;
		real r0 = 1 - r1;
		real s1 = y - xBelow;
		real s0 = 1 - s1;
		real t1 = z - xFront;
		real t0 = 1 - t1;

		int direction = 0, resVX = resX + 1, resVY = resY, resVZ = resZ;
		int w000 = isValidVelocity(xLeft, xBelow, xFront, direction);
		int w001 = isValidVelocity(xLeft, xBelow, xBack, direction);
		int w010 = isValidVelocity(xLeft, xAbove, xFront, direction);
		int w011 = isValidVelocity(xLeft, xAbove, xBack, direction);
		int w100 = isValidVelocity(xRight, xBelow, xFront, direction);
		int w101 = isValidVelocity(xRight, xBelow, xBack, direction);
		int w110 = isValidVelocity(xRight, xAbove, xFront, direction);
		int w111 = isValidVelocity(xRight, xAbove, xBack, direction);
		real vx000 = velocityX[getIndex(xLeft, xBelow, xFront, resVX, resVY, resVZ)];
		real vx100 = velocityX[getIndex(xRight, xBelow, xFront, resVX, resVY, resVZ)];
		real vx010 = velocityX[getIndex(xLeft, xAbove, xFront, resVX, resVY, resVZ)];
		real vx110 = velocityX[getIndex(xRight, xAbove, xFront, resVX, resVY, resVZ)];
		real vx001 = velocityX[getIndex(xLeft, xBelow, xBack, resVX, resVY, resVZ)];
		real vx101 = velocityX[getIndex(xRight, xBelow, xBack, resVX, resVY, resVZ)];
		real vx011 = velocityX[getIndex(xLeft, xAbove, xBack, resVX, resVY, resVZ)];
		real vx111 = velocityX[getIndex(xRight, xAbove, xBack, resVX, resVY, resVZ)];
		real weightSum =
			w000 * r0 * s0 * t0
			+ w001 * r0 * s0 * t1
			+ w010 * r0 * s1 * t0
			+ w011 * r0 * s1 * t1
			+ w100 * r1 * s0 * t0
			+ w101 * r1 * s0 * t1
			+ w110 * r1 * s1 * t0
			+ w111 * r1 * s1 * t1;
		real vSum =
			w000 * r0 * s0 * t0 * vx000
			+ w001 * r0 * s0 * t1 * vx001
			+ w010 * r0 * s1 * t0 * vx010
			+ w011 * r0 * s1 * t1 * vx011
			+ w100 * r1 * s0 * t0 * vx100
			+ w101 * r1 * s0 * t1 * vx101
			+ w110 * r1 * s1 * t0 * vx110
			+ w111 * r1 * s1 * t1 * vx111;
		if (weightSum == 0) {
			// cout << "interpolation failed" << endl;
			vx = 0;
		}
		else {
			vx = vSum / weightSum;
		}
	}
	{
		// compute vy;
		// indices
		int left = (int)x;
		int right = left + 1;
		int below = (int)yInY;
		int above = below + 1;
		int front = (int)z;
		int back = front + 1;
		// weights
		real r1 = x - left;
		real r0 = 1 - r1;
		real s1 = yInY - below;
		real s0 = 1 - s1;
		real t1 = z - front;
		real t0 = 1 - t1;

		int direction = 1, resVX = resX, resVY = resY + 1, resVZ = resZ;
		int w000 = isValidVelocity(left, below, front, direction);
		int w001 = isValidVelocity(left, below, back, direction);
		int w010 = isValidVelocity(left, above, front, direction);
		int w011 = isValidVelocity(left, above, back, direction);
		int w100 = isValidVelocity(right, below, front, direction);
		int w101 = isValidVelocity(right, below, back, direction);
		int w110 = isValidVelocity(right, above, front, direction);
		int w111 = isValidVelocity(right, above, back, direction);
		real vy000 = velocityY[getIndex(left, below, front, resVX, resVY, resVZ)];
		real vy100 = velocityY[getIndex(right, below, front, resVX, resVY, resVZ)];
		real vy010 = velocityY[getIndex(left, above, front, resVX, resVY, resVZ)];
		real vy110 = velocityY[getIndex(right, above, front, resVX, resVY, resVZ)];
		real vy001 = velocityY[getIndex(left, below, back, resVX, resVY, resVZ)];
		real vy101 = velocityY[getIndex(right, below, back, resVX, resVY, resVZ)];
		real vy011 = velocityY[getIndex(left, above, back, resVX, resVY, resVZ)];
		real vy111 = velocityY[getIndex(right, above, back, resVX, resVY, resVZ)];
		real weightSum =
			w000 * r0 * s0 * t0
			+ w001 * r0 * s0 * t1
			+ w010 * r0 * s1 * t0
			+ w011 * r0 * s1 * t1
			+ w100 * r1 * s0 * t0
			+ w101 * r1 * s0 * t1
			+ w110 * r1 * s1 * t0
			+ w111 * r1 * s1 * t1;
		real vSum =
			w000 * r0 * s0 * t0 * vy000
			+ w001 * r0 * s0 * t1 * vy001
			+ w010 * r0 * s1 * t0 * vy010
			+ w011 * r0 * s1 * t1 * vy011
			+ w100 * r1 * s0 * t0 * vy100
			+ w101 * r1 * s0 * t1 * vy101
			+ w110 * r1 * s1 * t0 * vy110
			+ w111 * r1 * s1 * t1 * vy111;
		if (weightSum == 0) {
			//cout << "interpolation failed" << endl;
			vy = 0;
		}
		else {
			vy = vSum / weightSum;
		}
	}
	{
		// compute vz;
		// indices
		int left = (int)x;
		int right = left + 1;
		int below = (int)y;
		int above = below + 1;
		int front = (int)zInZ;
		int back = front + 1;
		// weights
		real r1 = x - left;
		real r0 = 1 - r1;
		real s1 = y - below;
		real s0 = 1 - s1;
		real t1 = zInZ - front;
		real t0 = 1 - t1;

		int direction = 2, resVX = resX, resVY = resY, resVZ = resZ + 1;
		int w000 = isValidVelocity(left, below, front, direction);
		int w001 = isValidVelocity(left, below, back, direction);
		int w010 = isValidVelocity(left, above, front, direction);
		int w011 = isValidVelocity(left, above, back, direction);
		int w100 = isValidVelocity(right, below, front, direction);
		int w101 = isValidVelocity(right, below, back, direction);
		int w110 = isValidVelocity(right, above, front, direction);
		int w111 = isValidVelocity(right, above, back, direction);
		real vz000 = velocityZ[getIndex(left, below, front, resVX, resVY, resVZ)];
		real vz100 = velocityZ[getIndex(right, below, front, resVX, resVY, resVZ)];
		real vz010 = velocityZ[getIndex(left, above, front, resVX, resVY, resVZ)];
		real vz110 = velocityZ[getIndex(right, above, front, resVX, resVY, resVZ)];
		real vz001 = velocityZ[getIndex(left, below, back, resVX, resVY, resVZ)];
		real vz101 = velocityZ[getIndex(right, below, back, resVX, resVY, resVZ)];
		real vz011 = velocityZ[getIndex(left, above, back, resVX, resVY, resVZ)];
		real vz111 = velocityZ[getIndex(right, above, back, resVX, resVY, resVZ)];
		real weightSum =
			w000 * r0 * s0 * t0
			+ w001 * r0 * s0 * t1
			+ w010 * r0 * s1 * t0
			+ w011 * r0 * s1 * t1
			+ w100 * r1 * s0 * t0
			+ w101 * r1 * s0 * t1
			+ w110 * r1 * s1 * t0
			+ w111 * r1 * s1 * t1;
		real vSum =
			w000 * r0 * s0 * t0 * vz000
			+ w001 * r0 * s0 * t1 * vz001
			+ w010 * r0 * s1 * t0 * vz010
			+ w011 * r0 * s1 * t1 * vz011
			+ w100 * r1 * s0 * t0 * vz100
			+ w101 * r1 * s0 * t1 * vz101
			+ w110 * r1 * s1 * t0 * vz110
			+ w111 * r1 * s1 * t1 * vz111;
		if (weightSum == 0) {
			//cout << "interpolation failed" << endl;
			vz = 0;
		}
		else {
			vz = vSum / weightSum;
		}
	}
}

void StaggeredGrid::extrapolateVelocity() {
	real *alignedVelocityX = new real[totalCells];
	real *alignedVelocityY = new real[totalCells];
	real *alignedVelocityZ = new real[totalCells];
	for (int i = 0; i < totalCells; i++) {
		alignedVelocityX[i] = -3.4567f;
	}
	for (int index = 0; index < totalCells; index++) {
		if (cellType[index] == SOLID || cellType[index] == FLUID) {
			// in this case the velocity is well defined.
			continue;
		}
		Vec3f *p = closestPoint[index];
		real valueFromX = p->v[0], valueFromY = p->v[1], valueFromZ = p->v[2];
		if (closestPointIndex[index] != -1) {
			real debugx = closestPoint[closestPointIndex[index]]->v[0];
			assert(debugx == valueFromX);
			assert(closestPointIndex[closestPointIndex[index]] == -1);
		}
		real vx, vy, vz;
		getAveragedVelocityAt(valueFromX, valueFromY, valueFromZ, vx, vy, vz);
		alignedVelocityX[index] = vx;
		alignedVelocityY[index] = vy;
		alignedVelocityZ[index] = vz;
		assert(abs(vx) < 10000);
	}
	for (int index = 0; index < totalCells; index++) {
		if (cellType[index] == SOLID || cellType[index] == FLUID) {
			continue;
		}
		if (cellType[index + 1] == AIR) {
			int i, j, k;
			getPos(index, resX, resY, resZ, i, j, k);
			velocityX[getIndex(i + 1, j, k, resX + 1, resY, resZ)] = (alignedVelocityX[index] + alignedVelocityX[index + 1]) * 0.5;
		}
		if (cellType[index + resX] == AIR) {
			int i, j, k;
			getPos(index, resX, resY, resZ, i, j, k);
			velocityY[getIndex(i, j + 1, k, resX, resY + 1, resZ)] = (alignedVelocityY[index] + alignedVelocityY[index + resX]) * 0.5;
		}
		if (cellType[index + slabSize] == AIR) {
			int i, j, k;
			getPos(index, resX, resY, resZ, i, j, k);
			velocityZ[getIndex(i, j, k + 1, resX, resY, resZ + 1)] = (alignedVelocityZ[index] + alignedVelocityZ[index + slabSize]) * 0.5;
		}
	}
}

StaggeredGrid::StaggeredGrid(bool test) {
	resX = resY = 100; resZ = 5;
	totalCells = resX * resY * resZ;
	slabSize = resX * resY;
	cellType = (CELL_TYPE*)calloc(totalCells, sizeof(CELL_TYPE));
	levelSetInfo = (LEVEL_SET_INFO*)calloc(totalCells, sizeof(LEVEL_SET_INFO));
	closestPoint = (Vec3f**)calloc(totalCells, sizeof(Vec3f*));
	totalVX = (resX + 1) * resY * resZ;
	totalVY = resX * (resY + 1) * resZ;
	totalVZ = resX * resY * (resZ + 1);
	velocityX = (real*)calloc(totalVX, sizeof(real));
	velocityY = (real*)calloc(totalVY, sizeof(real));
	velocityZ = (real*)calloc(totalVZ, sizeof(real));
	signedDistanceField = (real*)calloc(totalCells, sizeof(real));
	closestPointIndex = (int*)calloc(totalCells, sizeof(int));
	for (int i = 0; i < totalCells; i++) {
		closestPoint[i] = new Vec3f(-1, -2, -3); // value set for debugging only.
		closestPointIndex[i] = -2;
	}
}
