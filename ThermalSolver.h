#pragma once

#include "meta.h"

class Field;
class Config;
class TimeStepController;

// Maybe non-isothermal solver should be derived from this class.
class ThermalSolver
{
protected:
	Field* rhoField;
	Field* vxField;
	Field* vyField;
	Field* vzField;
	Field* thetaField;

	Config* config;
	// Staggered Grid Size
	// For easier access and better readability
	// resX, Y, Z should be equal to size of rhoField
	int resX;
	int resY;
	int resZ;

	real h;
	real vdwA;
	real vdwB;
	real vdwInvWe;
	real vdwCv;

	real velConvergeTol;
	real rhoConvergeTol;
	real rhoRelaxCoef;

	real envGravity;
	real heatDiffuseSpeed;
	real targetTheta;
	real heatSpeed;

	// Positions of heaters.
	int heaterCount;
	int* heatersX;
	int* heatersY;
	int* heatersZ;

	void initializeHeaters();

	// Steps for solving SIMPLE.
	void computeVelocityStar(real dt, Field* vxGuessField, Field* vyGuessField,
		Field* vzGuessField, Field* rhoGuessField, Field* rhsRhoField,
		Field* vxStarField, Field* vyStarField, Field* vzStarField);
	void computeRhoPrime(real dt, Field* vxStarField, Field* vyStarField, Field* vzStarField,
		Field* rhoGuessField, Field* rhoPrimeField);
	void computeVelocityPrime(real dt,
		Field* rhoGuessField, Field* rhsRhoStarField, Field* rhsRhoStarStarField, 
		Field* vxPrimeField, Field* vyPrimeField, Field* vzPrimeField);
	void updateThetaField(real dt, Field* vxBackgroundField, Field* vyBackgroundField,
		Field* vzBackgroundField);

	// Other util functions.
	void advectVelocitySemiLagrange(real dt,
		Field * vxBackgroundField, Field * vyBackgroundField, Field * vzBackgroundField,
		Field * vxInterimField, Field *vyInterimField, Field *vzInterimField,
		Field * vxNewField, Field * vyNewField, Field * vzNewField);
	void advectFieldSemiLagrange(real dt, Field *vxField, Field *vyField, Field *vzField,
		Field *oldField, Field *newField);
	void fillVelocityFieldBorderZero(Field* xF, Field* yF, Field* zF);

	real updateRhoField(Field* rhoGuess, Field* rhoPrime);
	real updateVelocityField(Field* vxStarField, Field* vyStarField, Field* vzStarField,
		Field* vxPrimeField, Field* vyPrimeField, Field* vzPrimeField, 
		Field* vxGuessField, Field* vyGuessField, Field* vzGuessField);
	// Grid-dependent computation
	void laplacianFieldOnAlignedGrid(Field* f, Field* lapF);
	// void wdRhoOnAlignedGrid(Field* rho, Field* wdRho);
	// For fast and memory-efficient isothermal manipulation
	void rhsRhoOnAlignedGrid(Field* rho, Field* rhsRho);
	inline real thermalWd(real r, real t) const {
		static const real R = 1.0f;
		return -2 * vdwA * r + R * t * log(r / (vdwB - r)) +
			R * t * vdwB / (vdwB - r) - vdwCv * t * (log(t) - 1);
	}
	inline int ifloor(real x) {
		return (x > 0) ? (int)x : (int)(x - 1);
	}
public:
	void run(TimeStepController* step);
	ThermalSolver(Config* cfg, Field* initRhoField, Field* initVxField,
		Field* initVyField, Field* initVzField, Field* initThetaField);
	void stepSimple(real dt);
	virtual ~ThermalSolver();
};

