#pragma once

#include "meta.h"

class Field;
class Config;
class TimeStepController;

// Maybe non-isothermal solver should be derived from this class.
class SimpleSolver
{
protected:
	Field* rhoField;
	Field* vxField;
	Field* vyField;
	Field* vzField;

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
	real vdwTheta;
	real vdwInvWe;

	real velConvergeTol;
	real rhoConvergeTol;
	real rhoRelaxCoef;

	// Steps for solving SIMPLE.
	void computeVelocityStar(real dt, Field* vxGuessField, Field* vyGuessField,
		Field* vzGuessField, Field* rhoGuessField, Field* rhsRhoField,
		Field* vxStarField, Field* vyStarField, Field* vzStarField);
	void computeRhoPrime(real dt, Field* vxStarField, Field* vyStarField, Field* vzStarField,
		Field* rhoGuessField, Field* rhoPrimeField);
	void computeVelocityPrime(real dt,
		Field* rhoGuessField, Field* rhsRhoStarField, Field* rhsRhoStarStarField, 
		Field* vxPrimeField, Field* vyPrimeField, Field* vzPrimeField);

	// Other util functions.
	void advectVelocitySemiLagrange(real dt,
		Field * vxBackgroundField, Field * vyBackgroundField, Field * vzBackgroundField,
		Field * vxInterimField, Field *vyInterimField, Field *vzInterimField,
		Field * vxNewField, Field * vyNewField, Field * vzNewField);

	void fillVelocityFieldLoopBoundary(Field* xF, Field* yF, Field* zF);

	real updateRhoField(Field* rhoGuess, Field* rhoPrime);
	real updateVelocityField(Field* vxStarField, Field* vyStarField, Field* vzStarField,
		Field* vxPrimeField, Field* vyPrimeField, Field* vzPrimeField, 
		Field* vxGuessField, Field* vyGuessField, Field* vzGuessField);
	// Grid-dependent computation
	void laplaceRhoOnAlignedGrid(Field* rho, Field* lapRho);
	// void wdRhoOnAlignedGrid(Field* rho, Field* wdRho);
	// For fast and memory-efficient isothermal manipulation
	void rhsRhoOnAlignedGrid(Field* rho, Field* rhsRho);
	inline real isothermalWd(real r) const {
		return -2 * vdwA * r + vdwTheta * log(r / (vdwB - r)) +
			vdwTheta * vdwB / (vdwB - r);
	}
	inline int ifloor(real x) {
		return (x > 0) ? (int)x : (int)(x - 1);
	}
	inline int loopIndex(int pos, int dim) {
		int modular = pos % dim;
		return (modular >= 0) ? modular : modular + dim;
	}
public:
	void run(TimeStepController* step);
	SimpleSolver(Config* cfg, Field* initRhoField, Field* initVxField,
		Field* initVyField, Field* initVzField);
	void stepSimple(real dt);
	virtual ~SimpleSolver();
};

