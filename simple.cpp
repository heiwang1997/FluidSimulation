#include "stdafx.h"
#include "StaggeredGrid.h"
#include "FieldManipulation.h"
#include "config.h"

using std::cout;
using std::endl;

void StaggeredGrid::setInitGuess(real * out_vxGuess, real * out_vyGuess, real * out_vzGuess, real * out_rhoGuess) {
	memcpy(out_vxGuess, velocityX, totalVX * sizeof(real));
	memcpy(out_vyGuess, velocityY, totalVY * sizeof(real));
	memcpy(out_vzGuess, velocityZ, totalVZ * sizeof(real));
	memcpy(out_rhoGuess, rho, totalCells * sizeof(real));
}


void StaggeredGrid::computeVelocityStar(real dt, real * vxGuess, real * vyGuess, real * vzGuess, real * rhoGuess, real * out_vxStar, real * out_vyStar, real * out_vzStar) {
	if (debugOutput) {
		cout << "computeVelocityStar()" << endl;
		// semilagrange u1 = advect u_n in u_guess, dt
		checkFieldStatus();
		cout << " input:" << endl;
		cout << "  vxguess: " << fieldMax(vxGuess, totalVX) << endl;
		cout << "  out_vxStar: " << fieldMax(out_vxStar, totalVX) << endl;
		cout << "  v " << fieldMax(velocityX, totalVX) << endl;
		cout << "  rhoGuess: " << fieldMax(rhoGuess, totalCells, true) << endl;
		cout << " =====" << endl;
	}
	advectVelocitySemiLagrange(dt, vxGuess, vyGuess, vzGuess, velocityX, velocityY, velocityZ, out_vxStar, out_vyStar, out_vzStar);
	if (debugOutput) cout << "  vxStar after advection: " << fieldMax(out_vxStar, totalVX) << endl;
	// eular step u* = u1 + f(rho_guess) * dt
	real *laplaceRho = new real[totalCells];
	if (debugOutput) cout << "  rho: " << fieldMax(rho, totalCells, true) << endl;
	laplaceRhoOnAlignedGrid(rho, laplaceRho);
	if (debugOutput) cout << "  laplaceRho: " << fieldMax(laplaceRho, totalCells) << endl;
	// W'(rho)
	real *WdRho = new real[totalCells];
	for (int i = 0; i < totalCells; i++) {
		WdRho[i] = funcWd(rho[i]);

	}
	// res+1 is copied from 0
	if (debugOutput) cout << "  WdRho: " << fieldMax(WdRho, totalCells) << endl;

	// cout << "dx, dt " << dx << " " << dt << endl;
	// !!! laplace rho can be very large due to the discontinuity of rho. (1)
	// cout << "laplace rho " << fieldMax(laplaceRho, totalCells) << endl;
	// cout << "func(rho)" << fieldMax(WdRho, totalCells) << endl;
	// update vx
	// !!! due to (1), v* can be very large. (2)
	for (int k = 0; k < resZ; k++) {
		for (int j = 0; j < resY; j++) {
			for (int i = 0; i < resX; i++) {
				int vIndex = getIndex(i, j, k, resX + 1, resY, resZ);
				int centerRight = getIndex(i, j, k, resX, resY, resZ);
				int centerLeft = getIndex((i - 1 + resX) % resX, j, k, resX, resY, resZ);
				out_vxStar[vIndex] += -(WdRho[centerRight] - WdRho[centerLeft]) / dx * dt;
				out_vxStar[vIndex] += V_INVWE * (laplaceRho[centerRight] - laplaceRho[centerLeft]) / dx * dt;
			}
			out_vxStar[getIndex(0, j, k, resX + 1, resY, resZ)] = out_vxStar[getIndex(resX, j, k, resX + 1, resY, resZ)];
		}
	}
	// update vy
	for (int i 
		= 0; i < resX; i++) {
		for (int k = 0; k < resZ; k++) {
			for (int j = 0; j < resY; j++) {
				int vIndex = getIndex(i, j, k, resX, resY + 1, resZ);
				int centerAbove = getIndex(i, j, k, resX, resY, resZ);
				int centerBelow = getIndex(i, (j - 1 + resY) % resY, k, resX, resY, resZ);
				out_vyStar[vIndex] += -(WdRho[centerAbove] - WdRho[centerBelow]) / dx * dt;
				out_vyStar[vIndex] += V_INVWE * (laplaceRho[centerAbove] - laplaceRho[centerBelow]) / dx * dt;
			}
			out_vyStar[getIndex(i, 0, k, resX, resY + 1, resZ)] = out_vyStar[getIndex(i, resY, k, resX, resY + 1, resZ)];
		}
	}
	// update vz
	for (int i = 0; i < resX; i++) {
		for (int j = 0; j < resY; j++) {
			for (int k = 0; k < resZ; k++) {
				int vIndex = getIndex(i, j, k, resX, resY, resZ + 1);
				int centerBack = getIndex(i, j, k, resX, resY, resZ);
				int centerFront = getIndex(i, j, (k - 1 + resZ) % resZ, resX, resY, resZ);
				out_vzStar[vIndex] += -(WdRho[centerBack] - WdRho[centerFront]) / dx * dt;
				out_vzStar[vIndex] += V_INVWE * (laplaceRho[centerBack] - laplaceRho[centerFront]) / dx * dt;
			}
			out_vzStar[getIndex(i, j, 0, resX, resY, resZ + 1)] = out_vzStar[getIndex(i, j, resZ, resX, resY, resZ + 1)];
		}
	}
	if (debugOutput) cout << "  out_vxStar after euler: " << fieldMax(out_vxStar, totalVX) << endl;
}

void StaggeredGrid::computeRhoPrime(real dt, real * vxStar, real * vyStar, real * vzStar, real * rhoStar, real * out_rhoPrime) {
	if (debugOutput) {
		cout << "computeRhoPrime()" << endl;
		checkFieldStatus();
		cout << " input:" << endl;
		cout << "  rhoStar: " << fieldMax(rhoStar, totalCells, true) << endl;
		cout << " ======" << endl;
	}
	// here we first calculate rho' + rho*
	// (rho' - rho) / dt + u . grad(rho) + rho * div(u) = 0

	// !!! due to (2), v* can have very large value. (3)
	if (debugOutput) cout << "  rho (that is advected): " << fieldMax(rho, totalCells, true) << endl;
	advectFieldSemiLagrange(dt, vxStar, vyStar, vzStar, rho, out_rhoPrime, true);
	if (debugOutput) cout << "  out_rhoPrime (after advection): " << fieldMax(out_rhoPrime, totalCells, true) << endl;
	for (int i = 0; i < resX; i++) {
		for (int j = 0; j < resY; j++) {
			for (int k = 0; k < resZ; k++) {
				int index = getIndex(i, j, k, resX, resY, resZ);
				real vRight = vxStar[getIndex(i + 1, j, k, resX + 1, resY, resZ)];
				real vLeft = vxStar[getIndex(i, j, k, resX + 1, resY, resZ)];
				real vTop = vyStar[getIndex(i, j + 1, k, resX, resY + 1, resZ)];
				real vBottom = vyStar[getIndex(i, j, k, resX, resY + 1, resZ)];
				real vBack = vzStar[getIndex(i, j, k + 1, resX, resY, resZ + 1)];
				real vFront = vzStar[getIndex(i, j, k, resX, resY, resZ + 1)];
				real divergence = (vRight - vLeft + vTop - vBottom + vBack - vFront) / dx;
				// eular
				out_rhoPrime[index] /= (1 + divergence * dt);
			}
		}
	}
	if (debugOutput) cout << "  out_rhoPrime (after eular): " << fieldMax(out_rhoPrime, totalCells, true) << endl;
	for (int i = 0; i < totalCells; i++) {
		out_rhoPrime[i] -= rhoStar[i];
	}
	if (debugOutput) cout << "  out_rhoPrime (after subtraction): " << fieldMax(out_rhoPrime, totalCells, true) << endl;
}

void StaggeredGrid::computeVelocityPrime(real dt, real * vxStar, real * vyStar, real * vzStar, real * rhoStar, real * rhoPrime, real * out_vxPrime, real * out_vyPrime, real * out_vzPrime) {
	return;
}

bool StaggeredGrid::updateGuesses(real dt, real * io_vxGuess, real * io_vyGuess, real * io_vzGuess, real * in_vxPrime, real * in_vyPrime, real * in_vzPrime, real * io_rhoGuess, real * in_rhoPrime) {
	if (debugOutput) {
		cout << "updateGuesses()" << endl;
		checkFieldStatus();
		cout << " input:" << endl;
		cout << "  i_vxGuess: " << fieldMax(io_vxGuess, totalVX) << endl;
		cout << "  i_rhoPrime: " << fieldMax(in_rhoPrime, totalCells, true) << endl;
		cout << "  i_rhoGuess: " << fieldMax(io_rhoGuess, totalCells, true) << endl;
		cout << " =======" << endl;
	}
	real lambda_u = 0.0, lambda_rho = 1.0;
	double squaredNormVPrime = 0.0, squaredNormRhoPrime = 0.0;
	for (int i = 0; i < totalVX; i++) {
		io_vxGuess[i] += lambda_u * in_vxPrime[i];
		squaredNormVPrime += in_vxPrime[i] * in_vxPrime[i];
	}
	for (int i = 0; i < totalVY; i++) {
		io_vyGuess[i] += lambda_u * in_vyPrime[i];
		squaredNormVPrime += in_vyPrime[i] * in_vyPrime[i];
	}
	for (int i = 0; i < totalVZ; i++) {
		io_vzGuess[i] += lambda_u * in_vzPrime[i];
		squaredNormVPrime += in_vzPrime[i] * in_vzPrime[i];
	}
	for (int i = 0; i < totalCells; i++) {
		io_rhoGuess[i] += lambda_rho * in_rhoPrime[i];
		squaredNormRhoPrime += in_rhoPrime[i] * in_rhoPrime[i];
	}
	real tol = 1e-8;
	if (squaredNormRhoPrime < tol && squaredNormVPrime < tol) { // norm(rho')^2 < tol
		return true;
	}
	else {
		return false;
	}
}

void StaggeredGrid::stepSIMPLE(real dt) {
	real *vxGuess = new real[totalVX];
	real *vyGuess = new real[totalVY];
	real *vzGuess = new real[totalVZ];
	real *rhoGuess = new real[totalCells];
	if (debugOutput) cout << "set initial guess" << endl;
	checkFieldStatus(true);
	setInitGuess(vxGuess, vyGuess, vzGuess, rhoGuess);
	int loopcnt = 0;
	while (true) {
		loopcnt++;
		if (debugOutput) cout << loopcnt << "-th iteration" << endl;
		if (loopcnt > 100) {
			cout << "not converging" << endl;
			exit(0);
		}

		// 1. assemble and solve momentum equation for v*
		// cout << "compute v*" << endl;
		real *vxStar = new real[totalVX];
		real *vyStar = new real[totalVY];
		real *vzStar = new real[totalVZ];
		computeVelocityStar(dt, vxGuess, vyGuess, vzGuess, rhoGuess, vxStar, vyStar, vzStar);
		//cout << "vx' " << fieldMax(vxStar, totalVX) << endl;


		// 2. compute rho* using the equation of state ...
		// skipped step
		real *rhoStar = new real[totalCells];
		memcpy(rhoStar, rhoGuess, totalCells * sizeof(real));

		// 3. assemble and solve *continuity equation for rho'
		// cout << "compute rho'" << endl;
		real *rhoPrime = new real[totalCells];
		computeRhoPrime(dt, vxStar, vyStar, vzStar, rhoStar, rhoPrime);

		// 4. skipped step
		// cout << "compute v'" << endl;
		real *vxPrime = new real[totalVX];
		real *vyPrime = new real[totalVY];
		real *vzPrime = new real[totalVZ];
		computeVelocityPrime(dt, vxStar, vyStar, vzStar, rhoStar, rhoPrime, vxPrime, vyPrime, vzPrime);

		memcpy(vxGuess, vxStar, totalVX * sizeof(real));
		memcpy(vyGuess, vyStar, totalVY * sizeof(real));
		memcpy(vzGuess, vzStar, totalVZ * sizeof(real));

		// cout << "update guess" << endl;
		bool converged = updateGuesses(dt, vxGuess, vyGuess, vzGuess, vxPrime, vyPrime, vzPrime, rhoGuess, rhoPrime);
		delete[] rhoStar;
		delete[] vxStar;
		delete[] vyStar;
		delete[] vzStar;
		delete[] rhoPrime;
		delete[] vxPrime;
		delete[] vyPrime;
		delete[] vzPrime;
		if (converged) {
			cout << "converged in " << loopcnt << " iterations" << endl;
			break;
		}
		//cout << fieldMax(velocityX, totalVX) << "---" << fieldMax(vxGuess, totalVX) << endl;
		if (debugOutput) {
			checkFieldStatus();
			cout << "rhoGuess: " << fieldMax(rhoGuess, totalCells, true) << endl;
			cout << "end of iteration\n" << endl;
		}
	}
	//cout << fieldMax(velocityX, totalVX) << " " << fieldMax(vxGuess, totalVX) << endl;
	memcpy(velocityX, vxGuess, totalVX * sizeof(real));
	memcpy(velocityY, vyGuess, totalVY * sizeof(real));
	memcpy(velocityZ, vzGuess, totalVZ * sizeof(real));
	memcpy(rho, rhoGuess, totalCells * sizeof(real));
	delete[] vxGuess;
	delete[] vyGuess;
	delete[] vzGuess;
	delete[] rhoGuess;
	cout << endl;
}

void StaggeredGrid::addBubble() {
	float xTotal = dx*resX;
	float yTotal = dx*resY;
	float zTotal = dx*resZ;

	Vec3f bubble1 = Vec3f(0.41, 0.50, 0.50) * xTotal;// *_xRes;
	Vec3f bubble2 = Vec3f(0.67, 0.50, 0.50) * yTotal;// *_yRes;
	float Rb1 = 0.16 * xTotal;// *_xRes;
	float Rb2 = 0.08 * yTotal;// *_xRes;

	double totalmass = 0;

	for (int z = 0; z < resZ; z++)
		for (int y = 0; y < resY; y++)
			for (int x = 0; x < resX; x++)
			{
				int index = x + y * resX + z * slabSize;
				rho[index] = 500.0f; //liquid dens
				
				Vec3f gc = Vec3f(x + 0.5, y + 0.5, z + 0.5) * dx;

				//bubble1
				Vec3f dis = gc - bubble1;

				if (mag2(dis)<Rb1*Rb1) {
					rho[index] = 0.5f;  //vapor dens
				}

				//bubble2
				dis = gc - bubble2;

				if (mag2(dis)<Rb2*Rb2) {
					rho[index] = 0.5f;  //vapor dens
				}
				totalmass += rho[index];
			}
	printf("initial total mass: %f\n", totalmass);
	return;
}

void StaggeredGrid::runSIMPLE() {
	timeStep = new TimeStepController(config->totalFrame(), config->gridFPS(), config->gridDt());

	addBubble();

	while (!timeStep->isFinished()) {
		stepSIMPLE(timeStep->getStepDt());
		int fmCnt;
		if (timeStep->isFrameTime(fmCnt)) {
#ifndef _DEBUG
			dumpSlicePreview(fmCnt, rho, 50, 0.2);
#endif // !_DEBUG
		}
		std::cout << "This step elapsed *** seconds." << std::endl;
	}
	delete timeStep;
	std::cout << "Total elapsed *** seconds." << std::endl;
}

void StaggeredGrid::laplaceRhoOnAlignedGrid(real * rho, real * out_laplacianRho) {
	// var for debugging purpose, remove when optimize.
	bool allZero = true;

	for (int i = 0; i < resX; i++) {
		for (int j = 0; j < resY; j++) {
			for (int k = 0; k < resZ; k++) {
				int center	= getIndex(i, j, k,						resX, resY, resZ);
				int left	= getIndex((i - 1 + resX) % resX, j, k, resX, resY, resZ);
				int right	= getIndex((i + 1) % resX, j, k,		resX, resY, resZ);
				int up		= getIndex(i, (j + 1) % resY, k,		resX, resY, resZ);
				int bottom	= getIndex(i, (j - 1 + resY) % resY, k, resX, resY, resZ);
				int front	= getIndex(i, j, (k - 1 + resZ) % resZ, resX, resY, resZ);
				int back	= getIndex(i, j, (k + 1) % resZ,		resX, resY, resZ);
				out_laplacianRho[center] = (rho[left] + rho[right] + rho[up] + rho[bottom] + rho[front] + rho[back] - 6 * rho[center]) / dx / dx;
				if (out_laplacianRho[center] != 0) {
					//cout << out_laplacianRho[center] << endl;
					allZero = false;
				}
			}
		}
	}
	if (allZero) {
		cout << "laplacian all zero" << endl;
		exit(0);
	}
}


real StaggeredGrid::funcWd(real r) {
	// if for debugging purpose, remove when optimize.
	if (r - V_PB == 0) {
		cout << "illegal value for funcWd, error 1" << endl;
		exit(0);
	}
	if (r == 0) {
		cout << "illegal value for funcWd, error 2" << endl;
		exit(0);
	}
	real ans = -2 * V_PA * r + V_RTM * log(r / (V_PB - r)) + V_RTM * V_PB / (r - V_PB);
	if (ans != ans) {
		cout << "funcWd(): " << r << endl;
	}
	return ans;
}

// debug function, does not modify the arrays, remove usage when optimize
void StaggeredGrid::checkFieldStatus(bool summary) {
	real vxm = fieldMax(velocityX, totalVX);
	real vym = fieldMax(velocityY, totalVY);
	real vzm = fieldMax(velocityZ, totalVZ);
	cout << "velocity values checked; ";
	real rhom = fieldMax(rho, totalCells, true);
	cout << "rho checked." << endl;
	if (summary) {
		cout << "summary: " << vxm << ", " << vym << ", " << vzm << "; " << rhom << endl;
	}
}
