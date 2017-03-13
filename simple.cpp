#include "stdafx.h"
#include <Eigen/Sparse>
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


/* This function compute V* by solving the Momentum Equation.
 * We ignore the shear tensor and gravity.
 * The PDE is solved in a splitting manner.
 * Du/Dt = - grad W'(\rho) + 1 / We * grad laplacian rho
 */
void StaggeredGrid::computeVelocityStar(real dt, real * vxGuess, real * vyGuess, real * vzGuess, 
	real * rhoGuess, real * out_vxStar, real * out_vyStar, real * out_vzStar) {
	// Optimize Suggestion: 1. high-speed cache usage: loop seq.; 2. allocate once.
	if (debugOutput) {
		cout << "++ computeVelocityStar()" << endl;
		// semilagrange u1 = advect u_n in u_guess, dt
		// checkFieldStatus();
		cout << " input:" << endl;
		cout << "  vxguess: " << fieldMax(vxGuess, totalVX) << endl;
		cout << "  v " << fieldMax(velocityX, totalVX) << endl;
		cout << "  rhoGuess: " << fieldMax(rhoGuess, totalCells, true) << endl;
		cout << " =====" << endl;
	}
	// First Step:
	// Solve: Du/Dt = 0
	advectVelocitySemiLagrange(dt, vxGuess, vyGuess, vzGuess, 
		velocityX, velocityY, velocityZ,
		out_vxStar, out_vyStar, out_vzStar);

	//cout << "----------- After Advection.-------------------" << endl;
	//cout << "  out_vxStar: " << fieldMax(out_vxStar, totalVX) << endl;

	// Second Step:
	// Euler forward u* = u1 + f(rho_guess) * dt
	real *laplaceRho = new real[totalCells];
	laplaceRhoOnAlignedGrid(rhoGuess, laplaceRho);

	//real *WdRho = new real[totalCells];
	//for (int i = 0; i < totalCells; i++) {
	//	WdRho[i] = funcWd(rhoGuess[i]);
	//}

	//if (debugOutput) cout << "  WdRho: " << fieldMax(WdRho, totalCells) << endl;

	// !!! laplace rho can be very large due to the discontinuity of rho. (1)
	// !!! due to (1), v* can be very large. (2)

	// We employ a loop boundary, so res is compied from 0
	// update vx
	for (int k = 0; k < resZ; k++) {
		for (int j = 0; j < resY; j++) {
			for (int i = 0; i < resX; i++) {
				int vIndex = getIndex(i, j, k, resX + 1, resY, resZ);
				int centerRight = getIndex(i, j, k, resX, resY, resZ);
				int centerLeft = getIndex((i - 1 + resX) % resX, j, k, resX, resY, resZ);
				// out_vxStar[vIndex] += -(WdRho[centerRight] - WdRho[centerLeft]) / dx * dt;
				out_vxStar[vIndex] += -(rhoGuess[centerRight] - rhoGuess[centerLeft]) * 
					(funcWd(0.5f * (rhoGuess[centerLeft] + rhoGuess[centerRight]))) / dx * dt;
				out_vxStar[vIndex] += V_INVWE * (laplaceRho[centerRight] - laplaceRho[centerLeft]) / dx * dt;
			}
			out_vxStar[getIndex(resX, j, k, resX + 1, resY, resZ)] = out_vxStar[getIndex(0, j, k, resX + 1, resY, resZ)];
		}
	}
	// update vy
	for (int i = 0; i < resX; i++) {
		for (int k = 0; k < resZ; k++) {
			for (int j = 0; j < resY; j++) {
				int vIndex = getIndex(i, j, k, resX, resY + 1, resZ);
				int centerAbove = getIndex(i, j, k, resX, resY, resZ);
				int centerBelow = getIndex(i, (j - 1 + resY) % resY, k, resX, resY, resZ);
				// out_vyStar[vIndex] += -(WdRho[centerAbove] - WdRho[centerBelow]) / dx * dt;
				out_vyStar[vIndex] += -(rhoGuess[centerAbove] - rhoGuess[centerBelow]) *
					(funcWd(0.5f * (rhoGuess[centerAbove] + rhoGuess[centerBelow]))) / dx * dt;
				out_vyStar[vIndex] += V_INVWE * (laplaceRho[centerAbove] - laplaceRho[centerBelow]) / dx * dt;
			}
			out_vyStar[getIndex(i, resY, k, resX, resY + 1, resZ)] = out_vyStar[getIndex(i, 0, k, resX, resY + 1, resZ)];
		}
	}
	// update vz
	for (int i = 0; i < resX; i++) {
		for (int j = 0; j < resY; j++) {
			for (int k = 0; k < resZ; k++) {
				int vIndex = getIndex(i, j, k, resX, resY, resZ + 1);
				int centerBack = getIndex(i, j, k, resX, resY, resZ);
				int centerFront = getIndex(i, j, (k - 1 + resZ) % resZ, resX, resY, resZ);
				//out_vzStar[vIndex] += -(WdRho[centerBack] - WdRho[centerFront]) / dx * dt;
				out_vzStar[vIndex] += -(rhoGuess[centerBack] - rhoGuess[centerFront]) *
					(funcWd(0.5f * (rhoGuess[centerBack] + rhoGuess[centerFront]))) / dx * dt;
				out_vzStar[vIndex] += V_INVWE * (laplaceRho[centerBack] - laplaceRho[centerFront]) / dx * dt;
			}
			out_vzStar[getIndex(i, j, resZ, resX, resY, resZ + 1)] = out_vzStar[getIndex(i, j, 0, resX, resY, resZ + 1)];
		}
	}
	if (debugOutput) cout << "  out_vxStar after euler: " << fieldMax(out_vxStar, totalVX) << endl;

	delete[] laplaceRho;
	// delete[] WdRho;
}

void printSlicePreview(real* arr, int resX, int resY, int resZ, int z) {
	cout << "Slice Preview of slab " << z << endl;
	for (int i = 0; i < resX; ++i) {
		for (int j = 0; j < resY; ++j) {
			cout << arr[getIndex(i, j, z, resX, resY, resZ)] << '\t';
		}
		cout << endl << endl << endl;
	}
	cout << "End of dump" << endl;
}

// Rho Prime is calculated using C.E.
// See rhoprime.docx for more info.
void StaggeredGrid::computeRhoPrime(real dt, real * vxStar, real * vyStar, real * vzStar, 
	real * rhoStar, real * out_rhoPrime) {
	typedef Eigen::SparseMatrix<real> SpMat;
	typedef Eigen::Triplet<real> Triplet;

	if (debugOutput) {
		cout << "++ computeRhoPrime()" << endl;
		// semilagrange u1 = advect u_n in u_guess, dt
		// checkFieldStatus();
		cout << " input:" << endl;
		cout << "  vxStar: " << fieldMax(vxStar, totalVX) << endl;
		cout << "  rhoStar " << fieldMax(rhoStar, totalCells, true) << endl;
		cout << " =====" << endl;
	}

	std::vector<Triplet> coefficients;
	coefficients.reserve(totalCells * 7);

	Eigen::VectorXf b(totalCells);

	for (int k = 0; k < resZ; ++k) {
		for (int j = 0; j < resY; ++j) {
			for (int i = 0; i < resX; ++i) {
				int center = getIndex(i, j, k, resX, resY, resZ);
				int left = getIndex((i - 1 + resX) % resX, j, k, resX, resY, resZ);
				int right = getIndex((i + 1) % resX, j, k, resX, resY, resZ);
				int up = getIndex(i, (j + 1) % resY, k, resX, resY, resZ);
				int bottom = getIndex(i, (j - 1 + resY) % resY, k, resX, resY, resZ);
				int front = getIndex(i, j, (k - 1 + resZ) % resZ, resX, resY, resZ);
				int back = getIndex(i, j, (k + 1) % resZ, resX, resY, resZ);
				real vRight = vxStar[getIndex(i + 1, j, k, resX + 1, resY, resZ)];
				real vLeft = vxStar[getIndex(i, j, k, resX + 1, resY, resZ)];
				real vTop = vyStar[getIndex(i, j + 1, k, resX, resY + 1, resZ)];
				real vBottom = vyStar[getIndex(i, j, k, resX, resY + 1, resZ)];
				real vBack = vzStar[getIndex(i, j, k + 1, resX, resY, resZ + 1)];
				real vFront = vzStar[getIndex(i, j, k, resX, resY, resZ + 1)];

				coefficients.push_back(Triplet(center, center, 2 * dx / dt + vRight - vLeft + vTop - vBottom + vBack - vFront));
				coefficients.push_back(Triplet(center, right, vRight));
				coefficients.push_back(Triplet(center, left, -vLeft));
				coefficients.push_back(Triplet(center, up, vTop));
				coefficients.push_back(Triplet(center, bottom, -vBottom));
				coefficients.push_back(Triplet(center, back, vBack));
				coefficients.push_back(Triplet(center, front, -vFront));

				real rhs = -((rhoStar[center] + rhoStar[right]) * vRight
					- (rhoStar[center] + rhoStar[left]) * vLeft
					+ (rhoStar[center] + rhoStar[up]) * vTop
					- (rhoStar[center] + rhoStar[bottom]) * vBottom
					+ (rhoStar[center] + rhoStar[back]) * vBack
					- (rhoStar[center] + rhoStar[front]) * vFront
					+ 2 * dx / dt * (rhoStar[center] - rho[center]));
				// Should be equivalent to (b << rhs;)
				b[center] = rhs;
			}
		}
	}

	SpMat A(totalCells, totalCells);
	A.setFromTriplets(coefficients.begin(), coefficients.end());

	// Eigen::ConjugateGradient;
	// Eigen::LeastSquaresConjugateGradient;
	Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<real> > solver;
	Eigen::VectorXf x;
	solver.compute(A);
	x = solver.solve(b);
	if (solver.info() != Eigen::Success) {
		if (debugOutput) cout << "Solver Failed." << endl;
		exit(0);
	}

	memcpy(out_rhoPrime, x.data(), sizeof(real) * totalCells);

	if (debugOutput) {
		cout << " output:" << endl;
		//cout << "  vxStar: " << fieldMax(vxStar, totalVX) << endl;
		cout << "  rhoPrime " << fieldMax(out_rhoPrime, totalCells, true) << endl;
		printSlicePreview(out_rhoPrime, resX, resY, resZ, resZ / 2);
	}
}

void checkFuncWrong(real* arr, int pos1, int pos2) {
	if (arr[pos1] + arr[pos2] == 0) {
		cout << "Check point matched!" << endl;
		cout << arr[pos1] << ' ' << arr[pos2] << endl;
		cout << pos1 << ' ' << pos2 << endl;
		cout << endl;
	}
}
 
/* This function compute V' using rho* and rho'
 * The intuition behind this is that M.E. needs to be enhanced.
 * See vprime.docx for more info.
 */
void StaggeredGrid::computeVelocityPrime(real dt, real * vxStar, real * vyStar, real * vzStar,
	real * rhoStar, real * rhoPrime, real * out_vxPrime, real * out_vyPrime, real * out_vzPrime) {
	
	real* rhoStarStar = new real[totalCells];
	for (int i = 0; i < totalCells; ++i) {
		rhoStarStar[i] = rhoStar[i] + rhoPrime[i];
		if (rhoStar[i] == 0) {
			cout << "rhoStar" << endl;
		}
		if (rhoStarStar[i] == 0) {
			cout << "rhoStarStar" << endl;
		}
	}
	real *laplaceRhoStarStar = new real[totalCells];
	laplaceRhoOnAlignedGrid(rhoStarStar, laplaceRhoStarStar);
	real *laplaceRhoStar = new real[totalCells];
	laplaceRhoOnAlignedGrid(rhoStar, laplaceRhoStar);

	// We employ a loop boundary, so res is compied from 0
	// update vx
	for (int k = 0; k < resZ; k++) {
		for (int j = 0; j < resY; j++) {
			for (int i = 0; i < resX; i++) {
				int vIndex = getIndex(i, j, k, resX + 1, resY, resZ);
				int centerRight = getIndex(i, j, k, resX, resY, resZ);
				int centerLeft = getIndex((i - 1 + resX) % resX, j, k, resX, resY, resZ);

				checkFuncWrong(rhoStarStar, centerLeft, centerRight);
				checkFuncWrong(rhoStar, centerLeft, centerRight);

				out_vxPrime[vIndex] += -(rhoStarStar[centerRight] - rhoStarStar[centerLeft]) *
					(funcWd(0.5f * (rhoStarStar[centerLeft] + rhoStarStar[centerRight]))) / dx * dt;
				out_vxPrime[vIndex] += V_INVWE * (laplaceRhoStarStar[centerRight] - laplaceRhoStarStar[centerLeft]) / dx * dt;

				out_vxPrime[vIndex] -= -(rhoStar[centerRight] - rhoStar[centerLeft]) *
					(funcWd(0.5f * (rhoStar[centerLeft] + rhoStar[centerRight]))) / dx * dt;
				out_vxPrime[vIndex] -= V_INVWE * (laplaceRhoStar[centerRight] - laplaceRhoStar[centerLeft]) / dx * dt;
			}
			out_vxPrime[getIndex(resX, j, k, resX + 1, resY, resZ)] = out_vxPrime[getIndex(0, j, k, resX + 1, resY, resZ)];
		}
	}
	// update vy
	for (int i = 0; i < resX; i++) {
		for (int k = 0; k < resZ; k++) {
			for (int j = 0; j < resY; j++) {
				int vIndex = getIndex(i, j, k, resX, resY + 1, resZ);
				int centerAbove = getIndex(i, j, k, resX, resY, resZ);
				int centerBelow = getIndex(i, (j - 1 + resY) % resY, k, resX, resY, resZ);

				out_vyPrime[vIndex] += -(rhoStarStar[centerAbove] - rhoStarStar[centerBelow]) *
					(funcWd(0.5f * (rhoStarStar[centerAbove] + rhoStarStar[centerBelow]))) / dx * dt;
				out_vyPrime[vIndex] += V_INVWE * (laplaceRhoStarStar[centerAbove] - laplaceRhoStarStar[centerBelow]) / dx * dt;

				out_vyPrime[vIndex] -= -(rhoStar[centerAbove] - rhoStar[centerBelow]) *
					(funcWd(0.5f * (rhoStar[centerAbove] + rhoStar[centerBelow]))) / dx * dt;
				out_vyPrime[vIndex] -= V_INVWE * (laplaceRhoStar[centerAbove] - laplaceRhoStar[centerBelow]) / dx * dt;
			}
			out_vyPrime[getIndex(i, resY, k, resX, resY + 1, resZ)] = out_vyPrime[getIndex(i, 0, k, resX, resY + 1, resZ)];
		}
	}
	// update vz
	for (int i = 0; i < resX; i++) {
		for (int j = 0; j < resY; j++) {
			for (int k = 0; k < resZ; k++) {
				int vIndex = getIndex(i, j, k, resX, resY, resZ + 1);
				int centerBack = getIndex(i, j, k, resX, resY, resZ);
				int centerFront = getIndex(i, j, (k - 1 + resZ) % resZ, resX, resY, resZ);

				out_vzPrime[vIndex] += -(rhoStarStar[centerBack] - rhoStarStar[centerFront]) *
					(funcWd(0.5f * (rhoStarStar[centerBack] + rhoStarStar[centerFront]))) / dx * dt;
				out_vzPrime[vIndex] += V_INVWE * (laplaceRhoStarStar[centerBack] - laplaceRhoStarStar[centerFront]) / dx * dt;

				out_vzPrime[vIndex] -= -(rhoStar[centerBack] - rhoStar[centerFront]) *
					(funcWd(0.5f * (rhoStar[centerBack] + rhoStar[centerFront]))) / dx * dt;
				out_vzPrime[vIndex] -= V_INVWE * (laplaceRhoStar[centerBack] - laplaceRhoStar[centerFront]) / dx * dt;
			}
			out_vzPrime[getIndex(i, j, resZ, resX, resY, resZ + 1)] = out_vzPrime[getIndex(i, j, 0, resX, resY, resZ + 1)];
		}
	}

	delete[] rhoStarStar;
	delete[] laplaceRhoStarStar;
	delete[] laplaceRhoStar;

}

// I wonder why writing in C is like writing verilog... inout wire[31:0] io_vxGuess
bool StaggeredGrid::updateGuesses(real * io_vxGuess, real * io_vyGuess, 
	real * io_vzGuess, real * in_vxStar, real * in_vyStar, real * in_vzStar,
	real * in_vxPrime, real * in_vyPrime, real * in_vzPrime, 
	real * io_rhoGuess, real * in_rhoPrime) {
	if (debugOutput) {
		cout << "updateGuesses()" << endl;
		checkFieldStatus();
		cout << " input:" << endl;
		cout << "  i_vxGuess: " << fieldMax(io_vxGuess, totalVX) << endl;
		cout << "  i_rhoPrime: " << fieldMax(in_rhoPrime, totalCells, true) << endl;
		cout << "  i_rhoGuess: " << fieldMax(io_rhoGuess, totalCells, true) << endl;
		cout << " =======" << endl;
	}

	static const real lambda_rho = 1.0;
	static const real eps = 1e-6;

	double squaredNormVPrime = 0.0, squaredNormRhoPrime = 0.0;
	for (int i = 0; i < totalVX; i++) {
		real finalVx = in_vxStar[i] + in_vxPrime[i];
		real vxDelta = finalVx - io_vxGuess[i];
		io_vxGuess[i] = finalVx;
		squaredNormVPrime += vxDelta * vxDelta;
	}
	for (int i = 0; i < totalVY; i++) {
		real finalVy = in_vyStar[i] + in_vyPrime[i];
		real vyDelta = finalVy - io_vyGuess[i];
		io_vyGuess[i] = finalVy;
		squaredNormVPrime += vyDelta * vyDelta;
	}
	for (int i = 0; i < totalVZ; i++) {
		real finalVz = in_vzStar[i] + in_vzPrime[i];
		real vzDelta = finalVz - io_vzGuess[i];
		io_vzGuess[i] = finalVz;
		squaredNormVPrime += vzDelta * vzDelta;
	}
	for (int i = 0; i < totalCells; i++) {
		io_rhoGuess[i] += lambda_rho * in_rhoPrime[i];
		squaredNormRhoPrime += in_rhoPrime[i] * in_rhoPrime[i];
	}
	if (squaredNormRhoPrime < eps && squaredNormVPrime < eps) {
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
	// checkFieldStatus(true);
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
		real *vxStar = new real[totalVX];
		real *vyStar = new real[totalVY];
		real *vzStar = new real[totalVZ];
		computeVelocityStar(dt, vxGuess, vyGuess, vzGuess, rhoGuess, vxStar, vyStar, vzStar);

		// 2. compute rho* using the equation of state - skipped.
		// real *rhoStar = new real[totalCells];
		// memcpy(rhoStar, rhoGuess, totalCells * sizeof(real));
		real* rhoStar = rhoGuess;

		// 3. assemble and solve *continuity equation for rho'
		real *rhoPrime = new real[totalCells];
		computeRhoPrime(dt, vxStar, vyStar, vzStar, rhoStar, rhoPrime);

		// 4. Compute v prime.
		// cout << "compute v'" << endl;
		real *vxPrime = new real[totalVX];
		real *vyPrime = new real[totalVY];
		real *vzPrime = new real[totalVZ];
		computeVelocityPrime(dt, vxStar, vyStar, vzStar, rhoStar, rhoPrime, vxPrime, vyPrime, vzPrime);

		// cout << "update guess" << endl;
		bool converged = updateGuesses(vxGuess, vyGuess, vzGuess,
			vxStar, vyStar, vzStar, vxPrime, vyPrime, vzPrime, rhoGuess, rhoPrime);
		// delete[] rhoStar;
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
	float xTotal = dx * resX;
	float yTotal = dx * resY;
	float zTotal = dx * resZ;

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

				// Grid center position
				Vec3f gc = Vec3f(x + 0.5, y + 0.5, z + 0.5) * dx;

				//bubble1
				Vec3f dis = gc - bubble1;
				if (mag2(dis) < Rb1 * Rb1) {
					rho[index] = 0.5f;  //vapor dens
				}

				//bubble2
				dis = gc - bubble2;
				if (mag2(dis) < Rb2 * Rb2) {
					rho[index] = 0.5f;  //vapor dens
				}
				totalmass += rho[index];
			}
	printf("initial total mass: %f\n", totalmass);
	return;
}

void StaggeredGrid::runSIMPLE() {
	timeStep = new TimeStepController(config->totalFrame(), config->gridFPS(), config->gridDt());
	
	bool oldLoopBoundary = loopBoundary;
	loopBoundary = true;
	
	addBubble();

	while (!timeStep->isFinished()) {
		stepSIMPLE(timeStep->getStepDt());
		int fmCnt;
		if (timeStep->isFrameTime(fmCnt)) {
#ifndef _DEBUG
			dumpSlicePreview(fmCnt, rho, 50, 0.2f);
#endif // !_DEBUG
		}
		std::cout << "This step elapsed *** seconds." << std::endl;
	}
	delete timeStep;
	
	loopBoundary = oldLoopBoundary;
	std::cout << "Total elapsed *** seconds." << std::endl;
}

void StaggeredGrid::laplaceRhoOnAlignedGrid(real * rho, real * out_laplacianRho) {
	// var for debugging purpose, remove when optimize.
	bool allZero = true;
	// Question: How to discretize a Laplacian?
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
	static const real theta = 1.0f;

	//real ans = -2 * V_PA * r + V_RTM * log(r / (V_PB - r)) + V_RTM * V_PB / (r - V_PB);
	real ans = -2 * V_PA + (V_PB * V_PB * V_RTM * theta) / (r * (V_PB - r) * (V_PB - r));
	
	if (ans != ans) { // NaN. && r < 0 ???
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
