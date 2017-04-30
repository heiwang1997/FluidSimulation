#include "stdafx.h"
#include <ctime>
#include <climits>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "io.h"
#include "ThermalSolver.h"
#include "TimeStepController.h"
#include "config.h"
#include "Field.h"
#include <direct.h>

//#define USE_ENERGY_MODEL

// const real ThermalSolver::extraoplateLength = 4.0f;
const real ThermalSolver::R = 1.0f;

template<class T>
inline T scalarClamp(T target, T lower, T upper) {
	if (target > upper) return upper;
	if (target < lower) return lower;
	return target;
}

void ThermalSolver::initializeHeaters()
{
	heaterCount = resX * resZ;
	heatersX = new int[heaterCount];
	heatersY = new int[heaterCount];
	heatersZ = new int[heaterCount];

	int idx = 0;
	for (int x = 0; x < resX; ++x) {
		for (int z = 0; z < resZ; ++z) {
			heatersX[idx] = x;
			heatersY[idx] = 0;
			heatersZ[idx] = z;
			++idx;
		}
	}
	/*

	heaterCount = 0;
	// Initialize using intervals
	int hintvX = config->heaterIntervalX();
	int hintvZ = config->heaterIntervalZ();

	int hCountX = resX / hintvX - 1;
	int hCountZ = resZ / hintvZ - 1;

	heaterCount = hCountX * hCountZ;
	if (hCountX <= 0 || hCountZ <= 0) {
	heaterCount = 0;
	LOG(WARNING) << "No heaters detected.";
	return;
	}
	heatersX = new int[heaterCount];
	heatersY = new int[heaterCount];
	heatersZ = new int[heaterCount];

	int idx = 0;
	for (int x = 0; x < hCountX; ++x) {
	for (int z = 0; z < hCountZ; ++z) {
	heatersX[idx] = (x + 1) * hintvX;
	heatersY[idx] = 0;
	heatersZ[idx] = (z + 1) * hintvZ;
	++idx;
	}
	}
	
	*/
	CHECK_EQ(idx, heaterCount);
	LOG(INFO) << heaterCount << " heaters added.";
}

void ThermalSolver::computeVelocityStar(real dt,
	Field * vxGuessField, Field * vyGuessField, Field * vzGuessField, 
	Field * rhoGuessField, Field* rhsRhoField,
	Field * vxStarField, Field * vyStarField, Field * vzStarField)
{
	// First Step:
	// Solve: Du/Dt = 0
	advectVelocitySemiLagrange(dt, vxGuessField, vyGuessField, vzGuessField,
		vxField, vyField, vzField,
		vxStarField, vyStarField, vzStarField);

	// Second Step:
	// Euler forward u* = u1 + grad(rhs(rho_guess)) * dt
	rhsRhoOnAlignedGrid(rhoGuessField, rhsRhoField, environmentRho, environmentTheta);

	Field* vStarField[3] = {vxStarField, vyStarField, vzStarField};
	real* vStarContent[3] = { vxStarField->content, vyStarField->content, vzStarField->content };
	for (int d = 0; d < 3; ++d) {
		for (int k = (d == 2); k < resZ; ++k) {
			for (int j = (d == 1); j < resY + (d == 1); ++j) {
				for (int i = (d == 0); i < resX; ++i) {
					int vIndex = vStarField[d]->getIndex(i, j, k);
					real rhsPlus = (j == resY) ? rhsRhoField->topSlice[rhsRhoField->getTopSliceIndex(i, k)] :
						rhsRhoField->content[rhsRhoField->getIndex(i, j, k)];
					int centerMinus = rhsRhoField->getIndex(i - (d == 0), j - (d == 1), k - (d == 2));
					vStarContent[d][vIndex] += (dt / h) * (rhsPlus - rhsRhoField->content[centerMinus]);
					// Add gravity at y axis.
					if (d == 1) {
						vStarContent[d][vIndex] -= envGravity * dt;
					}
				}
			}
		}
	}
}

// Rho Prime is calculated using C.E.
// See rhoprime.docx for more info.
void ThermalSolver::computeRhoPrime(real dt, Field* vxStarField, Field* vyStarField, Field* vzStarField,
	Field* rhoGuessField, Field* rhoPrimeField) {
	typedef Eigen::SparseMatrix<real> SpMat;
	typedef Eigen::Triplet<real> Triplet;
	typedef Eigen::VectorXf EVector;

	std::vector<Triplet> coefficients;
	int totalCells = rhoGuessField->totalSize;
	coefficients.reserve(totalCells * 7);

	EVector b(totalCells);

	real* vxStar = vxStarField->content;
	real* vyStar = vyStarField->content;
	real* vzStar = vzStarField->content;
	real* rhoStar = rhoGuessField->content;

	for (int k = 0; k < resZ; ++k) {
		for (int j = 0; j < resY; ++j) {
			for (int i = 0; i < resX; ++i) {
				int center = rhoGuessField->getIndex(i, j, k);
				int left = rhoGuessField->getIndexLoopBoundary(i - 1, j, k);
				int right = rhoGuessField->getIndexLoopBoundary(i + 1, j, k);
				int up = rhoGuessField->getIndexLoopBoundary(i, j + 1, k);
				int bottom = rhoGuessField->getIndexLoopBoundary(i, j - 1, k);
				int front = rhoGuessField->getIndexLoopBoundary(i, j, k - 1);
				int back = rhoGuessField->getIndexLoopBoundary(i, j, k + 1);
				real vRight = vxStar[vxStarField->getIndex(i + 1, j, k)];
				real vLeft = vxStar[vxStarField->getIndex(i, j, k)];
				real vTop = vyStar[vyStarField->getIndex(i, j + 1, k)];
				real vBottom = vyStar[vyStarField->getIndex(i, j, k)];
				real vBack = vzStar[vzStarField->getIndex(i, j, k + 1)];
				real vFront = vzStar[vzStarField->getIndex(i, j, k)];

				coefficients.push_back(Triplet(center, center, 2 * h / dt + vRight - vLeft + vTop - vBottom + vBack - vFront));
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
					+ 2 * h / dt * (rhoStar[center] - rhoField->content[center]));
				// Should be equivalent to (b << rhs;)
				b[center] = rhs;
			}
		}
	}

	SpMat A(totalCells, totalCells);
	A.setFromTriplets(coefficients.begin(), coefficients.end());

	Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<real> > solver;
	EVector x;
	solver.compute(A);
	x = solver.solve(b);
	CHECK_EQ(solver.info(), Eigen::Success) << "Solver Failed";

	memcpy(rhoPrimeField->content, x.data(), sizeof(real) * totalCells);
}

void ThermalSolver::computeVelocityPrime(real dt, 
	Field * rhoGuessField, Field * rhsRhoStarField, Field * rhsRhoStarStarField, 
	Field * vxPrimeField, Field * vyPrimeField, Field * vzPrimeField)
{
	// Fill v' border with zero.
	// Top open boundary will be reset.
	fillVelocityFieldBorderZero(vxPrimeField, vyPrimeField, vzPrimeField);

	rhsRhoOnAlignedGrid(rhoGuessField, rhsRhoStarStarField, environmentRho, environmentTheta);
	Field* vPrimeField[3] = { vxPrimeField, vyPrimeField, vzPrimeField };
	real* vPrimeContent[3] = { vxPrimeField->content, vyPrimeField->content, vzPrimeField->content };
	for (int d = 0; d < 3; ++d) {
		for (int k = (d == 2); k < resZ; ++k) {
			for (int j = (d == 1); j < resY + (d == 1); ++j) {
				for (int i = (d == 0); i < resX; ++i) {
					int vIndex = vPrimeField[d]->getIndex(i, j, k);
					int topSliceIdx = rhsRhoStarField->getTopSliceIndex(i, k);
					int centerPlus = (j == resY) ? -1 : rhsRhoStarField->getIndex(i, j, k);

					real rhsRhoStarStarPlus = (j == resY) ? rhsRhoStarStarField->topSlice[topSliceIdx] :
						rhsRhoStarStarField->content[centerPlus];
					real rhsRhoStarPlus = (j == resY) ? rhsRhoStarField->topSlice[topSliceIdx] :
						rhsRhoStarField->content[centerPlus];

					int centerMinus = rhoGuessField->getIndex(i - (d == 0), j - (d == 1), k - (d == 2));
					vPrimeContent[d][vIndex] = (dt / h) * (rhsRhoStarStarPlus -
						rhsRhoStarStarField->content[centerMinus]);
					vPrimeContent[d][vIndex] -= (dt / h) * (rhsRhoStarPlus -
						rhsRhoStarField->content[centerMinus]);
					// The gravity term cancels.
				}
			}
		}
	}
}

void ThermalSolver::updateThetaField(real dt, Field * vxGuessField, 
	Field * vyGuessField, Field * vzGuessField, Field* rhoGuessField)
{
	// Solve Dtheta/Dt = 0
	advectFieldSemiLagrange(dt, vxGuessField, vyGuessField, vzGuessField, 
		thetaField, thetaGuessField, environmentTheta);

	Field* lapThetaField = new Field(thetaGuessField, false);
	// This step computes lap_Theta with theta_bar.
	laplacianFieldOnAlignedGrid(thetaGuessField, lapThetaField, environmentTheta);

	real* theta = thetaGuessField->content;
	real* lapTheta = lapThetaField->content;

#ifdef USE_ENERGY_MODEL
	Field* lapRhoField = new Field(rhoGuessField, false);
	laplacianFieldOnAlignedGrid(rhoGuessField, lapRhoField, environmentRho);
	`
	real* vx = vxGuessField->content;
	real* vy = vyGuessField->content;
	real* vz = vzGuessField->content;
	real* rho = rhoGuessField->content;
	// Simple aliasing.
	Field* vxF = vxGuessField, *vyF = vyGuessField, *vzF = vzGuessField;
	// Eigen::Vector3f* piArr = new Eigen::Vector3f[6];
	// Get pressure, this enables the 'microPi' optimization, which may save 2 * field space.
	Field* tempPressure = new Field(rhoGuessField, false);
	real* pressure = tempPressure->content;
	for (int k = 0; k < resZ; ++k) {
		for (int j = 0; j < resY; ++j) {
			for (int i = 0; i < resX; ++i) {
				int index = tempPressure->getIndex(i, j, k);
				pressure[index] = vdwEqState(rho[index], theta[index]);
			}
		}
	}

	for (int k = 0; k < resZ; ++k) {
		for (int j = 0; j < resY; ++j) {
			for (int i = 0; i < resX; ++i) {
				int index = thetaGuessField->getIndex(i, j, k);
				int left = thetaGuessField->getIndexClampBoundary(i - 1, j, k);
				int right = thetaGuessField->getIndexClampBoundary(i + 1, j, k);
				int up = thetaGuessField->getIndexClampBoundary(i, j + 1, k);
				int down = thetaGuessField->getIndexClampBoundary(i, j - 1, k);
				int front = thetaGuessField->getIndexClampBoundary(i, j, k - 1);
				int back = thetaGuessField->getIndexClampBoundary(i, j, k + 1);
				real rhoUp = (j == resY - 1) ? environmentRho : rho[up];

				Eigen::Vector3f gradRho;
				gradRho << (rho[right] - rho[left]) / (2 * h),
					(rhoUp - rho[down]) / (2 * h),
					(rho[back] - rho[front]) / (2 * h);
				// TODO!
				Eigen::Matrix3f kortewegTensor = vdwInvWe * (
					(rho[index] * lapRhoField->content[index] +
						0.5 * gradRho.squaredNorm()) * Eigen::MatrixXf::Identity(3, 3) -
					gradRho * gradRho.transpose());

				Eigen::Matrix3f cauchyTensor =
					-pressure[index] * Eigen::MatrixXf::Identity(3, 3) +
					kortewegTensor;

				// grad u is a tensor
				Eigen::Matrix3f gradVelocity;
				gradVelocity << (vx[vxF->getIndex(i + 1, j, k)] - vx[vxF->getIndex(i, j, k)]) / h,
					(vx[vxF->getIndexClampBoundary(i, j + 1, k)] + vx[vxF->getIndexClampBoundary(i + 1, j + 1, k)]
						- vx[vxF->getIndexClampBoundary(i, j - 1, k)] - vx[vxF->getIndexClampBoundary(i + 1, j - 1, k)]) / (4 * h),
						(vx[vxF->getIndexClampBoundary(i, j, k + 1)] + vx[vxF->getIndexClampBoundary(i + 1, j, k + 1)]
							- vx[vxF->getIndexClampBoundary(i, j, k - 1)] - vx[vxF->getIndexClampBoundary(i + 1, j, k - 1)]) / (4 * h),
							(vy[vyF->getIndexClampBoundary(i + 1, j, k)] + vy[vyF->getIndexClampBoundary(i + 1, j + 1, k)]
								- vy[vyF->getIndexClampBoundary(i - 1, j, k)] - vy[vyF->getIndexClampBoundary(i - 1, j + 1, k)]) / (4 * h),
								(vy[vyF->getIndex(i, j + 1, k)] - vy[vyF->getIndex(i, j, k)]) / h,
					(vy[vyF->getIndexClampBoundary(i, j, k + 1)] + vy[vyF->getIndexClampBoundary(i, j + 1, k + 1)]
						- vy[vyF->getIndexClampBoundary(i, j, k - 1)] - vy[vyF->getIndexClampBoundary(i, j + 1, k - 1)]) / (4 * h),
						(vz[vzF->getIndexClampBoundary(i + 1, j, k)] + vz[vzF->getIndexClampBoundary(i + 1, j, k + 1)]
							- vz[vzF->getIndexClampBoundary(i - 1, j, k)] - vz[vzF->getIndexClampBoundary(i - 1, j, k + 1)]) / (4 * h),
							(vz[vzF->getIndexClampBoundary(i, j + 1, k)] + vz[vzF->getIndexClampBoundary(i, j + 1, k + 1)]
								- vz[vzF->getIndexClampBoundary(i, j - 1, k)] - vz[vzF->getIndexClampBoundary(i, j - 1, k + 1)]) / (4 * h),
								(vz[vzF->getIndex(i, j, k + 1)] - vz[vzF->getIndex(i, j, k)]) / h;

				// div u is a scalar
				real divVelocity = gradVelocity.trace();

				Eigen::Vector3f microPi = vdwInvWe * rho[index] * divVelocity * gradRho;
				//if (i == 113 && j == 109 && k == 8) {
				//	piArr[0] = microPi;
				//}
				//if (i == 115 && j == 109 && k == 8) piArr[1] = microPi;
				//if (i == 114 && j == 108 && k == 8) piArr[2] = microPi;
				//if (i == 114 && j == 110 && k == 8) piArr[3] = microPi;
				//if (i == 114 && j == 109 && k == 7) piArr[4] = microPi;
				//if (i == 114 && j == 109 && k == 9) piArr[5] = microPi;
				// MicroPi is added to element needed it.
				if (i >= 1) theta[left] += 0.5f * (dt / h) * (-microPi[0]) / (rho[left] * vdwCv);
				if (i < resX - 1) theta[right] -= 0.5f * (dt / h) * (-microPi[0]) / (rho[right] * vdwCv);
				if (j >= 1) theta[down] += 0.5f * (dt / h) * (-microPi[1]) / (rho[down] * vdwCv);
				if (j < resY - 1) theta[up] -= 0.5f * (dt / h) * (-microPi[1]) / (rho[up] * vdwCv);
				if (k >= 1) theta[front] += 0.5f * (dt / h) * (-microPi[2]) / (rho[front] * vdwCv);
				if (k < resZ - 1) theta[back] -= 0.5f * (dt / h) * (-microPi[2]) / (rho[back] * vdwCv);
				// For boundary elements, updates are missed, get them back here.
				real centerCoef = 0.5f * (dt / h) / (rho[index] * vdwCv);
				if (i == 0) theta[index] -= centerCoef * (-microPi[0]);
				if (i == resX - 1) theta[index] += centerCoef * (-microPi[0]);
				if (j == 0) theta[index] -= centerCoef * (-microPi[1]);
				if (j == resY - 1) {
					// MicroPi at (i, resY, k), divVelocity is same due to extrapolation.
					real microPi1 = vdwInvWe * environmentRho * divVelocity * 
						(environmentRho - rho[index]) / (2 * h);
					theta[index] += centerCoef * (-microPi1);
					//theta[index] += centerCoef * (-microPi[1]);
				}
				if (k == 0) theta[index] -= centerCoef * (-microPi[2]);
				if (k == resZ - 1) theta[index] += centerCoef * (-microPi[2]);

				//real motherFucker = dt * ((-vdwA / vdwCv * rho[index] * divVelocity));
				//real fatherFucker = dt * (
				//	((cauchyTensor.cwiseProduct(gradVelocity.transpose()).sum())) / (rho[index] * vdwCv));
				//real diffuse = dt * (
				//	(heatDiffuseSpeed * lapThetaField->content[i]) / (rho[index] * vdwCv));

				//if (j == 109) {
				//	std::cout << i << ' ' << j << ' ' << k << ' ' << rho[index] << std::endl;
				//	std::cout << motherFucker << ' ' << fatherFucker << ' ' << diffuse << std::endl;
				//}

				//divVelocity << ' ' << (vy[vyF->getIndex(i, j + 1, k)] - vy[vyF->getIndex(i, j, k)]) << std::endl;
				//if (motherFucker < -10) {
				//	std::cout << "MF: " << motherFucker << std::endl;
				//}
				//if (fatherFucker < -10) {
				//	std::cout << "FF: " << fatherFucker << std::endl;
				//}

				theta[index] += dt * ((-vdwA / vdwCv * rho[index] * divVelocity) +
					((cauchyTensor.cwiseProduct(gradVelocity.transpose()).sum()) +
						heatDiffuseSpeed * lapThetaField->content[i]) / (rho[index] * vdwCv));
			}
		}
	}
	//for (int i = 0; i < 6; ++i) {
	//	std::cout << piArr[i] << std::endl;
	//}
	//std::cout << (piArr[1][0] - piArr[0][0]) / (2 * h) + (piArr[3][1] - piArr[2][1]) / (2 * h) + (piArr[5][2] - piArr[4][2]) / (2 * h) << std::endl;
	// Additional heat term is added using splitting.
	for (int i = 0; i < heaterCount; ++i) {
		int heaterIndex = thetaField->getIndex(heatersX[i], heatersY[i], heatersZ[i]);
		theta[heaterIndex] += (1 - exp(- heatSpeed * dt)) * (targetTheta - theta[heaterIndex]);
		// theta[heaterIndex] = targetTheta;
	}
	/*
	------------ BELOW: NOT APPLICABLE ANY MORE -------------------
	for (int i = 0; i < thetaField->totalSize; ++i) {
	theta[i] += dt * heatDiffuseSpeed * lapThetaField->content[i];
	}
	------------ ABOVE: NOT APPLICABLE ANY MORE -------------------
	*/
	delete tempPressure;
	delete lapRhoField;
#else
	for (int i = 0; i < thetaField->totalSize; ++i) {
		theta[i] += dt * heatDiffuseSpeed * lapThetaField->content[i];
	}
	// Additional heat term is added using splitting.
	for (int i = 0; i < heaterCount; ++i) {
		int heaterIndex = thetaField->getIndex(heatersX[i], heatersY[i], heatersZ[i]);
		theta[heaterIndex] += (1 - exp(- heatSpeed * dt)) * (targetTheta - theta[heaterIndex]);
		// theta[heaterIndex] = targetTheta;
	}
#endif // USE_ENERGY_MODEL
	LOG(INFO) << "Field after update: MAX = " << thetaGuessField->getMax() <<
		", MIN = " << thetaGuessField->getMin();

	delete lapThetaField;
}

real ThermalSolver::updateRhoField(Field * rhoGuess, Field * rhoPrime)
{
	real delta = 0.0f;
	for (int i = 0; i < rhoGuess->totalSize; ++i) {
		rhoGuess->content[i] += rhoRelaxCoef * rhoPrime->content[i];
		delta += (rhoPrime->content[i] * rhoPrime->content[i]);
	}
	return delta;
}

real ThermalSolver::updateVelocityField(Field * vxStarField, Field * vyStarField, Field * vzStarField,
	Field * vxPrimeField, Field * vyPrimeField, Field * vzPrimeField, 
	Field * vxGuessField, Field * vyGuessField, Field * vzGuessField)
{
	real delta = 0.0f;
	for (int i = 0; i < vxStarField->totalSize; i++) {
		real finalVx = vxStarField->content[i] + vxPrimeField->content[i];
		real vxDelta = finalVx - vxGuessField->content[i];
		vxGuessField->content[i] = finalVx;
		delta += vxDelta * vxDelta;
	}
	for (int i = 0; i < vyStarField->totalSize; i++) {
		real finalVy = vyStarField->content[i] + vyPrimeField->content[i];
		real vyDelta = finalVy - vyGuessField->content[i];
		vyGuessField->content[i] = finalVy;
		delta += vyDelta * vyDelta;
	}
	for (int i = 0; i < vzStarField->totalSize; i++) {
		real finalVz = vzStarField->content[i] + vzPrimeField->content[i];
		real vzDelta = finalVz - vzGuessField->content[i];
		vzGuessField->content[i] = finalVz;
		delta += vzDelta * vzDelta;
	}
	return delta;
}

void ThermalSolver::laplacianFieldOnAlignedGrid(Field* f, Field* lapF, real envVal)
{
	real* fc = f->content;
	real* lapFc = lapF->content;
	bool useEnvVal = true;
	for (int k = 0; k < resZ; ++ k) {
		for (int j = 0; j < resY; ++ j) {
			for (int i = 0; i < resX; ++ i) {
				int center = f->getIndex(i, j, k);
				int left = f->getIndexClampBoundary(i - 1, j, k);
				int right = f->getIndexClampBoundary(i + 1, j, k);
				int up = f->getIndexClampBoundary(i, j + 1, k);
				int bottom = f->getIndexClampBoundary(i, j - 1, k);
				int front = f->getIndexClampBoundary(i, j, k - 1);
				int back = f->getIndexClampBoundary(i, j, k + 1);
				real fcUp = (useEnvVal && (j == resY - 1)) ? envVal : fc[up];
				lapFc[center] = (fc[left] + fc[right] +
					fcUp + fc[bottom] + fc[front] + fc[back] - 
					6 * fc[center]) / h / h;
			}
		}
	}
}

void ThermalSolver::rhsRhoOnAlignedGrid(Field * rF, Field * rhsRF, real envRho, real envTheta)
{
	real* rho = rF->content;
	real* rhsRho = rhsRF->content;
	bool useEnvRho = true;

	for (int k = 0; k < resZ; ++k) {
		for (int j = 0; j < resY; ++j) {
			for (int i = 0; i < resX; ++i) {
				int center = rF->getIndex(i, j, k);
				int left = rF->getIndexClampBoundary(i - 1, j, k);
				int right = rF->getIndexClampBoundary(i + 1, j, k);
				int up = rF->getIndexClampBoundary(i, j + 1, k);
				int bottom = rF->getIndexClampBoundary(i, j - 1, k);
				int front = rF->getIndexClampBoundary(i, j, k - 1);
				int back = rF->getIndexClampBoundary(i, j, k + 1);
				real rhoUp = (useEnvRho && (j == resY - 1)) ? envRho : rho[up];
				rhsRho[center] = - thermalWd(rho[center], thetaGuessField->content[center]);
				rhsRho[center] += (vdwInvWe / h / h) * (rho[left] + rho[right] +
					rhoUp + rho[bottom] + rho[front] + rho[back] -
					6 * rho[center]);
			}
		}
	}

	real* rhsRhoTopSlice = rhsRF->topSlice;
	for (int k = 0; k < resZ; ++k) {
		for (int i = 0; i < resX; ++i) {
			int topSliceCenter = rhsRF->getTopSliceIndex(i, k);
			int topSliceBottom = rF->getIndex(i, resY - 1, k);
			rhsRhoTopSlice[topSliceCenter] = -thermalWd(envRho, envTheta);
			rhsRhoTopSlice[topSliceCenter] += (vdwInvWe / h / h) * (rho[topSliceBottom] - envRho);
		}
	}

	// std::cout << rhsRho[rhsRF->getIndex(64, 127, 16)] << ' ' << rhsRhoTopSlice[rhsRF->getTopSliceIndex(64, 16)] << std::endl;
}

void ThermalSolver::advectVelocitySemiLagrange(real dt,
	Field * vxBackgroundField, Field * vyBackgroundField, Field * vzBackgroundField,
	Field * vxInterimField, Field *vyInterimField, Field *vzInterimField,
	Field * vxNewField, Field * vyNewField, Field * vzNewField) {
	// Get the content of the fields.
	real* vxBackground = vxBackgroundField->content;
	real* vyBackground = vyBackgroundField->content; 
	real* vzBackground = vzBackgroundField->content;
	real* vxInterim = vxInterimField->content; 
	real* vyInterim = vyInterimField->content; 
	real* vzInterim = vzInterimField->content;
	real* vxNew = vxNewField->content; 
	real* vyNew = vyNewField->content; 
	real* vzNew = vzNewField->content;
	// Fill the velocity at solid still boundaries
	// The open boundary at the top will be fixed below.
	fillVelocityFieldBorderZero(vxNewField, vyNewField, vzNewField);
	// velocity X
	for (int z = 0; z < resZ; z++) {
		for (int y = 0; y < resY; y++) {
			for (int x = 1; x < resX; x++) {
				int index = vxBackgroundField->getIndex(x, y, z);
				real velx = vxBackground[index];
				real vely = (
					vyBackground[vyBackgroundField->getIndex(x - 1, y,     z)] +
					vyBackground[vyBackgroundField->getIndex(x - 1, y + 1, z)] +
					vyBackground[vyBackgroundField->getIndex(x    , y    , z)] +
					vyBackground[vyBackgroundField->getIndex(x    , y + 1, z)]) * 0.25f;
				real velz = (
					vzBackground[vzBackgroundField->getIndex(x    , y, z)] +
					vzBackground[vzBackgroundField->getIndex(x - 1, y, z)] +
					vzBackground[vzBackgroundField->getIndex(x    , y, z + 1)] +
					vzBackground[vzBackgroundField->getIndex(x - 1, y, z + 1)]) * 0.25f;

				// backtrace
				real xTrace = x - (dt / h) * velx;
				real yTrace = y - (dt / h) * vely;
				real zTrace = z - (dt / h) * velz;

				// Detect y coord is at open boundary, true if obo > 0.
				// real openBoundaryOverflow = yTrace - (real)(resY - 1.0f);

				// clamp backtrace to grid boundaries
				xTrace = scalarClamp(xTrace, 0.0f, (real)resX);
				yTrace = scalarClamp(yTrace, -0.5f, (real)(resY - 0.5f));
				zTrace = scalarClamp(zTrace, -0.5f, (real)(resZ - 0.5f));

				// locate neighbors to interpolate
				const int x0 = ifloor(xTrace);
				const int x1 = x0 + 1;
				const int y0 = ifloor(yTrace);
				const int y1 = y0 + 1;
				const int z0 = ifloor(zTrace);
				const int z1 = z0 + 1;

				// get interpolation weights
				const real s1 = xTrace - floor(xTrace);
				const real s0 = 1.0f - s1;
				const real t1 = yTrace - floor(yTrace);
				const real t0 = 1.0f - t1;
				const real u1 = zTrace - floor(zTrace);
				const real u0 = 1.0f - u1;

				const int i000 = vxBackgroundField->getIndexClampBoundary(x0, y0, z0);
				const int i010 = vxBackgroundField->getIndexClampBoundary(x0, y1, z0);
				const int i100 = vxBackgroundField->getIndexClampBoundary(x1, y0, z0);
				const int i110 = vxBackgroundField->getIndexClampBoundary(x1, y1, z0);
				const int i001 = vxBackgroundField->getIndexClampBoundary(x0, y0, z1);
				const int i011 = vxBackgroundField->getIndexClampBoundary(x0, y1, z1);
				const int i101 = vxBackgroundField->getIndexClampBoundary(x1, y0, z1);
				const int i111 = vxBackgroundField->getIndexClampBoundary(x1, y1, z1);

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
	// Note: modified for open boundary.
	for (int z = 0; z < resZ; z++) {
		for (int y = 1; y /* < */ <= resY; y++) {
			for (int x = 0; x < resX; x++) {
				int index = vyBackgroundField->getIndex(x, y, z);
				real velx = (
					vxBackground[vxBackgroundField->getIndex(x, y - 1, z)] +
					vxBackground[vxBackgroundField->getIndexClampBoundary(x, y, z)] +
					vxBackground[vxBackgroundField->getIndex(x + 1, y - 1, z)] +
					vxBackground[vxBackgroundField->getIndexClampBoundary(x + 1, y, z)]) * 0.25f;
				real vely = vyBackground[index];
				real velz = (
					vzBackground[vzBackgroundField->getIndexClampBoundary(x, y, z)] +
					vzBackground[vzBackgroundField->getIndexClampBoundary(x, y, z + 1)] +
					vzBackground[vzBackgroundField->getIndex(x, y - 1, z)] +
					vzBackground[vzBackgroundField->getIndex(x, y - 1, z + 1)]) * 0.25f;

				// backtrace
				real xTrace = x - (dt / h) * velx;
				real yTrace = y - (dt / h) * vely;
				real zTrace = z - (dt / h) * velz;

				// clamp backtrace to grid boundaries
				xTrace = scalarClamp(xTrace, -0.5f, (real)(resX - 0.5f));
				yTrace = scalarClamp(yTrace, 0.0f, (real)resY);
				zTrace = scalarClamp(zTrace, -0.5f, (real)(resZ - 0.5f));

				// locate neighbors to interpolate
				const int x0 = ifloor(xTrace);
				const int x1 = x0 + 1;
				const int y0 = ifloor(yTrace);
				const int y1 = y0 + 1;
				const int z0 = ifloor(zTrace);
				const int z1 = z0 + 1;

				// get interpolation weights
				const real s1 = xTrace - floor(xTrace);
				const real s0 = 1.0f - s1;
				const real t1 = yTrace - floor(yTrace);
				const real t0 = 1.0f - t1;
				const real u1 = zTrace - floor(zTrace);
				const real u0 = 1.0f - u1;

				const int i000 = vyBackgroundField->getIndexClampBoundary(x0, y0, z0);
				const int i010 = vyBackgroundField->getIndexClampBoundary(x0, y1, z0);
				const int i100 = vyBackgroundField->getIndexClampBoundary(x1, y0, z0);
				const int i110 = vyBackgroundField->getIndexClampBoundary(x1, y1, z0);
				const int i001 = vyBackgroundField->getIndexClampBoundary(x0, y0, z1);
				const int i011 = vyBackgroundField->getIndexClampBoundary(x0, y1, z1);
				const int i101 = vyBackgroundField->getIndexClampBoundary(x1, y0, z1);
				const int i111 = vyBackgroundField->getIndexClampBoundary(x1, y1, z1);

				// interpolate
				vyNew[index] = u0 * (s0 * (t0 * vyInterim[i000] +
					t1 * vyInterim[i010]) +
					s1 * (t0 * vyInterim[i100] +
						t1 * vyInterim[i110])) +
					u1 * (s0 * (t0 * vyInterim[i001] +
						t1 * vyInterim[i011]) +
						s1 * (t0 * vyInterim[i101] +
							t1 * vyInterim[i111]));
				//if (x == 64 && y == 128 && z == 16) {
				//	std::cout << "Advect: " << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1 <<
				//		' ' << z0 << ' ' << z1 << std::endl;
				//	std::cout << "Weight: " << s0 << ' ' << s1 << ' ' << t0 << ' ' << t1 <<
				//		' ' << u0 << ' ' << u1 << std::endl;
				//	std::cout << "Result: " << vyNew[index] << std::endl;
				//}
			}
		}
	}
	// velocity Z
	for (int z = 1; z < resZ; z++) {
		for (int y = 0; y < resY; y++) {
			for (int x = 0; x < resX; x++) {
				int index = vzBackgroundField->getIndex(x, y, z);
				real velx = (
					vxBackground[vxBackgroundField->getIndex(x, y, z - 1)] +
					vxBackground[vxBackgroundField->getIndex(x, y, z)] +
					vxBackground[vxBackgroundField->getIndex(x + 1, y, z - 1)] +
					vxBackground[vxBackgroundField->getIndex(x + 1, y, z)]) * 0.25f;
				real vely = (
					vyBackground[vyBackgroundField->getIndex(x, y, z - 1)] +
					vyBackground[vyBackgroundField->getIndex(x, y + 1, z - 1)] +
					vyBackground[vyBackgroundField->getIndex(x, y, z)] +
					vyBackground[vyBackgroundField->getIndex(x, y + 1, z)]) * 0.25f;
				real velz = vzBackground[index];

				// backtrace
				real xTrace = x - (dt / h) * velx;
				real yTrace = y - (dt / h) * vely;
				real zTrace = z - (dt / h) * velz;

				// clamp backtrace to grid boundaries
				xTrace = scalarClamp(xTrace, -0.5f, (real)(resX - 0.5f));
				yTrace = scalarClamp(yTrace, -0.5f, (real)(resY - 0.5f));
				zTrace = scalarClamp(zTrace, 0.0f, (real)resZ);

				// locate neighbors to interpolate
				const int x0 = ifloor(xTrace);
				const int x1 = x0 + 1;
				const int y0 = ifloor(yTrace);
				const int y1 = y0 + 1;
				const int z0 = ifloor(zTrace);
				const int z1 = z0 + 1;

				// get interpolation weights
				const real s1 = xTrace - floor(xTrace);
				const real s0 = 1.0f - s1;
				const real t1 = yTrace - floor(yTrace);
				const real t0 = 1.0f - t1;
				const real u1 = zTrace - floor(zTrace);
				const real u0 = 1.0f - u1;

				const int i000 = vzBackgroundField->getIndexClampBoundary(x0, y0, z0);
				const int i010 = vzBackgroundField->getIndexClampBoundary(x0, y1, z0);
				const int i100 = vzBackgroundField->getIndexClampBoundary(x1, y0, z0);
				const int i110 = vzBackgroundField->getIndexClampBoundary(x1, y1, z0);
				const int i001 = vzBackgroundField->getIndexClampBoundary(x0, y0, z1);
				const int i011 = vzBackgroundField->getIndexClampBoundary(x0, y1, z1);
				const int i101 = vzBackgroundField->getIndexClampBoundary(x1, y0, z1);
				const int i111 = vzBackgroundField->getIndexClampBoundary(x1, y1, z1);

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


void ThermalSolver::advectFieldSemiLagrange(real dt, Field *vxBackgroundField, 
	Field *vyBackgroundField, Field *vzBackgroundField,
	Field *oldField, Field *newField, real envValue) {

	// Get the content of the fields.
	real* vx = vxBackgroundField->content;
	real* vy = vyBackgroundField->content;
	real* vz = vzBackgroundField->content;
	real* oldFieldContent = oldField->content;
	real* newFieldContent = newField->content;
	
	// If envValue != nan, USE ENV.
	bool useEnvValue = true;

	for (int z = 0; z < resZ; z++) {
		for (int y = 0; y < resY; y++) {
			for (int x = 0; x < resX; x++) {
				int index = newField->getIndex(x, y, z);

				int vxLeftIndex = vxBackgroundField->getIndex(x, y, z);
				int vxRightIndex = vxBackgroundField->getIndex(x + 1, y, z);
				real xTrace = x - (dt / h) * (vx[vxLeftIndex] + vx[vxRightIndex]) * 0.5f;

				int vyBelowIndex = vyBackgroundField->getIndex(x, y, z);
				int vyAboveIndex = vyBackgroundField->getIndex(x, y + 1, z);
				real yTrace = y - (dt / h) * (vy[vyBelowIndex] + vy[vyAboveIndex]) * 0.5f;

				int vzFrontIndex = vzBackgroundField->getIndex(x, y, z);
				int vzBackIndex = vzBackgroundField->getIndex(x, y, z + 1);
				real zTrace = z - (dt / h) * (vz[vzFrontIndex] + vz[vzBackIndex]) * 0.5f;

				if (useEnvValue && yTrace > (real)(resY - 0.5f)) {
					// Clamp x, z and give y.
					newFieldContent[index] = envValue;
					continue;
				}

				// clamp backtrace to grid boundaries
				xTrace = scalarClamp(xTrace, -0.5f, (real)(resX - 0.5f));
				yTrace = scalarClamp(yTrace, -0.5f, (real)(resY - 0.5f));
				zTrace = scalarClamp(zTrace, -0.5f, (real)(resZ - 0.5f));

				// locate neighbors to interpolate
				const int x0 = ifloor(xTrace);
				const int x1 = x0 + 1;
				const int y0 = ifloor(yTrace);
				const int y1 = y0 + 1;
				const int z0 = ifloor(zTrace);
				const int z1 = z0 + 1;

				// get interpolation weights
				const real s1 = xTrace - floor(xTrace);
				const real s0 = 1.0f - s1;
				const real t1 = yTrace - floor(yTrace);
				const real t0 = 1.0f - t1;
				const real u1 = zTrace - floor(zTrace);
				const real u0 = 1.0f - u1;

				const int i000 = oldField->getIndexClampBoundary(x0, y0, z0);
				const int i010 = oldField->getIndexClampBoundary(x0, y1, z0);
				const int i100 = oldField->getIndexClampBoundary(x1, y0, z0);
				const int i110 = oldField->getIndexClampBoundary(x1, y1, z0);
				const int i001 = oldField->getIndexClampBoundary(x0, y0, z1);
				const int i011 = oldField->getIndexClampBoundary(x0, y1, z1);
				const int i101 = oldField->getIndexClampBoundary(x1, y0, z1);
				const int i111 = oldField->getIndexClampBoundary(x1, y1, z1);

				// interpolates
				newFieldContent[index] = u0 * (s0 * (t0 * oldFieldContent[i000] +
					t1 * oldFieldContent[i010]) +
					s1 * (t0 * oldFieldContent[i100] +
						t1 * oldFieldContent[i110])) +
					u1 * (s0 * (t0 * oldFieldContent[i001] +
						t1 * oldFieldContent[i011]) +
						s1 * (t0 * oldFieldContent[i101] +
							t1 * oldFieldContent[i111]));

			}
		}
	}
}

void ThermalSolver::fillVelocityFieldBorderZero(Field * xF, Field * yF, Field * zF)
{
	real* xFc = xF->content;
	real* yFc = yF->content;
	real* zFc = zF->content;
	for (int k = 0; k < resZ; k++) {
		for (int j = 0; j < resY; j++) {
			xFc[xF->getIndex(0, j, k)] = 0.0f;
			xFc[xF->getIndex(resX, j, k)] = 0.0f;
		}
	}
	for (int i = 0; i < resX; i++) {
		for (int k = 0; k < resZ; k++) {
			yFc[yF->getIndex(i, 0, k)] = 0.0f;
			yFc[yF->getIndex(i, resY, k)] = 0.0f;
		}
	}
	for (int i = 0; i < resX; i++) {
		for (int j = 0; j < resY; j++) {
			zFc[zF->getIndex(i, j, 0)] = 0.0f;
			zFc[zF->getIndex(i, j, resZ)] = 0.0f;
		}
	}
}

void ThermalSolver::run(TimeStepController* timeStep)
{
	int& stepCount = timeStep->_stepCount;
	int snapshotInterval = config->snapshotInterval();
	while (!timeStep->isFinished()) {
		std::clock_t startTime = std::clock();
		stepSimple(timeStep->getStepDt());
		int fmCnt;
		if (timeStep->isFrameTime(fmCnt)) {
			LOG(INFO) << "Frame " << fmCnt << " generated";
		}
		std::clock_t endTime = std::clock();
		stepCount += 1;
		LOG(INFO) << "Time step " << stepCount << " finished, this step took " <<
			(real)((endTime - startTime) / CLOCKS_PER_SEC) << "s";
		// Write preview to destination folder.
		{
			std::string baseFolder = config->fieldOutputDir() + "/" + config->runName() + "/";
			_mkdir(baseFolder.c_str());
			auto getFieldOutputFilename = [&](std::string pfx) {
				return baseFolder + pfx + std::to_string(stepCount);
			};
			rhoField->writeSlabPreviewToFile(getFieldOutputFilename("rho"));
			vxField->writeSlabPreviewToFile(getFieldOutputFilename("vx"));
			vyField->writeSlabPreviewToFile(getFieldOutputFilename("vy"));
			thetaField->writeSlabPreviewToFile(getFieldOutputFilename("theta"));
		}
		// Dump to file every snapshot interval
		{
			if (stepCount % snapshotInterval == 0) {
				std::string baseFolder = config->snapshotOutputDir() + "/" + config->runName() + "/";
				_mkdir(baseFolder.c_str());
				std::string snapshotFilename = baseFolder + std::to_string(stepCount);
				LOG(INFO) << "Dumping solver state to " << snapshotFilename;
				io::dumpSolverToFile(snapshotFilename, rhoField, vxField, vyField, vzField, timeStep);
			}
		}
	}
}

ThermalSolver::ThermalSolver(Config * cfg, Field * initRhoField, Field * initVxField, 
	Field * initVyField, Field * initVzField, Field* initThetaField)
{
	resX = cfg->resX();
	resY = cfg->resY();
	resZ = cfg->resZ();
	h = cfg->h();

	vdwA = cfg->vdwPA();
	vdwB = cfg->vdwPB();
	vdwInvWe = 1.0f / cfg->vdwPWE();
	vdwCv = cfg->vdwCv();

	velConvergeTol = cfg->velConvergeTol();
	rhoConvergeTol = cfg->rhoConvergeTol();
	rhoRelaxCoef = cfg->rhoRelaxCoef();

	envGravity = cfg->gravity();
	heatDiffuseSpeed = cfg->heatDiffuse();
	targetTheta = cfg->targetTheta();
	heatSpeed = cfg->heatSpeed();

	environmentRho = cfg->envRho();
	environmentTheta = cfg->startTheta();

	if (initRhoField) rhoField = new Field(initRhoField);
	else rhoField = new Field(resX, resY, resZ, true);

	if (initVxField) vxField = new Field(initVxField);
	else vxField = new Field(resX + 1, resY, resZ, true);

	if (initVyField) vyField = new Field(initVyField);
	else vyField = new Field(resX, resY + 1, resZ, true);

	if (initVzField) vzField = new Field(initVzField);
	else vzField = new Field(resX, resY, resZ + 1, true);

	if (initThetaField) thetaField = new Field(initThetaField);
	else thetaField = new Field(resX, resY, resZ, true);

	config = cfg;

	initializeHeaters();
}

void ThermalSolver::stepSimple(real dt)
{
	Field* vxGuessField = new Field(vxField);
	Field* vyGuessField = new Field(vyField);
	Field* vzGuessField = new Field(vzField);
	Field* rhoGuessField = new Field(rhoField);
	thetaGuessField = new Field(thetaField);
	// For vStarField.
	Field* vxStarField = new Field(vxField, false);
	Field* vyStarField = new Field(vyField, false);
	Field* vzStarField = new Field(vzField, false);
	Field* rhoPrimeField = new Field(rhoField, false);
	// For rhoStarField.
	Field* rhsRhoStarField = new Field(rhoField, false);
	rhsRhoStarField->initTopSlice();
	Field* rhsRhoStarStarField = new Field(rhoField, false);
	rhsRhoStarStarField->initTopSlice();

	int loopCount = 0;
	real lastRhoDelta = std::numeric_limits<real>::max();
	real lastVelDelta = std::numeric_limits<real>::max();
	while (true) {
		loopCount += 1;
		LOG(INFO) << "Inner loop: # " << loopCount;

		// 1. Assemble and solve momentum equation for v*
		computeVelocityStar(dt, vxGuessField, vyGuessField, vzGuessField, rhoGuessField,
			rhsRhoStarField, vxStarField, vyStarField, vzStarField);
		// 2. Let rho* = rho Guess
		// 3. Assemble and solve continuity equation for rho'
		computeRhoPrime(dt, vxStarField, vyStarField, vzStarField, rhoGuessField, rhoPrimeField);
		// Correct rho now to be rhoStarStar, 
		// Notice: only for memory efficiency.
		real rhoDelta = updateRhoField(rhoGuessField, rhoPrimeField);
		// 4. Compute V prime. (move the news to outer for better performance)
		Field* vxPrimeField = new Field(vxField, false);
		Field* vyPrimeField = new Field(vyField, false);
		Field* vzPrimeField = new Field(vzField, false);
		// rhsRhoStarStarField is computed through rhoGuessField
		computeVelocityPrime(dt, rhoGuessField,
			rhsRhoStarField, rhsRhoStarStarField, vxPrimeField, vyPrimeField, vzPrimeField);
		// Correct v field.
		real velDelta = updateVelocityField(vxStarField, vyStarField, vzStarField,
			vxPrimeField, vyPrimeField, vzPrimeField, vxGuessField, vyGuessField, vzGuessField);
		delete vxPrimeField;
		delete vyPrimeField;
		delete vzPrimeField;
		// 5. Assemble T* using u**, rho** (Guess fields here)
		updateThetaField(dt, vxGuessField, vyGuessField, vzGuessField, rhoGuessField);
		// 6. Converge judgement
		LOG(INFO) << "Rho diff = " << rhoDelta << ", Velocity diff = " << velDelta;
		if (rhoDelta < rhoConvergeTol && velDelta < velConvergeTol) {
			break;
		}
		else {
			if (rhoDelta > lastRhoDelta || velDelta > lastVelDelta) {
				LOG(WARNING) << "Converge Delta is larger than previous iteration.";
			}
			lastRhoDelta = rhoDelta;
			lastVelDelta = velDelta;
		}
	}
	// Update current fields.
	rhoField->copyFrom(rhoGuessField);
	vxField->copyFrom(vxGuessField);
	vyField->copyFrom(vyGuessField);
	vzField->copyFrom(vzGuessField);
	thetaField->copyFrom(thetaGuessField);
	delete vxGuessField;
	delete vyGuessField;
	delete vzGuessField;
	delete rhoGuessField;
	delete thetaGuessField;
	delete vxStarField;
	delete vyStarField;
	delete vzStarField;
	delete rhoPrimeField;
	delete rhsRhoStarField;
	delete rhsRhoStarStarField;
	// std::cout << "VY last update = " << vyField->content[vyField->getIndex(64, 128, 16)] << std::endl;
	// Last step: update temperature.
	// updateThetaField(dt, vxField, vyField, vzField);
}

ThermalSolver::~ThermalSolver()
{
	delete rhoField;
	delete vxField;
	delete vyField;
	delete vzField;
	delete thetaField;
	delete[] heatersX;
	delete[] heatersY;
	delete[] heatersZ;
}
