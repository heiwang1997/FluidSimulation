#include "stdafx.h"
#include <ctime>
#include <climits>
#include <Eigen/Sparse>
#include "io.h"
#include "ThermalSolver.h"
#include "TimeStepController.h"
#include "config.h"
#include "Field.h"
#include <direct.h>

template<class T>
inline T scalarClamp(T target, T lower, T upper) {
	if (target > upper) return upper;
	if (target < lower) return lower;
	return target;
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
	rhsRhoOnAlignedGrid(rhoGuessField, rhsRhoField);
	Field* vStarField[3] = {vxStarField, vyStarField, vzStarField};
	real* vStarContent[3] = { vxStarField->content, vyStarField->content, vzStarField->content };
	for (int d = 0; d < 3; ++d) {
		for (int k = (d == 2); k < resZ; ++k) {
			for (int j = (d == 1); j < resY; ++j) {
				for (int i = (d == 0); i < resX; ++i) {
					int vIndex = vStarField[d]->getIndex(i, j, k);
					int centerPlus = rhsRhoField->getIndex(i, j, k);
					int centerMinus = rhsRhoField->getIndex(i - (d == 0), j - (d == 1), k - (d == 2));
					vStarContent[d][vIndex] += (dt / h) * (rhsRhoField->content[centerPlus] -
						rhsRhoField->content[centerMinus]);
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
	fillVelocityFieldBorderZero(vxPrimeField, vyPrimeField, vzPrimeField);

	rhsRhoOnAlignedGrid(rhoGuessField, rhsRhoStarStarField);
	Field* vPrimeField[3] = { vxPrimeField, vyPrimeField, vzPrimeField };
	real* vPrimeContent[3] = { vxPrimeField->content, vyPrimeField->content, vzPrimeField->content };
	for (int d = 0; d < 3; ++d) {
		for (int k = (d == 2); k < resZ; ++k) {
			for (int j = (d == 1); j < resY; ++j) {
				for (int i = (d == 0); i < resX; ++i) {
					int vIndex = vPrimeField[d]->getIndex(i, j, k);
					int centerPlus = rhoGuessField->getIndex(i, j, k);
					int centerMinus = rhoGuessField->getIndex(i - (d == 0), j - (d == 1), k - (d == 2));
					vPrimeContent[d][vIndex] = (dt / h) * (rhsRhoStarStarField->content[centerPlus] -
						rhsRhoStarStarField->content[centerMinus]);
					vPrimeContent[d][vIndex] -= (dt / h) * (rhsRhoStarField->content[centerPlus] -
						rhsRhoStarField->content[centerMinus]);
					// The gravity term cancels.
				}
			}
		}
	}
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

void ThermalSolver::laplacianFieldOnAlignedGrid(Field* f, Field* lapF)
{
	real* fc = f->content;
	real* lapFc = lapF->content;
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
				lapFc[center] = (fc[left] + fc[right] +
					fc[up] + fc[bottom] + fc[front] + fc[back] - 
					6 * fc[center]) / h / h;
			}
		}
	}
}

void ThermalSolver::rhsRhoOnAlignedGrid(Field * rF, Field * rhsRF)
{
	real* rho = rF->content;
	real* rhsRho = rhsRF->content;
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
				rhsRho[center] = - isothermalWd(rho[center]);
				rhsRho[center] += (vdwInvWe / h / h) * (rho[left] + rho[right] +
					rho[up] + rho[bottom] + rho[front] + rho[back] -
					6 * rho[center]);
			}
		}
	}
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
	for (int z = 0; z < resZ; z++) {
		for (int y = 1; y < resY; y++) {
			for (int x = 0; x < resX; x++) {
				int index = vyBackgroundField->getIndex(x, y, z);
				real velx = (
					vxBackground[vxBackgroundField->getIndex(x, y - 1, z)] +
					vxBackground[vxBackgroundField->getIndex(x, y, z)] +
					vxBackground[vxBackgroundField->getIndex(x + 1, y - 1, z)] +
					vxBackground[vxBackgroundField->getIndex(x + 1, y, z)]) * 0.25f;
				real vely = vyBackground[index];
				real velz = (
					vzBackground[vzBackgroundField->getIndex(x, y, z)] +
					vzBackground[vzBackgroundField->getIndex(x, y, z + 1)] +
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
	Field * initVyField, Field * initVzField)
{
	resX = cfg->resX();
	resY = cfg->resY();
	resZ = cfg->resZ();
	h = cfg->h();

	vdwA = cfg->vdwPA();
	vdwB = cfg->vdwPB();
	vdwTheta = cfg->vdwTheta();
	vdwInvWe = 1.0f / cfg->vdwPWE();

	velConvergeTol = cfg->velConvergeTol();
	rhoConvergeTol = cfg->rhoConvergeTol();
	rhoRelaxCoef = cfg->rhoRelaxCoef();

	envGravity = cfg->gravity();

	if (initRhoField) rhoField = new Field(initRhoField);
	else rhoField = new Field(resX, resY, resZ, true);

	if (initVxField) vxField = new Field(initVxField);
	else vxField = new Field(resX + 1, resY, resZ, true);

	if (initVyField) vyField = new Field(initVyField);
	else vyField = new Field(resX, resY + 1, resZ, true);

	if (initVzField) vzField = new Field(initVzField);
	else vzField = new Field(resX, resY, resZ + 1, true);

	config = cfg;
}

void ThermalSolver::stepSimple(real dt)
{
	Field* vxGuessField = new Field(vxField);
	Field* vyGuessField = new Field(vyField);
	Field* vzGuessField = new Field(vzField);
	Field* rhoGuessField = new Field(rhoField);
	// For vStarField.
	Field* vxStarField = new Field(vxField, false);
	Field* vyStarField = new Field(vyField, false);
	Field* vzStarField = new Field(vzField, false);
	Field* rhoPrimeField = new Field(rhoField, false);
	// For rhoStarField.
	Field* rhsRhoStarField = new Field(rhoField, false);
	Field* rhsRhoStarStarField = new Field(rhoField, false);

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

		LOG(INFO) << "Rho diff = " << rhoDelta << ", Velocity diff = " << velDelta;
		if (rhoDelta < rhoConvergeTol && velDelta < velConvergeTol) {
			break;
		}
		else {
			if (rhoDelta > lastRhoDelta || velDelta > lastVelDelta) {
				LOG(WARNING) << "Converge Delta reaches its minimum. Stop iteration";
				break;
			}
			else {
				lastRhoDelta = rhoDelta;
				lastVelDelta = velDelta;
			}
		}
	}
	// Update current fields.
	rhoField->copyFrom(rhoGuessField);
	vxField->copyFrom(vxGuessField);
	vyField->copyFrom(vyGuessField);
	vzField->copyFrom(vzGuessField);
	delete vxGuessField;
	delete vyGuessField;
	delete vzGuessField;
	delete rhoGuessField;
	delete vxStarField;
	delete vyStarField;
	delete vzStarField;
	delete rhoPrimeField;
	delete rhsRhoStarField;
	delete rhsRhoStarStarField;
}

ThermalSolver::~ThermalSolver()
{
	delete rhoField;
	delete vxField;
	delete vyField;
	delete vzField;
}
