#include "BvpOde.h"
#include<iostream>
#include<fstream>
#include<cassert>
#include<Eigen/IterativeLinearSolvers>
using namespace Eigen;
typedef Eigen::SparseMatrix<double> SparseMatrixXd;
BvpOde::BvpOde(SecondOrderOde* pOde, BoundaryConditions* pBcs, int numNodes)
{
	mpOde = pOde;
	mpBconds = pBcs;
	mNumNodes = numNodes;
	mpGrid = new FiniteDifferenceGrid(mNumNodes, pOde->mXmin, pOde->mXmax);
	mpSolVec = Eigen::VectorXd(mNumNodes);
	mpRhsVec = Eigen::VectorXd(mNumNodes);
	//mpLhsMat = new SparseMatrixXd(mNumNodes, mNumNodes);
	mpLhsMat = SparseMatrixXd(mNumNodes, mNumNodes);
	mFileName = "ode_output.dat";
}

BvpOde::~BvpOde()
{
	//delete mpSolVec;
	//delete mpRhsVec;
	//delete mpLhsMat;
	delete mpGrid;
}

void BvpOde::Solve()
{
	PopulateMatrix();
	PopulateVector();
	ApplyBoundaryConditions();
	BiCGSTAB<SparseMatrix<double> > solver;
	solver.compute(mpLhsMat);
	mpSolVec = solver.solve(mpRhsVec);
	std::cout << "#iterations:     " << solver.iterations() << std::endl;
	std::cout << "estimated error: " << solver.error() << std::endl;
	/* ... update b ... */
	mpSolVec = solver.solve(mpRhsVec); // solve again
	//std::cout << mpSolVec << std::endl;
	WriteSolutionFile();
}

void BvpOde::PopulateMatrix()
{
	mpLhsMat.reserve(VectorXd::Constant(mNumNodes, 3));
	for (int i = 1; i < mNumNodes - 1; i++) 
	{
		double xm = mpGrid->mNodes(i - 1);
		double x = mpGrid->mNodes(i);
		double xp = mpGrid->mNodes(i + 1);
		double alpha = 2.0 / (xp - xm) / (x - xm);
		double beta = -2.0 / (xp - x) / (x - xm);
		double gamma = 2.0 / (xp - xm) / (xp - x);
		mpLhsMat.coeffRef(i , i-1) = (mpOde->mCoeffUxx) * alpha - (mpOde->mCoeffUx) / (xp - xm);
		mpLhsMat.coeffRef(i , i ) = (mpOde->mCoeffUxx) * beta + (mpOde->mCoeffU);
		mpLhsMat.coeffRef(i , i+1) = (mpOde->mCoeffUxx) * gamma + (mpOde->mCoeffUx) / (xp - xm);
	}
	
}

void BvpOde::PopulateVector()
{
	for (int i = 1; i < mNumNodes - 1; i++)
	{
		double x = mpGrid->mNodes(i);
		 mpRhsVec(i) = mpOde->mpRhsFunc(x);

	}
}

void BvpOde::ApplyBoundaryConditions()
{
	bool left_bc_applied = false;
	bool right_bc_applied = false;
	if (mpBconds->mLhsBcIsDirichlet) {
		mpLhsMat.coeffRef(0, 0) = 1.0;
		mpRhsVec(0) = mpBconds->mLhsBcValue;
		left_bc_applied = true;

	}
	if (mpBconds->mRhsBcIsDirichlet) {
		mpLhsMat.coeffRef(mNumNodes-1, mNumNodes-1) = 1.0;
		mpRhsVec(mNumNodes - 1) = mpBconds->mRhsBcValue;
		right_bc_applied = true;
		std::cout << mpBconds->mRhsBcValue << std::endl;
	}
	if (mpBconds->mLhsBcIsNeumann) {
		assert(left_bc_applied == false);
		double h = mpGrid->mNodes(1) - mpGrid->mNodes(0);
		mpLhsMat.coeffRef(0, 0) = -1.0 / h;
		mpLhsMat.coeffRef(0, 1) = 1.0 / h;
		mpRhsVec(0) = mpBconds->mLhsBcValue;
		left_bc_applied = true;
	}

	if (mpBconds->mRhsBcIsNeumann) {
		assert(right_bc_applied == false);
		double h = mpGrid->mNodes(mNumNodes-2) - mpGrid->mNodes(mNumNodes-1);
		mpLhsMat.coeffRef(mNumNodes - 2, mNumNodes - 1) = -1.0 / h;
		mpLhsMat.coeffRef(mNumNodes - 1, mNumNodes - 1) = 1.0 / h;
		mpRhsVec(mNumNodes - 1) = mpBconds->mRhsBcValue;
		left_bc_applied = true;
	}
	//std::cout << mpLhsMat << std::endl;
	//std::cout << mpRhsVec << std::endl;
}

void BvpOde::WriteSolutionFile()
{
	std::ofstream output_file(mFileName.c_str());
	assert(output_file.is_open());
	for (int i = 0; i < mNumNodes; i++)
	{
		output_file << mpGrid->mNodes(i) << " " << mpSolVec(i) << "\n";
	}
	output_file.close();
}