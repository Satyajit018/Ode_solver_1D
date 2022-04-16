#ifndef BVPODEHEADERDEF
#define BVPODEHEADERDEF
#include <Eigen/Core>
#include<Eigen/Sparse>
#include "FiniteDifferenceGrid.h"
#include "SecondOrderOde.h"
#include "BoundaryConditions.h"
#include <string>
#include<Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>
typedef Eigen::SparseMatrix<double> SparseMatrixXd;
class BvpOde
{
private:
	BvpOde(const BvpOde& otherBvpOde) {}

	int mNumNodes;
	FiniteDifferenceGrid* mpGrid;
	//pointer to an instance of an ODE
	SecondOrderOde* mpOde;
    // pointer to an instance of BC
	BoundaryConditions* mpBconds;

	// pointer to solution vector
	Eigen::VectorXd mpSolVec;

	//Right hand side vector

	Eigen::VectorXd mpRhsVec;

	//matrix assembly
	SparseMatrixXd mpLhsMat;

	//Output file name
	std::string mFileName;

	//methods for solving

	void PopulateMatrix();
	void PopulateVector();
	void ApplyBoundaryConditions();

public:
	//constructor

	BvpOde(SecondOrderOde* pOde, BoundaryConditions* pBcs, int numNodes);

	//destructor

	~BvpOde();
	void SetFileName(const std::string& name)
	{
		mFileName = name;
	}
	void Solve();
	void WriteSolutionFile();
};

#endif
