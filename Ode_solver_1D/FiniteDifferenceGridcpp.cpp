#include<cassert>
#include<iostream>
#include "FiniteDifferenceGrid.h"
#include <Eigen/Core>
FiniteDifferenceGrid::FiniteDifferenceGrid(int numNodes, double xMin, double xMax)
{
	double stepsize = (xMax-xMin) / double(numNodes-1);
	Eigen::VectorXd nodes(numNodes); //location of the FV center
	for (int i = 0; i < numNodes; i++)
	{
		nodes(i) = xMin + i * stepsize;
		mNodes = nodes;
	}
	assert(mNodes.size() == numNodes);
	//std::cout << nodes << std::endl;
}