#ifndef FINITEDIFFERENCEHEADERDEF
#define FINITEDIFFERENCEHEADERDEF
#include <Eigen/Core>

class FiniteDifferenceGrid
{
public:
	friend class BvpOde;
private: 
	Eigen::VectorXd mNodes;
public: 
	FiniteDifferenceGrid(int numNodes, double xMin, double xMax);
};
#endif
