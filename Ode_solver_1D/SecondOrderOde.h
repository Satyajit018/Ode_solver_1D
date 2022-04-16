#ifndef SECONDORDERODEHEADERDEF
#define SECONDORDERODEHEADERDEF
class SecondOrderOde
{
public:
	friend class BvpOde;
private: 
	//coeff of the ODE
	double mCoeffUxx;
	double mCoeffUx;
	double mCoeffU;

	//RHS function of the ODE
	double (*mpRhsFunc) (double x);

	//Interval for domain

	double mXmin;
	double mXmax;

public:
 SecondOrderOde(double coeffUxx, double coeffUx, double coeffU, double(*rhsFunc) (double), double xMin, double xMax);

};

#endif