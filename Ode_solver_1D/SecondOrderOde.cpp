#include "SecondOrderOde.h"
 SecondOrderOde::SecondOrderOde(double coeffUxx, double coeffUx, double coeffU, double(*rhsFunc) (double), double xMin, double xMax)
{
	 mCoeffUxx = coeffUxx;
	 mCoeffUx = coeffUx;
	 mCoeffU = coeffU;
	 mpRhsFunc = rhsFunc;
	 mXmin = xMin;
	 mXmax = xMax;

}