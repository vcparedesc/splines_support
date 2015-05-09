#include "base_spline.hpp"

SplineSupport::SplineSupport(int n_splines, int n_outputs): nSplines(n_splines), nOutputs(n_outputs)
{
	
}

void SplineSupport::addEqualityConstraint()
{

}

void SplineSupport::addInequalityConstraint()
{

}

void SplineSupport::setNormalizer()
{

}

VectorXd SplineSupport::computeOutput(double time_parameter)
{
	double timeInterval = this->tConvergence;
	double indexSpline = floor(time_parameter / timeInterval);
}
