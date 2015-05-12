#include <iostream>
#include "base_spline.hpp"

int main(int argc, char ** argv)
{
	VectorXd P1 = VectorXd::Zero(5);
	VectorXd P2 = VectorXd::Zero(5);
	VectorXd DP1 = VectorXd::Zero(5);
	VectorXd DP2 = VectorXd::Zero(5);
	VectorXd DDP1 = VectorXd::Zero(5);
	VectorXd DDP2 =  VectorXd::Zero(5);
	P1 << 0.1008, 1.1933, 0.1542, 0.0488, -0.5269;
	P2 << 0.2992, 0.8425, 0.7477, -0.2467, -0.0007;
	DP1 << -0.0007, 0.4600, -0.7760, -0.110, -0.1654;
	DP2 << -0.1891, -5.0513, -0.2383, 0.1498, -0.1654;
	DDP1 << 0.0194, -22.1782, -1.2521, -0.0465, 19.6795;
	DDP2 << 5.2354, -4.9979, -14.4895, 1.2320, -19.5341;

	SplineSupport TrajConnector(4,5,0.001,0.4,0.2);
	TrajConnector.addInequalityConstraint(6);
	TrajConnector.addBoundaryConditions(P1,P2,DP1,DP2,DDP1,DDP2);
	TrajConnector.buildReferenceVector();

	return 0;
}