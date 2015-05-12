#include <eigen3/Eigen/Eigen>
#include <math.h>
#include <iostream>

using namespace Eigen;

class SplineSupport
{
private:
	MatrixXd nBlockDiag(MatrixXd basal_matrix, int n_repetitions);
	void buildNormalVector();


	int typeOptimization; // 0: Equality Constraint, 1: Equality and Inequality Constraint
	int readyFlag; // Indicates that all conditions are set
	int nSplines;
	int nOutputs;
	/* Represents Objective Function: || Ax - b || */
	MatrixXd A;
	VectorXd b;
	/* Represents Equality Constraint: Cx = d */
	MatrixXd C;
	VectorXd d;
	/* Represents Inequality Constraint: Fx < g */
	MatrixXd F;
	VectorXd g;
	/* Slack variable normalizer: u * || sx || */
	double uNormal;

	/* Velocity Limit */
	double vLimit;

	/* Augmented Matrices considering slack variables */
	MatrixXd aA;
	VectorXd ab;
	MatrixXd aC;
	VectorXd ad;

	/* First Optimality Condition: N x* = X --> x* = N^-1 * X */
	MatrixXd N;
	VectorXd X;

	/* Time-parameter at which we use the splines */
	double tSwitch;

	/* Additional Time-parameter at it converges */
	double tConvergence;

	/* Vector of Cubic-Spline parameters */
	VectorXd lambda;

	/* Useful Interprocess Matrices (to save computation time) */
	MatrixXd timeBlock;
	VectorXd splineOutput;
	VectorXd reducedLambda;
public:
	SplineSupport(int n_splines, int n_outputs, double normalizer, double time_switch, double time_convergence);
	~SplineSupport(void);
	void addBoundaryConditions(VectorXd P1, VectorXd P2, VectorXd DP1, VectorXd DP2, VectorXd DDP1, VectorXd DDP2);
	void addInequalityConstraint(double limit);
	void setNormalizer(double normalizer);
	VectorXd referenceTraj(int time_parameter);
	VectorXd computeOutput(double time_parameter);
	void buildReferenceVector();

};
