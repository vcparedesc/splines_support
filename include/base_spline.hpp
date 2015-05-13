#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <math.h>
#include <iostream>

using namespace Eigen;

class SplineSupport
{
private:
	MatrixXd nBlockDiag(MatrixXd basal_matrix, int n_repetitions);
	void buildReferenceVector();
	void buildNormalVector();

public:


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
	FullPivLU<MatrixXd> solver;
	//JacobiSVD<MatrixXd> solver;

	/* Bounday Condition Vectors */
	VectorXd _P1, _P2, _DP1, _DP2, _DDP1, _DDP2;

	/* Time-parameter at which we use the splines */
	double tSwitch;

	/* Additional Time-parameter at it converges */
	double tConvergence;

	/* Vector of Cubic-Spline parameters */
	VectorXd lambda;

	/* Useful Interprocess Matrices (to save computation time) */
	MatrixXd timeBlock;
	VectorXd reducedLambda;

	SplineSupport(int n_splines, int n_outputs, double normalizer, double time_switch, double time_convergence);
	~SplineSupport(void);
	void addBoundaryConditions(VectorXd P1, VectorXd P2, VectorXd DP1, VectorXd DP2, VectorXd DDP1, VectorXd DDP2);
	void addInequalityConstraint(double limit);
	void setNormalizer(double normalizer);
	void solveSplines();
	virtual VectorXd referenceTraj(double time_parameter);
	VectorXd computeOutput(double time_parameter);


};
