#include <eigen3/Eigen/Eigen>
#include <math.h>

using namespace Eigen;

class SplineSupport
{
private:
	SplineSupport(int n_splines, int n_outputs);
	~SplineSupport(void);

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
	void addEqualityConstraint();
	void addInequalityConstraint();
	void setNormalizer();
	VectorXd computeOutput(double time_parameter);

};
