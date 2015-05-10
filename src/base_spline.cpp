#include "base_spline.hpp"

SplineSupport::SplineSupport(int n_splines, int n_outputs, int optimization_type): nSplines(n_splines), nOutputs(n_outputs), typeOptimization(optimization_type)
{
	double Dti;
	uNormal = 0.001;
	tConvergence = 0.2;
	MatrixXd tA = MatrixXd::Zero(nSplines + 1, nSplines * 4);

	VectorXd tb = VectorXd::Zero(nSplines + 1);

	MatrixXd tPos = MatrixXd::Zero((nSplines + 1) , nSplines * 4);

	MatrixXd tVel = MatrixXd::Zero((nSplines + 1) , nSplines * 4);

	MatrixXd tAcc = MatrixXd::Zero((nSplines + 1) , nSplines * 4);

	MatrixXd tC = MatrixXd::Zero((nSplines + 1) * 3, nSplines * 4);

	VectorXd td = VectorXd::Zero((nSplines + 1) * 3);

	MatrixXd tF = MatrixXd::Zero(nSplines + 1, nSplines * 4);

	VectorXd tg = VectorXd::Zero(nSplines + 1);

	// Allocates space for the matrices
	this->lambda.resize(nSplines * nOutputs * 4);

	this->A = MatrixXd::Zero((nSplines + 1) * nOutputs,nSplines * nOutputs * 4);
	this->b = VectorXd::Zero((nSplines + 1) * nOutputs);
	this->C = MatrixXd::Zero((nSplines + 1) * nOutputs * 3, nSplines * nOutputs * 4);
	this->d = VectorXd::Zero((nSplines + 1) * nOutputs * 3);
	this->F = MatrixXd::Zero((nSplines + 1) * nOutputs, nSplines * nOutputs * 4);
	this->g = VectorXd::Zero((nSplines + 1) * nOutputs);
	this->aA = MatrixXd::Zero((nSplines + 1) * nOutputs * 2, nSplines * nOutputs * 4 + (nSplines + 1) * nOutputs);
	this->ab = VectorXd::Zero(nSplines * nOutputs * 4 + (nSplines + 1) * nOutputs);
	this->aC = MatrixXd::Zero((nSplines + 1) * nOutputs * 4, nSplines *nOutputs * 4 + (nSplines + 1) * nOutputs);
	this->ad = VectorXd::Zero((nSplines + 1) * nOutputs * 3 + (nSplines + 1)*nOutputs);

	this->N = MatrixXd::Zero(nSplines * nOutputs * 4 + (nSplines + 1) * nOutputs * 5, nSplines * nOutputs * 4 + (nSplines + 1) * nOutputs * 5);
	this->X = VectorXd::Zero(nSplines * nOutputs * 4 + (nSplines + 1) * nOutputs);

	std::cout<<"Building Basal Spline Matrices"<<std::endl;

	// Build the matrices and store them in memory for future use!
	for (int i = 0; i < (n_splines+1); i++)
	{
		// Time used for this point
		Dti = (double)(tConvergence * (double)i / nSplines);

		// Matrix tA
		if(i == 0)
		{
			tA.block<1,4>(i,i) << 0,0,0,0;
		}else if(i == nSplines)
		{
			tA.block<1,4>(i, 4*(i-1)) << 0,0,0,0;
		}else
		{
			tA.block<1,4>(i, 4*(i-1)) << 1,Dti , pow(Dti,2), pow(Dti,3);
		}

		// Matrix tPos
		if(i == 0)
		{
			tPos.block<1,4>(i,i) << 1,0,0,0;
		}else if(i == nSplines)
		{
			tPos.block<1,4>(i, 4*(i-1)) << 1, tConvergence, pow(tConvergence,2), pow(tConvergence,3);
		}else
		{
			tPos.block<1,8>(i, 4*(i-1)) << 1, Dti, pow(Dti,2), pow(Dti,3), -1, -Dti, -pow(Dti,2), -pow(Dti,3);
		}

		// Matrix tVel
		if(i == 0)
		{
			tVel.block<1,4>(i,i) << 0,1,0,0;
		}else if(i == nSplines)
		{
			tVel.block<1,4>(i, 4*(i-1)) << 0, 1, 2 * tConvergence, 3 * pow(tConvergence,2);
		}else
		{
			tVel.block<1,8>(i, 4*(i-1)) << 0, 1, 2 * Dti, 3 * pow(Dti,2), 0, -1, -2 * Dti, -3 * pow(Dti,2);
		}


		// Matrix tAcc
		if(i == 0)
		{
			tAcc.block<1,4>(i,i) << 0,0,2,0;
		}else if(i == nSplines)
		{
			tAcc.block<1,4>(i, 4*(i-1)) << 0, 0, 2, 6 * Dti;
		}else
		{
			tAcc.block<1,8>(i, 4*(i-1)) << 0, 0, 2, 6 * Dti, 0, 0, -2, -6 * Dti;
		}

		// Velocity as Inequality constraint tFx < tg
		if(i == 0)
		{
			tF.block<1,4>(i,i) << 0,1,2 * Dti, 3 * pow(Dti,2);
		}else
		{
			tF.block<1,4>(i, 4*(i-1)) << 0, 1, 2 * Dti, 3 * pow(Dti,2);
		}
	}

	std::cout<<"Obtaining Block Diagonal Matrices"<<std::endl;
	A = nBlockDiag(tA,nOutputs);
	C.block((nSplines + 1) * nOutputs * 0,0, (nSplines + 1) * nOutputs, nSplines * nOutputs * 4 ) = nBlockDiag(tPos,5);
	C.block((nSplines + 1) * nOutputs * 1,0, (nSplines + 1) * nOutputs, nSplines * nOutputs * 4 ) = nBlockDiag(tVel,5);
	C.block((nSplines + 1) * nOutputs * 2,0, (nSplines + 1) * nOutputs, nSplines * nOutputs * 4 ) = nBlockDiag(tAcc,5);

	F = nBlockDiag(tF, nOutputs);

	std::cout<<"Augmented Matrix with slack variables"<<std::endl;
	/* Augmented Matrix due slack variables */
	aA.block(0,0,(nSplines + 1) * nOutputs, nSplines * nOutputs * 4) = A;
	aA.block(0, nSplines * nOutputs * 4, (nSplines + 1) * nOutputs, (nSplines + 1) * nOutputs) = MatrixXd::Zero((nSplines + 1) * nOutputs,(nSplines + 1) * nOutputs);
	aA.block( (nSplines + 1) * nOutputs, 0, (nSplines + 1) * nOutputs, nSplines * nOutputs * 4) = MatrixXd::Zero((nSplines + 1) * nOutputs, nSplines * nOutputs * 4);
	aA.block( (nSplines + 1) * nOutputs, nSplines * nOutputs * 4, (nSplines + 1) * nOutputs, (nSplines + 1) * nOutputs) = uNormal * MatrixXd::Identity((nSplines + 1) * nOutputs, (nSplines + 1) * nOutputs);
	std::cout<<"Matrix C"<<std::endl;
	aC.block(0,0, (nSplines + 1) * nOutputs * 3, nSplines * nOutputs * 4) = C;
	aC.block(0, nSplines * nOutputs * 4, (nSplines + 1) * nOutputs, (nSplines + 1) * nOutputs) = MatrixXd::Zero((nSplines + 1) * nOutputs, (nSplines + 1) * nOutputs);
	aC.block((nSplines + 1) * nOutputs * 3, 0, (nSplines + 1) * nOutputs, nSplines * nOutputs * 4) = F;
	aC.block((nSplines + 1) * nOutputs * 3, nSplines * nOutputs * 4, (nSplines + 1) * nOutputs, (nSplines + 1) * nOutputs) = MatrixXd::Identity((nSplines + 1) * nOutputs, (nSplines + 1) * nOutputs);
	/* Now we have || aA X - ab || st aC X = ad */
	std::cout<<"Building Normal Matrix"<<std::endl;
	/* Building Normal Matrix */
	N.block(0,0, aA.cols(), aA.cols()) = 2 * aA.transpose() * aA;
	N.block(0, aA.cols(), aC.cols(), aC.rows()) = aC.transpose();
	std::cout<<"here"<<std::endl;
	N.block(aA.cols(), 0, aC.rows(), aC.cols()) = aC;
	N.block(aA.cols(), aC.rows(), aC.rows(), aC.rows()) = MatrixXd::Zero(aC.rows(), aC.rows());

	std::cout<<N<<std::endl;

}

SplineSupport::~SplineSupport(void)
{
	/* Freeing Memory! */
	A.resize(0,0);
	b.resize(0);
	C.resize(0,0);
	d.resize(0);
	F.resize(0,0);
	g.resize(0);
	aA.resize(0,0);
	ab.resize(0);
	aC.resize(0,0);
	ad.resize(0);
	N.resize(0,0);
	X.resize(0);
}

MatrixXd SplineSupport::nBlockDiag(MatrixXd basal_matrix, int n_repetitions)
{
	MatrixXd nBlockMatrix = MatrixXd::Zero(basal_matrix.rows() * n_repetitions, basal_matrix.cols() * n_repetitions);

	for (int i = 0; i < n_repetitions; i++)
	{
		nBlockMatrix.block(basal_matrix.rows() * i, basal_matrix.cols() * i, basal_matrix.rows(), basal_matrix.cols()) = basal_matrix;
	}

	return nBlockMatrix;
}

void SplineSupport::addBoundaryConditions(VectorXd P1, VectorXd P2, VectorXd DP1, VectorXd DP2, VectorXd DDP1, VectorXd DDP2)
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

	// Make Sure the index spline is in the correct range
	if(indexSpline >= nSplines)
		indexSpline = nSplines;


	for (int i = 0; i < nOutputs; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			// Calculate reduced lambda: Lambda set at time_parameter
			this->splineOutput(i) =  pow(time_parameter,j) * this->lambda(4*nSplines*i + 4*indexSpline + j);
		}
	}

	return this->splineOutput;
}
