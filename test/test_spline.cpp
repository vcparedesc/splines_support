#include <iostream>
#include <stdio.h>
#include <gtest/gtest.h>
#include <yaml_eigen_utilities/yaml_eigen_utilities.hpp>
#include <roslib_utilities/ros_package_utilities.hpp>
#include "base_spline.hpp"
#include <vector>

using namespace std;
using namespace YAML;
using boost::filesystem::path;
using namespace yaml_utilities;
using namespace Eigen;
using namespace common;

TEST(SplineSupportTest, SplineConstructionTimeTest)
{
	VectorXd P1 = VectorXd::Zero(5);
	VectorXd P2 = VectorXd::Zero(5);
	VectorXd DP1 = VectorXd::Zero(5);
	VectorXd DP2 = VectorXd::Zero(5);
	VectorXd DDP1 = VectorXd::Zero(5);
	VectorXd DDP2 =  VectorXd::Zero(5);
	VectorXd outputs = VectorXd::Zero(5);

	P1 << 0.1008, 1.1933, 0.1542, 0.0488, -0.5269;
	P2 << 0.2992, 0.8425, 0.7477, -0.2467, -0.0007;
	DP1 << -0.0007, 0.4600, -0.7760, -0.0110, -0.1654;
	DP2 << -0.1891, -5.0513, -0.2383, 0.1498, -0.1654;
	DDP1 << 0.0194, -22.1782, -1.2521, -0.0465, 19.6795;
	DDP2 << 5.2354, -4.9979, -14.4895, 1.2320, -19.5341;

	SplineSupport TrajConnector(4,5,0.001,0.4,0.2);
	TrajConnector.addInequalityConstraint(3);
	TrajConnector.addBoundaryConditions(P1,P2,DP1,DP2,DDP1,DDP2);
	TrajConnector.solveSplines();
}

TEST(SplineSupportTest, BuildingTest)
{
	Node node;
	string filepath = roslib_utilities::resolve_local_url("package://splines_support/test/YamlTest/SplineModel_1.yaml").string();
	yaml_utilities::yaml_read_file(filepath, node);

	const Node &aA = node["aA"];
 	const Node &ab = node["ab"];
 	const Node &aC = node["aC"];
 	const Node &ad = node["ad"];
 	const Node &N = node["N"];
 	const Node &X = node["X"];
 	const Node &lambda = node["lambda"];

	VectorXd P1 = VectorXd::Zero(5);
	VectorXd P2 = VectorXd::Zero(5);
	VectorXd DP1 = VectorXd::Zero(5);
	VectorXd DP2 = VectorXd::Zero(5);
	VectorXd DDP1 = VectorXd::Zero(5);
	VectorXd DDP2 =  VectorXd::Zero(5);
	VectorXd outputs = VectorXd::Zero(5);

	P1 << 0.1008, 1.1933, 0.1542, 0.0488, -0.5269;
	P2 << 0.2992, 0.8425, 0.7477, -0.2467, -0.0007;
	DP1 << -0.0007, 0.4600, -0.7760, -0.0110, -0.1654;
	DP2 << -0.1891, -5.0513, -0.2383, 0.1498, -0.1654;
	DDP1 << 0.0194, -22.1782, -1.2521, -0.0465, 19.6795;
	DDP2 << 5.2354, -4.9979, -14.4895, 1.2320, -19.5341;

	SplineSupport TrajConnector(4,5,0.001,0.4,0.2);
	TrajConnector.addInequalityConstraint(3);
	TrajConnector.addBoundaryConditions(P1,P2,DP1,DP2,DDP1,DDP2);
	TrajConnector.solveSplines();

	MatrixXd mat;
	aA >> mat;
	mat = mat - TrajConnector.aA;
	std::cout<<"aA max:"<<mat.maxCoeff()<<std::endl;
	EXPECT_LE(abs(mat.mean()), 0.001);

	aC >> mat;
	mat = mat - TrajConnector.aC;
	std::cout<<"aC max:"<<mat.maxCoeff()<<std::endl;
	EXPECT_LE(abs(mat.mean()), 0.001);

	N >> mat;
    mat = mat - TrajConnector.N;
	std::cout<<"N max:"<<mat.maxCoeff()<<std::endl;
    EXPECT_LE(abs(mat.mean()), 0.001);

	MatrixXd vec;

	ab >> vec;
	vec = vec - TrajConnector.ab;
	std::cout<<"ab max:"<<vec.maxCoeff()<<std::endl;
	EXPECT_LE(abs(vec.mean()), 0.001);

	ad >> vec;
	vec = vec - TrajConnector.ad;
	std::cout<<"ad max:"<<vec.maxCoeff()<<std::endl;
	EXPECT_LE(abs(vec.mean()), 0.001);

	lambda >> vec;
	vec = vec - TrajConnector.lambda;
	std::cout<<"lambda max:"<<vec.maxCoeff()<<std::endl;
	EXPECT_LE(abs(vec.mean()), 0.001);

}

TEST(SplineSupportTest, TrajectoryIntegrity)
{
    VectorXd P1 = VectorXd::Zero(1);
    VectorXd P2 = VectorXd::Zero(1);
    VectorXd DP1 = VectorXd::Zero(1);
    VectorXd DP2 = VectorXd::Zero(1);
    VectorXd DDP1 = VectorXd::Zero(1);
    VectorXd DDP2 =  VectorXd::Zero(1);
    VectorXd outputs = VectorXd::Zero(1);

    P1 << 0.6109;
    P2 << 0.1008;
    DP1 << 0;
    DP2 << 0;
    DDP1 << 0;
    DDP2 << 0;

    SplineSupport TrajConnector(4,1,0.001,0,0.4212);
    TrajConnector.addInequalityConstraint(3);
    TrajConnector.addBoundaryConditions(P1,P2,DP1,DP2,DDP1,DDP2);
    TrajConnector.solveSplines();

    std::cout<<"OUTPUT: "<<std::endl;

    for(int i = 1; i < 100; i ++)
    {
        outputs << TrajConnector.computeOutput(i / 99.0 * 0.4212);
        std::cout<<outputs<<std::endl;
    }
}


int main(int argc, char ** argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
