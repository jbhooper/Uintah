#include "testmatrix3.h"
#include <Uintah/Components/MPM/Util/Matrix3.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <iostream>

using namespace std;

void doMatrixSolvingTests(Suite& suite);
void doEigenTests(Suite& suite);
void doEigenPlaneTests(Suite& suite);

void addSolveTests(Suite& suite, const string& test_name,
	       const Matrix3& M, const Vector& rhs, bool exp_return,
	       int exp_xg_basis_size);

//bool testEigenValue(const Matrix3& M, double eigen_value);
bool testEigenValue(const Matrix3& M, double eigen_value, double max_eigen_value);

void randomlyMixup(Matrix3& M, Vector& rhs);

bool equal_enough(double x1, double x2, double rel_scale);
bool equal_enough(Vector v1, Vector v2, double rel_scale);

double getRandom(); // returns a random number between -100.000 and 100.000

// returns 0 with the given probability or a random number
// between -100.000 and 100.000 otherwise.
double getRandomOrZero(float probability_for_zero); 

Vector randomVector();
Matrix3 randomMatrix();

#define NEAR_ZERO 1e-7

void displayEigen(Matrix3 M);

SuiteTree* matrix3TestTree()
{
  srand(getpid());  
  SuiteTreeNode* matrix3Tests = new SuiteTreeNode("Matrix3");
    
  Suite* solvingTests = new Suite("Solving Ax=b");
  Suite* eigenTests = new Suite("Eigen values/vectors");
  Suite* eigenPlaneTests = new Suite("Eigen plane values");
  doMatrixSolvingTests(*solvingTests);

  doEigenTests(*eigenTests);
  doEigenPlaneTests(*eigenPlaneTests);
  matrix3Tests->addSuite(solvingTests);
  matrix3Tests->addSuite(eigenTests);
  matrix3Tests->addSuite(eigenPlaneTests);

  //Matrix3 matrix(1e-10, 0, 0, 0, 1e-10, 0, 0, 0, 1e-10);
  //Matrix3 matrix(1e6, 0, 0, 0, 1e6, 0, 0, 0, 1e6);

  /* Was doing some testing here
  Matrix3 M = randomMatrix();
  displayEigen(M);
  displayEigen(M*.0001);
  */

  return matrix3Tests;

}

void displayEigen(Matrix3 M)
{
  double e[3];
  int n = M.getEigenValues(e[0], e[1], e[2]);
  for (int i = 0; i < n; i++) {
    cout << "Eigen value: " << e[i] << endl;
    std::vector<Vector> eigenVectors = M.getEigenVectors(e[i], M.MaxAbsElem());
    cout << "Eigen vectors:\n";
    for (int j = 0; j < eigenVectors.size(); j++)
      cout << eigenVectors[j] << endl;
  }
  cout << endl;
}

void doMatrixSolvingTests(Suite& suite)
{
  Vector rhs;
  Vector xp;
  std::vector<Vector> xg_basis;
  Matrix3 M;

  int i = 0;
  string test_name_base = "Point: {{1, 0, 0} {0, 1, 0} {0, 0, 1}}";
  bool has_solution = true;
  double a, b, c, e;

  // run 10 random tests of matrix equations of one point solution
  for (i = 0; i < 10; i++) {
    M.Identity();
    rhs = randomVector();
    randomlyMixup(M, rhs); // this is what will really test it
    addSolveTests(suite, test_name_base, M, rhs, has_solution, 0);
  }
  
  // run 10 random tests for matrices that reduce to the following type
  test_name_base = "Line: {{1, 0, c} {0, 1, e} {0, 0, 0}}";
  xg_basis.resize(1);
  for (i = 0; i < 10; i++) {
    M.Identity();
    c = getRandom();
    e = getRandom();
    M(3, 3) = 0; M(1, 3) = c; M(2, 3) = e;
    rhs = randomVector();
    rhs.z(getRandomOrZero(0.75)); // make rhs.z zero 75% of the time
    has_solution = (rhs.z() == 0);

    randomlyMixup(M, rhs); // this is what will really test it
    addSolveTests(suite, test_name_base, M, rhs, has_solution, 1);
  }

   // run 10 random tests for matrices that reduce to the following type
  test_name_base = "Line: {{1, b, 0} {0, 0, 1} {0, 0, 0}}";
  for (i = 0; i < 10; i++) {
    M.Identity();
    b = getRandom();
    M(3, 3) = 0; M(2, 2) = 0; M(2, 3) = 1; M(1, 2) = b;
    rhs = randomVector();
    rhs.z(getRandomOrZero(0.75)); // make rhs.z zero 75% of the time
    has_solution = (rhs.z() == 0);

    randomlyMixup(M, rhs); // this is what will really test it
    addSolveTests(suite, test_name_base, M, rhs, has_solution, 1);
  }

  test_name_base = "Line: {{0, 1, 0} {0, 0, 1} {0, 0, 0}}";
  for (i = 0; i < 10; i++) {
    M.Identity();
    b = getRandom();
    M(1, 1) = 0; M(1, 2) = 1; M(2, 2) = 0; M(2, 3) = 1; M(3, 3) = 0;
    rhs = randomVector();
    rhs.z(getRandomOrZero(0.75)); // make rhs.z zero 75% of the time
    has_solution = (rhs.z() == 0);

    randomlyMixup(M, rhs); // this is what will really test it
    addSolveTests(suite, test_name_base, M, rhs, has_solution, 1);
  }
  test_name_base = "Plane: {{a, b, c} {0, 0, 0} {0, 0, 0}}";
  for (i = 0; i < 100; i++) {
    M.Identity();
    a = getRandomOrZero(0.5);
    b = getRandomOrZero(0.5);
    c = (a == 0 && b == 0) ? 1 : getRandomOrZero(0.5);
    M(1, 1) = a; M(1, 2) = b; M(1, 3) = c; M(2, 2) = 0; M(3, 3) = 0;
    rhs = randomVector();
    rhs.y(getRandomOrZero(0.80)); // make rhs.y zero 80% of the time
    rhs.z(getRandomOrZero(0.80)); // make rhs.z zero 80% of the time
    has_solution = ((rhs.z() == 0) && (rhs.y() == 0));

    randomlyMixup(M, rhs); // this is what will really test it
    addSolveTests(suite, test_name_base, M, rhs, has_solution, 2);
  }
  test_name_base = "All Space: {{0, 0, 0} {0, 0, 0} {0, 0, 0}}";
  M = Matrix3();
  rhs = Vector(0, 0, 0);
  addSolveTests(suite, test_name_base, M, rhs, true, 3);

  for (i = 0; i < 10; i++) {
    rhs = randomVector();
    addSolveTests(suite, test_name_base, M, rhs, (rhs == Vector(0, 0, 0)), 3);
  }
}

void doEigenTests(Suite& suite)
{  
  Test* eigenValueOrderTest = suite.addTest("Eigenvalue e1, e2, e3 order");
  Test* threeEigenTestA = suite.addTest("Three eigenvalues/vectors, e1");
  Test* threeEigenTestB = suite.addTest("Three eigenvalues/vectors, e2");
  Test* threeEigenTestC = suite.addTest("Three eigenvalues/vectors, e3");
  Test* oneEigenTest = suite.addTest("One eigenvalue/vector");
  Matrix3 M;
  double e1, e2, e3;
  int num_eigen_values;
  
  for (int i = 0; i < 100; i++) {
    // run a bunch of random tests
    M = randomMatrix();
    num_eigen_values = M.getEigenValues(e1, e2, e3);
    if (num_eigen_values == 1) {
      oneEigenTest->setResults(testEigenValue(M, e1, e1));
    }
    else if (num_eigen_values == 3) {
      threeEigenTestA->setResults(testEigenValue(M, e1, e1));
      threeEigenTestB->setResults(testEigenValue(M, e2, e1));
      threeEigenTestC->setResults(testEigenValue(M, e3, e1));
      eigenValueOrderTest->setResults((e1 > e2) && (e2 > e3));
    }
    else if (num_eigen_values == 2) {
      // two eigen values are the same
      // just treat this under the three eigen test case since
      // this is rare on random data
      threeEigenTestA->setResults(testEigenValue(M, e1, e1));
      threeEigenTestB->setResults(testEigenValue(M, e2, e1));
      eigenValueOrderTest->setResults(e1 > e2);    
    }
    else {
      suite.addTest("Bad number of eigen values", false);
    }
  }
  
}

void doEigenPlaneTests(Suite& suite)
{  
  Test* eigenValueOrderTest = suite.addTest("e1 > e2");
  Test* eigenTestA[3];
  Test* eigenTestB[3];

  eigenTestA[0] = suite.addTest("YZ e1");
  eigenTestB[0] = suite.addTest("YZ e2");
  eigenTestA[1] = suite.addTest("XZ e1");
  eigenTestB[1] = suite.addTest("XZ e2");
  eigenTestA[2] = suite.addTest("XY e1");
  eigenTestB[2] = suite.addTest("XY e2");
  
  Matrix3 M;
  Matrix3 planeM; // 3d version of sub-matrix (with zeroes at bottom and right)
  double e1, e2;
  int num_eigen_values;
  
  for (int i = 0; i < 5; i++) {
    // run a bunch of random tests
    M = randomMatrix();
    for (int plane = 1; plane <= 3; plane++) {
      if (plane == 1) {
	num_eigen_values = M.getYZEigenValues(e1, e2);
	planeM = Matrix3(M(2, 2), M(2, 3), 0, M(3, 2), M(3, 3), 0, 0, 0, 0);
      }
      else if (plane == 2) {
	num_eigen_values = M.getXZEigenValues(e1, e2);
	planeM = Matrix3(M(1, 1), M(1, 3), 0, M(3, 1), M(3, 3), 0, 0, 0, 0);
      }
      else {
 	num_eigen_values = M.getXYEigenValues(e1, e2);
	planeM = Matrix3(M(1, 1), M(1, 2), 0, M(2, 1), M(2, 2), 0, 0, 0, 0);
      }

      // Use the 3x3 test with planeM (which is the sub-matrix
      // with 0's on the bottom and right).  This works because
      // this every eigenvalue of the 2x2 sub-matrix should be
      // an eigenvalue of planeM.
      
      if (num_eigen_values == 1) {
	// two eigen values are the same
	eigenTestA[plane-1]->setResults(testEigenValue(planeM, e1, e1));
      }
      else if (num_eigen_values == 2) {
	eigenTestA[plane-1]->setResults(testEigenValue(planeM, e1, e1));
	eigenTestB[plane-1]->setResults(testEigenValue(planeM, e2, e1));
	eigenValueOrderTest->setResults(e2 < e1);    
      }
      else if (num_eigen_values != 0) {
	suite.addTest("Bad number of eigen values", false);
      }
    }
  }

}

bool testEigenValue(const Matrix3& M, double eigen_value, double max_eigen_value)
{
  double rel_scale = fabs(max_eigen_value); //M.MaxAbsElem(); //M.Norm();
  std::vector<Vector> eigenVectors = M.getEigenVectors(eigen_value, rel_scale);

  if (eigenVectors.size() == 0) {
    cout << "oh shit!\n";
    return false; // should have at least one Vector for an eigen_value
  }
  else if (eigenVectors.size() == 1) {
    // only a line of eigenVectors -- just test the one vector
    return equal_enough(M*eigenVectors[0], eigenVectors[0]*eigen_value,
			rel_scale);
  }
  else {
    // try different random linear combinations of eigenVectors
    bool success = true;
    for (int trials = 0; trials < 10; trials++) {
      Vector x(0, 0, 0);
      for (int i = 0; i < eigenVectors.size(); i++) {
	x = x + eigenVectors[i] * getRandom();
      }
      if (!equal_enough(M*x, x*eigen_value, rel_scale))
	success = false;
    }
    return success;
  }
}

void addSolveTests(Suite& suite, const string& test_name, const Matrix3& M,
		   const Vector& rhs, bool exp_return, int exp_xg_basis_size)
{
  Vector xp;
  std::vector<Vector> xg_basis;
  double rel_scale = M.MaxAbsElem();
  bool result = M.solve(rhs, xp, xg_basis, rel_scale);

  suite.findOrAddTest(test_name + ", existence", exp_return == result);

  if (result == true) {
    // this other stuff is only relevent if result == true
    suite.findOrAddTest(test_name + ", xp", equal_enough(M * xp, rhs,
							 rel_scale));
    /*    if (!equal_enough(M * xp, rhs, rel_scale)) {
      cout << rel_scale << endl;
      cout << xp << " != " << rhs << endl;
    }*/
	

    if (exp_xg_basis_size != xg_basis.size()) {
      suite.findOrAddTest(test_name + ", xg size", false);
      suite.findOrAddTest(test_name + ", xg_basis", false);
    }
    else {
      // test the basis by making arbitrary linear combinations
      // of basis vectors and testing them out
      bool success = true;
      for (int trial = 0; trial < 50; trial++) {
	Vector x(0, 0, 0);
	for (int i = 0; i < xg_basis.size(); i++) {
	  x = x + xg_basis[i] * getRandom();
	}
	if (!equal_enough(M*x, Vector(0, 0, 0), rel_scale))
	  success = false;
      }
      suite.findOrAddTest(test_name + ", xg_basis", success);
    }
  }
}

// returns a random number between -100.000 to 100.000
double getRandom()
{
  return (rand() % 200000 - 100000)/1000.0;
}

double getRandomOrZero(float probability_for_zero)
{
  if (rand() < probability_for_zero * RAND_MAX)
    return 0;
  else
    return getRandom();
}


// randomly mix up a matrix
void randomlyMixup(Matrix3& M, Vector& rhs)
{
  Matrix3 new_M(M);
  Vector new_rhs(rhs);
  int i, j, k;
  
  for (i = 1; i <= 3; i++) {
    for (j = 1; j <= 3; j++) {
      double mult = getRandom();
      if (i != j) {
	for (k = 1; k <= 3; k++)
	  new_M(i, k) += M(j, k) * mult;
	new_rhs(i-1) += rhs(j-1) * mult;
      }
      else {
	if (fabs(mult) > 1) {
	  for (k = 1; k <= 3; k++)
	    new_M(i, k) *= mult;
	  new_rhs(i-1) *= mult;
	}
      }
    }
  }

  // now randomly swap rows
  int row_orders[6][3] = {{1, 2, 3}, {1, 3, 2}, {2, 1, 3}, {2, 3, 1},
			  {3, 1, 2}, {3, 2, 1}};
  int* row_order = row_orders[rand() % 6];
  for (i = 1; i <= 3; i++) {
    for (j = 1; j <= 3; j++) {
      M(i, j) = new_M(row_order[i-1], j);
    }
    rhs(i-1) = new_rhs(row_order[i-1]-1);
  }
}

Vector randomVector()
{
  return Vector(getRandom(), getRandom(), getRandom());
}

Matrix3 randomMatrix()
{
  Matrix3 M;
  for (int i = 1; i <= 3; i++)
    for (int j = 1; j <= 3; j++)
      M(i, j) = getRandom();
  return M * pow(10.0, rand() % 10 - 3);
}

bool equal_enough(double x1, double x2, double rel_scale)
{
  return fabs(x2 - x1) <= NEAR_ZERO * (rel_scale > 1 ? rel_scale : 1);
}

bool equal_enough(Vector v1, Vector v2, double rel_scale)
{
  for (int i = 0; i < 3; i++) {
    if (!equal_enough(v1(i), v2(i), rel_scale))
      return false;
  }
  return true;
}



