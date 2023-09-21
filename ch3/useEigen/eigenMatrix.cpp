#include <iostream>

using namespace std;

#include <ctime>

// Eigen core part
// This header file includes the core functionalities of the Eigen library. It provides the basic matrix and array classes along with the necessary methods to manipulate them.
// It contains definitions for objects like MatrixXd (a matrix of doubles), VectorXd (a column vector of doubles), etc.
// It is essential to include this header file when you are using Eigen to work with matrices and vectors in C++.
#include <Eigen/Core>

// Algebraic operations (inverse, eigenvalues, etc.) for dense matrices
// This header file includes dense matrix and vector operations like matrix inverses, determinants, eigenvalues, etc.
// Basically, it extends the core functionalities provided by <Eigen/Core> to include more advanced linear algebra operations.
// This header file is necessary when you want to perform operations that go beyond simple matrix and vector manipulations.
#include <Eigen/Dense>

// Geometry module
// This header file includes the necessary functionalities for working with 2D and 3D rotations and transformations.
// It provides the necessary classes and methods to represent rotations and transformations in 2D and 3D.
#include <Eigen/Geometry>

using namespace Eigen;

#define MATRIX_SIZE 50

/****************************
 * This program demonstrates the basic usage of Eigen types
 ****************************/

int main(int argc, char **argv)
{
     // In Eigen, all vectors and matrices are Eigen::Matrix, which is a template class. Its first three parameters are: data type, rows, columns
     // Declare a 2x3 float matrix
     Matrix<float, 2, 3> matrix_23;
     // Declare a 4x4 float matrix
     Eigen::Matrix<float, 4, 4> matrix_44 = Eigen::Matrix<float, 4, 4>::Zero();

     // Also, Eigen provides many built-in types through typedef, but the underlying layer is still Eigen::Matrix
     // For example, Vector3d is essentially Eigen::Matrix<double, 3, 1>, which is a three-dimensional vector
     Vector3d v_3d;
     // This is the same
     Matrix<float, 3, 1> vd_3d;
     // It is also the same as
     Eigen::VectorXd v3d = Eigen::VectorXd::Zero(3);
     // Can even set it as random
     Eigen::VectorXd v3d_random = Eigen::VectorXd::Random(3);
     Eigen::VectorXf vd(3);
     // Here because the size of vd is 3, but 5 values are assigned, so the last two values are not assigned
     vd << 1, 2, 3, 4, 5;

     // Matrix3d is essentially Eigen::Matrix<double, 3, 3>
     Matrix3d matrix_33 = Matrix3d::Zero(); // Initialize to zero
     // If the size of the matrix is uncertain, you can use a dynamically sized matrix
     Matrix<double, Dynamic, Dynamic> matrix_dynamic;
     // To use this:
     matrix_dynamic = MatrixXd::Random(3, 3); // Random number matrix

     // Simpler
     MatrixXd matrix_x;
     // There are many types like this, we will not list them one by one

     // Below are operations on Eigen matrices
     // Input data (initialization)
     matrix_23 << 1, 2, 3, 4, 5, 6;
     // Output
     cout << "matrix 2x3 from 1 to 6: \n"
          << matrix_23 << endl;

     // Use () to access elements in the matrix
     cout << "print matrix 2x3: " << endl;
     for (int i = 0; i < 2; i++)
     {
          for (int j = 0; j < 3; j++)
               cout << matrix_23(i, j) << "\t";
          cout << endl;
     }

     // Matrix and vector multiplication (actually still matrix and matrix)
     v_3d << 3, 2, 1;
     vd_3d << 4, 5, 6;

     // But in Eigen you can't mix two different types of matrices, this is wrong
     // Matrix<double, 2, 1> result_wrong_type = matrix_23 * v_3d;
     // Should be explicitly converted
     Matrix<double, 2, 1> result = matrix_23.cast<double>() * v_3d;
     cout << "[1,2,3;4,5,6]*[3,2,1]=" << result.transpose() << endl;

     Matrix<float, 2, 1> result2 = matrix_23 * vd_3d;
     cout << "[1,2,3;4,5,6]*[4,5,6]: " << result2.transpose() << endl;

     // Similarly, you can't get the dimensions of the matrix wrong
     // Try to uncomment the following line and see what error Eigen reports
     // Eigen::Matrix<double, 2, 3> result_wrong_dimension = matrix_23.cast<double>() * v_3d;

     // Some matrix operations
     // Basic arithmetic operations are not demonstrated, you can directly use +-*/.
     matrix_33 = Matrix3d::Random(); // Random number matrix
     cout << "random matrix: \n"
          << matrix_33 << endl;
     cout << "transpose: \n"
          << matrix_33.transpose() << endl;          // Transpose
     cout << "sum: " << matrix_33.sum() << endl;     // Sum of elements
     cout << "trace: " << matrix_33.trace() << endl; // Trace
     cout << "times 10: \n"
          << 10 * matrix_33 << endl; // Scalar multiplication
     cout << "inverse: \n"
          << matrix_33.inverse() << endl;                // Inverse
     cout << "det: " << matrix_33.determinant() << endl; // Determinant

     // Some vector operations
     float dot_result = v_3d.dot(vd_3d.cast<double>()); // Dot product
     cout << "dot product: " << dot_result << endl;
     Eigen::Vector3d cross_result = v_3d.cross(vd_3d.cast<double>()); // Cross product
     cout << "cross product: " << cross_result.transpose() << endl;


     // Eigenvalues
     // Real symmetric matrices can guarantee successful diagonalization
     SelfAdjointEigenSolver<Matrix3d> eigen_solver(matrix_33.transpose() * matrix_33);
     cout << "Eigen values = \n"
          << eigen_solver.eigenvalues() << endl;
     cout << "Eigen vectors = \n"
          << eigen_solver.eigenvectors() << endl;

     // Quaternions
     Eigen::Quaterniond q_from_matrix_33 = Eigen::Quaterniond(matrix_33);                         // Initialize from a matrix
     Eigen::Quaterniond q_from_angle_axis(Eigen::AngleAxisd(M_PI / 4, Eigen::Vector3d::UnitZ())); // Initialize from angle axis
     // Quaterniond can be converted to rotation matrix through .toRotationMatrix()
     Eigen::Matrix3d rot = q_from_matrix_33.normalized().toRotationMatrix();

     // Angle axis
     Eigen::AngleAxisd angle_axis = Eigen::AngleAxisd(M_PI / 4, Eigen::Vector3d::UnitZ());
     // Convert rotation matrix to angle axis
     Eigen::AngleAxisd angle_axis_from_rotation_matrix = Eigen::AngleAxisd(matrix_33);
     Eigen::AngleAxisd angle_axis_from_rotation_matrix2;
     angle_axis_from_rotation_matrix2.fromRotationMatrix(matrix_33);
     // Convert quaternion to angle axis
     Eigen::AngleAxisd angle_axis_from_quaternion = Eigen::AngleAxisd(q_from_matrix_33);

     // Euler angles
     // The numbers 0, 1, and 2 specify the rotation axes order for Euler angles decomposition,
     // representing rotations around the x, y, and z axes, respectively.
     Eigen::Vector3d euler_angles = rot.eulerAngles(0, 1, 2);
     // Quaternion can also be converted to Euler angles
     Eigen::Vector3d euler_angles_from_quaternion = q_from_matrix_33.normalized().toRotationMatrix().eulerAngles(0, 1, 2);
     // Convert angle axis to Euler angles
     Eigen::Vector3d euler_angles_from_angle_axis = angle_axis.angle() * angle_axis.axis();
     // Convert Euler angles to rotation matrix
     Eigen::Matrix3d rot_matrix_from_euler_angles;
     rot_matrix_from_euler_angles = Eigen::AngleAxisd(euler_angles[0], Eigen::Vector3d::UnitX()) *
                                    Eigen::AngleAxisd(euler_angles[1], Eigen::Vector3d::UnitY()) *
                                    Eigen::AngleAxisd(euler_angles[2], Eigen::Vector3d::UnitZ());

     // Solve equations
     // We solve the equation matrix_NN * x = v_Nd
     // The size of N is defined in the macro above, it is generated by random numbers
     // Direct inversion is naturally the most direct, but inversion is computationally intensive

     Matrix<double, MATRIX_SIZE, MATRIX_SIZE> matrix_NN = MatrixXd::Random(MATRIX_SIZE, MATRIX_SIZE);
     matrix_NN = matrix_NN * matrix_NN.transpose(); // Ensure semi-positive definite
     Matrix<double, MATRIX_SIZE, 1> v_Nd = MatrixXd::Random(MATRIX_SIZE, 1);

     clock_t time_stt = clock(); // Timing
     // Direct inversion
     Matrix<double, MATRIX_SIZE, 1> x = matrix_NN.inverse() * v_Nd;
     cout << "time of normal inverse is "
          << 1000 * (clock() - time_stt) / (double)CLOCKS_PER_SEC << "ms" << endl;
     cout << "x = " << x.transpose() << endl;

     // Usually we use matrix decomposition to solve, such as QR decomposition, which will be much faster
     time_stt = clock();
     x = matrix_NN.colPivHouseholderQr().solve(v_Nd);
     cout << "time of Qr decomposition is "
          << 1000 * (clock() - time_stt) / (double)CLOCKS_PER_SEC << "ms" << endl;
     cout << "x = " << x.transpose() << endl;

     // For positive definite matrices, you can also use cholesky decomposition to solve equations
     time_stt = clock();
     x = matrix_NN.ldlt().solve(v_Nd);
     cout << "time of ldlt decomposition is "
          << 1000 * (clock() - time_stt) / (double)CLOCKS_PER_SEC << "ms" << endl;
     cout << "x = " << x.transpose() << endl;

     return 0;
}
