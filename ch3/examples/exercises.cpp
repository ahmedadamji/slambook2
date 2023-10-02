#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <iostream>
#include <vector>


using namespace std;
using namespace Eigen;

class Exercises{
public:

    Exercises()
    {

    }

    void editMatrixBlock(MatrixXd &matrix, int startRow, int startCol, int numRows, int numCols, MatrixXd &block)
    {
        
        std::cout << "current block: " << matrix.block(startRow, startCol, numRows, numCols) << std::endl;
        // Set the block
        matrix.block(startRow, startCol, numRows, numCols) = block;
        std::cout << "matrix: " << matrix << std::endl;

    }

    // Function to solve Ax = b for x
    // Conditions for Unique Solution:
    // 1. A should be a square matrix (n x n)
    // 2. A should be of full rank (all rows/columns are linearly independent)
    // 3. A should be non-singular (det(A) != 0)
    //
    // Numerical Method:
    // Uses Eigen's colPivHouseholderQr() method for solving Ax = b.
    Eigen::VectorXd solveLinearEquation(const Eigen::MatrixXd& A, const Eigen::VectorXd& b)
    {

        // Do not use A.inverse() due to reasons:
        // 1. High computational cost (O(n^3)).
        // 2. Numerical stability issues.
        // 3. Not suitable for singular/nearly singular matrices.

        Eigen::VectorXd x;

        // Check if A is non-singular
        if (A.determinant() == 0.0)
        {
            std::cout << "A is singular. Cannot solve." << std::endl;
            exit(1); // or handle error as appropriate
        }

        // Check if A is not square
        if (A.rows() != A.cols())
        {
            // Either least-norm solution for underdetermined (m < n) 
            // or least-squares solution for overdetermined (m > n).
            // Both can be calculated using Moore-Penrose pseudoinverse via SVD.
            x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
        }

        // Solve Ax = b for x
        int n = A.rows();
        if (n < 100) { // For small matrices
            // Use FullPivLU for general-purpose, dense matrices
            x = A.fullPivLu().solve(b);
        } else if (A.isApprox(A.transpose())) { // For symmetric matrices
            // Use LLT for symmetric, positive definite matrices
            x = A.llt().solve(b);
        } else if (n < 500) { // For medium-size matrices
            // Use QR decomposition for more stable solutions
            x = A.colPivHouseholderQr().solve(b);
        } else { // For large matrices
            // Use Conjugate Gradient for large, sparse, symmetric, positive-definite matrices
            // Assuming A is sparse and of type SparseMatrix<double>
            // Eigen::SparseMatrix<double> spA = A.sparseView();
            // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
            // cg.compute(spA);
            // x = cg.solve(b);
        }


        return x;
    }
    

};

int main(int argc, char **argv)
{
    Exercises exercises;
    MatrixXd matrix = MatrixXd::Random(5, 5);
    MatrixXd block = Matrix3d::Identity();
    std::cout << "matrix: " << matrix << std::endl;
    std::cout << "block: " << block << std::endl;
    exercises.editMatrixBlock(matrix, 0, 0, 3, 3, block);

    // Solve Ax = b for x
    MatrixXd A = MatrixXd::Random(3, 3);
    VectorXd b = VectorXd::Random(3);
    std::cout << "A: " << A << std::endl;
    std::cout << "b: " << b << std::endl;
    VectorXd x = exercises.solveLinearEquation(A, b);
    std::cout << "x: " << x << std::endl;

    return 0;
}