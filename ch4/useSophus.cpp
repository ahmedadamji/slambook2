#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "sophus/se3.hpp" // In Sophus, when you include the header file for SE(3),
                          // you are also inherently bringing in the functionality for SO(3) as it's a part of the SE(3) representation.
                          // The library is structured in such a way that the necessary SO(3) functionality is available when you include the SE(3) header.

using namespace std;
using namespace Eigen;

/// This program demonstrates the basic usage of Sophus

int main(int argc, char **argv) {

  // Rotation matrix of 90 degrees about the Z axis
  Matrix3d R = AngleAxisd(M_PI / 2, Vector3d(0, 0, 1)).toRotationMatrix();
  // Or quaternion
  Quaterniond q(R);
  Sophus::SO3d SO3_R(R);              // Sophus::SO3d can be directly constructed from a rotation matrix
  Sophus::SO3d SO3_q(q);              // It can also be constructed through quaternion
  // Both are equivalent
  cout << "SO(3) from matrix:\n" << SO3_R.matrix() << endl;     // The matrix() function converts a Sophus::SO3d object to a matrix which is of type Eigen::Matrix3d
                                                                // as Sophus is developed directly on top of Eigen.
  cout << "SO(3) from quaternion:\n" << SO3_q.matrix() << endl;
  cout << "they are equal" << endl;

  // Use logarithmic mapping to get its Lie algebra
  Vector3d so3 = SO3_R.log();
  cout << "so3 = " << so3.transpose() << endl;
  // hat is the operator for vector to anti-symmetric matrix
  cout << "so3 hat=\n" << Sophus::SO3d::hat(so3) << endl;
  // Conversely, vee is the operator for anti-symmetric matrix to vector
  cout << "so3 hat vee= " << Sophus::SO3d::vee(Sophus::SO3d::hat(so3)).transpose() << endl;

  // Update model of incremental disturbance
  Vector3d update_so3(1e-4, 0, 0); // Assume the update amount is this much
  Sophus::SO3d SO3_updated = Sophus::SO3d::exp(update_so3) * SO3_R;
  cout << "SO3 updated = \n" << SO3_updated.matrix() << endl;

  cout << "*******************************" << endl;
  // Operations on SE(3) are very similar
  Vector3d t(1, 0, 0);           // Translation along X axis by 1
  Sophus::SE3d SE3_Rt(R, t);           // Construct SE(3) from R,t
  Sophus::SE3d SE3_qt(q, t);            // Construct SE(3) from q,t
  cout << "SE3 from R,t= \n" << SE3_Rt.matrix() << endl;
  cout << "SE3 from q,t= \n" << SE3_qt.matrix() << endl;
  // Lie algebra se(3) is a six-dimensional vector, for convenience, typedef it first
  typedef Eigen::Matrix<double, 6, 1> Vector6d;
  Vector6d se3 = SE3_Rt.log();
  cout << "se3 = " << se3.transpose() << endl;
  // Observe the output, you will find that in Sophus, se(3) has translation in front, rotation at the end.
  // Similarly, there are hat and vee two operators
  cout << "se3 hat = \n" << Sophus::SE3d::hat(se3) << endl;
  cout << "se3 hat vee = " << Sophus::SE3d::vee(Sophus::SE3d::hat(se3)).transpose() << endl;

  // Finally, demonstrate an update
  Vector6d update_se3; // Update amount
  update_se3.setZero();
  update_se3(0, 0) = 1e-4;
  Sophus::SE3d SE3_updated = Sophus::SE3d::exp(update_se3) * SE3_Rt;
  cout << "SE3 updated = " << endl << SE3_updated.matrix() << endl;

  return 0;
}
