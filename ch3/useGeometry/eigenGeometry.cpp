#include <iostream>
#include <cmath>

using namespace std;

#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

// This program demonstrates how to use the Eigen Geometry module

int main(int argc, char **argv) {

  // The Eigen/Geometry module provides various representations for rotation and translation
  // For 3D rotation matrices, directly use Matrix3d or Matrix3f
  Matrix3d rotation_matrix = Matrix3d::Identity();
  // For rotation vectors, use AngleAxis. It's not directly a matrix, but can be used like one (because operators are overloaded)
  AngleAxisd rotation_vector(M_PI / 4, Vector3d(0, 0, 1));  // 45 degrees rotation along Z axis
  cout.precision(3);
  cout << "rotation matrix =\n" << rotation_vector.matrix() << endl;  // Convert to matrix using .matrix()
  // Direct assignment is also possible
  rotation_matrix = rotation_vector.toRotationMatrix();
  // Coordinate transformation can be done using AngleAxis
  Vector3d v(1, 0, 0);
  Vector3d v_rotated = rotation_vector * v;
  cout << "(1,0,0) after rotation (by angle axis) = " << v_rotated.transpose() << endl;
  // Or use the rotation matrix
  v_rotated = rotation_matrix * v;
  cout << "(1,0,0) after rotation (by matrix) = " << v_rotated.transpose() << endl;

  // Euler Angles: The rotation matrix can be directly converted to Euler angles
  Vector3d euler_angles = rotation_matrix.eulerAngles(2, 1, 0);  // ZYX order, i.e., yaw-pitch-roll
  cout << "yaw pitch roll = " << euler_angles.transpose() << endl;

  // For Euclidean transformation matrices, use Eigen::Isometry
  Isometry3d T = Isometry3d::Identity();  // Despite being called 3d, it's actually a 4x4 matrix
  T.rotate(rotation_vector);              // Rotate according to rotation_vector
  T.pretranslate(Vector3d(1, 3, 4));      // Set the translation vector to (1,3,4)
  cout << "Transform matrix = \n" << T.matrix() << endl;

  // Coordinate transformation using the transformation matrix
  Vector3d v_transformed = T * v;  // Equivalent to R*v+t
  cout << "v transformed = " << v_transformed.transpose() << endl;

  // For affine and projective transformations, use Eigen::Affine3d and Eigen::Projective3d, omitted here

  // Quaternions
  // You can directly assign AngleAxis to quaternion, and vice versa
  Quaterniond q = Quaterniond(rotation_vector);
  cout << "quaternion from rotation vector = " << q.coeffs().transpose()
       << endl;  // Note that the order of coeffs is (x,y,z,w), w is the real part, and the first three are the imaginary parts
  // You can also assign it from a rotation matrix
  q = Quaterniond(rotation_matrix);
  cout << "quaternion from rotation matrix = " << q.coeffs().transpose() << endl;
  // To rotate a vector using a quaternion, use overloaded multiplication
  v_rotated = q * v;  // Mathematically, this is qvq^{-1}
  cout << "(1,0,0) after rotation = " << v_rotated.transpose() << endl;
  // If using normal vector multiplication, it should be calculated as follows
  cout << "should be equal to " << (q * Quaterniond(0, 1, 0, 0) * q.inverse()).coeffs().transpose() << endl;

  return 0;
}
