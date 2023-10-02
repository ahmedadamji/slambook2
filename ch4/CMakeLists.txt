cmake_minimum_required(VERSION 3.0)
project(useSophus)

# Searches for the Sophus library in the system. The REQUIRED flag indicates that this is a mandatory dependency.
# As it is a required dependency, if Sophus is not found, CMake will stop processing the project and report an error.
# If it is not a required dependency, CMake will print a warning and continue processing.
find_package(Sophus REQUIRED)

# Eigen
include_directories("/usr/include/eigen3")
add_executable(useSophus useSophus.cpp)
target_link_libraries(useSophus Sophus::Sophus)

add_subdirectory(example)
