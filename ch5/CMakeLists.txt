cmake_minimum_required(VERSION 3.0)  # Updated minimum version for better compatibility

set(CMAKE_BUILD_TYPE "Release")
# Enable C++11
set(CMAKE_CXX_FLAGS "-std=c++11 -O2")

# Eigen
include_directories("/usr/include/eigen3")

# OpenCV
set(OpenCV_DIR "/usr/local/include/opencv4")
# Find OpenCV 4.x
find_package(OpenCV 4 REQUIRED)
# Include OpenCV headers
include_directories(${OpenCV_INCLUDE_DIRS})

add_subdirectory(imageBasics)
add_subdirectory(stereo)
add_subdirectory(rgbd)
