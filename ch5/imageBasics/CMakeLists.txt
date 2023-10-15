project(imageBasics)

add_executable(imageBasics imageBasics.cpp)

# OpenCV
set(OpenCV_DIR "/usr/local/lib/cmake/opencv4")
# Find OpenCV 4.x
find_package(OpenCV 4 REQUIRED)
# Include OpenCV headers
include_directories(${OpenCV_INCLUDE_DIRS})

target_link_libraries(imageBasics ${OpenCV_LIBS})

add_executable(undistortImage undistortImage.cpp)

target_link_libraries(undistortImage ${OpenCV_LIBS})
