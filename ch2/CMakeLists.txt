# ---------------------------------------------------------------------------------------------
# Specifies the minimum required version of CMake and project name

cmake_minimum_required(VERSION 3.10)
project(HelloSLAM)

# ---------------------------------------------------------------------------------------------
# Sets the C++ standard to C++17 and makes it a requirement

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# ---------------------------------------------------------------------------------------------
# Sets the source and include directories relative to the project root

set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

# ---------------------------------------------------------------------------------------------
# Adding Executables

# Adds executable 'helloSLAM'
add_executable(helloSLAM "${SOURCE_DIR}/helloSLAM.cpp")
target_include_directories(helloSLAM PUBLIC "${INCLUDE_DIR}")

# Adds executable 'useHello'
add_executable(useHello "${SOURCE_DIR}/useHello.cpp")
target_include_directories(useHello PUBLIC "${INCLUDE_DIR}")

# ---------------------------------------------------------------------------------------------
# Adding Libraries and Specifying include directories

# Adds static library 'Hello'
add_library(Hello "${SOURCE_DIR}/libHelloSLAM.cpp")
target_include_directories(Hello PUBLIC "${INCLUDE_DIR}")

# Adds shared library 'HelloSharedLib'
add_library(HelloSharedLib SHARED "${SOURCE_DIR}/libHelloSLAM.cpp")
target_include_directories(HelloSharedLib PUBLIC "${INCLUDE_DIR}")

# ---------------------------------------------------------------------------------------------
# Specifying Instructions for Installation on local machine for use by other projects in any directory

# Specifies that 'libHelloSLAM.h' is the public header for 'HelloSharedLib'
set_target_properties(HelloSharedLib PROPERTIES PUBLIC_HEADER "${INCLUDE_DIR}/libHelloSLAM.h")

# Sets installation rules for 'HelloSharedLib'
install(TARGETS HelloSharedLib
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include
)

# ---------------------------------------------------------------------------------------------
# Linking Libraries with Executables

# Links 'useHello' executable with 'HelloSharedLib'
target_link_libraries(useHello HelloSharedLib)

# ---------------------------------------------------------------------------------------------