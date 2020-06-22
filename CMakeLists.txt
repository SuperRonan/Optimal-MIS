cmake_minimum_required(VERSION 3.10)

project("Optimal-MIS" VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

OPTION(ONE_CONTIGUOUS_ARRAY "Store Optimal estimator buffer in one contiguous array" ON)

IF (ONE_CONTIGUOUS_ARRAY)
  ADD_DEFINITIONS ( -D ONE_CONTIGUOUS_ARRAY )
ENDIF()

SET(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
SET(EXAMPLE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/example)

SET(OPTIMIS_FILES ${INCLUDE_DIR}/Estimators.h ${INCLUDE_DIR}/ImageEstimators.h)
SET(EXAMPLE_FILES ${EXAMPLE_DIR}/Example.cpp)

# Eigen
FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)
INCLUDE_DIRECTORIES(ext/Eigen)

# Actual Optimal-MIS

add_library(OptiMIS INTERFACE)

target_sources(OptiMIS INTERFACE ${OPTIMIS_FILES})
target_include_directories(OptiMIS INTERFACE ${INCLUDE_DIR})
target_link_libraries(OptiMIS INTERFACE Eigen3::Eigen)

# Example 
add_executable(Example ${EXAMPLE_FILES})

target_link_libraries(Example OptiMIS)