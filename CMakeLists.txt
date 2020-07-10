cmake_minimum_required(VERSION 3.10)

project("Optimal-MIS" VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

IF (ONE_CONTIGUOUS_ARRAY)
  ADD_DEFINITIONS ( -D ONE_CONTIGUOUS_ARRAY )
ENDIF()

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/MIS)
SET(EXAMPLE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/example)

SET(OPTIMIS_FILES 
${SRC_DIR}/Estimator.h ${SRC_DIR}/SimpleEstimator.h ${SRC_DIR}/BalanceEstimator.h ${SRC_DIR}/DirectEstimator.h ${SRC_DIR}/PowerEstimator.h ${SRC_DIR}/NaiveEstimator.h ${SRC_DIR}/CutOffEstimator.h ${SRC_DIR}/MaximumEstimator.h
${SRC_DIR}/ImageEstimator.h ${SRC_DIR}/ImageSimpleEstimator.h ${SRC_DIR}/ImageBalanceEstimator.h ${SRC_DIR}/ImageDirectEstimator.h ${SRC_DIR}/ImagePowerEstimator.h ${SRC_DIR}/ImageNaiveEstimator.h ${SRC_DIR}/ImageCutOffEstimator.h ${SRC_DIR}/ImageMaximumEstimator.h
)

SET(OPTIMIS_UTILS_FILES ${SRC_DIR}/utils/Parallel.h ${SRC_DIR}/utils/SpectrumWrapper.h ${SRC_DIR}/utils/settings.h)

SET(EXAMPLE_FILES ${EXAMPLE_DIR}/Example.cpp ${EXAMPLE_DIR}/Spectrum.h)

# Eigen
FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)
INCLUDE_DIRECTORIES(ext/Eigen)

# Actual Optimal-MIS

add_library(OptiMIS INTERFACE)

SOURCE_GROUP("MIS" OPTIMIS_FILES)
SOURCE_GROUP("MIS/utils" OPTIMIS_UTILS_FILES)

target_sources(OptiMIS INTERFACE ${OPTIMIS_FILES} ${OPTIMIS_UTILS_FILES})
target_include_directories(OptiMIS INTERFACE ${SRC_DIR})
target_link_libraries(OptiMIS INTERFACE Eigen3::Eigen)

# Example 
add_executable(Example ${EXAMPLE_FILES})
INCLUDE_DIRECTORIES(src)
target_link_libraries(Example OptiMIS)