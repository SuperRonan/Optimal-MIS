# Optimal-MIS
Practical implementation of Optimal MIS for Light Transport.

This repository provides a C++ header only library of Practical Optimal Multiple Importance Sampling, as presented in the 
[TODO: Practical Optimal MIS technical report](https://github.com/SuperRonan/Optimal-MIS), based on the original 
[Optimal MIS paper](https://cgg.mff.cuni.cz/~jaroslav/papers/2019-optimal-mis/).
We also made a (https://github.com/SuperRonan/PBRT-Optimal-MIS) using this library for an example in a renderer. 


The library provides two abstract classes, `Estimator` (in `src/MIS/Estimator.h`) and `ImageEstimator` (in `src/MIS/ImageEstimator.h`).
`Estimator` is entended to be used by algorithms computing one pixel at once (like a path tracer), whereas `ImageEstimator` 
is intended to be used by algorithms for which the any pixel can be updated (like a bidirectional path tracer, for the contribution of the light tracer). 
The implementation of `ImageEstimator` is optimized over types like `Image<Estimator>`. The library uses the namespace `MIS`.

There are several implementation of Estimators for different weigthing functions
- Eric Veach's balance, power, cutoff, maximum heuristics in `MIS::BalanceEstimator`, `MIS::PowerEstimator`, `MIS::CutOffEstimator`,`MIS::MaximumEstimator`
and `MIS::ImageBalanceEstimator`, `MIS::ImagePowerEstimator`, `MIS::ImageCutOffEstimator`,`MIS::ImageMaximumEstimator` respectively.
- A naive heuristic (1 / N) in `MIS::NaiveEstimator` and `MIS::ImageNaiveEstimator`.
- The direct estimator of the Optimal MIS weights in `MIS::DirectEstimator` and `MIS::ImageDirectEstimator` (slightly biased but efficient).
- TODO The progressive estimator of the Optimal MIS weights in `MIS::ProgressiveEstimator` and `MIS::ImageProgressiveEstimator` 
(unbiased, but far less efficient, both in runtime and variance reduction).

The files `src/MIS/Estimators.h` includes all the estimators, and `src/MIS/ImageEstimators.h` includes all the image estimators.

The estimators are templated by `Spectrum` and `Float`. 
- `Spectrum` is intended to provide the functions described in `Example/Spectrum.h`, but can also be set to `float` or `double`. 
- `Float` is intended to be `float` or `double`.

Using the library
-----------------

To clone the library: `git clone https://github.com/SuperRonan/Optimal-MIS.git --recurse-submodules`. 
Do not forget the `--recurse-submodules` to also get Eigen.

Since this is a header only library, you need to add the `src` directory to the include directories of your project. 
For example with CMake, use: `INCLUDE_DIRECTORIES(PathToOptimal-MIS/src)`. 
Optimal-MIS uses [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), another header only library for linear algebra.
If Eigen is not already included, it is a submodule of this repository, under `ext/Eigen`. 
To include it with CMake, use: `INCLUDE_DIRECTORIES(PathToOptimal-MIS/ext/Eigen)`. 
Make sure it is `ext/Eigen` that is added to your include directories, not `ext/Eigen/Eigen`.
You can also do some CMake black magic I don't really know about (like `FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)` or something).
Finaly, OptiMIS can benefit from OpenMP. 
It can easily be activated in CMake with:  
`FIND_PACKAGE(OpenMP)`  
`SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")`  
`SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")`  
