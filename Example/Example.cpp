#include <iostream>
#include <Estimators.h>

int main(int argc, char ** argv)
{
    MIS::BalanceEstimator<double, double> estimator(1);
    estimator.addEstimate(1.0, nullptr, 0);
    std::cout<<"Hello World!"<<std::endl;
    std::cout << "Estimator result: " << estimator.solve(1) << std::endl;
}