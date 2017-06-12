#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  //cout << "CalculateRMSE" << endl;

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  if (estimations.size() == 0)
    cout << "Estimation vector size is zero" << endl;

  //  * the estimation vector size should equal ground truth vector size'
  if (estimations.size() != ground_truth.size())
    cout << "Estimation vector size is not equal to ground truth vector size" << endl;

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    VectorXd e = estimations[i] - ground_truth[i];
    e = e.array() * e.array();
    rmse += e;
  }

  //calculate the mean
  // ... your code here
  rmse = rmse/estimations.size();

  //calculate the squared root
  // ... your code here
  rmse = rmse.array().sqrt();

  //cout << "end CalculateRMSE" << endl;

  //return the result
  return rmse;
}
