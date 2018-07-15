// tools.cpp- function to calculate RMSE and the Jacobian matrix

#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

// #################
// CALCULATE RMSE
// #################
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // Prepare RMSE variable and initialize with zeros
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  // Script abortion: estimation vector too small or size unequal to ground truth vector
  if (estimations.size() < 1) {
    return rmse;
  }
  if (estimations.size() != ground_truth.size()) {
    return rmse;
  }
  
  // Calculate squared residuals
  for (int j=1; j < estimations.size(); ++j) {
    // Calculate 
    VectorXd resid = estimations[j] - ground_truth[j];
    
    // Coefficient-wise multiplication
    resid = resid.array() * resid.array();
    rmse += resid;
  }
  
  // Calculate the mean
  rmse = rmse / estimations.size();
  
  // Calculate square root of RMSE
  rmse = rmse.array().sqrt();
  
  // Return RMSE result, exit function
  return rmse;
}


// #################
// CALCULATE JACOBIAN
// #################
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  // Define variable to hold Jacobian matrix
  MatrixXd Hj(3,4);
  
  // Recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  // Check for division by zero, exit function
  float pxy2 = px*px + py*py;
  if (pxy2 < 0.00001) {
    cout << "Division by zero detected!" << "\n";
    px += 0.001;
    py += 0.001;
    pxy2 = px*px + py*py;
  }
  
  // Compute the Jacobian matrix
  float pxy2_sqrt = sqrt(pxy2);
  float vxpy_vypx = vx*py - vy*px;
  float l3_term = vxpy_vypx / (pxy2 * pxy2_sqrt);
  
  Hj <<   px/pxy2_sqrt,   py/pxy2_sqrt,  0,               0,
          -py/pxy2,       px/pxy2,       0,               0,
          py*l3_term,     -px*l3_term,   px/pxy2_sqrt,    py/pxy2_sqrt;
  
  return Hj;
}
