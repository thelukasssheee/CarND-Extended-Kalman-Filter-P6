// kalman_filter.cpp- defines the predict function, the update function for lidar, and the update function for radar
//

#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

// Predict position and Covariance matrix
void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

// Update position (laser measurements)
void KalmanFilter::Update(const VectorXd &z) {
//  cout << "[DEBUG] EKF z \n";
//  cout << z << endl;
//  cout << "[DEBUG] EKF H \n";
//  cout << H_ << endl;
  VectorXd y = z - H_ * x_;
//  cout << "[DEBUG] EKF y \n";
//  cout << y << endl;
  
  // Calculate innovation covariance S
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  
  // Calculcate optimal Kalman gain
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  
  // Calculate new position
  x_ = x_ + K * y;
  
  // Calculate new covariance
  MatrixXd I = MatrixXd::Identity(x_.size(),x_.size());
  P_ = (I - K * H_) * P_;
}

// Update position (radar measurements)
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  // Calculate h-function (1st element)
  VectorXd h_func(3);
  h_func(0) = sqrt(px*px + py*py);
  
  // Check for division by zero
  if (h_func(0) < 0.000001) {
    px += 0.001;
    py += 0.001;
    h_func(0) = sqrt(px*px + py*py);
  }

  // Calculate h-function (2nd & 3rd element)
  h_func(1) = atan2(py,px);
  h_func(2) = (px*vx + py*vy) / h_func(0);

  // Calculate predicted position vector y in measurement space
  VectorXd y = z - h_func;

  // Normalize angle phi in y to be in the range between -pi and +pi
  while (y(1) < -M_PI) {
    y(1) += 2*M_PI;
  }
  while (y(1) > M_PI) {
    y(1) -= 2*M_PI;
  }
  
  // Calculate innovation covariance S
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  
  // Calculcate optimal Kalman gain
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  
  // Calculate new position
  x_ = x_ + K * y;
   
  // Calculate new covariance
  MatrixXd I = MatrixXd::Identity(x_.size(),x_.size());
  P_ = (I - K * H_) * P_;
}
