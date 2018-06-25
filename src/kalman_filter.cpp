#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
    VectorXd y = z - H_ * x_;
    UpdateY(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);
    if (px != 0 && py != 0) {
//    if(px < 0.0001) px = 0.0001;
//    if(py < 0.0001) py = 0.0001;
    
    double rho = sqrt(px * px + py * py);
    double phi = atan2(py, px);
    double rhoD = (px*vx + py*vy) / rho;
    
    VectorXd h = VectorXd(3);
    h << rho, phi, rhoD;
    
    VectorXd y = z - h;
    while ( y(1) > M_PI || y(1) < -M_PI ) {
        if ( y(1) > M_PI ) {
            y(1) -= M_PI;
        } else {
            y(1) += M_PI;
        }
    }
    UpdateY(y);
    }
}

void KalmanFilter::UpdateY(const VectorXd &y) {
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = P_ * Ht * Si;
    
    x_ = x_ + K * y;
//    MatrixXd I = MatrixXd::Identity(4, 4);
    int x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}