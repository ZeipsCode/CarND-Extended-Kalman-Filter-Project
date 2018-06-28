#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  

  //measurement covariance - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

    //state covariance matrix P

    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;
    

    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
                0, 1, 0, 1,
                0, 0, 1, 0,
                0, 0, 0, 1;


    //measurement matrix
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);      
    ekf_.H_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;      
    ekf_.H_ = H_laser_;
    

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /*
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
            
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      
      double vx = 0;
      double vy = 0;

      ekf_.x_ << px, py, vx, vy;  


      return;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

      /**
      Initialize state.
      */

      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      
    }

    // done initializing, no need to predict or update

    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    
    return;

  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

    double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
  

    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;            
    
    double noise_ax = 9.0;
    double noise_ay = 9.0;
    
    double ax_2 = pow(noise_ax, 2);
    double ay_2 = pow(noise_ay, 2);
    
    double dt_2 = pow(dt, 2);
    double dt_3 = pow(dt, 3);
    double dt_4 = pow(dt, 4);
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << (dt_4 * ax_2)/4, 0, (dt_3 * ax_2)/2, 0,
                0, (dt_4 * ay_2)/4, 0, (dt_3 * ay_2)/2,
                (dt_3 * ax_2)/2, 0, dt_2 * ax_2, 0,
                0, (dt_3 * ay_2)/2, 0, dt_2 * ay_2;
    

    ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
  {
    // Radar updates
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } 
  
  else 

  {
    // Lidar updates
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}