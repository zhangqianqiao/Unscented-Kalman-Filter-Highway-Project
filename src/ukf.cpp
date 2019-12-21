#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 4;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  is_initialized_ = false;

  //State dimension
  n_x_ = 5;

  // Augemented state dimension
  n_aug_ = 7;

  //sigma point spreading parameter;
  lambda_ = 3 - n_aug_;

  //initial sigma predict points matrix;
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  //initial sigma points matrix;
  Xsig_out = MatrixXd(n_aug_, 2*n_aug_+1);

  //initial weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);

  //initial time
  time_us_ = 0;

  //initial radar measurement noise
  R_radar = MatrixXd(3, 3);

  //initial ladar measurement noise
  R_lidar = MatrixXd(2, 2);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /************************************************************
   * Initialization the state for the first time
   ************************************************************/

  if(!is_initialized_)
  { 
    //Initial state
    x_ << 1, 1, 0, 0, 0;

    //Initial state covariance matrix P
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, M_PI, 0,
          0, 0, 0, 0, 1;
    
    // Initialize weights.
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < 2*n_aug_+1; ++i)
	  {
		  weights_(i) = 0.5 / (lambda_ + n_aug_);
	  }

    //initialize lidar measurement noise;
	  R_lidar << std_laspx_ * std_laspx_, 0,
		            0, std_laspy_*std_laspy_;
    
    //initiallize radar measurement noise
    R_radar << std_radr_ * std_radr_, 0, 0,
		        0, std_radphi_*std_radphi_, 0,
		        0, 0, std_radrd_*std_radrd_;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    { 
      //first measurement is lidar
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);

      time_us_ = meas_package.timestamp_;
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      //first measurement is Radar
      float ro = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float ro_dot = meas_package.raw_measurements_(2);

      x_(0) = ro * cos(phi);
      x_(1) = ro * sin(phi);
    
      time_us_ = meas_package.timestamp_;
    }
    
    is_initialized_ = true;
    return;
  }

  /************************************************************
   * Prediction
   ************************************************************/
  double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;

  //Precict state and covriance
  Prediction(dt);

  /************************************************************
   * Update
   ************************************************************/
  
  if(meas_package.sensor_type_==MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
    //std::cout << "This Radar nothing to do"<<"\n";
  }else{
    UpdateLidar(meas_package);
  }
}

void UKF::GenerateSigmaPoints()
{
  /*
   * Generate the sigma points
   */

  //create augmented mean matrix;
  VectorXd x_aug_ = VectorXd(7);

  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5, 5) = std_a_*std_a_;
  P_aug(6, 6) = std_yawdd_*std_yawdd_;

  //Calculate square root of P_aug
  MatrixXd A_aug_ = P_aug.llt().matrixL();

  //Set first column of sigma points matrix
  Xsig_out.col(0) = x_aug_;

  //Set remaining sigma points
  for(int i=0; i < n_aug_; ++i)
  {
    Xsig_out.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_)*A_aug_.col(i);
    Xsig_out.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_)*A_aug_.col(i);
  }

}

void UKF::SigmaPointPrediction(double delta_t)
{
  
  //predict sigma points
  for(int i=0; i<2*n_aug_+1; ++i)
  {
    double p_x = Xsig_out(0, i);
    double p_y = Xsig_out(1, i);
    double v = Xsig_out(2, i);
    double yaw = Xsig_out(3, i);
    double yawd = Xsig_out(4, i);

    double nu_a = Xsig_out(5,i);
    double nu_yawdd = Xsig_out(6,i);

    //predict state value
    double px_p, py_p;

    //avoid division by zero
    if(fabs(yawd)>0.001)
    {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

     // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;

  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //Generater the sigma points in Sig_out
  GenerateSigmaPoints();

  //Predict sigma points in Sig_pred_
  SigmaPointPrediction(delta_t);


  //predicted state mean
  x.fill(0.0);
  for(int i=0; i<2*n_aug_+1; ++i)
  {
     x = x + weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //Update state and covariance
  x_ = x;
  P_ = P;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  //Set measurement dimension
  int n_lidar_ = 2;

  //create martix for sigma point in measurement space
  MatrixXd Zsig = MatrixXd(n_lidar_, 2*n_aug_+1);

  //Mean predicted meansurement
  VectorXd z_pred = VectorXd(n_lidar_);

  //Measurement covariance matrix S;
  MatrixXd S_ = MatrixXd(n_lidar_, n_lidar_);

  for (int i = 0; i < 2 * n_aug_+1; ++i)
  {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }

  //mean predicted measurement
  z_pred.fill(0.0);
  for(int i = 0; i<2*n_aug_+1; ++i)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation convariance matrix s
  S_.fill(0.0);

  for(int i = 0; i < 2*n_aug_+1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S_ = S_ + weights_(i)*z_diff*z_diff.transpose();
  }

  //add measurement noise
  S_ = S_ + R_lidar;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_lidar_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; ++i)
  { 
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // create vector for incoming radar measurement
  VectorXd z = VectorXd(n_lidar_);
  z <<  meas_package.raw_measurements_(0),
        meas_package.raw_measurements_(1);
  
  //kalman gain k;
  MatrixXd K = Tc * S_.inverse();

  VectorXd z_diff = z - z_pred;

  // angle normalization
  //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  double nis_ = z_diff.transpose() * S_.inverse() * z_diff;

  //cout <<"The Lidar meansurement NIS: " << nis_ << "\n";

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_*K.transpose();

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  //Set measurement dimension
  int n_radar_ = 3;

  //create martix for sigma point in measurement space
  MatrixXd Zsig = MatrixXd(n_radar_, 2*n_aug_+1);

  //Mean predicted meansurement
  VectorXd z_pred = VectorXd(n_radar_);

  //Measurement covariance matrix S;
  MatrixXd S = MatrixXd(n_radar_, n_radar_);


  for (int i = 0; i < 2 * n_aug_+1; ++i)
  {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    //measurement model;
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y); // r
    Zsig(1, i) = atan2(p_y, p_x);         //phi
    Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y); //r_dot
  }

  //mean predicted measurement
  z_pred.fill(0.0);
  for(int i = 0; i<2*n_aug_+1; ++i)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation convariance matrix s
  S.fill(0.0);

  for(int i = 0; i < 2*n_aug_+1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i)*z_diff*z_diff.transpose();
  }

  //add measurement noise
  
  S = S + R_radar;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_radar_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; ++i)
  { 
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // create vector for incoming radar measurement
  VectorXd z = VectorXd(n_radar_);
  z <<  meas_package.raw_measurements_(0),
        meas_package.raw_measurements_(1),
        meas_package.raw_measurements_(2);


  //kalman gain k;
  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z - z_pred;

  // angle normalization
  //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  double nis_ = z_diff.transpose() * S.inverse() * z_diff;
  cout << nis_ << " " <<"\n";

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}