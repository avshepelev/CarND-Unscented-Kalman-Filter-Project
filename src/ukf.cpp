#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  is_initialized_ = false;

  previous_timestamp_ = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2; //30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2; //30;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;

  for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
    weights_(i) = 0.5 / (n_aug_ + lambda_);
  }

  ///* the current NIS for radar
  NIS_radar_;

  ///* the current NIS for laser
  NIS_laser_;

  //MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  //Xsig_aug_.fill(0.0);
  MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  //cout << "ProcessMeasurement: " << endl;
  /*****************************************************************************
  *  Initialization
  ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    //cout << "UKF: " << endl;
    //x_ = VectorXd(5);
    x_ << 0, 0, 0, 0, 0;

    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = meas_package.raw_measurements_(0); //Range
      float phi = meas_package.raw_measurements_(1); //Bearing (relative heading)
      float p_x = rho * cos(phi); //x coordinates along vehicle longitudinal axis
      float p_y = rho * sin(phi);

      x_ << p_x, p_y, 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*if (use_radar_ == false && meas_package.sensor_type_ == MeasurementPackage::RADAR) 
    return;
  if (use_laser_ == false && meas_package.sensor_type_ == MeasurementPackage::LASER) 
    return;*/

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // Compute elapsed time delta_t
  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  //if ( dt > 0.001 ) {
    Prediction(dt);
  //}

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_) {
    UpdateRadar(meas_package);
  } else if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_) {
    UpdateLidar(meas_package);
  } else {
    return;
  }

  // print the output
  //cout << "x_ = " << x_ << endl << endl;
  //cout << "P_ = " << P_ << endl << endl;

  //cout << "end ProcessMeasurement: " << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //cout << "Prediction" << endl;

  //cout << n_aug_ << endl;

  //Predict sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);
  AugmentedSigmaPoints(Xsig_aug, P_, x_);
  //std::cout << "In Prediction Xsig_aug = " << std::endl << Xsig_aug << std::endl;
  SigmaPointPrediction(&Xsig_pred_, Xsig_aug, delta_t);
  //std::cout << "In Prediction Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl;

  //Predict the state, and the state covariance matrix
  PredictMeanAndCovariance(x_, P_);
  //cout << "end Prediction" << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //cout << "UpdateLidar: " << endl;

  int n_z = 2;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //create vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //create matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z, n_z);

  PredictLidarMeasurement(z_pred, S, Zsig, Xsig_pred_);
  UpdateState(x_, P_, meas_package.raw_measurements_, z_pred, Zsig, S, Xsig_pred_);

  //cout << "end UpdateLidar: " << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //cout << "UpdateRadar: " << endl;

  int n_z = 3;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //create vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //create matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z, n_z);

  PredictRadarMeasurement(z_pred, S, Zsig, Xsig_pred_);
  UpdateState(x_, P_, meas_package.raw_measurements_, z_pred, Zsig, S, Xsig_pred_);

  //cout << "end UpdateRadar: " << endl;
}

void UKF::PredictLidarMeasurement(VectorXd &z_pred, MatrixXd &S, MatrixXd &Zsig, MatrixXd &Xsig_pred) {

  //cout << "PredictLidarMeasurement: " << endl;

  int n_z = 2;

  //cout << "Xsig_pred = " << Xsig_pred << endl << endl;
  //cout << "Zsig = " << Zsig << endl << endl;

  //create matrix for sigma points in measurement space
  //MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    Zsig(0, i) = p_x;  
    Zsig(1, i) = p_y;
  }

  //cout << "PredictLidarMeasurement 2" << endl;

  //mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //cout << "PredictLidarMeasurement 3" << endl;

  //measurement covariance matrix S
  //MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //cout << "PredictLidarMeasurement 4" << endl;
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,                
          0,                 std_laspy_*std_laspy_;
  S = S + R;

  //print result
  /*std::cout << "Zsig: " << std::endl << Zsig << std::endl;
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;*/

  //write result
  //*z_out = z_pred;
  //*S_out = S;

  //cout << "end PredictLidarMeasurement: " << endl;
}

void UKF::PredictRadarMeasurement(VectorXd &z_pred, MatrixXd &S, MatrixXd &Zsig, MatrixXd &Xsig_pred) {

  //cout << "PredictRadarMeasurement: " << endl;

  int n_z = 3;

  //create matrix for sigma points in measurement space
  //MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);          //r
    Zsig(1, i) = atan2(p_y, p_x);                      //phi

    if (Zsig(0, i) > 0.001) {
      Zsig(2, i) = (p_x * v1 + p_y * v2) / Zsig(0, i);      //r_dot
    } else {
      Zsig(2, i) = 0;
    } 
  }

  //mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  //MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0,                     0,
          0,                 std_radphi_*std_radphi_, 0,
          0,                 0,                     std_radrd_*std_radrd_;
  S = S + R;

  //print result
  /*std::cout << "Zsig: " << std::endl << Zsig << std::endl;
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;*/

  //write result
  //*z_out = z_pred;
  //*S_out = S;

  //cout << "end PredictRadarMeasurement: " << endl;
}

void UKF::GenerateSigmaPoints(MatrixXd &Xsig) {

  //cout << "GenerateSigmaPoints: " << endl;

  //create sigma point matrix
  //MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
  Xsig.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda_ + n_x_) * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
  }

  //print result
  //std::cout << "Xsig = " << std::endl << Xsig << std::endl;

  //write result
  //*Xsig_out = Xsig;

  //cout << "GenerateSigmaPoints: " << endl;
}

void UKF::AugmentedSigmaPoints(MatrixXd &Xsig_aug, MatrixXd &P, VectorXd &x) {

  //cout << "AugmentedSigmaPoints" << endl;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  //MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //cout << L << endl << endl;
  ////cout << x_aug << endl << endl;
  //cout << Xsig_aug << endl << endl;

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i + 1)          = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  
  //print result
  //std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  //write result
  //*Xsig_out = Xsig_aug;

  //cout << "AugmentedSigmaPoints end" << endl;
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, MatrixXd &Xsig_aug, double delta_t) {

  //cout << "SigmaPointPrediction" << endl;
  //std::cout << "In SigmaPointPrediction delta_t = " << std::endl << delta_t << std::endl;
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred.fill(0.0);

  //std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
  //std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  //create example sigma point matrix
  //MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  //create matrix with predicted sigma points as columns
  //MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  //double delta_t = 0.1; //time diff in sec

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd * delta_t) );
    }
    else {
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0, i) = px_p;
    Xsig_pred(1, i) = py_p;
    Xsig_pred(2, i) = v_p;
    Xsig_pred(3, i) = yaw_p;
    Xsig_pred(4, i) = yawd_p;

    //*Xsig_out = Xsig_pred;

  }

  //print result
  //std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  //write result
  *Xsig_out = Xsig_pred;
  //cout << "end SigmaPointPrediction" << endl;
}

void UKF::PredictMeanAndCovariance(VectorXd &x, MatrixXd &P) {

  //cout << "PredictMeanAndCovariance" << endl;

  /*std::cout << "State" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Covariance matrix" << std::endl;
  std::cout << P << std::endl;
  std::cout << "Weights" << std::endl;
  std::cout << weights_ << std::endl;*/

  //create vector for weights
  //VectorXd weights = VectorXd(2 * n_aug_ + 1);
  
  //create vector for predicted state
  //VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  //MatrixXd P = MatrixXd(n_x_, n_x_);

  // set weights
  /*double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }*/

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2. * M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //print result
  /*std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;*/

  //write result
  //*x_out = x;
  //*P_out = P;

  //cout << "end PredictMeanAndCovariance" << endl;
}

void UKF::UpdateState(VectorXd &x, MatrixXd &P, VectorXd &z, VectorXd &z_pred, MatrixXd &Zsig, MatrixXd &S, MatrixXd &Xsig_pred) {

  //cout << "UpdateState" << endl;

  int n_x = x.rows();
  int n_z = z_pred.rows();

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x = x + K * z_diff;
  P = P - K*S*K.transpose();

  //print result
  /*std::cout << "Updated state x: " << std::endl << x << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P << std::endl;*/

  //write result
  //*x_out = x;
  //*P_out = P;

  //cout << "end UpdateState" << endl;
}
