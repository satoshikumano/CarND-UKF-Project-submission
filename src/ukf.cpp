#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 10;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 10;

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

  // Parameters above this line are scaffolding, do not modify
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;
  Xsig_pred_ = MatrixXd(this->n_x_, 2 * this->n_aug_ + 1);
  weights_ = VectorXd(2 * this->n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!this->is_initialized_) {
    this->x_.setZero();
    this->P_ << 1, 0, 0, 0, 0,
               0, 1, 0, 0, 0,
               0, 0, 1, 0, 0,
               0, 0, 0, 1, 0,
               0, 0, 0, 0, 1;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "Initialize RADAR:" << endl;
      double ro = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double ro_dot = meas_package.raw_measurements_(2);
      double px = ro * cos(phi);
      double py = ro * sin(phi);
      this->x_(0) = px;
      this->x_(1) = py;
    } else if  (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      cout << "Initialize LASER:" << endl;
      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);
      this->x_(0) = px;
      this->x_(1) = py;
    }
    this->is_initialized_ = true;
    this->previous_timestamp_ = meas_package.timestamp_;
    return;
  }
  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  this->Prediction(dt);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    this->UpdateRadar(meas_package);
  } else if  (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    this->UpdateLidar(meas_package);
  }
  this->previous_timestamp_ = meas_package.timestamp_;
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
}

void UKF::PredictSigmaPoints(double delta_t) {
  // Predict sigma points
  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7, 7);
  MatrixXd Xsig_aug = MatrixXd(this->n_aug_, 2 * this->n_aug_ + 1);

  x_aug.head(5) = this->x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  MatrixXd Q = MatrixXd(2,2);
  Q << this->std_a_*this->std_a_, 0,
       0, this->std_yawdd_*this->std_yawdd_;
  P_aug.setZero();
  P_aug.topLeftCorner(5,5) = this->P_;
  P_aug.bottomRightCorner(2,2) = Q;

  MatrixXd P_aug_sqrt = P_aug.llt().matrixL();

  Xsig_aug.col(0) = this->x_;

  MatrixXd cal = sqrt(this->lambda_ + this->n_aug_) * P_aug_sqrt;

  for (int i=0; i<14; ++i) {
    if (i < 7) {
      Xsig_aug.col(i+1) = x_aug + cal.col(i);
    } else {
      Xsig_aug.col(i+1) = x_aug - cal.col(i-7);
    }
  }

  double delta_t_2 = delta_t * delta_t;
  double delta_t_2_div2 = delta_t_2 / 2.0;
  for (int i=0; i<15; ++i) {
    VectorXd col = Xsig_aug.col(i);
    double nu_a = col(5);
    double nu_si_dot_dot = col(6);
    VectorXd xk = col.head(5);
    double v = xk(2);
    double si = xk(3);
    double si_dot = xk(4);
    VectorXd nu_vector_t = VectorXd(5);
    nu_vector_t(0) = delta_t_2_div2 * cos(si) * nu_a;
    nu_vector_t(1) = delta_t_2_div2 * sin(si) * nu_a;
    nu_vector_t(2) = delta_t * nu_a;
    nu_vector_t(3) = delta_t_2_div2 * nu_si_dot_dot;
    nu_vector_t(4) = delta_t * nu_si_dot_dot;

    VectorXd xk_plus1 = VectorXd(5);
    xk_plus1 << 0,0,0,0,0;
    VectorXd pred = VectorXd(5);
    if (si_dot != 0) {
      pred(0) = (v / si_dot) * (sin(si + si_dot * delta_t) - sin(si));
      pred(1) = (v / si_dot) * (-cos(si + si_dot * delta_t) + cos(si));
      pred(2) = 0;
      pred(3) = si_dot * delta_t;
      pred(4) = 0;
    } else {
      pred(0) = v * cos(si) * delta_t;
      pred(1) = v * sin(si) * delta_t;
      pred(2) = 0;
      pred(3) = 0;
      pred(4) = 0;
    }
    xk_plus1 = xk + pred + nu_vector_t;
    this->Xsig_pred_.col(i) = xk_plus1;
   }
}

void UKF::PredictMeanAndCovariance() {
  //set weights
  this->weights_(0) = this->lambda_ / (this->lambda_ + this->n_aug_);
  this->weights_.segment(1, 2 * this->n_aug_).fill(0.5 / (this->lambda_ + this->n_aug_));

  for (int i=0; i<2*this->n_aug_+1; ++i) {
    this->x_ = this->x_ + (this->weights_(i) * this->Xsig_pred_.col(i));
  }
  for (int i=0; i<2*this->n_aug_+1; ++i) {
    VectorXd cal = this->Xsig_pred_.col(i) - this->x_;
    this->P_ = this->P_ + this->weights_(i) * (cal * cal.transpose());
  }
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {
  // TODO: implement it.
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
}
