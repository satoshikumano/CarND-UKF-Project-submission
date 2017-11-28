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
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  weights_ = VectorXd(2 * n_aug_ + 1);
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
  if (!is_initialized_) {
    x_.setZero();
    P_ << 1, 0, 0, 0, 0,
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
      x_(0) = px;
      x_(1) = py;
    } else if  (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      cout << "Initialize LASER:" << endl;
      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);
      x_(0) = px;
      x_(1) = py;
    }
    is_initialized_ = true;
    previous_timestamp_ = meas_package.timestamp_;
    return;
  }
  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  Prediction(dt);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if  (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
  previous_timestamp_ = meas_package.timestamp_;
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
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  MatrixXd Q = MatrixXd(2,2);
  Q << std_a_*std_a_, 0,
       0, std_yawdd_*std_yawdd_;
  P_aug.setZero();
  P_aug.topLeftCorner(5,5) = P_;
  P_aug.bottomRightCorner(2,2) = Q;

  MatrixXd P_aug_sqrt = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_;

  MatrixXd cal = sqrt(lambda_ + n_aug_) * P_aug_sqrt;

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
    Xsig_pred_.col(i) = xk_plus1;
   }
}

void UKF::PredictMeanAndCovariance() {
  //set weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  weights_.segment(1, 2 * n_aug_).fill(0.5 / (lambda_ + n_aug_));

  for (int i=0; i<2*n_aug_+1; ++i) {
    x_ = x_ + (weights_(i) * Xsig_pred_.col(i));
  }
  for (int i=0; i<2*n_aug_+1; ++i) {
    VectorXd cal = Xsig_pred_.col(i) - x_;
    P_ = P_ + weights_(i) * (cal * cal.transpose());
  }
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  z_pred_radar_ = VectorXd(n_z);
  S_radar_ = MatrixXd(n_z, n_z);
  for (int i=0; i<=n_aug_*2; ++i) {
    VectorXd col = Xsig_pred_.col(i);
    double px = col(0);
    double py = col(1);
    double v = col(2);
    double psi = col(3);
    double psi_dot = col(4);
   
    double eps = 0.000001;
    if (fabs(px) < eps) {
      px = eps;
    }
    if (fabs(py)< eps) {
      py = eps;
    } 
    double ro = sqrt(px*px + py*py);
    double phi = atan2(py, px);
    double ro_dot = (px * cos(psi) * v + py * sin(psi) * v) / ro;
    
    Zsig.col(i) << ro, phi, ro_dot;
  }
  //calculate mean predicted measurement
  z_pred_radar_ << 0, 0, 0;
  for (int i=0; i<= n_aug_*2; ++i) {
    z_pred_radar_ += weights_(i) * Zsig.col(i);
  }
  //calculate measurement covariance matrix S
  S_radar_.setZero();
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  for (int i=0; i<= n_aug_*2; ++i) {
    VectorXd cal = Zsig.col(i) - z_pred_radar_;
    S_radar_ += weights_(i) * (cal * cal.transpose());
  }
  S_radar_ += R;
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
