#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd::Zero(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/10;

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
    
  is_initialized_ = false;
  Xsig_pred_ = MatrixXd::Zero(5,15);
  time_us_ = 0.0;
    
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

    if (is_initialized_ == false) {
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_ << meas_package.raw_measurements_(0),
                  meas_package.raw_measurements_(1),
                  0,
                  M_PI/10,
                  M_PI/6;
            
            P_ << (std_laspx_ * std_laspx_), 0, 0, 0, 0,
                  0, (std_laspy_ * std_laspy_), 0, 0, 0,
                  0, 0, 16.0, 0, 0,
                  0, 0, 0, 1, 0,
                  0, 0, 0, 0, 1;
            std::cout << "Initialization with Laser. \n";
            std::cout << "x_ = \n" << x_ << "\n\n";
            std::cout << "P_ = \n" << P_ << "\n\n";
        }
        
        else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            double theta = meas_package.raw_measurements_[1];
            
            theta = fmod(theta, 2*M_PI);
            
            x_ << meas_package.raw_measurements_(0)*cos(theta),
                  meas_package.raw_measurements_(0)*sin(theta),
                  4.0,
                  M_PI/6,
                  M_PI/10;
            
            P_ << 1, 0, 0, 0, 0,
                  0, 1, 0, 0, 0,
                  0, 0, 16.0, 0, 0,
                  0, 0, 0, (M_PI * M_PI / 36), 0,
                  0, 0, 0, 0, (M_PI * M_PI / 49);
            std::cout << "Initialization with Radar. \n";
            std::cout << "x_ = \n" << x_ << "\n\n";
            std::cout << "P_ = \n" << P_ << "\n\n";
        }
        
        else {
            x_ << 0,0,0,0,0;
            
            P_ << 1,0,0,0,0,
                  0,1,0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1;
            std::cout << "Error - the data file has some weird starting value. \n";
        }
        
        // Initialize R once and for all
        R_radar = MatrixXd(3,3);
        R_radar << 0,0,0,0,0,0,0,0,0;
        R_radar(0,0) = std_radr_ * std_radr_;
        R_radar(1,1) = std_radphi_ * std_radphi_;
        R_radar(2,2) = std_radrd_ * std_radrd_;
        
        R_laser = MatrixXd(2,2);
        R_laser << 0,0,0,0;
        R_laser(0,0) = std_laspx_ * std_laspx_;
        R_laser(1,1) = std_laspy_ * std_laspy_;
        
        
        is_initialized_ = true;
        
        time_us_ = meas_package.timestamp_;
        
        return;
    };
    
    if (is_initialized_ == true) {
        double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
        
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            if (use_laser_ == true) {
                Prediction(delta_t);
                UpdateLidar(meas_package);
                time_us_ = meas_package.timestamp_;
            }
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            if (use_radar_ == true) {
                Prediction(delta_t);
                UpdateRadar(meas_package);
                time_us_ = meas_package.timestamp_;
            }
        }
        
        else {
            std::cout << "Error - the data file has some weird starting value. \n";
            return;
        }
        
        return;
        
    }
    
    
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
    
    // First sigma points
    n_x = 5;
    lambda = 3 - n_x;
    
    MatrixXd Xsig = MatrixXd::Zero(n_x, 2*n_x+1);
    MatrixXd A = P_.llt().matrixL();
    VectorXd xsig_0 = x_;
    MatrixXd xsig_1(5,5);
    MatrixXd xsig_2(5,5);
    
    for (int i=0; i<A.cols(); i++) {
        xsig_1.col(i) = xsig_0 + (sqrt(3) * A.col(i));
    };
    
    for (int i=0; i<A.cols(); i++) {
        xsig_2.col(i) = xsig_0 - (sqrt(3) * A.col(i));
    };
    
    
    Xsig.col(0) = xsig_0;
    Xsig.col(1) = xsig_1.col(0);
    
    for (int i=0; i<xsig_1.cols(); i=i+1) {
        Xsig.col(i+1) = xsig_1.col(i);
    };
    
    for (int k=0; k<5; k=k+1) {
        Xsig.col(k+6) = xsig_2.col(k);
    };
    
    
    // Augmentation so that we can account for process noise.
    n_aug = 7;
    lambda = 3 - n_aug;
    
    VectorXd x_aug;
    x_aug = VectorXd(7);
    x_aug = VectorXd::Zero(7);
    MatrixXd P_aug;
    P_aug = MatrixXd(7, 7);
    P_aug = MatrixXd::Zero(7,7);
    MatrixXd Xsig_aug;
    Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
    
    x_aug.head(5) = x_;
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_ * std_a_;
    P_aug(6,6) = std_yawdd_ * std_yawdd_;
    
    
    MatrixXd L = P_aug.llt().matrixL();
    
    
    xsig_0 = x_aug;
    xsig_1 = MatrixXd(7,7);
    xsig_2 = MatrixXd(7,7);
    
    for (int i=0; i<L.cols(); i++) {
        xsig_1.col(i) = xsig_0 + (sqrt(3) * L.col(i));
    };
    
    for (int i=0; i<L.cols(); i++) {
        xsig_2.col(i) = xsig_0 - (sqrt(3) * L.col(i));
    };
    
    
    Xsig_aug.col(0) = xsig_0;
    
    for (int i=0; i<xsig_1.cols(); i=i+1) {
        Xsig_aug.col(i+1) = xsig_1.col(i);
    };
    
    for (int k=0; k<xsig_2.cols(); k=k+1) {
        Xsig_aug.col(k+8) = xsig_2.col(k);
    };
    
    // Prediction using Xsig_aug
    double dt2 = dt * dt;
    
    for (int i=0; i<Xsig_aug.cols(); i=i+1) {
        VectorXd x(7);
        VectorXd xp(5);
        
        x = Xsig_aug.col(i);
        double px = x(0);
        double py = x(1);
        double v = x(2);
        double phi = x(3);
        
        phi = fmod(phi, 2*M_PI);
        
        double phi_dot = x(4);
        double ua = x(5);
        double ub = x(6);
        
        if ((phi_dot < 0.004) && (phi_dot > -0.004)) {
            xp(0) = px + (v * cos(phi) * dt) + (0.5 * dt2 * cos(phi) * ua);
            xp(1) = py + (v * sin(phi) * dt) + (0.5 * dt2 * sin(phi) * ua);
            xp(2) = v + (dt * ua);
            xp(3) = phi + (0.5 * dt2 * ub);
            xp(4) = phi_dot + (dt * ub);
        }
        else {
            xp(0) = px + (v*(sin(phi + phi_dot*dt) - sin(phi))/phi_dot) + (0.5 * dt2 * cos(phi) * ua);
            xp(1) = py + (v*(-cos(phi + phi_dot*dt) + cos(phi))/phi_dot) + (0.5 * dt2 * sin(phi) * ua);
            xp(2) = v + (dt * ua);
            xp(3) = phi + (phi_dot * dt) + (0.5 * dt2 * ub);
            xp(4) = phi_dot + (dt * ub);
        };
        Xsig_pred_.col(i) = xp;
    };
    
    
    //Getting new values of x_ and P_ from Xsig_pred!
    weights_ = VectorXd(2*n_aug+1);
    weights_(0) = lambda / (lambda + n_aug);
    
    for (int i=1; i<2*n_aug+1; i++) {
        weights_(i) = 0.5 / (lambda + n_aug);
    };
    
    int cols = Xsig_pred_.cols();
    x_ = VectorXd::Zero(5);
    P_ = MatrixXd::Zero(5,5);
    for (int i=0; i<cols; i=i+1) {
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    };
    
    for (int i=0; i<cols; i=i+1) {
        P_ = P_ + weights_(i) * (Xsig_pred_.col(i) - x_) * (Xsig_pred_.col(i) - x_).transpose();
    };
    
    //Prediction step finished!
    
    return;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

    n_x = 5;
    n_aug = 7;
    int n_z = 2;
    lambda = 3 - n_aug;
    
    weights_ = VectorXd(2*n_aug+1);
    double weight_0 = lambda/(lambda+n_aug);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug+1; i++) {
        double weight = 0.5/(n_aug+lambda);
        weights_(i) = weight;
    }
    
    // Create sigma point matrix of current X in Z-space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
    VectorXd z_pred = VectorXd::Zero(n_z);
    MatrixXd S = MatrixXd(n_z,n_z);
    
    S << 0,0,0,0;
    
    int cols = Xsig_pred_.cols();
    
    for (int i=0; i<cols; i++) {
        double px = Xsig_pred_.col(i)(0);
        double py = Xsig_pred_.col(i)(1);
        Zsig.col(i) << px, py;
    };
    
    
    for (int j=0; j<cols; j++) {
        z_pred = z_pred + weights_(j)*Zsig.col(j);
    };
    std::cout << "z pred_ laser = \n" << z_pred << "\n\n";
    
    
    for (int k=0; k<cols; k++) {
        S = S + weights_(k) * (Zsig.col(k) - z_pred) * (Zsig.col(k) - z_pred).transpose();
    };
    
    S = S + R_laser;
    
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);
    std::cout << "z laser = \n" << z << "\n\n";
    
    MatrixXd Tc = MatrixXd::Zero(n_x, n_z);
    
    cols = Zsig.cols();
    
    for (int i=0; i<cols; i++) {
        Tc = Tc + (weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose());
    };
    
    
    MatrixXd K = MatrixXd::Zero(5,2);
    K = Tc * S.inverse();
    
    NIS_laser_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
    std::cout << "nis laser = \n" << NIS_laser_ << "\n\n";
    
    //update state mean and covariance matrix
    
    x_ = x_ + (K * (z - z_pred));
    P_ = P_ - (K * S * K.transpose());
    
    return;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    
    n_x = 5;
    n_aug = 7;
    int n_z = 3;
    lambda = 3 - n_aug;
    
    weights_ = VectorXd(2*n_aug+1);
    double weight_0 = lambda/(lambda+n_aug);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug+1; i++) {
        double weight = 0.5/(n_aug+lambda);
        weights_(i) = weight;
    }
    
    // Create sigma point matrix of current X in Z-space
    MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug + 1);
    VectorXd z_pred = VectorXd::Zero(n_z);
    MatrixXd S = MatrixXd::Zero(n_z,n_z);

    int cols = Xsig_pred_.cols();
    
    for (int i=0; i<cols; i++) {
        double px = Xsig_pred_.col(i)(0);
        double py = Xsig_pred_.col(i)(1);
        double v = Xsig_pred_.col(i)(2);
        double phi = Xsig_pred_.col(i)(3);
        
        phi = fmod(phi, 2*M_PI);
        
        double rho = sqrt(px*px + py*py);
        double omega = atan2(py, px);
        
        omega = fmod(omega, 2*M_PI);
        
        if (rho < 0.005) {
            std::cout << "Zero rho Xsig_pred radar = \n" << Xsig_pred_.col(i) << "\n\n";
        }
        
        double rho_dot = v*(px*cos(phi) + py*sin(phi))/rho;
        
        Zsig.col(i) << rho, omega, rho_dot;
    };
    
    for (int j=0; j<cols; j++) {
        z_pred = z_pred + weights_(j)*Zsig.col(j);
    };
    
    z_pred(1) = fmod(z_pred(1), 2*M_PI);
    
    
    // Create S (covariance matrix)
    
    for (int k=0; k<cols; k++) {
        S = S + (weights_(k) * (Zsig.col(k) - z_pred)) * ((Zsig.col(k) - z_pred).transpose());
    };
    S = S + R_radar;
    
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_(0),
         meas_package.raw_measurements_(1),
         meas_package.raw_measurements_(2);
    
    z(1) = fmod(z(1), 2*M_PI);
    
    MatrixXd Tc = MatrixXd::Zero(5,3);
    
    cols = Zsig.cols();
    
    for (int i=0; i<cols; i++) {
        Tc = Tc + (weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose());
    };
    
    // Calculating NIS
    NIS_radar_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
    std::cout << "nis radar = \n" << NIS_radar_ << "\n\n";
    
    //calculate Kalman gain K;
    MatrixXd K;
    K = MatrixXd::Zero(5,3);
    K = Tc * S.inverse();
    
    //update state mean and covariance matrix
    
    x_ = x_ + (K * (z - z_pred));
    P_ = P_ - (K * S * K.transpose());
    
    return;
}
