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
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/12;

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
    
    is_initialized_ = false;
    Xsig_pred_ = MatrixXd(5,15);
    time_us_ = 0.0;
    
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
    if (is_initialized_ == false) {
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 4.0, M_PI/10, M_PI/6;
            
            P_ << (std_laspx_ * std_laspx_), 0, 0, 0, 0,
                  0, (std_laspy_ * std_laspy_), 0, 0, 0,
                  0, 0, 16.0, 0, 0,
                  0, 0, 0, (M_PI * M_PI / 36), 0,
                  0, 0, 0, 0, (M_PI * M_PI / 49);
            //std::cout << "Initialization with Laser. \n";
            //std::cout << "x_ = \n" << x_ << "\n\n";
            //std::cout << "P_ = \n" << P_ << "\n\n";
        }
        
        else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            double theta = meas_package.raw_measurements_[1];
            x_ << meas_package.raw_measurements_(0)*cos(theta), meas_package.raw_measurements_(0)*sin(theta), 4.0, M_PI/10, M_PI/6;
            
            P_ << (std_laspx_ * std_laspx_), 0, 0, 0, 0,
                  0, (std_laspy_ * std_laspy_), 0, 0, 0,
                  0, 0, 16.0, 0, 0,
                  0, 0, 0, (M_PI * M_PI / 36), 0,
                  0, 0, 0, 0, (M_PI * M_PI / 49);
            std::cout << "Initialization with Radar. \n";
            std::cout << "x_ = \n" << x_ << "\n\n";
            std::cout << "P_ = \n" << P_ << "\n\n";
        }
        
        else {
            x_ << 0,0,0,0,0;
            
            P_ << 0,0,0,0,0,
                  0,0,0,0,0,
                  0,0,0,0,0,
                  0,0,0,0,0,
                  0,0,0,0,0;
            std::cout << "Error - the data file has some weird starting value. \n";
        }
        
        is_initialized_ = true;
        
        time_us_ = meas_package.timestamp_;
        
        return;
    };
    
    if (is_initialized_ == true) {
        
        double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
        
        Prediction(delta_t);
        
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            if (use_laser_ == true) {
                UpdateLidar(meas_package);
            }
            else {
                return;
            }
            
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            if (use_radar_ == true) {
                UpdateRadar(meas_package);
            }
            else {
                return;
            }
        }
        else std::cout << "Error - the data file has some weird starting value. \n";
        
        // Debugging purposes
        //std::cout << "x_ = \n" << x_ << "\n\n";
        //std::cout << "P_ = \n" << P_ << "\n\n";
        
        time_us_ = meas_package.timestamp_;
        
        return;
    }
    
    
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    // First sigma points
    int n_x = 5;
    int lambda = 3 - n_x;
    
    MatrixXd Xsig(n_x, 2*n_x+1);
    MatrixXd A = P_.llt().matrixL();
    
    Xsig.col(0)  = x_;
    
    //set remaining sigma points
    for (int i = 0; i < n_x; i++)
    {
        Xsig.col(i+1)     = x_ + sqrt(lambda+n_x) * A.col(i);
        Xsig.col(i+1+n_x) = x_ - sqrt(lambda+n_x) * A.col(i);
    }
    std::cout << "xsig_ = \n" << Xsig << "\n\n";
    // Augmentation so that we can account for process noise.
    int n_aug = 7;
    lambda = 3 - n_aug;
    
    VectorXd x_aug = VectorXd(7);
    MatrixXd P_aug = MatrixXd(7, 7);
    MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
    
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_ * std_a_;
    P_aug(6,6) = std_yawdd_ * std_yawdd_;
    
    MatrixXd L = P_aug.llt().matrixL();
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug) * L.col(i);
        Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * L.col(i);
    }
    std::cout << "xsig aug_ = \n" << Xsig_aug << "\n\n";
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
        double phi_dot = x(4);
        double ua = x(5);
        double ub = x(6);
        
        if (phi < -M_PI) {
            phi = phi + 2*M_PI;
        }
        
        if (phi > M_PI) {
            phi = phi - 2*M_PI;
        }
        
        if (phi_dot < 0.004) {
            xp(0) = px + (v * cos(phi) * dt) + (0.5 * dt2 * cos(phi) * ua);
            xp(1) = py + (v * sin(phi) * dt) + (0.5 * dt2 * sin(phi) * ua);
            xp(2) = v + (dt * ua);
            xp(3) = phi + (0.5 * dt2 * ub);
            xp(4) = (dt * ub);
        }
        else {
            xp(0) = px + v*(sin(phi + phi_dot*dt) - sin(phi))/phi_dot + (0.5 * dt2 * cos(phi) * ua);
            xp(1) = py + v*(-cos(phi + phi_dot*dt) + cos(phi))/phi_dot + (0.5 * dt2 * sin(phi) * ua);
            xp(2) = v + (dt * ua);
            xp(3) = phi + (phi_dot * dt) + (0.5 * dt2 * ub);
            xp(4) = phi_dot + (dt * ub);
        }
        
        Xsig_pred_.col(i)(0) = xp(0);
        Xsig_pred_.col(i)(1) = xp(1);
        Xsig_pred_.col(i)(2) = xp(2);
        Xsig_pred_.col(i)(3) = xp(3);
        Xsig_pred_.col(i)(4) = xp(4);
    };
    std::cout << "xsig pred_ = \n" << Xsig_pred_ << "\n\n";
    
    //Getting new values of x_ and P_ from Xsig_pred!
    VectorXd weights = VectorXd(2*n_aug+1);
    weights(0) = lambda / (lambda + n_aug);
    
    for (int i=1; i<2*n_aug+1; i++) {
        weights(i) = 0.5 / (lambda + n_aug);
    };
    
    int cols = Xsig_pred_.cols();
    for (int i=0; i<cols; i=i+1) {
        x_ = x_ + weights(i) * Xsig_pred_.col(i);
    };
    
    for (int i=0; i<cols; i=i+1) {
        P_ = P_ + weights(i) * (Xsig_pred_.col(i) - x_) * (Xsig_pred_.col(i) - x_).transpose();
    };
    
    std::cout << "x_ = \n" << x_ << "\n\n";
    //Prediction step finished!
    
    return;
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
    int n_x = 5;
    int n_aug = 7;
    int n_z = 2;
    double lambda = 3 - n_aug;
    
    VectorXd weights = VectorXd(2*n_aug+1);
    double weight_0 = lambda/(lambda+n_aug);
    weights(0) = weight_0;
    for (int i=1; i<2*n_aug+1; i++) {
        double weight = 0.5/(n_aug+lambda);
        weights(i) = weight;
    }
    
    // Create sigma point matrix of current X in Z-space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
    VectorXd z_pred = VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z,n_z);
    
    int cols = Xsig_pred_.cols();
    
    for (int i=0; i<cols; i++) {
        double px = Xsig_pred_.col(i)(0);
        double py = Xsig_pred_.col(i)(1);
        Zsig.col(i) << px, py;
    };
    
    std::cout << "zsig laser = \n" << Zsig << "\n\n";
    
    for (int j=0; j<cols; j++) {
        z_pred = z_pred + weights(j)*Zsig.col(j);
    };
    std::cout << "z pred_ laser = \n" << z_pred << "\n\n";
    MatrixXd R(n_z,n_z);
    R(0,0) = std_laspx_ * std_laspx_;
    R(1,1) = std_laspy_ * std_laspy_;
    
    for (int k=0; k<cols; k++) {
        S = S + weights(k) * (Zsig.col(k) - z_pred) * (Zsig.col(k) - z_pred).transpose();
    };
    
    S = S + R;
    std::cout << "S laser = \n" << S << "\n\n";
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);
    
    MatrixXd Tc = MatrixXd(n_x, n_z);
    
    cols = Zsig.cols();
    
    for (int i=0; i<cols; i++) {
        Tc = Tc + (weights(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose());
    };
    
    MatrixXd K(5,2);
    K = Tc * S.inverse();
    
    NIS_laser_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
    
    //update state mean and covariance matrix
    
    x_ = x_ + (K * (z - z_pred));
    P_ = P_ - (K * S * K.transpose());
    std::cout << "x after laser = \n" << x_ << "\n\n";
    
    return;
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
    
    int n_x = 5;
    int n_aug = 7;
    int n_z = 3;
    double lambda = 3 - n_aug;
    
    VectorXd weights = VectorXd(2*n_aug+1);
    double weight_0 = lambda/(lambda+n_aug);
    weights(0) = weight_0;
    for (int i=1; i<2*n_aug+1; i++) {
        double weight = 0.5/(n_aug+lambda);
        weights(i) = weight;
    }
    
    // Create sigma point matrix of current X in Z-space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
    VectorXd z_pred = VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z,n_z);

    int cols = Xsig_pred_.cols();
    
    for (int i=0; i<cols; i++) {
        double px = Xsig_pred_.col(i)(0);
        double py = Xsig_pred_.col(i)(1);
        double v = Xsig_pred_.col(i)(2);
        double phi = Xsig_pred_.col(i)(3);
        /*
        if (phi < -M_PI) {
            phi = phi + 2*M_PI;
        }
        
        if (phi > M_PI) {
            phi = phi - 2*M_PI;
        }
        */
        double rho = sqrt(px*px + py*py);
        double omega = atan2(py, px);
        /*
        if (omega < -M_PI) {
            omega = omega + 2*M_PI;
        }
        
        if (omega > M_PI) {
            omega = omega - 2*M_PI;
        }
        */
        double rho_dot = v*(px*cos(phi) + py*sin(phi))/rho;
        
        Zsig.col(i) << rho, omega, rho_dot;
    };
    
    std::cout << "zsig radar = \n" << Zsig << "\n\n";
    
    for (int j=0; j<cols; j++) {
        z_pred = z_pred + weights(j)*Zsig.col(j);
    };
    
    std::cout << "z pred_ radar = \n" << z_pred << "\n\n";
    
    // Create S (covariance matrix)
    MatrixXd R(3,3);
    R(0,0) = std_radr_ * std_radr_;
    R(1,1) = std_radphi_ * std_radphi_;
    R(2,2) = std_radrd_ * std_radrd_;
    
    for (int k=0; k<cols; k++) {
        S = S + weights(k) * (Zsig.col(k) - z_pred) * ((Zsig.col(k) - z_pred).transpose());
    };
    S = S + R;
    
    std::cout << "S radar = \n" << S << "\n\n";
    
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), meas_package.raw_measurements_(2);
    
    MatrixXd Tc = MatrixXd(n_x, n_z);
    
    cols = Zsig.cols();
    
    for (int i=0; i<cols; i++) {
        Tc = Tc + (weights(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose());
    };
    
    // Calculating NIS
    NIS_radar_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
    
    
    //calculate Kalman gain K;
    MatrixXd K(5,3);
    K = Tc * S.inverse();
    
    //update state mean and covariance matrix
    
    x_ = x_ + (K * (z - z_pred));
    P_ = P_ - (K * S * K.transpose());
    
    std::cout << "x_ radar = \n" << x_ << "\n\n";
    //std::cout << "P_ = \n" << P_ << "\n\n";
    
    return;
}
