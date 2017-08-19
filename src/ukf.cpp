#include "ukf.h"
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
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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
  n_x_ = 5;
  n_aug_ = 7;
P_ = MatrixXd::Identity(n_x_,n_x_);
  lambda_ = 3 - n_aug_;
  weights_ = VectorXd(2*n_aug_+1);
   double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
  iterations = 0;
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
	if(!is_initialized_)
	{

	time_us_ = meas_package.timestamp_;
		if(meas_package.sensor_type_ == MeasurementPackage::RADAR)			{
			x_ = VectorXd(5);
float ro,theta, ro_dot;

      ro = meas_package.raw_measurements_[0];
      theta = meas_package.raw_measurements_[1];
      ro_dot = meas_package.raw_measurements_[2];
      float px,py;
      px = ro*cos(theta);
      py = ro*sin(theta);
			lastpx = px;
		
 
				//x_  << px, py, 1/sin(atan2(py,px))*py, atan2(py,px), 0; 
x_  << px, py, 2, 2, 0; 

		}
		else
		{
			lastpx = meas_package.raw_measurements_[0];
			x_ = VectorXd(5);
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0; 
		}
	
		is_initialized_ = true;
		return;
	}
	float dt = ((float)meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;


	Prediction(dt); 
	if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
	{
		

		UpdateRadar(meas_package); 
	}
	else
	{
		UpdateLidar(meas_package); 	
	}
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
 //set state dimension
  

  //define spreading parameter

  //set example state
 VectorXd x = x_;
  
  MatrixXd P = P_;
  

  
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
 Xsig_aug.fill(0.0);
/*******************************************************************************
 * Student part begin
 ******************************************************************************/
  x_aug.fill(0.0); 
  //create augmented mean state
  x_aug.head(5) = x;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  
  //predict sigma points
  Xsig_pred.fill(0.0);
  Xsig_pred = Xsig_aug.block(0,0,n_x_, 2 * n_aug_ + 1);
  //avoid division by zero
  //write predicted sigma points into right column
  int i = 0;
  for(i =0 ; i<(2*n_aug_+1); i++)
  {
      MatrixXd influence = MatrixXd(n_x_, 1);
      MatrixXd integral = MatrixXd(n_x_, 1);
      float yawraterate = Xsig_aug( n_aug_-1, i);
      float yawk= Xsig_aug( n_aug_-3, i);
      float vak = Xsig_aug(n_aug_-2, i);
      float vk = Xsig_aug(2, i);
      float yaw = Xsig_aug(3, i);
      integral(0, 0) = 0.5*delta_t*delta_t*cos(yaw)*vak;
      integral(1, 0) = 0.5*delta_t*delta_t*sin(yaw)*vak;
      integral(2, 0) = delta_t * vak;
      integral(3, 0) = 0.5*delta_t*delta_t*yawraterate;
      integral(4, 0) = delta_t*yawraterate;
      if(yawk==0)
      {
          
          influence(0,0) = vk*cos(yaw)*delta_t;
          influence(1,0) = vk*sin(yaw)*delta_t;
          influence(2,0) = 0;
          influence(3,0) = yawk*delta_t;
          influence(4,0) = 0;
      }
      else
      {
          
          influence(0,0) = vk/yawk*(sin(yaw+delta_t*yawk)-sin(yaw));
          influence(1,0) = vk/yawk*(-1*cos(yaw+delta_t*yawk)+cos(yaw));
          influence(2,0) = 0;
          influence(3,0) = yawk*delta_t;
          influence(4,0) = 0;
          
           
            
      }
       Xsig_pred.col(i) += influence+integral;
  }
  //create vector for weights
  
  //create vector for predicted state


  //create covariance matrix for prediction



/*******************************************************************************
 * Student part begin
 ******************************************************************************/
for(int i = 0; i < 2 * n_aug_ + 1; i++)
{
while(Xsig_pred(3,i)>M_PI)
{
Xsig_pred(3,i) -= 2*M_PI;
}
while(Xsig_pred(3,i)<-M_PI)
{
Xsig_pred(3,i) += 2*M_PI;
}

}
x.fill(0.0);
  //predicted state mean
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights_(i) * Xsig_pred.col(i);
  }
  //predicted state covariance matrix
 while(x(3) > M_PI)
  {
 x(3) -= 2*M_PI;
  }
  while(x(3) < -M_PI)
  {
  x(3) += 2*M_PI;
  }


  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
  while(x_diff(3) > M_PI)
  {
 x_diff(3) -= 2*M_PI;
  }
  while(x_diff(3) < -M_PI)
  {
  x_diff(3) += 2*M_PI;
  }


    P = P + weights_(i) * x_diff * x_diff.transpose();

cout << "p at " <<i <<" "<< P<<endl;
  }
x_ = x;
while(P(3)>M_PI)
{
 P(3) -= 2*M_PI;

}
while(P(3)<-M_PI)
{
 P(3) += 2*M_PI;
}

P_ = P;
Xsig_pred_ = Xsig_pred;
cout <<"xsig is " << Xsig_pred << endl;
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
int n_z = n_x_ - 2;
VectorXd z = VectorXd(n_z);
z << meas_package.raw_measurements_[0],  meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  Zsig.fill(0.0);
/*******************************************************************************
 * Student part begin
 ******************************************************************************/
cout << Xsig_pred_ << endl;
  //transform sigma points into measurement space
  MatrixXd Noise = MatrixXd(n_z,n_z);
  int i = 0;
  for(i = 0; i<2*n_aug_+1; i++)
  {
  float theta;
  
  theta = atan2(Xsig_pred_(1,i),Xsig_pred_(0,i));
  while(theta > M_PI)
  {
  theta -= 2*M_PI;
  }
  while(theta < -M_PI)
  {
  theta += 2*M_PI;
  }
  float p = sqrt( pow(Xsig_pred_(0,i), 2)+ pow(Xsig_pred_(1,i), 2));    
  
  
  Zsig(0,i)= p;
  Zsig(1,i)= theta;
  Zsig(2,i)= ((Xsig_pred_(0,i)*cos(Xsig_pred_(3,i))*Xsig_pred_(2,i)+Xsig_pred_(1,i)*sin(Xsig_pred_(3,i))*Xsig_pred_(2,i))/p);
  }

z_pred.fill(0.0);
  //calculate mean predicted measurement
    for(i = 0; i<2*n_aug_+1; i++)
  {
z_pred+= weights_(i)*Zsig.col(i);
}
Noise.fill(0.0);


  //radar measurement noise standard deviation radius change in m/s

Noise(0,0) = std_radr_*std_radr_;
Noise(1,1) = std_radphi_*std_radphi_;
Noise(2,2) = std_radrd_*std_radrd_;
  //calculate measurement covariance matrix S
S.fill(0.0);
    for(i = 0; i<2*n_aug_+1; i++)
  {
      S += weights_(i)*(Zsig.col(i).colwise()-z_pred)*(Zsig.col(i).colwise()-z_pred).transpose();
  }
  S += Noise;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  Tc.fill(0.0);
  i = 0;
  for(i = 0; i<2 * n_aug_ + 1; i++)
  {
  //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd Kgain = Tc*S.inverse();
  //update state mean and covariance matrix
  MatrixXd Kz = (Kgain*((-Zsig).colwise()+z_pred));
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
      x_+=Kgain*z_diff;
  
  P_ = (P_ - Kgain*S*Kgain.transpose());
}
