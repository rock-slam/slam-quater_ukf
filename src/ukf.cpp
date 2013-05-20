/**\file ukf.cpp
 *
 * This class has the primitive methods for an Unscented Kalman Filter implementation
 * for an Attitude and Heading Reference System - AHRS. The filter is Quaternion
 * based using accelerometers, gyroscopes and magnetometers (only when set to true).
 * The filter performs the prediction step based on the gyroscopes and therefore quaternion integration.
 * The measurement is formed in a way tht multiple (N) vector measurements can be concatenated.
 * 
 * This Unscented Kalman filter is based on the paper:  John L. Crassidis and F. Landis Markley
 * "Unscented Filtering for Spacecraft Attitude Estimation"
 * A copy if the manuscript can be found in the /doc folder of the library.
 * 
 * @author Javier Hidalgo Carrio | DFKI RIC Bremen | javier.hidalgo_carrio@dfki.de
 * @date March 2012.
 * @version 1.0.
 */

#include <iostream> /**< IO C++ Standard library */
#include <algorithm> /**< Algorithm C++ Standard library */
#include <Eigen/LU> /**< Lineal algebra of Eigen */
#include <Eigen/SVD> /**< Singular Value Decomposition (SVD) of Eigen */
#include <Eigen/Cholesky> /**< Cholesky module **/
#include "ukf.hpp" /**< Unscented Kalman Filter */

/** WGS-84 ellipsoid constants (Nominal Gravity Model and Earth angular velocity) **/
#ifndef Re
#define Re	6378137 /**< Equatorial radius in meters **/
#endif
#ifndef Rp
#define Rp	6378137 /**< Polar radius in meters **/
#endif
#ifndef ECC
#define ECC  0.0818191908426 /**< First eccentricity **/
#endif
#ifndef GRAVITY
#define GRAVITY 9.79766542 /**< Mean value of gravity value in m/s^2 **/
#endif
#ifndef GWGS0
#define GWGS0 9.7803267714 /**< Gravity value at the equator in m/s^2 **/
#endif
#ifndef GWGS1
#define GWGS1 0.00193185138639 /**< Gravity formula constant **/
#endif
#ifndef EARTHW
#define EARTHW  7.292115e-05 /**< Earth angular velocity in rad/s **/
#endif



namespace filter
{
    /** Namespaces to use **/
    using namespace Eigen;
    using namespace std;
    
    /** Indirect Kalman Filter methods **/
 
    /**
    * @brief This function Initilize the vectors and matrix of the UKF   
    */
    void ukf::Init(Matrix <double,ukf::UKFSTATEVECTORSIZE,1> *x_0, Eigen::Matrix< double, ukf::UKFSTATEVECTORSIZE , ukf::UKFSTATEVECTORSIZE  >* P_0, Eigen::Matrix< double, ukf::UKFSTATEVECTORSIZE , ukf::UKFSTATEVECTORSIZE  > *Q, Eigen::Matrix< double, ukf::NUMAXIS , ukf::NUMAXIS  > *R, Eigen::Quaternion <double> *at_q, double a, double f, double lambda, double g)
    {
	
      /** Gravitation acceleration **/
      gtilde << 0, 0, g;
      
      /** Set the parameters **/
      this->a = a;
      this->f = f;
      this->lambda = lambda;
      
      /** Set the state vector **/
      this->x = *x_0;
      
      /** Set the matrices **/
      this->Px = *P_0;
      this->Q = *Q;
      this->R = *R;
      
      /** Set the initial attitude **/
      this->at_q = *at_q;
      
      /** Print values **/
//       std::cout<<"a: "<<a<<"\n";
//       std::cout<<"f: "<<f<<"\n";
//       std::cout<<"lambda: "<<lambda<<"\n";
//       std::cout<<"x:\n"<<this->x<<"\n";
//       std::cout<<"at_q:\n"<<this->at_q.x()<<" "<<this->at_q.y()<<" "<<this->at_q.z()<<" "<<this->at_q.w()<<"\n";
//       std::cout<<"Px:\n"<<Px<<"\n";
//       std::cout<<"Q:\n"<<this->Q<<"\n";
//       std::cout<<"R:\n"<<this->R<<"\n";
      
      return;
    }
    
    /**
    * @brief This function Initilize Attitude
    */
    bool ukf::setAttitude(Eigen::Quaternion< double > *initq)
    {
      if (initq != NULL)
      {
	/** Initial orientation **/
	 this->at_q = (*initq);
	
	return true;
      }
      
      return false;
    }
    
    /**
    * @brief Gets the current orientation in Euler angles (rad)
    */
    Eigen::Matrix< double, ukf::NUMAXIS , 1  > ukf::getEuler()
    {
      Eigen::Matrix <double, ukf::NUMAXIS, 1> euler;
      
      Vector3d e = Eigen::Matrix3d(at_q).eulerAngles(2,1,0);
       euler(0) = e[2]; 
       euler(1) = e[1]; 
       euler(2) = e[0]; 
      
      return euler;
    }
    
    /**
    * @brief Gets the current orientation in Quaternion
    */
    Eigen::Quaternion< double > ukf::getAttitude()
    {
      return at_q;
    }

    /**
    * @brief Gets the current state vector of the filter
    */
    Eigen::Matrix< double, ukf::UKFSTATEVECTORSIZE , 1  > ukf::getState()
    {
      return x;

    }
    
    /**
    * @brief Gets Noise covariance matrix
    */
    Eigen::Matrix< double, ukf::UKFSTATEVECTORSIZE , ukf::UKFSTATEVECTORSIZE> ukf::getCovariance()
    {
	return Px;
    }
    
    /**
    * @brief This computes the theoretical gravity value according to the WGS-84 ellipsoid earth model.
    */
    double ukf::GravityModel(double latitude, double altitude)
    {
      double g; /**< g magnitude at zero altitude **/

      /** Nominal Gravity model **/
      g = GWGS0*((1+GWGS1*pow(sin(latitude),2))/sqrt(1-pow(ECC,2)*pow(sin(latitude),2)));

      /** Gravity affects by the altitude (aprox the value r = Re **/
      g = g*pow(Re/(Re+altitude), 2);

      std::cout<<"Theoretical gravity for this location (WGS-84 ellipsoid model): "<< g<<" [m/s^2]\n";

      return g;

    }
    
    /**
    * @brief Substract the Earth rotation from the gyroscopes readout
    */
    void ukf::SubstractEarthRotation(Eigen::Matrix <double, ukf::NUMAXIS, 1> *u, Eigen::Quaternion <double> *qb_g, double latitude)
    {
      Eigen::Matrix <double, ukf::NUMAXIS, 1> v (EARTHW*cos(latitude), 0, EARTHW*sin(latitude)); /**< vector of earth rotation components expressed in the geografic frame according to the latitude **/

      /** Compute the v vector expressed in the body frame **/
      v = (*qb_g) * v;

      /** Subtract the earth rotation to the vector of inputs (u = u-v**/
      (*u)  = (*u) - v;
      
      return;
    }
    
    /**
     * @brief Discrete-time quaternion kinematic equation
     */
    void ukf::Omega(Eigen::Quaternion< double >* quat, Eigen::Matrix< double, ukf::NUMAXIS , 1  >* angvelo, double dt)
    {

	register int j, l;
	double auxvar;
	Eigen::Matrix <double, ukf::NUMAXIS, ukf::NUMAXIS> Crossproduct;
	Eigen::Matrix <double, ukf::QUATERSIZE, ukf::QUATERSIZE> Omega;
	Eigen::Matrix< double, ukf::QUATERSIZE , 1  > q;
	Eigen::Matrix <double, ukf::NUMAXIS, ukf::NUMAXIS> Upper;
	Eigen::Matrix <double, ukf::NUMAXIS, 1> psi;
	
	Crossproduct = Eigen::Matrix <double, ukf::NUMAXIS, ukf::NUMAXIS>::Zero();
	Omega = Eigen::Matrix <double, ukf::QUATERSIZE, ukf::QUATERSIZE>::Zero();

	/** If angular velocity is not zero **/
	if ((double)(*angvelo).norm() != 0.00)
	{
	    /** Copy the quaternion **/
	    q(3) = quat->w();
	    q(0) = quat->x();
	    q(1) = quat->y();
	    q(2) = quat->z();

	    /** Psi vector calculation **/
	    auxvar = sin (0.5*(*angvelo).norm()*dt);
	    psi(0) = (double) auxvar * ((*angvelo)(0)/(*angvelo).norm());
	    psi(1) = (double) auxvar * ((*angvelo)(1)/(*angvelo).norm());
	    psi(2) = (double) auxvar * ((*angvelo)(2)/(*angvelo).norm());

	    /** Create the cross-product matrix from the angular velocity **/
	    Crossproduct (0,1) = -psi(2);
	    Crossproduct (1,0) = psi(2);
	    Crossproduct (0,2) = psi(1);
	    Crossproduct (2,0) = -psi(1);
	    Crossproduct (1,2) = -psi(0);
	    Crossproduct (2,1) = psi(0);

	    /** Upper matrix **/
	    Upper = (double)cos(0.5*(*angvelo).norm()*(dt)) * Eigen::Matrix <double, ukf::NUMAXIS, ukf::NUMAXIS>::Identity() - Crossproduct;

	    /** Create the omega transition matrix **/
	    for (j=0; j<ukf::QUATERSIZE; j++) /** Columns **/
	    {
		for (l=0; l<ukf::QUATERSIZE; l++) /** Rows **/
		{
		if (l<3)
		{
		    if (j<3)
		    Omega (l,j) = Upper (l,j);
		    else
		    Omega (l,j) = psi (l);
		}
		else
		{
		    if (j<3)
		    Omega (l,j) = -psi (j);
		    else
		    Omega (l,j) = (double)cos(0.5*(*angvelo).norm()*(dt));
		}
		}
	    }

	    /** Propagate forward in time, y = (A) x + y **/
	    q = Omega*q;

	    /** Store the update quaternion in the argument quaternion **/
	    quat->w() = q(3);
	    quat->x() = q(0);
	    quat->y() = q(1);
	    quat->z() = q(2);

	}

	return;
    }
    
    /**
    * @brief Performs the prediction step of the filter.
    */
    void ukf::predict(Eigen::Matrix< double, ukf::NUMAXIS , 1  >* u, double dt)
    {
	register int i;
	double q4; /**< scalar part of a quaternion **/
	Eigen::Matrix <double, ukf::QUATERSIZE-1, 1> vectorq; /**< vectorial part of a quaternion **/
	Eigen::Matrix <double, ukf::UKFSTATEVECTORSIZE, ukf::UKFSTATEVECTORSIZE> M;
	Eigen::Quaternion <double> auxq; /**< Auxiliar quaternion for operations **/
	Eigen::Matrix <double, ukf::NUMAXIS, 1> p_sig_point; /**< Vector containing the rodrigues parameters , upper part of the state vector and the sigma points  **/
	Eigen::Matrix <double, ukf::NUMAXIS, 1> u_plus; /**< Vector of corrected angular velocity **/
	Eigen::Matrix <double, ukf::UKFSTATEVECTORSIZE, 1> sumvar = Eigen::Matrix <double, ukf::UKFSTATEVECTORSIZE, 1>::Zero(); /**< Summation variable **/
	Eigen::Matrix <double, ukf::UKFSTATEVECTORSIZE, ukf::UKFSTATEVECTORSIZE> sumM = Eigen::Matrix <double, ukf::UKFSTATEVECTORSIZE, ukf::UKFSTATEVECTORSIZE>::Zero(); /**< Summation matrix **/
	
	/** Compute the sigma points **/
// 	std::cout<<"(Px + Q):\n"<<(Px + Q)<<"\n (n+lambda)"<<(ukf::UKFSTATEVECTORSIZE+lambda)<<"\n";
	M = (ukf::UKFSTATEVECTORSIZE+lambda)*(Px + Q);
// 	std::cout<<"M:\n"<<M<<"\n";
	Eigen::LLT<Eigen::MatrixXd> lltOfM(M);
	M = lltOfM.matrixL(); //square root of a semi-definited positive matrix 
	
//  	std::cout<<"Sqrt(M):\n"<<M<<"\n";

	sig_point.col(0) = x;
	
	for (i=1; i<=ukf::UKFSTATEVECTORSIZE; i++)
	{
	    sig_point.col(i) = x + M.col(i-1);
	}
	
	for (i=(ukf::UKFSTATEVECTORSIZE+1); i<ukf::SIGPOINTSIZE; i++)
	{
	    sig_point.col(i) = x - M.col(i-(ukf::UKFSTATEVECTORSIZE+1));
	}
	
	/** Calculate Error Quaternion **/
	auxq.setIdentity();
	
	
	e_q[0] = auxq; //Error quaternion (0) is the identity quaternion.
	
	for (i=1; i<ukf::SIGPOINTSIZE; i++)
	{
	    /** Calculate **/
	    p_sig_point = (sig_point.col(i)).block<ukf::NUMAXIS,1>(0,0);
	    
	    q4 = ((-a * p_sig_point.squaredNorm())+ f * sqrt(pow(f,2)+(1-pow(a,2))*p_sig_point.squaredNorm()))/(pow(f,2)+p_sig_point.squaredNorm());
	    vectorq = (1/f)*(a + q4) * p_sig_point;
	    
	    
	    /** Store in the quaternion **/
	    /** Note the order of the arguments in a quaternion:
	    * the real w coefficient first, while internally the
	    * coefficients are stored in the following order: [x, y, z, w] **/
	    e_q[i] = Eigen::Quaternion <double> (q4, vectorq[0], vectorq[1], vectorq[2]);
	    
//  	    std::cout<<"Error Quaternion["<<i<<"]:\n"<<e_q[i].x()<<" "<<e_q[i].y()<<" "<<e_q[i].z()<<" "<<e_q[i].w()<<"\n";
//  	    e_q.block<ukf::NUMAXIS,ukf::NUMAXIS>(0,0).col(i) = vectorq;
//  	    e_q.col(i)[3] = q4;
	}
	
	/** Compute sigma point quaternions from error quaternions **/
 	sig_q[0] = at_q;
	
	for (i=1;i<ukf::SIGPOINTSIZE;i++)
	{
	    sig_q[i] = e_q[i]*at_q;
	}
	

	/** Propagate quaternions forward (sigma point quaternions) u is the vector with the angular velocities **/
	for (i=0; i<ukf::SIGPOINTSIZE;i++)
	{
	    /** The estimated angular velocities are given by substracting the bias **/
 	    u_plus = *u - (sig_point.col(i)).block<ukf::NUMAXIS,1>(3,0);
	    
	    /** Attitude quaternion dynamic matrix Omega **/
	    Omega(&(sig_q[i]), &(u_plus), dt);
	    
//   	    std::cout<<"Sigma Quaternion (after)["<<i<<"]:\n"<<sig_q[i].x()<<" "<<sig_q[i].y()<<" "<<sig_q[i].z()<<" "<<sig_q[i].w()<<"\n";
	    
	    /** Store the inverse value of the quaternion when i = 0 **/
	    if (i==0)
	    {
		auxq = sig_q[i].inverse(); 
	    }
	    
	    /** Propagate error quaternions **/
	    e_q[i] = sig_q[i] * auxq;
	    
//  	    std::cout<<"Error Quaternion (after)["<<i<<"]:\n"<<e_q[i].x()<<" "<<e_q[i].y()<<" "<<e_q[i].z()<<" "<<e_q[i].w()<<"\n";
// 	    std::cout<<"************************\n";
	}
	
	
	
	/** Compute the sigma points **/
	sig_point.col(0).block<ukf::NUMAXIS,1>(0,0) << 0.00, 0.00, 0.00;
	
	for (i=1; i<ukf::SIGPOINTSIZE; i++)
	{
	    /** Vectorial part of a quaternion **/
 	    vectorq << e_q[i].x(), e_q[i].y(), e_q[i].z();
	    sig_point.col(i).block<ukf::NUMAXIS,1>(0,0) = f*(vectorq/(a+e_q[i].w()));
	}
	
	/** Predicted mean **/
	for (i=1; i<ukf::SIGPOINTSIZE; i++)
	{
	    sumvar += sig_point.col(i);
	}
	sumvar = (lambda * sig_point.col(0))+(0.5*sumvar);
	x = (1/(ukf::UKFSTATEVECTORSIZE + lambda)) * sumvar;
	
	
	/** Predicted covariance **/
	for (i=1; i<ukf::SIGPOINTSIZE; i++)
	{
	    sumM += ((sig_point.col(i)-x) * (sig_point.col(i)-x).transpose());
	}
	
	Px = (1/(ukf::UKFSTATEVECTORSIZE + lambda))* (lambda*((sig_point.col(0)-x) * (sig_point.col(0)-x).transpose()) + (0.5 * sumM)) + Q;
	
	
	return;
    }
    
    /**
    * @brief Performs the update (measurement and correction) step of the filter.
    */
    void ukf::update(Eigen::Matrix< double, ukf::NUMAXIS , 1  >* acc, Eigen::Matrix< double, ukf::NUMAXIS , 1  >* mag)
    {
	register int i;
	double q4; /**< scalar part of a quaternion **/
	Eigen::Matrix <double, ukf::NUMAXIS, 1> r1, euler;
	Eigen::Matrix <double, ukf::QUATERSIZE-1, 1> vectorq; /**< vectorial part of a quaternion **/
	Eigen::Matrix <double, ukf::NUMAXIS, 1> sumvar = Eigen::Matrix <double, ukf::NUMAXIS, 1>::Zero(); /**< Summation variable **/
	Eigen::Matrix <double, ukf::NUMAXIS, ukf::NUMAXIS> sumM = Eigen::Matrix <double, ukf::NUMAXIS, ukf::NUMAXIS>::Zero(); /**< Summation matrix **/
	Eigen::Matrix <double, ukf::UKFSTATEVECTORSIZE, ukf::NUMAXIS> sumPxz = Eigen::Matrix <double, ukf::UKFSTATEVECTORSIZE, ukf::NUMAXIS>::Zero(); /**< Summation matrix **/
	Eigen::Quaternion <double> rotation;
	
	
	/** Reference vector **/
	r1 << 1,1,1;
	
	/** Compute observation from the state prediction (from the propagated sigma point quaternion)) **/
	for (i=0;i<ukf::SIGPOINTSIZE;i++)
	{
	    gamma.col(i) = sig_q[i]*r1;
	    
	}
	
	/** Predicted mean of the observation **/
	for (i=1; i<ukf::SIGPOINTSIZE; i++)
	{
	    sumvar += gamma.col(i);
	}
	sumvar = (lambda * gamma.col(0))+(0.5*sumvar);
	z_e = (1/(ukf::UKFSTATEVECTORSIZE + lambda)) * sumvar;
	
// 	std::cout<<"z_e\n"<<z_e<<"\n";
	
	/** Predicted covariance of the observation **/
	for (i=1; i<ukf::SIGPOINTSIZE; i++)
	{
	    sumM += ((gamma.col(i)-z_e) * (gamma.col(i)-z_e).transpose());
	}
	
	Pzz = (1/(ukf::UKFSTATEVECTORSIZE + lambda)) * (lambda*((gamma.col(0)-z_e) * (gamma.col(0)-z_e).transpose()) + (0.5 * sumM));
	
	
	/** Update using Accelerometers **/
	euler[0] = (double) asin((double)(*acc)[1]/ (double)acc->norm()); // Roll
	euler[1] = (double)-atan((*acc)[0]/(*acc)[2]); //Pitch
	euler[2] = (double) sig_q[0].toRotationMatrix().eulerAngles(2,1,0)[0];//Yaw
	
// 	std::cout<<"Measurement euler angles: "<<euler(0)*R2D<<" "<<euler(1)*R2D<<" "<<euler(2)*R2D<<"\n";
	
	/** Compute quaternion from euler angles **/
	rotation = Eigen::Quaternion <double> (Eigen::AngleAxisd(euler[2], Eigen::Vector3d::UnitZ())*
 			    Eigen::AngleAxisd(euler[1], Eigen::Vector3d::UnitY()) *
 			    Eigen::AngleAxisd(euler[0], Eigen::Vector3d::UnitX()));

	/** Measurement vector **/
	z_r = rotation * r1;
	
	/** Innovation covariance **/
	Pnu = Pzz + R;
	
	/** Cross-correlation matrix **/	
	for (i=1; i<ukf::SIGPOINTSIZE; i++)
	{
	    sumPxz += ((sig_point.col(i)-x) * (gamma.col(i)-z_e).transpose());
	}
	Pxz = (1/(ukf::UKFSTATEVECTORSIZE + lambda))*(((sig_point.col(0)-x) * (gamma.col(0)-z_e).transpose()) + (0.5 * sumPxz));
	
	/** Compute the Kalman Gain **/
	K = Pxz * Pnu.inverse();
	
	
	/** Innovation in the measurement **/
	Nu = z_r - z_e;
	
	
	/** Correction of the mean value (state)**/
	x = x + K*Nu;
	
	/** Covariance matrix of the process **/
	Px = Px - K*Pnu*K.transpose();
	
	/** Update the attitude quaternion using the state vector (rodrigues parameter)**/
	q4 = ((-a * (x.block<ukf::NUMAXIS,1>(0,0)).squaredNorm())+ f * sqrt(pow(f,2)+(1-pow(a,2))*(x.block<ukf::NUMAXIS,1>(0,0)).squaredNorm()))/(pow(f,2)+(x.block<ukf::NUMAXIS,1>(0,0)).squaredNorm());
	vectorq = (1/f)*(a + q4) * x.block<ukf::NUMAXIS,1>(0,0);
	
	rotation = Eigen::Quaternion <double> (q4, vectorq[0], vectorq[1], vectorq[2]);

	at_q = rotation * sig_q[0];
	
	/** Set to zero for the next propagation **/
	x.block<ukf::NUMAXIS,1>(0,0) = Eigen::Matrix <double, ukf::NUMAXIS, 1>::Zero();
	
    }


    /**
    * @brief Performs the update step of the filter. Only update the attitude quaternion from a given rotation in rodrigues parameter form
    */
    void ukf::attitudeUpdate()
    {
	
	double q4; /**< scalar part of a quaternion **/
	Eigen::Matrix <double, ukf::QUATERSIZE-1, 1> vectorq; /**< vectorial part of a quaternion **/
	Eigen::Quaternion <double> rotation;

	/** Update the attitude quaternion using the state vector (rodrigues parameter)**/
	q4 = ((-a * (x.block<ukf::NUMAXIS,1>(0,0)).squaredNorm())+ f * sqrt(pow(f,2)+(1-pow(a,2))*(x.block<ukf::NUMAXIS,1>(0,0)).squaredNorm()))/(pow(f,2)+(x.block<ukf::NUMAXIS,1>(0,0)).squaredNorm());
	vectorq = (1/f)*(a + q4) * x.block<ukf::NUMAXIS,1>(0,0);

	rotation = Eigen::Quaternion <double> (q4, vectorq[0], vectorq[1], vectorq[2]);

	at_q = rotation * sig_q[0];

// 	std::cout<<"at_q:\n"<<at_q.x()<<" "<<at_q.y()<<" "<<at_q.z()<<" "<<at_q.w()<<"\n";

	/** Set to zero for the next propagation **/
	x.block<ukf::NUMAXIS,1>(0,0) = Eigen::Matrix <double, ukf::NUMAXIS, 1>::Zero();
    }





   
    

}
