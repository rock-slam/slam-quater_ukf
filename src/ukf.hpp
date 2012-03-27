/**\file ukf.hpp
 * Header function file and defines
 */
#ifndef _FILTER_UKF_HPP_
#define _FILTER_UKF_HPP_

#include <iostream>
#include <Eigen/Geometry> /**< Eigen data type for Matrix, Quaternion, etc... */

namespace filter
{
    
    using namespace Eigen;

    /** General defines **/
    #ifndef OK
    #define OK	0  /**< Integer value in order to return when everything is all right. */
    #endif
    #ifndef ERROR
    #define ERROR	-1  /**< Integer value in order to return when an error occured. */
    #endif

    /** IKF constant parameters **/
    #define STATEVECTORSIZE 6 /**< Number of variables of the vector state-space representation **/
    #define QUATERSIZE 4 /**< Number of parameters of a quaternion **/
    #define SIGPOINTSIZE (2*STATEVECTORSIZE) + 1 /**< Number of Sigma Points **/

    #ifndef PI
    #define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286 /**< Pi Number */
    #endif
    #define EAST 1 /**< EAST is 1 and means positive magnetic declination **/
    #define WEST 2 /**< WEST is 2 and means negative magnetic declination **/

    #define D2R PI/180.00 /**< Convert degree to radian **/
    #define R2D 180.00/PI /**< Convert radian to degree **/

    /** Sensors constant parameters **/
    #ifndef NUMAXIS
    #define NUMAXIS 3 /**< Number of axis sensed by the sensors **/
    #endif

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


  
    class ukf
    {
	
	/**
	 * Filters members
	 **/
	private:
	    double f, a, lambda; /**< Parameters for the Unscented KF  */
	    Eigen::Matrix <double,STATEVECTORSIZE,1> x; /**< State vector */
	    Eigen::Matrix <double,NUMAXIS,1> gtilde; /**< gravitation acceleration */
	    Eigen::Matrix <double,STATEVECTORSIZE, STATEVECTORSIZE> Px; /**< State covariance matrix */
	    Eigen::Quaternion <double> at_q;  /**< Attitude quaternion. Note the order of the arguments: the real w coefficient first, while internally the coefficients are stored in the following order: [x, y, z, w] */
	    Eigen::Matrix <double,STATEVECTORSIZE, STATEVECTORSIZE> Q; /**< Process noise covariance matrix */
	    Eigen::Matrix <double,NUMAXIS, NUMAXIS>  R; /**< Measurements noise variance and covar matrix */
	    Eigen::Matrix <double,STATEVECTORSIZE, SIGPOINTSIZE> sig_point; /**< Sigma points for SPKF */
	    Eigen::Matrix <Eigen::Quaternion <double>, SIGPOINTSIZE, 1> e_q; /**< Error quaternions */
	    Eigen::Matrix <Eigen::Quaternion <double>, SIGPOINTSIZE, 1> sig_q; /**< Sigma point quaternions */
	    Eigen::Matrix <double,NUMAXIS, SIGPOINTSIZE> gamma; /**< Observation in the model */
	    Eigen::Matrix <double,NUMAXIS, 1> z_e; /**<Predictive observation */ 
	    Eigen::Matrix <double,NUMAXIS, 1> z_r; /**<Measurement of the observation */
	    Eigen::Matrix <double,NUMAXIS, 1> Nu; /**< Innovation */
	    Eigen::Matrix <double,NUMAXIS, NUMAXIS> Pzz; /** Output covariance matrix */
	    Eigen::Matrix <double,NUMAXIS, NUMAXIS> Pnu; /** Innovation covariance matrix */
	    Eigen::Matrix <double,STATEVECTORSIZE, NUMAXIS> Pxz; /** Cross-correlation matrix */
	    Eigen::Matrix <double,STATEVECTORSIZE, NUMAXIS> K; /**< Kalman Gain sometimes also called W */
	
	public: 
	    
	    
	    /**
	    * @brief Gets the current state vector of the filter
	    * 
	    * @author Javier Hidalgo Carrio.
	    *
	    * @return State Vector
	    *
	    */
	    Eigen::Matrix <double,STATEVECTORSIZE,1> getState();


	    /**
	    * @brief Gets the current orientation in Quaternion
	    * 
	    * @author Javier Hidalgo Carrio.
	    *
	    * @return Quaternion with the current orientation.
	    *
	    */
	    Eigen::Quaternion <double> getAttitude();

	    /**
	    * @brief Gets the current orientation in Euler angles
	    * 
	    * @author Javier Hidalgo Carrio.
	    *
	    * @return Current orientation in Euler angles.
	    *
	    */
	    Eigen::Matrix <double, NUMAXIS, 1> getEuler();

	    /**
	    * @brief Gets Noise covariance matrix
	    * 
	    * @author Javier Hidalgo Carrio.
	    *
	    * @return Matrix P of the covariance of the state vector
	    *
	    */
	    Eigen::Matrix <double,STATEVECTORSIZE,STATEVECTORSIZE> getCovariance();

	    /**
	    * @brief This function Initilize Attitude
	    * 
	    * Initial orientation value beforeestart the IKF 
	    *
	    * @author Javier Hidalgo Carrio.
	    *
	    * @param[in] *initq pointer to quaternion with the initial orientation
	    *
	    * @return OK is everything all right. ERROR on other cases.
	    *
	    */
	    int setAttitude (Eigen::Quaternion <double> *initq);

	    	    
	    /**
	    * @brief This function initilize the filter
	    * 
	    * This method receives the measurement noise matrix of the sensors
	    * and filter parameters
	    *
	    * @author Javier Hidalgo Carrio.
	    *
	    * @param[in] *x_0 initial state vector
	    * @param[in] *P_0 initial convariance matrix of the process.
	    * @param[in] *Q noise covariance matrix of the model
	    * @param[in] *R noise covariance matrix of the measurement
	    * @param[in] *a parameter for the UKF
	    * @param[in] *f parameter for the UKF\
	    * @param[in] *lambda parameter for the UKF, to define the distance of the sigma point with respect to the mean
	    *
	    * @return void
	    *
	    */
	    void Init (Matrix <double,STATEVECTORSIZE,1> *x_0, Eigen::Matrix <double, STATEVECTORSIZE, STATEVECTORSIZE> *P_0, Eigen::Matrix <double, STATEVECTORSIZE, STATEVECTORSIZE> *Q, Eigen::Matrix <double, NUMAXIS, NUMAXIS> *R,
		       Eigen::Quaternion <double> *at_q, double a, double f, double lambda, double g);
	    
	    /**
	    * @brief This computes the theoretical gravity value according to the WGS-84 ellipsoid earth model.
	    *
	    * @author Javier Hidalgo Carrio.
	    *
	    * @param[in] latitude double the latitude value in radian
	    * @param[in] altitude double with the altitude value in meters
	    *
	    * @return double. the theoretical value of the local gravity
	    *
	    */
	    double GravityModel (double latitude, double altitude);
	    
	    /**
	    * @brief Substract the Earth rotation from the gyroscopes readout
	    *
	    * This function computes the substraction of the rotation of the Earth (EARTHW)
	    * from the gyroscope values. This function uses quaternion of transformation from
	    * the body to the geographic frame and the latitude in radians.
	    *
	    * @author Javier Hidalgo Carrio.
	    *
	    * @param[in, out] *u pointer to angular velocity
	    * @param[in] *qb_g quaternion from body frame to geographic frame
	    * @param[in] latitude location latitude angle in radians
	    *
	    * @return void
	    *
	    */
	    void SubstractEarthRotation(Eigen::Matrix <double, NUMAXIS, 1> *u, Eigen::Quaternion <double> *qb_g, double latitude);
      
	     /**
	    * @brief Discrete-time quaternion kinematic equation
	    * 
	    * It computes the quaternion kinematics from angular velocity.
	    * It is the quaternions integration (discrete version)
	    *
	    * @author Javier Hidalgo Carrio.
	    *
	    * @param [in,out] *quat the quaternion to propagate
	    * @param[in] *angvelo pointer to vector with the angular velocity
	    * @param[in] dt delta time between samples
	    *
	    * @return void
	    *
	    */
	    void Omega(Eigen::Quaternion< double >* quat, Eigen::Matrix< double, NUMAXIS , 1  >* angvelo, double dt);

	    
	    /**
	    * @brief Performs the prediction step of the filter.
	    * 
	    * It computes the sigma point, the error quaternions, the sigma point quaternion form the
	    * previous error quaternions and the propagation step of the filter
	    *
	    * @author Javier Hidalgo Carrio.
	    *
	    * @param[in] *u pointer to vector with the angular velocity
	    * @param[in] dt delta time between samples
	    *
	    * @return void
	    *
	    */
	    void predict(Eigen::Matrix <double,NUMAXIS,1>  *u, double dt);
	    
	    
	    
	    /**
	    * @brief Performs the measurement and correction steps of the filter.
	    * 
	    * The UKf measurement step
	    *
	    * @author Javier Hidalgo Carrio.
	    *
	    * @param[in] *u pointer to vector with the angular velocity
	    * @param[in] dt delta time between samples
	    *
	    * @return void
	    *
	    */
	    void update(Eigen::Matrix <double,NUMAXIS,1>  *acc, Eigen::Matrix <double,NUMAXIS,1>  *mag);
	    
	    /**
	    * @brief Performs attitude quaternion update.
	    * 
	    * It only performs the update of the attitude quaternion from the state vector (x)
	    * After the ccomputation set x to zero for next propagation (only the Rodriguez parameters part of x)
	    *
	    * @author Javier Hidalgo Carrio.
	    *
	    * @param[in,out] *at_q pointer to qttitude quaternion
	    * @param[in,out] *x pointer to the KF state vector
	    *
	    * @return void
	    *
	    */
	    void attitudeUpdate();
	    
	

    };

} // end namespace dummy_project

#endif // _FILTER_UKF_HPP_
