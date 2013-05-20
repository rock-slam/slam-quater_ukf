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

    class ukf
    {
        public:
            static const int  UKFSTATEVECTORSIZE = 6; /** Number of variables of the vector state-space representation **/
            static const int  QUATERSIZE = 4; /** Number of parameters of a quaternion **/
            static const int  SIGPOINTSIZE =  (2*ukf::UKFSTATEVECTORSIZE) + 1; /** Number of Sigma Points **/

            /** Sensors constant parameters **/
            static const int  NUMAXIS = 3; /** Number of axis sensed by the sensors **/

            /** UKF constant parameters **/
            enum DECLINATION_CONSTS {
                EAST = 1, /** EAST is 1 and means positive magnetic declination **/
                WEST = 2 /** WEST is 2 and means negative magnetic declination **/
            };

	
	/**
	 * Filters members
	 **/
	private:
	    double f, a, lambda; /**< Parameters for the Unscented KF  */
	    Eigen::Matrix <double,ukf::UKFSTATEVECTORSIZE,1> x; /**< State vector */
	    Eigen::Matrix <double,ukf::NUMAXIS,1> gtilde; /**< gravitation acceleration */
	    Eigen::Matrix <double,ukf::UKFSTATEVECTORSIZE, ukf::UKFSTATEVECTORSIZE> Px; /**< State covariance matrix */
	    Eigen::Quaternion <double> at_q;  /**< Attitude quaternion. Note the order of the arguments: the real w coefficient first, while internally the coefficients are stored in the following order: [x, y, z, w] */
	    Eigen::Matrix <double,ukf::UKFSTATEVECTORSIZE, ukf::UKFSTATEVECTORSIZE> Q; /**< Process noise covariance matrix */
	    Eigen::Matrix <double,ukf::NUMAXIS, ukf::NUMAXIS>  R; /**< Measurements noise variance and covar matrix */
	    Eigen::Matrix <double,ukf::UKFSTATEVECTORSIZE, ukf::SIGPOINTSIZE> sig_point; /**< Sigma points for SPKF */
	    Eigen::Matrix <Eigen::Quaternion <double>, ukf::SIGPOINTSIZE, 1> e_q; /**< Error quaternions */
	    Eigen::Matrix <Eigen::Quaternion <double>, ukf::SIGPOINTSIZE, 1> sig_q; /**< Sigma point quaternions */
	    Eigen::Matrix <double,ukf::NUMAXIS, ukf::SIGPOINTSIZE> gamma; /**< Observation in the model */
	    Eigen::Matrix <double,ukf::NUMAXIS, 1> z_e; /**<Predictive observation */ 
	    Eigen::Matrix <double,ukf::NUMAXIS, 1> z_r; /**<Measurement of the observation */
	    Eigen::Matrix <double,ukf::NUMAXIS, 1> Nu; /**< Innovation */
	    Eigen::Matrix <double,ukf::NUMAXIS, ukf::NUMAXIS> Pzz; /** Output covariance matrix */
	    Eigen::Matrix <double,ukf::NUMAXIS, ukf::NUMAXIS> Pnu; /** Innovation covariance matrix */
	    Eigen::Matrix <double,ukf::UKFSTATEVECTORSIZE, ukf::NUMAXIS> Pxz; /** Cross-correlation matrix */
	    Eigen::Matrix <double,ukf::UKFSTATEVECTORSIZE, ukf::NUMAXIS> K; /**< Kalman Gain sometimes also called W */
	
	public: 
	    
	    
	    /**
	    * @brief Gets the current state vector of the filter
	    * 
	    * @author Javier Hidalgo Carrio.
	    *
	    * @return State Vector
	    *
	    */
	    Eigen::Matrix <double,ukf::UKFSTATEVECTORSIZE,1> getState();


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
	    Eigen::Matrix <double, ukf::NUMAXIS, 1> getEuler();

	    /**
	    * @brief Gets Noise covariance matrix
	    * 
	    * @author Javier Hidalgo Carrio.
	    *
	    * @return Matrix P of the covariance of the state vector
	    *
	    */
	    Eigen::Matrix <double,ukf::UKFSTATEVECTORSIZE,ukf::UKFSTATEVECTORSIZE> getCovariance();

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
	    bool setAttitude (Eigen::Quaternion <double> *initq);

	    	    
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
	    * @param[in] *at_q initial attitude quaternion
	    * @param[in] *a parameter for the UKF
	    * @param[in] *f parameter for the UKF\
	    * @param[in] *lambda parameter for the UKF, to define the distance of the sigma point with respect to the mean
	    * @param[in] g local gravity value
	    *
	    * @return void
	    *
	    */
	    void Init (Matrix <double,ukf::UKFSTATEVECTORSIZE,1> *x_0, Eigen::Matrix <double, ukf::UKFSTATEVECTORSIZE, ukf::UKFSTATEVECTORSIZE> *P_0, Eigen::Matrix <double, ukf::UKFSTATEVECTORSIZE, ukf::UKFSTATEVECTORSIZE> *Q, Eigen::Matrix <double, ukf::NUMAXIS, ukf::NUMAXIS> *R,
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
	    void SubstractEarthRotation(Eigen::Matrix <double, ukf::NUMAXIS, 1> *u, Eigen::Quaternion <double> *qb_g, double latitude);
      
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
	    void Omega(Eigen::Quaternion< double >* quat, Eigen::Matrix< double, ukf::NUMAXIS , 1  >* angvelo, double dt);

	    
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
	    void predict(Eigen::Matrix <double,ukf::NUMAXIS,1>  *u, double dt);
	    
	    
	    
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
	    void update(Eigen::Matrix <double,ukf::NUMAXIS,1>  *acc, Eigen::Matrix <double,ukf::NUMAXIS,1>  *mag);
	    
	    /**
	    * @brief Performs attitude quaternion update.
	    * 
	    * It only performs the update of the attitude quaternion from the state vector (x)
	    * After the ccomputation set x to zero for next propagation (only the Rodriguez parameters part of x)
	    *
	    * @author Javier Hidalgo Carrio.
	    *
	    *
	    * @return void
	    *
	    */
	    void attitudeUpdate();
	    
	

    };

} // end namespace dummy_project

#endif // _FILTER_UKF_HPP_
