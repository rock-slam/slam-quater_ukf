#ifndef _FILTER_UKF_HPP_
#define _FILTER_UKF_HPP_

#include <iostream>

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
    #define STATEVECTORSIZE 9 /**< Number of variables of the vector state-space representation **/
    #define QUATERSIZE 4 /**< Number of parameters of a quaternion **/

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
	    * @brief This function set the initial Omega matrix
	    * 
	    * Initial Omega matrix with angular velocity for 
	    * quaternion integration.
	    *
	    * @author Javier Hidalgo Carrio.
	    *
	    * @param[in] *u pointer to vector with the angular velocity
	    *
	    * @return OK is everything all right. ERROR on other cases.
	    *
	    */
	    int setOmega (Eigen::Matrix <double,NUMAXIS,1>  *u);
      
	    /**
	    * Print a welcome to stdout
	    * \return nothing
	    */

    };

} // end namespace dummy_project

#endif // _FILTER_UKF_HPP_
