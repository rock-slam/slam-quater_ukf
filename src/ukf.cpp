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
 * @date June 2011.
 * @version 1.0.
 */

#include <iostream> /**< IO C++ Standard library */
#include <algorithm> /**< Algorithm C++ Standard library */
#include <Eigen/LU> /**< Lineal algebra of Eigen */
#include <Eigen/SVD> /**< Singular Value Decomposition (SVD) of Eigen */
#include "ukf.hpp" /**< Unscented Kalman Filter */

namespace filter
{
    /** Namesapces to use **/
    using namespace Eigen;
    using namespace std;
    
    

}
