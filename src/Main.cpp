#include <iostream>
#include <quater_ukf/ukf.hpp>


void quater_mult (Eigen::Quaternion <double> *q1, Eigen::Quaternion <double> *q2, Eigen::Quaternion <double> *q3)
{

	q3->w() = q1->w()*q2->w() - q1->x()*q2->x() - q1->y()*q2->y() - q1->z()*q2->z(); //(pq1[3]*pq2[3] - pq1[0]*pq2[0]-pq1[1]*pq2[1] -pq1[2]*pq2[2]);
	
	q3->x() = q1->w()*q2->x() + q1->x()*q2->w() + q1->y()*q2->z() - q1->z()*q2->y(); //(pq1[3]*pq2[0] + pq1[0]*pq2[3] +pq1[1]*pq2[2] -pq1[2]*pq2[1]);
	
	q3->y() = q1->w()*q2->y() - q1->x()*q2->z() + q1->y()*q2->w() + q1->z()*q2->x(); //(pq1[3]*pq2[1] - pq1[0]*pq2[2] +pq1[1]*pq2[3] +pq1[2]*pq2[0]);
	
	q3->z() = q1->w()*q2->z() + q1->x()*q2->y() + q1->z()*q2->w() - q1->y()*q2->x(); //(pq1[3]*pq2[2] +pq1[0]*pq2[1]-pq1[1]*pq2[0] +pq1[2]*pq2[3]);

	return;
}
int main(int argc, char** argv)
{
	filter::ukf myukf;
	register int i;
	Eigen::Quaternion <double> quater;
	Eigen::Matrix <Eigen::Quaternion <double>, SIGPOINTSIZE, 1> e_q; /**< Error quaternions */
	Eigen::Matrix <double,NUMAXIS,NUMAXIS> L;
	Eigen::Matrix <double,NUMAXIS,1> mivector;
	Eigen::Matrix3d a;
	Eigen::Matrix <double, NUMAXIS, 1> euler;
	
	Eigen::Matrix <double,STATEVECTORSIZE,1> x_0; /**< Initial state vector */
	Eigen::Matrix <double,STATEVECTORSIZE,1> vector; /**< Initial state vector */
	Eigen::Matrix <double,NUMAXIS,1> gtilde; /**< gravitation acceleration */
	Eigen::Matrix <double,STATEVECTORSIZE, STATEVECTORSIZE> P_0; /**< Initial State covariance matrix */
	Eigen::Quaternion <double> at_q;  /**< Attitude quaternion. Note the order of the arguments: the real w coefficient first, while internally the coefficients are stored in the following order: [x, y, z, w] */
	Eigen::Matrix <double,STATEVECTORSIZE, STATEVECTORSIZE> Q; /**< Process noise covariance matrix */
	Eigen::Matrix <double,NUMAXIS, NUMAXIS>  R; /**< Measurements noise variance and covar matrix */
	Eigen::Matrix <double,NUMAXIS, 1>  u;
	Eigen::Matrix <double,NUMAXIS, 1>  v;
	
	
	mivector << 1, -4, 5;
	
	euler << 1.57, 0,0;
	a << 0.45, 0, 0,
		0, 0.56, 0,
		0, 0, 0.7;
	
	std::cout <<"a\n"<<a<<"\n";
	
// 	a.llt();
	Eigen::LLT<Eigen::Matrix3d> lltOfa(a);
	
	L = lltOfa.matrixL();
	
	std::cout <<"a\n"<<L<<"\n";
	
	std::cout <<"a\n"<<L*L<<"\n"; 
	
	e_q[0].w() = 1;
	e_q[0].x() = 0.6;
	e_q[0].y() = 0.3;
	e_q[0].z() = 0.56;
	
	e_q[1].w() = 0.4;
	e_q[1].x() = 0.2;
	e_q[1].y() = 0.5;
	e_q[1].z() = 0.6;
	
	e_q[10].w() = 0.9876;
	e_q[10].x() = 0.2;
	e_q[10].y() = 0.5;
	e_q[10].z() = 0.6;
	
	for (i=0; i<SIGPOINTSIZE; i++)
	{
	    std::cout << "Array("<<i<<") of quaternions"<< e_q[i].w() << "\n";
	}
	
	quater = e_q[0] * e_q[1];
	
	std::cout << "quater:\n"<< quater.x()<<" " << quater.y()<< " " << quater.z()<< " " << quater.w()<<"\n";
	
	quater.setIdentity();
	quater_mult(&(e_q[0]), &(e_q[1]), &(quater));
	std::cout << "quater:\n"<< quater.x()<<" " << quater.y()<< " " << quater.z()<< " " << quater.w()<<"\n";
	
	std::cout << "quater2RotationMatrix:\n"<< quater.toRotationMatrix() <<"\n";
	
	Eigen::Transform <double, NUMAXIS, NUMAXIS> trans(quater);
	
	std::cout << "transform\n"<< trans.matrix() <<"\n";
	
	std::cout << "transform mivector\n"<<mivector<<"\n"<< trans.matrix().block<NUMAXIS, NUMAXIS>(0,0)*mivector <<"\n";
	
	std::cout << "transform mivector\n"<<mivector<<"\n con quater\n"<< quater*mivector <<"\nsize = "<<(quater*mivector).size()<<"\n";
	
	e_q[4] = Eigen::Quaternion <double> (Eigen::AngleAxisd(euler[2], Eigen::Vector3d::UnitZ())*
 			    Eigen::AngleAxisd(euler[1], Eigen::Vector3d::UnitY()) *
 			    Eigen::AngleAxisd(euler[0], Eigen::Vector3d::UnitX()));
	
	
	
	e_q[5] = Eigen::Quaternion <double> (Eigen::AngleAxisd(euler[0], Eigen::Vector3d::UnitX())*
 			    Eigen::AngleAxisd(euler[1], Eigen::Vector3d::UnitY()) *
 			    Eigen::AngleAxisd(euler[2], Eigen::Vector3d::UnitZ()));
	
	std::cout << "e_q[4] is:\n"<< e_q[4].x()<<" " << e_q[4].y()<< " " << e_q[4].z()<< " " << e_q[4].w()<<"\n";
	std::cout << "e_q[5] is:\n"<< e_q[5].x()<<" " << e_q[5].y()<< " " << e_q[5].z()<< " " << e_q[5].w()<<"\n";
	
	
	/** INIT FIL:TER VALUES **/
	gtilde << 0,0, 9.81;
	at_q = Eigen::Quaternion<double>::Identity();
	x_0 << 0,0,0,0.1,0.1,0.1;
	vector << 0.25, 0.25, 0.25, 0.2, 0.2, 0.2;
	P_0 = vector.asDiagonal();
	vector << 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005;
	Q = vector.asDiagonal();
	R << 0.0002, 0.00, 0.00,
	    0.00, 0.0002, 0.00,
	    0.00, 0.00, 0.0002;
	u<< (10*D2R), (10*D2R), (10*D2R);
	
	myukf.Init(&x_0, &P_0, &Q, &R, &at_q,(double)1.00, (double)4.00, 1, gtilde[2]);
	
	/** LUNCH THE FILTER **/
	myukf.predict(&u, 0.01);
	
	std::cout << "Enter";
	std::cin >> i;
	
	mivector << 0.00, 0.00, 9.81;
	
	myukf.update(&v, &mivector);
	
	return 0;
}
