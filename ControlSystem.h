#pragma once
#include <Eigen/Dense>
#include "PID.h"
#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;

class ControlSystem
{

	private:
		double Kp = 0;
		double Kd = 0;
		PID pid = PID(0, 0, 0);
		PID pidX = PID(0, 0, 0);
		PID pidZ = PID(0, 0, 0);

		static double Sq(const double val)
		{
			return val * val;
		}

	public:
		ControlSystem(double Kp, double Kd, double KpA, double KdA, double KiA, double KpX, double KdX, double KiX, double KpZ, double KdZ, double KiZ)
		{
			this->Kp = Kp;
			this->Kd = Kd;
			this->pid = PID(KpA, KdA, KiA);
			this->pidX = PID(KpX, KdX, KiX);
			this->pidZ = PID(KpZ, KdZ, KiZ);
		}

		Vector3d getControlMomentStabilize(Quaterniond q_ref, Quaterniond q, Vector3d w)
		{
			Quaterniond q_err = q_ref * q.conjugate();
			Vector3d vector_part = q_err.vec();
			//Vector3d M = -Kp * Vector3d(vector_part(2), -vector_part(1), vector_part(0))  - Kd * w;
			Vector3d M = -Kp * Vector3d(-vector_part(0), -vector_part(1), -vector_part(2))  - Kd * w;
			return M;
		}

		Quaterniond getTargetQuaternion(Vector3d position, Vector3d target)
		{
			Vector3d targetRelPosition = target - position;

			if ( abs(targetRelPosition.x()+abs(targetRelPosition.z()) )<0.1)
				return Quaterniond(1, 0, 0, 0);

			double psi_asin = asin(targetRelPosition.z() / targetRelPosition.norm());
			double psi_acos = acos(targetRelPosition.x() / targetRelPosition.norm());
			double psi;
			if (psi_asin >= 0)
				psi = psi_acos;
			if (psi_asin <= 0)
				psi = -psi_acos;

			
			Quaterniond rq = AngleAxisd(0, Vector3d::UnitX())
						  * AngleAxisd(0, Vector3d::UnitY())
						  * AngleAxisd(psi, Vector3d::UnitZ());

			return rq;
		}

		double getAltitudeControl(double position, double target, double vy_g)
		{
			double y_err = target - position;
			return pid.GetControl(y_err, vy_g);
		}

		static double Constrain(const double val, double min, double max)
		{
			if (val > max) return max;
			if (val < min) return min;
			return val;
		}


		Quaterniond getPosControlQuaternionInvariantPsi(Quaterniond targetQuaternion, Vector3d position, Vector3d target, Vector3d velocity, double KP, double KD)
		{

			//DBG 
			velocity = velocity * 0;
			//DBG 


			Vector3d error = (target - position) - velocity * KD;
			error.y() = 0;	//(x, 0, z) 

			Vector3d axis = error.cross(Vector3d(0, 1, 0));
			axis.normalize();

			double angle = error.norm();

			angle = Constrain(angle, -1 / 10.0, +1 / 10.0);

			axis = axis * abs(sin(angle));
			Quaterniond q = Quaterniond(cos(angle), axis.x(), axis.y(), axis.z());

			//std::cout << q << std::endl;
			//std::cout << q * targetQuaternion.conjugate() << std::endl;

			//return q*targetQuaternion.conjugate();
			return targetQuaternion;
		}
			
		Quaterniond getPosControlQuaternion(Vector3d position, Vector3d target, Vector3d velocity, double KP, double KD)
		{
			Vector3d error = (target - position) - velocity*KD;
			error.y() = 0;	//(x, 0, z) 

			Vector3d axis = error.cross(Vector3d(0, 1, 0));
			axis.normalize();

			double angle = error.norm();

			angle = Constrain(angle, -1 / 10.0, +1 / 10.0);

			axis = axis * abs(sin(angle));
			Quaterniond q = Quaterniond(cos(angle), axis.x(), axis.y(), axis.z());

			return q;
		}

		//TODO: bad version
		Quaterniond keeperQuaternion(double PidX, double PidZ, Quaterniond q)
		{
			double th = asin(2*( q.w()*q.y() - q.z()*q.x() ));

			

			double th_star = atan(PidX / 9.81);
			double fi_star = -atan(PidZ * cos(th) / 9.81);


			Quaterniond rq =  AngleAxisd(fi_star, Vector3d::UnitX())
							* AngleAxisd(0, Vector3d::UnitY())
							* AngleAxisd(th_star, Vector3d::UnitZ());

			return rq;
		}
		double getZControl(double position, double target, double vz_g)
		{
			double z_err = target - position;
			return pidZ.GetControl(z_err, vz_g);
		}
		double getXControl(double position, double target, double vx_g)
		{
			double x_err = target - position;
			return pidX.GetControl(x_err, vx_g);
		}


		//Vector3d getControlAltitude(double h_ref, double h, double dh)
		//{
		//	double h_err = h_ref - h;
		//	re
		//}
};

