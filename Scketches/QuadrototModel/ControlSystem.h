#pragma once
#include <Eigen/Dense>

using namespace Eigen;

class ControlSystem
{
	private:
		double Kp = 0;
		double Kd = 0;
		Matrix3d I = Matrix3d::Zero();

	public:
		ControlSystem(double Kp, double Kd, Matrix3d I)
		{
			this->Kp = Kp;
			this->Kd = Kd;
			this->I  = I;
		}

		Vector3d getControl(double ro_d, double lambda_d, double mu_d, double nu_d,
							double ro,   double lambda,   double mu,   double nu,
							double wx, double wy, double wz)
		{
			Vector3d w = Vector3d(wx, wy, wz);
			Quaterniond q_ref = Quaterniond(lambda_d, mu_d, nu_d, ro_d);
			Quaterniond q	  = Quaterniond(lambda, mu, nu, ro);
			
			Quaterniond q_err = q_ref * q.conjugate();
			Vector3d vvv = q_err.vec();
			Vector3d M = -Kp * Vector3d(vvv(2), -vvv(1), vvv(0))  - Kd * w;
			return M;
		}

		/*Vector3d getControl(double ro_d, double lambda_d, double mu_d, double nu_d,
							double ro,   double lambda,   double mu,   double nu,
							double wx,   double wy,       double wz)
		{
			Vector4d qd	   = Vector4d(lambda_d, mu_d, nu_d, ro_d);
			Vector4d q_vec = Vector4d(lambda,   mu,   nu,   ro);
			Vector3d w     = Vector3d(wx, wy, wz);

			Vector4d tau_d_vec = Kp * (qd - qd.dot(q_vec) * q_vec);

			Quaterniond tau_d_2 = Quaterniond(2 * tau_d_vec);
			Quaterniond q = Quaterniond(q_vec);

			Quaterniond wd = (q.conjugate() * tau_d_2);
			Vector3d wd_vec = wd.vec();

			
			Vector3d Eps_d = Kd * (wd_vec - w);
			
			Vector3d ControlMoment = (I * Eps_d) + wd_vec.cross(I * wd_vec);
			return ControlMoment;

		}*/
};

