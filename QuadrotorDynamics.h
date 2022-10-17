#pragma once
#include <iostream>
#include "NDSolver.h"
#include "ControlSystem.h"
#include "PID.h"
#include <Eigen/Dense>
#include <math.h>
#include "Network.h"

#define _USE_MATH_DEFINES
#define g 9.81


using namespace Eigen;

class QuadrotorDynamics
{
	private:
        
        double time = 0;
        //Начальные условия задача Коши
        VectorXd X0 = VectorXd(13);

   		static double Sq(const double val)
		{
			return val * val;
		}

		static double Constrain(const double val, double min, double max)
		{
			if (val > max) return max;
			if (val < min) return min;
			return val;
		}

        static VectorXd quadrotorDE(VectorXd X, double t, Vector4d control)
        {
          double m = 4;
          double Ix = 0.1;
          double Iy = 0.2;
          double Iz = 0.1;
          double Ku = 2.3E-6;
          double L = 0.22;
          double b = L;

          double u1 = control(0);
          double u2 = control(1);
          double u3 = control(2);
          double u4 = control(3);

            VectorXd dxdt = VectorXd(X.size());

            double x = X[0];
            double y = X[1];
            double z = X[2];

            double vx = X[3];
            double vy = X[4];
            double vz = X[5];

            double wx = X[6];
            double wy = X[7];
            double wz = X[8];

            double ro = X[9];
            double lambda = X[10];
            double mu = X[11];
            double nu = X[12];

            double u = sqrt(Sq(u1) + Sq(u2) + Sq(u3) + Sq(u4));

            double th  = asin(2 * (ro * nu + lambda * mu));
            double fi  = atan(2 * (ro * lambda - nu * mu) / (Sq(ro) + Sq(mu) - Sq(nu) - Sq(lambda)));
            double psi = atan(2 * (ro * mu - lambda * nu) / (Sq(ro) + Sq(lambda) - Sq(mu) - Sq(nu)));

            /*1-x*/       dxdt[0] = vx;
            /*2-y*/       dxdt[1] = vy;
            /*3-z*/       dxdt[2] = vz;

            /*4-vx*/      dxdt[3] = Ku / m * (sin(psi) * sin(fi) + cos(psi) * cos(fi) * sin(th)) * Sq(u);
            /*5-vy*/      dxdt[4] = Ku / m * (cos(th) * cos(fi)) * Sq(u) - g;
            /*6-vz*/      dxdt[5] = Ku / m * (cos(fi) * sin(psi) * sin(th) - cos(psi) * sin(fi)) * Sq(u);

            /*7-wx*/      dxdt[6] = L * Ku / Ix * (Sq(u1) - Sq(u3)) - (Iz - Iy) / Ix * wy * wz;
            /*8-wy*/      dxdt[7] = b * Ku / Iy * (Sq(u1) - Sq(u2) + Sq(u3) - Sq(u4)) - (Ix - Iz) / Iy * wx * wz;
            /*9-wz*/      dxdt[8] = L * Ku / Iz * (Sq(u2) - Sq(u4)) - (Iy - Ix) / Iz * wx * wy;

            /*10-ro*/     dxdt[9] = -0.5 * (wx * lambda + wy * mu + wz * nu);
            /*11-lambda*/ dxdt[10] = 0.5 * (wx * ro - wy * nu + wz * mu);
            /*11-mu*/     dxdt[11] = 0.5 * (wx * nu + wy * ro - wz * lambda);
            /*12-nu*/     dxdt[12] = 0.5 * (-wx * mu + wy * lambda + wz * ro);

            return dxdt;
        }

        NDSolver quadrotor = NDSolver();


    public:
        QuadrotorDynamics(Vector3d Position, Vector3d Velocity, Vector3d AngelarVelocity, Quaterniond quaternion)
        {

            X0(0) = Position[0];  // x
            X0(1) = Position[1];  // y
            X0(2) = Position[2];  // z

            X0(3) = Velocity[0];  // vx
            X0(4) = Velocity[1];  // vy
            X0(5) = Velocity[2];  // vz

            X0(6) = AngelarVelocity[0];  // wx
            X0(7) = AngelarVelocity[1];  // wy
            X0(8) = AngelarVelocity[2];  // wz

            X0( 9) = quaternion.w();  // ro
            X0(10) = quaternion.x(); // lambda
            X0(11) = quaternion.y(); // mu
            X0(12) = quaternion.z(); // nu

            VectorXd(*dEquations)(VectorXd x, double t, Vector4d c);
            dEquations = quadrotorDE;
            quadrotor = NDSolver(dEquations, X0, 0);
        }


        void SolveStep(const Vector4d control, double h)
        {

            double ro = quadrotor.getStateVector()[9];
            double lambda = quadrotor.getStateVector()[10];
            double mu = quadrotor.getStateVector()[11];
            double nu = quadrotor.getStateVector()[12];

            double th = asin(2 * (ro * nu + lambda * mu));
            double fi = atan(2 * (ro * lambda + nu * mu) / (Sq(ro) + Sq(mu) - Sq(nu) - Sq(lambda)));
            double psi = atan(2 * (ro * mu + lambda * nu) / (Sq(ro) + Sq(lambda) - Sq(mu) - Sq(nu)));

            double wx = quadrotor.getStateVector()[6];
            double wy = quadrotor.getStateVector()[7];
            double wz = quadrotor.getStateVector()[8];

            quadrotor.SolveStep(h, control);   
            time += h;
        }

        double GetTime()
        {
            return time;
        }
        
        Vector3d GetPosition()
        {
            Vector3d pos = Vector3d(quadrotor.getStateVectorCoordinate(0), quadrotor.getStateVectorCoordinate(1), quadrotor.getStateVectorCoordinate(2));
            return pos;
        }
        Vector3d GetVelocity()
        {
            Vector3d vel = Vector3d(quadrotor.getStateVectorCoordinate(3), quadrotor.getStateVectorCoordinate(4), quadrotor.getStateVectorCoordinate(5));
            return vel;
        }
        Vector3d GetAngularVelocity()
        {
            Vector3d angVel = Vector3d(quadrotor.getStateVectorCoordinate(6), quadrotor.getStateVectorCoordinate(7), quadrotor.getStateVectorCoordinate(8));
            return angVel;
        }
        Quaterniond GetQuaternion()
        {
            Quaterniond quat = Quaterniond(
                quadrotor.getStateVectorCoordinate(9),
                quadrotor.getStateVectorCoordinate(10), 
                quadrotor.getStateVectorCoordinate(11),
                quadrotor.getStateVectorCoordinate(12)
            );
            return quat;
        }

        std::string GetData()
        {
            std::string data =  std::to_string(quadrotor.getStateVector()[9]) + " " +
                                std::to_string(quadrotor.getStateVector()[10]) + " " +
                                std::to_string(quadrotor.getStateVector()[11]) + " " +
                                std::to_string(quadrotor.getStateVector()[12]) + " " +
                                std::to_string(0+1*quadrotor.getStateVector()[ 0]) + " " +
                                std::to_string(0+1*quadrotor.getStateVector()[ 1]) + " " +
                                std::to_string(0+1*quadrotor.getStateVector()[2]) + " " ;

          return data;
        }

        Vector3d GetEulerAngles()
        {
            double ro = quadrotor.getStateVector()[9];
            double lambda = quadrotor.getStateVector()[10];
            double mu = quadrotor.getStateVector()[11];
            double nu = quadrotor.getStateVector()[12];

            double th = asin(2 * (ro * nu + lambda * mu));
            double fi = atan(2 * (ro * lambda + nu * mu) / (Sq(ro) + Sq(mu) - Sq(nu) - Sq(lambda)));
            double psi = atan(2 * (ro * mu + lambda * nu) / (Sq(ro) + Sq(lambda) - Sq(mu) - Sq(nu)));

            return Vector3d( th, psi, fi );
        }

        VectorXd GetStateVector()
        {
            return quadrotor.getStateVector();
        }

};

