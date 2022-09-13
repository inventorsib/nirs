#include <iostream>
#include "NDSolver.h"
#include "ControlSystem.h"
#include "PID.h"
#include <Eigen/Dense>
#include "Network.h"

#define PI 3.14159265358979323851

const double m = 0.1*4 + 2;
const double Ix = 0.1, Iy = 0.2, Iz = 0.1;

//Коэффициент тяги
const double Ku = 1;

//Длина луча
const double L = 0.22, b = L;

const double g = 9.81;

double Sq(double val)
{
    return val * val;
}

double u1 = 0, u2 = 0, u3 = 0, u4 = 0;

VectorXd quadrotorDE(VectorXd X, double t)
{
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

    double ro     = X[ 9];
    double lambda = X[10];
    double mu     = X[11];
    double nu     = X[12];

    double u = sqrt(Sq(u1) + Sq(u2) + Sq(u3) + Sq(u4));

    double th = asin(2 * (ro * nu + lambda * mu ) );
    double fi = atan(2 * (ro * lambda - nu * mu) / ( Sq(ro) + Sq(mu) - Sq(nu) - Sq(lambda) ) );
    double psi = atan(2 * (ro * mu - lambda * nu) / ( Sq(ro) + Sq(lambda) - Sq(mu) - Sq(nu) ) );

    /*1-x*/       dxdt[ 0] = vx;
    /*2-y*/       dxdt[ 1] = vy;
    /*3-z*/       dxdt[ 2] = vz;

	/*4-vx*/      dxdt[ 3] = Ku / m * ( sin(psi) * sin(fi) + cos(psi) * cos(fi) * sin(th) ) * Sq(u);
    /*5-vy*/      dxdt[ 4] = Ku / m * ( cos(th ) * cos(fi)                                ) * Sq(u) - g;
    /*6-vz*/      dxdt[ 5] = Ku / m * ( cos(fi ) * sin(psi) * sin(th) - cos(psi) * sin(fi)) * Sq(u);

	/*7-wx*/      dxdt[ 6] = L * Ku / Ix * (Sq(u1) - Sq(u3))                   - (Iz - Iy) / Ix * wy * wz;
    /*8-wy*/      dxdt[ 7] = b * Ku / Iy * (Sq(u1) - Sq(u2) + Sq(u3) - Sq(u4)) - (Ix - Iz) / Iy * wx * wz;
    /*9-wz*/      dxdt[ 8] = L * Ku / Iz * (Sq(u2) - Sq(u4))                   - (Iy - Ix) / Iz * wx * wy;

    /*10-ro*/     dxdt[ 9] = - 0.5 * (  wx * lambda + wy * mu + wz * nu    );
    /*11-lambda*/ dxdt[10] =   0.5 * (  wx * ro     - wy * nu + wz * mu    );
    /*11-mu*/     dxdt[11] =   0.5 * (  wx * nu     + wy * ro - wz * lambda);
    /*12-nu*/     dxdt[12] =   0.5 * (- wx * mu     + wy*lambda+wz * ro    );

    return dxdt;
}

int main()
{
    //UDP
    std::string IP = "127.0.0.1";
    int PORT = 8888;

    WSASession Session;
    UDPSocket Socket;
    char buffer[100];

    //Начальные условия задача Коши
    Eigen::VectorXd X0(13);
    
    X0(0) = 0;  // x
    X0(1) = 0;  // y
    X0(2) = 0;  // z

    X0(3) = 0;  // vx
    X0(4) = 0;  // vy
    X0(5) = 0;  // vz

    X0(6) = 0;  // wx
    X0(7) = 0;  // wy
    X0(8) = 0;  // wz

    X0(9) = 1;  // ro
    X0(10) = 0; // lambda
    X0(11) = 0; // mu
    X0(12) = 0; // nu

    VectorXd(*dEquations)(VectorXd x, double t);
    dEquations = quadrotorDE;
    NDSolver quadrotor = NDSolver(dEquations, X0, 0);


    Matrix3d I;
    I << Ix,  0,  0,
          0, Iy,  0,
          0,  0, Iz;

    ControlSystem controller = ControlSystem(0.01, 0.01, I);
    
    
   
    while (true)
    {
        double ro = quadrotor.getStateVector()[9];
        double lambda = quadrotor.getStateVector()[10];
        double mu = quadrotor.getStateVector()[11];
        double nu = quadrotor.getStateVector()[12];
       
        std::string  data = std::to_string(quadrotor.getStateVector()[9]) + " " +
                            std::to_string(quadrotor.getStateVector()[10]) + " " +
                            std::to_string(quadrotor.getStateVector()[11]) + " " +
                            std::to_string(quadrotor.getStateVector()[12]) + " ";

        Socket.SendTo(IP, PORT, data.c_str(), data.size());

        
         
        double th = asin(2 * (ro * nu + lambda * mu));
        double fi = atan(2 * (ro * lambda + nu * mu) / (Sq(ro) + Sq(mu) - Sq(nu) - Sq(lambda)));
        double psi = atan(2 * (ro * mu + lambda * nu) / (Sq(ro) + Sq(lambda) - Sq(mu) - Sq(nu)));

        double wx = quadrotor.getStateVector()[6];
        double wy = quadrotor.getStateVector()[7];
        double wz = quadrotor.getStateVector()[8];

        Vector3d M = controller.getControl(  0.7071068, -0.7071068, 0, 0,
                                             ro, lambda, mu, nu,
                                             wx, wy, wz                 );

        u1 = 10 + M(0) + M(1);
        u3 = 10 - M(0) + M(1);

        u2 = 10 + M(2) - M(1);
        u4 = 10 - M(2) - M(1);

        std::cout << data << "  " << psi*180/PI << std::endl;
        
        quadrotor.SolveStep(0.001);
    }
    
}


