#include <iostream>
#include "NDSolver.h"
#include <Eigen/Dense>

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

VectorXd quadrotorDE(VectorXd x, double t)
{
    VectorXd dxdt = VectorXd(x.size());
    
    dxdt[0] = sigma * (x[1] - x[0]);
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];

    return dxdt;
}

int main()
{
    //Начальные условия задача Коши
    Eigen::VectorXd X0(3);

    X0(0) = 10;
    X0(1) = 1;
    X0(2) = 1;


    VectorXd(*dEquations)(VectorXd x, double t);
    dEquations = quadrotorDE;

    NDSolver quadrotor = NDSolver(dEquations, X0, 0);

    for (int i = 0; i < 10000; i++)
    {
        quadrotor.SolveStep(0.01);
        if (i%100==0)
            std::cout << quadrotor.getStateVector()[0]<< " " << quadrotor.getStateVector()[1] << "\n";
    }


    /*Можно явно интегрировать с отрицательным шагом.
    Это позволяет проверить настройку шага метода.*/

    //quadrotor.SolveStep(0.001);
    //std::cout<<quadrotor.getStateVector()<<"\n\n";

    //quadrotor.SolveStep(0.001);
    //std::cout << quadrotor.getStateVector() << "\n\n";

    //quadrotor.SolveStep(-0.001);
    //std::cout << quadrotor.getStateVector() << "\n\n";

}


