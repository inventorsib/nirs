#pragma once
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

class NDSolver
{
	private:
		//вектор состояния
		VectorXd x;
		double t;

		VectorXd k1, k2, k3, k4;

		//функция вычисления правой части
		VectorXd(*rightPartDE)(VectorXd x, double t, Vector4d c);


	public:
		NDSolver(VectorXd dxdt(VectorXd x, double t, Vector4d c), VectorXd x0, double t0)
		{
			x = x0;
			t = t0;
			rightPartDE = (dxdt);

			k1 = VectorXd(x.size());
			k2 = VectorXd(x.size());
			k3 = VectorXd(x.size());
			k4 = VectorXd(x.size());
		}

		//NDSolver(VectorXd dxdt(VectorXd x, double t, Vector c), VectorXd x0, double t0)
		//{
		//	x = x0;
		//	t = t0;
		//	rightPartDE = (dxdt);

		//	k1 = VectorXd(x.size());
		//	k2 = VectorXd(x.size());
		//	k3 = VectorXd(x.size());
		//	k4 = VectorXd(x.size());
		//}

		NDSolver()
		{

		}
		void SolveStep(double h, Vector4d c)
		{
			k1 = rightPartDE(				x, t,         c);
			k2 = rightPartDE(0.5 * h * k1 + x, t + h / 2, c);
			k3 = rightPartDE(0.5 * h * k2 + x, t + h / 2, c);
			k4 = rightPartDE(      h * k3 + x, t + h ,    c);

			x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
			t = t + h;
		}

		VectorXd getStateVector()
		{
			return this->x;
		}
		double getTime()
		{
			return t; 
		}

		void setStateVector(VectorXd stateVector)
		{
			if(stateVector.size() == x.size())
			x = stateVector;
		}
		double getStateVectorCoordinate(int coordNum)
		{
			return	x[coordNum];
		}
		void setStateVectorCoordinate(int coordNum, double value)
		{
			if (x.size() > coordNum && coordNum>=0)
				x[coordNum] = value;
		}


};

