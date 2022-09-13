#pragma once
class PID
{
	private:
		double Kp;
		double Kd;
		double Ki;
		double errorIntegral = 0;

	public:
		PID(double Kp, double Kd, double Ki)
		{
			this->Kp = Kp;
			this->Kd = Kd;
			this->Ki = Ki;
		}

		double GetControl(double value, double dValue)
		{
			errorIntegral += value;
			return Kp * value + Kd * dValue + Ki * errorIntegral;
		}
};

