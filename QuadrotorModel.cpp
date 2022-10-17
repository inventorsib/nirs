#include <iostream>
#include "QuadrotorDynamics.h"

#include "ControlSystem.h"
#include "PID.h"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <math.h>
#include "Network.h"
#include <vector>
#include <fstream>

void tokenize(std::string const& str, const char delim,
    std::vector<std::string>& out)
{
    size_t start;
    size_t end = 0;

    while ((start = str.find_first_not_of(delim, end)) != std::string::npos)
    {
        end = str.find(delim, start);
        out.push_back(str.substr(start, end - start));
    }
}

double Constrain(const double val, double min, double max)
{
    if (val > max) return max;
    if (val < min) return min;
    return val;
}

int main()
{

    std::ofstream record;          // поток для записи
    record.open("C:\\Users\\Ярослав\\Desktop\\НИРС\\UDP_Examples\\NDSolver\\version-0.1.2\\QuadrotorModel\\state.txt"); // окрываем файл для записи

    std::ifstream in("C:\\Users\\Ярослав\\Desktop\\НИРС\\UDP_Examples\\NDSolver\\version-0.1.2\\QuadrotorModel\\init.txt"); // окрываем файл для чтения
    std::string line;
    if (in.is_open())
    {
        getline(in, line);
    }
    in.close();     // закрываем файл

    std::vector<std::string> state;

    const char delim = ' ';
    tokenize(line, delim, state);


    //Quaterniond quaternion = Quaterniond(0.243, 0.299, 0.634, 0.670);
    //Quaterniond quaternion = Quaterniond(0.98, 0, 0, 0.197);
    //Quaterniond quaternion = Quaterniond(1, 0, 0, 0);
    Quaterniond quaternion = Quaterniond(0.711, 0.194, 0.640, 0.216);

    Vector3d position = Vector3d(std::stod(state[0]), std::stod(state[1]), std::stod(state[2]));
    //Vector3d position = Vector3d(4, 4, 4); // Vector3d::Ones();
    Vector3d velocity = Vector3d::Zero();
    Vector3d angularVelocity = Vector3d::Zero();

    QuadrotorDynamics quadcopter = QuadrotorDynamics(position, velocity, angularVelocity, quaternion);

    ControlSystem controller = ControlSystem(15E6, 1.2E6, 
                                             100E5, 100E5, 0.00,
                                             15, 75, 0,
                                             15, 75, 0    );

    int counter = 0;

    std::string timeStates = "";

    while (quadcopter.GetTime()<10)
    {
        Quaterniond targetQuaternion = controller.getPosControlQuaternion(quadcopter.GetPosition(), Vector3d(1, 1, 1), quadcopter.GetVelocity(), 0.1, 1);
        Vector3d M = controller.getControlMomentStabilize(targetQuaternion, quadcopter.GetQuaternion(), quadcopter.GetAngularVelocity());

        double vy_g = quadcopter.GetVelocity().y();
        double AltControl = controller.getAltitudeControl(quadcopter.GetPosition().y(), 1, vy_g);

        double kAlt = 1;
        Vector4d control = Vector4d(
            Constrain(500 + M(0) + M(1) + kAlt * AltControl, 100, 5000), //1
            Constrain(500 + M(2) - M(1) + kAlt * AltControl, 100, 5000),   //2
            Constrain(500 - M(0) + M(1) + kAlt * AltControl, 100, 5000), //3
            Constrain(500 - M(2) - M(1) + kAlt * AltControl, 100, 5000)    //4
        );


        double trust = ( control.dot(Vector4d(1, 1, 1, 1)) ) * 2.3E-6;

        if (counter % 2 == 0)
        {
            std::string data = quadcopter.GetData();
           
            timeStates += std::to_string(quadcopter.GetTime()) + " " + quadcopter.GetData()+"\n";

//            std::cout << quadcopter.GetTime() << std::endl;
        }
 
       

        quadcopter.SolveStep(control, 0.001);

        counter++;
    }

    std::cout << "Complete" << std::endl;
    record << timeStates << std::endl;
    
}


