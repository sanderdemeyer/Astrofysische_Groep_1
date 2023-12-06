#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream> 
#include <ctime>
#include <vector>
#include <functional>
#include <string>
#include <unordered_map>
#include <map>
#include <cassert>
#include <filesystem>
#include <algorithm>

#include "classes.cpp"

//#define ADAPTIVE_TIME_STEP true


typedef void (*integ) (NSystem&, double);

int main(){

    // Each new integrator must be added to this map
    std::unordered_map<std::string, integ> functions ={
        {"Forward Euler", Forward_Euler}, // 1 driver evaluation
        {"RK2", RK2_step}, // 2 driver evaluations
        {"Heun", Heun}, // 2 driver evaluations
        {"Heun3", Heun3}, // 3 driver evaluations
        {"Ralston", Ralston}, // 2 driver evaluations
        {"Ralston3", Ralston3}, // 3 driver evaluations
        {"RK3", RK3_step}, // 3 driver evaluations
        {"RK4", RK4_step}, // 4 driver evaluations
        {"Forest Ruth", Forest_Ruth_friend}, // 3 driver evaluations
        {"PEFRL", PEFRL_friend}, // 4 driver evaluations
        {"Velocity Verlet", Velocity_Verlet_friend}, // 2 driver evaluations
        {"Position Verlet", Position_Verlet_friend}, // 1 driver evaluation
        {"Leapfrog", Leapfrog_friend}, // 1 driver evaluation
        {"Yoshida_4", Yoshida_4_friend} // 3 driver evaluations
    };

    std::vector<std::string> driver_evaluations_1 = {"Forward Euler", "Position Verlet", "Leapfrog"};

    std::string in_cond;
    std::string integrator;
    double h;
    // int iter;
    double tmax;
    double t = 0;
    bool ADAPTIVE_TIME_STEP
;
    
    /*
    // Getting the inputs
    std::cout << "Run Start" << std::endl;
    // get the initial conditions file
    std::cout << "Provide the initial conditions: ";
    std::cin >> in_cond;
    // get the integrator to be used: must be out of the ones defined above
    std::cout << "Which integrator would you like to use: ";
    std::cin >> integrator;
    // get the timestep
    std::cout << "Provide the initial timestep: ";
    std::cin >> h;
    // get the number of iterations
    std::cout << "Provide the number of iterations: ";
    std::cin >> iter;
    */

    h = 0.001;
    // iter = 50000;
    tmax = 100;
    integrator = "RK4";
    in_cond = "rings.txt";
    ADAPTIVE_TIME_STEP = false;

    double Delta_max = pow(10, -10);
    double Delta_min = pow(10, -15);

    std::string SystemName = in_cond.substr(0, in_cond.size()-4);

    // start the execution
    int time_start = time(NULL);

    NSystem z = getvalues("Initial_conditions/" + in_cond);
    int N = z.n();

    //std::string filename= std::to_string(N)+ "_body_" + integrator + ".txt";
    std::string filename = SystemName + "_" + integrator + "_" + std::to_string(tmax) + "_" + std::to_string(h);
    if (ADAPTIVE_TIME_STEP){
        filename += "_adaptive";
    }

    std::ofstream outfile("traj/" + filename + ".txt");
    outfile << std::setprecision(8);
    outfile << t;
    for (int body_number = 0; body_number < N; body_number++){
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z() << ' ';
        }
    outfile << '\n';


    std::ofstream outfile_energy("energy/" + filename + ".txt");
    outfile_energy << std::setprecision(8);
    outfile_energy << t << ' ' << z.get_energy() << '\n';

    std::ofstream outfile_distances("temperatures/" + filename + ".txt");

    integ integrator_function = functions[integrator];
    // General_integrator integrator_function = General_integrator(integrator);

    int number_of_iterations;
    int driver_evaluations = get_driver_evaluations(integrator);
    std::cout << "Number of driver evaluations per time step is " << driver_evaluations << std::endl;

    while(t < tmax){
        t+= h;
        number_of_iterations++;

        if (ADAPTIVE_TIME_STEP){
            NSystem y = z;

            integrator_function(z, h);
            integrator_function(y, h/2);
            integrator_function(y, h/2);

            double error = compare_solutions(z, y);
            //std::cout << "For h = " << h << ", the error is " << error << std::endl;
            if (error > Delta_max) {
                h /= 2;
            } else if (error < Delta_min) {
                h *= 2;
            }
        } else {
            integrator_function(z, h);
        }

        outfile << t;
        for (int body_number = 0; body_number < N; body_number++){
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z() << ' ';
        }
        outfile << '\n';
        
        outfile_energy << t << ' ' << z.get_energy() << '\n';

        double temperature = 0;
        for (int body_number = 0; body_number < N-1; body_number++){
            temperature += 1/(z.positions()[body_number] - z.positions()[3]).norm2();
        }
        outfile_distances << t << ' ' << pow(temperature,0.25) << '\n';

    };
    
    /*
    for (int i = 0; i <= iter; i++){
        t+= h;
        number_of_iterations++;

        if (ADAPTIVE_TIME_STEP){
            NSystem y = z;

            integrator_function(z, h);
            integrator_function(y, h/2);
            integrator_function(y, h/2);

            double error = compare_solutions(z, y);
            //std::cout << "For h = " << h << ", the error is " << error << std::endl;
            if (error > Delta_max) {
                h /= 2;
            } else if (error < Delta_min) {
                h *= 2;
            }
        } else {
            integrator_function(z, h);
        }

        if (i % 10000 == 0){
            std::cout << " i = " << i << std::endl;
        }

        outfile << t;
        for (int body_number = 0; body_number < N; body_number++){
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z() << ' ';
        }
        outfile << '\n';
        
        outfile_energy << t << ' ' << z.get_energy() << '\n';

        double temperature = 0;
        for (int body_number = 0; body_number < N-1; body_number++){
            temperature += 1/(z.positions()[body_number] - z.positions()[3]).norm2();
        }
        outfile_distances << t << ' ' << pow(temperature,0.25) << '\n';

    }
    */

    int time_stop = time(NULL);
    std::cout << "Run Done" << std::endl;
    int time_diff = time_stop - time_start;
    std::cout << "Run Time = " << time_diff << " seconds" << std::endl;
    std::cout << "Total number of driver evaluations per time step is \n" << number_of_iterations << " iterations with " << driver_evaluations << " driver evaluations per time step gives " << driver_evaluations*number_of_iterations << " driver evaluations in total." << std::endl;

    return 1;
}