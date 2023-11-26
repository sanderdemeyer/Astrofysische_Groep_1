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

#include "classes.cpp"

#define ADAPTIVE_TIME_STEP true


typedef void (*integ) (NSystem&, double);

int main(){

    // Each new integrator must be added to this map
    std::unordered_map<std::string, integ> functions ={
        {"Forward Euler", Forward_Euler},
        {"RK2", RK2_step},
        {"Heun", Heun},
        {"Heun3", Heun3},
        {"Ralston", Ralston},
        {"Ralston3", Ralston3},
        {"RK3", RK3_step},
        {"RK4", RK4_step},
        {"Forest Ruth", Forest_Ruth_friend},
        {"PEFRL", PEFRL_friend},
        {"Velocity Verlet", Velocity_Verlet_friend},
        {"Position Verlet", Position_Verlet_friend},
        {"Leapfrog", Leapfrog_friend},
        {"Yoshida_4", Yoshida_4_friend}
    };

    std::string in_cond;
    std::string integrator;
    double h;
    double t = 0;
    int iter;
    
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
    iter = 100000;
    integrator = "RK4";
    // in_cond = "perturbed_criss_cross.txt";
    in_cond = "Burrau.txt";
    //in_cond = "100gauss.txt";

    double Delta_max = pow(10, -10);
    double Delta_min = pow(10, -15);

    std::string SystemName = in_cond.substr(0, in_cond.size()-4);

    // start the execution
    int time_start = time(NULL);

    NSystem z = getvalues("Initial_conditions/" + in_cond);
    int N = z.n();

    //std::string filename= std::to_string(N)+ "_body_" + integrator + ".txt";
    std::string filename = SystemName + "_" + integrator + "_" + std::to_string(iter) + "_" + std::to_string(h);
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

    for (int i = 0; i <= iter; i++){
        t+= h;

        if (ADAPTIVE_TIME_STEP){
            NSystem y = z;

            functions[integrator](z, h);
            functions[integrator](y, h/2);
            functions[integrator](y, h/2);

            double error = compare_solutions(z, y);
            std::cout << "For h = " << h << ", the error is " << error << std::endl;
            if (error > Delta_max) {
                h /= 2;
            } else if (error < Delta_min) {
                h *= 2;
            }
        } else {
            functions[integrator](z, h);
        }


        outfile << t;
        for (int body_number = 0; body_number < N; body_number++){
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z() << ' ';
        }
        outfile << '\n';
        
        outfile_energy << t << ' ' << z.get_energy() << '\n';

    }

    int time_stop = time(NULL);
    std::cout << "Run Done" << std::endl;
    int time_diff = time_stop - time_start;
    std::cout << "Run Time = " << time_diff << " seconds" << std::endl;
    return 1;
}