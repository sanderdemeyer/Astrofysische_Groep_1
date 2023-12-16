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

#include "classes_regularization.cpp"

namespace fs = std::filesystem;


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
    double tmax;
    bool ADAPTIVE_TIME_STEP;

    // **** Parameters to be chosen - Start **** //
    h = 0.001; // The (initial) timestep to be used for integration
    tmax = 300; // The total time to be simulated
    integrator = "PEFRL"; // the integrator used
    in_cond = "two-body-plane.txt"; // The initial conditions
    ADAPTIVE_TIME_STEP = false; // whether or not to use dapative timestep
    int i = 0; // the number of iterations
    double transform_distance = 0.5; // The distance below which the system is regularized

    double Delta_max = pow(10, -10); // Parameter for when ADAPTIVE_TIME_STEP = true
    double Delta_min = pow(10, -15); // Parameter for when ADAPTIVE_TIME_STEP = true
    // **** Parameters to be chosen - End **** //


    bool regularized = false; // whether or not the system should be regularized
    double dtau = h/transform_distance; // dtau of the regularized coordinates

    std::string SystemName = in_cond.substr(0, in_cond.size()-4);

    // start the execution
    int time_start = time(NULL);

    NSystem z_help = getvalues("Initial_conditions/" + in_cond); // Define the system based on the initial conditions
    int N = z_help.n(); // System size
    std::vector<double> initial_masses; // Initialize the masses of the regularized system to random values
    Regularized_coo_2D regul_coo; // Initialize the regularized coordinates of the regularized system to random values
    NSystem_reg_2D z = NSystem_reg_2D(z_help, false, regul_coo, initial_masses);

    int NoB = z.nsystem().positions().size(); // NoB = Number of Bodies

    // Define the relevant files for the output
    std::string filename = "Reg_2D_" + SystemName + "_" + integrator + "_" + std::to_string(tmax) + "_" + std::to_string(h) + "_" + std::to_string(transform_distance);
    if (ADAPTIVE_TIME_STEP){
        filename += "_adaptive";
    }

    fs::create_directories("traj_reg");

    std::ofstream outfile("traj_reg/" + filename + ".txt");
    outfile << std::setprecision(8);
    outfile << t;
    for (int body_number = 0; body_number < N; body_number++){
            outfile << ' ' << z.nsystem().positions()[body_number].x() << ' ' << z.nsystem().positions()[body_number].y() << ' ' << z.nsystem().positions()[body_number].z() << ' ';
        }
        outfile << '\n';

    fs::create_directories("energy_reg");
    std::ofstream outfile_energy("energy_reg/" + filename + ".txt");
    outfile_energy << std::setprecision(8);
    outfile_energy << 0 << ' ' << z.nsystem().get_energy() << '\n';

    Vec u, r;
    std::vector<int> should_be_regularized; // This is a vector of 3 elements. The first element denotes whether or not the system should be regularized. If this is 1, the 2 other values indicate which bodies should be regularized.

    while(t < tmax){
        t += h;
        i++;

        // Check whether the system should be regularized
        should_be_regularized = z.check_separation(transform_distance);

        // Based on the above, either transform the system forward or backward, or do nothing.
        if (should_be_regularized[0] && (!regularized)){
            z = z.transform_forward(should_be_regularized[1], should_be_regularized[2]);
            regularized = !regularized;
        }
        else if ((!should_be_regularized[0]) && (regularized)){
            // Because of this construction, the labels of the different bodies can change upon regularization.
            z = z.transform_backward();
            regularized = !regularized;
        }

        // Perform a time step, based on whether an adaptive timestep is used or not.
        if ((ADAPTIVE_TIME_STEP) && (!regularized)) {
            NSystem_reg_2D y = z;

            z.timestep(h, dtau);
            y.timestep(h/2, dtau);
            y.timestep(h/2, dtau);

            double error = compare_solutions(z, y);
            if (error > Delta_max) {
                h /= 2;
            } else if (error < Delta_min) {
                h *= 2;
            }
        } else {
            z.timestep(h, dtau);
        }

        // Output the energy if the system is not regularized.
        if (!regularized){
            outfile << t;
            for (int body_number = 0; body_number < NoB; body_number++){
                outfile << ' ' << z.nsystem().positions()[body_number].x() << ' ' << z.nsystem().positions()[body_number].y() << ' ' << z.nsystem().positions()[body_number].z();
            }
            outfile << '\n';
            
            outfile_energy << t << ' ' << z.nsystem().get_energy() << '\n';
        } 
        
        // Every 1000 iterations, the number of iterations, together with the current value of t and h are displayed. This serves as an indication of how far the simulation has progressed.
        if (i % 1000 == 0){
            std::cout << "iterations = " << i << ", t = " << t << ", h = " << h << std::endl;
        }

    }

    int time_stop = time(NULL);
    std::cout << "Run Done" << std::endl;
    int time_diff = time_stop - time_start;
    std::cout << "Run Time = " << time_diff << " seconds" << std::endl;
    return 1;
}