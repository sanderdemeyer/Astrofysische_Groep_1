#define _USE_MATH_DEFINES
#include <cmath>
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
#include <fstream>
#include <chrono>

#include "classes.cpp"

using namespace std::chrono;

typedef void (*integ) (NSystem&, double);

int main(){

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
        {"Yoshida 4", Yoshida_4_friend}
    };

    double h = 0.01;
    int iter = 1000;
    bool ADAPTIVE = true;

    double Delta_max = pow(10, -10); // Parameter for when ADAPTIVE_TIME_STEP = true
    double Delta_min = pow(10, -11); // Parameter for when ADAPTIVE_TIME_STEP = true

    std::vector<std::string> files;

    for (auto& file : std::filesystem::directory_iterator{ "random_gauss" }){
            std::string name = file.path().string();
            files.push_back(name);
    };

    for ( const auto &[integrator, function]: functions) {
        if (ADAPTIVE){
            // make the output files if they don't exist
            std::ofstream outfile_scaling("scaling_adaptive_friend/" + integrator + ".txt");
            // clear the files from the previous run if they exist
            std::ofstream outfile_clear;
            outfile_clear.open("scaling_adaptive_friend/" + integrator + ".txt", std::ios::out | std::ios::trunc);
        } else {
            // make the output files if they don't exist
            std::ofstream outfile_scaling("scaling_friend/" + integrator + ".txt");
            // clear the files from the previous run if they exist
            std::ofstream outfile_clear;
            outfile_clear.open("scaling_friend/" + integrator + ".txt", std::ios::out | std::ios::trunc);
        }
    };

    std::string adaptive = "friend/";
    if (ADAPTIVE){
        adaptive = "adaptive_friend/";
    };

    for (auto file: files){
        std::cout << "Integrating: " << file << '\n' ;
        NSystem z = getvalues(file);
        for ( const auto &[integrator, function]: functions) {
            std::ofstream outfile;
            outfile.open("scaling_" + adaptive + integrator + ".txt", std::ofstream::app);

            auto start = high_resolution_clock::now();

            for (int i = 0; i <= iter; i++){
                if (ADAPTIVE){
                    NSystem y = z; // Make a copy of the system, such that the system can be updated twice and independently. This is necessary to compute the error.

                    function(z, h); // Update the system once with time step h
                    function(y, h/2); // Update the system twice with time step h/2
                    function(y, h/2);

                    double error = compare_solutions(z, y); // Compute the error between the two solutions
                    if (error > Delta_max) {
                        h /= 2; // If the error is too large, the timestep is halved.
                    } else if (error < Delta_min) {
                        h *= 2; // If the error is too small, the timestep is doubled.
                    } // If neither is true, the timestep remains fixed.
                } else {
                    function(z, h);
                }
            };
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        outfile << duration.count() << '\n';
        };
    };
}