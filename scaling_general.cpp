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

    std::vector <std::string> functions = {"RK4", "RK6", "Wray3", "SSPRK3", "3_over_8", "Ralston4", "RK8", "RK5"};


    double h = 0.01;
    int iter = 1000;

    std::vector<std::string> files;

    for (auto& file : std::filesystem::directory_iterator{ "random_gauss" }){
            std::string name = file.path().string();
            files.push_back(name);
    };

    for ( const auto &integrator: functions) {
        // make the output files if they don't exist
        std::ofstream outfile_scaling("scaling_general/" + integrator + ".txt");
        // clear the files from the previous run if they exist
        std::ofstream outfile_clear;
        outfile_clear.open("scaling_general/" + integrator + ".txt", std::ios::out | std::ios::trunc);
    };

    for (auto file: files){
        std::cout << "Integrating: " << file << '\n' ;
        NSystem z = getvalues(file);
        for ( const auto &integrator: functions) {
            std::ofstream outfile;
            outfile.open("scaling_general/" + integrator + ".txt", std::ofstream::app);
            General_integrator integrator_function = General_integrator(integrator);

            auto start = high_resolution_clock::now();
            for (int i = 0; i <= iter; i++){
                integrator_function(z, h);
            };
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(stop - start);
            outfile << duration.count() << '\n';
        };
    };
}