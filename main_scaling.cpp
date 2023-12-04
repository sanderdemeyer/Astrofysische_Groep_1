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
        {"Yoshida_4", Yoshida_4_friend}
    };

    std::string integrator;
    double h;
    double tmax;
    double t = 0;

    h = 0.01;
    tmax = 100;
    integrator = "RK4";

    std::ofstream outfile_scaling("scaling/" + integrator + ".txt");
    for (int i=2; i<100; i++){
        std::string help = std::to_string(i);
        std::string file = help + ".txt";
        NSystem z = getvalues("random_gauss/" + file);
        int iterataions = tmax/h;
        int time_start = time(NULL);
        while(t < tmax){
            t += h;
            functions[integrator](z, h);
        };
        int time_stop = time(NULL);
        int time_diff = time_stop - time_start;
        int t_per_step = time_diff/iterataions;
        outfile_scaling << i << ' ' << t_per_step << '\n';
    }
}