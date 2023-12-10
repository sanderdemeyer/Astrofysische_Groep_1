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

    std::vector<std::string> function_list = {"Forward Euler", "RK2", "Heun", "Heun3", "Ralston", "Ralston3", "RK3", "RK4", "Forest Ruth", "PEFRL", "Velocity Verlet", "Position Verlet", "Leapfrog", "Yoshida_4"};
    // std::vector<std::string> function_list = {"RK6", "Wray3", "SSPRK3", "3_over_8", "RK8", "RK5"};

    std::vector<std::string> driver_evaluations_1 = {"Forward Euler", "Position Verlet", "Leapfrog"};

    std::string in_cond;
    std::string integrator;
    double h;
    // int iter;
    double tmax;
    double t = 0;
    bool ADAPTIVE_TIME_STEP;
    
    std::ofstream outfile("accuracy_vs_cost_adaptive_more_data.txt");
    outfile << std::setprecision(8);


    for (int j = 0; j < 14; j++) { //14
        integrator = function_list[j];
        outfile << integrator << ' ';

        //for (double exponent = 3.0; exponent < 15.5; exponent++) { // 8.5
        for (double exponent = 1.0; exponent < 7.5; exponent++) { // 8.5
            h = pow(10, -exponent/3);

            t = 0;
            tmax = 100;
            in_cond = "perturbed-criss-cross.txt";
            ADAPTIVE_TIME_STEP = true;

            double Delta_max = pow(10, -10.0+exponent);
            double Delta_min = pow(10, -12.0+exponent);

            std::string SystemName = in_cond.substr(0, in_cond.size()-4);

            // start the execution
            int time_start = time(NULL);

            NSystem z = getvalues("Initial_conditions/" + in_cond);
            int N = z.n();

            integ integrator_function = functions[integrator];
            //General_integrator integrator_function = General_integrator(integrator);

            int number_of_iterations = 0;
            int driver_evaluations = get_driver_evaluations(integrator);
            std::cout << "Number of driver evaluations per time step is " << driver_evaluations << std::endl;

            double energy_start;
            double energy_current;

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

                if (number_of_iterations == 3) {
                    energy_start = z.get_energy();
                }
                energy_current = z.get_energy();

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

            // outfile << h << ' ' << driver_evaluations*number_of_iterations << ' ' << (energy_current-energy_start)/energy_start << ' ';
            outfile << exponent << ' ' << driver_evaluations*number_of_iterations << ' ' << (energy_current-energy_start)/energy_start << ' ';

            int time_stop = time(NULL);
            std::cout << "Run Done" << std::endl;
            int time_diff = time_stop - time_start;
            std::cout << "Run Time = " << time_diff << " seconds" << std::endl;
            std::cout << "Total number of driver evaluations per time step is \n" << number_of_iterations << " iterations with " << driver_evaluations << " driver evaluations per time step gives " << driver_evaluations*number_of_iterations << " driver evaluations in total." << std::endl;
            //std::cout << "This was integrator = " << integrator << ", with h = " << h << std::endl;
            std::cout << "This was integrator = " << integrator << ", with exp = " << exponent << std::endl;
        }
        outfile << '\n';
    }

    return 1;
}