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
#include <list>

#include "classes.cpp"

namespace fs = std::filesystem;


typedef void (*integ) (NSystem&, double);

/**
 * This function takes some given N-body initial conditions and calculates their trajectories and total energy for a given maximum time using a given integrator and initial timestep. All integrators can be used with an adaptive timestep and if ADAPTIVE_RK45 is ture RK45 will be used.
 *
 * @param in_cond A string with the filename of the N-body initial condition to be read from the `initial_conditions` folder. The format of the file has been explained in the README file.
 * @param integrator A string with the name of the integrator to be used. This should be one of the friend void integrators defined in `classes.cpp`. The available integrators are listed in the README file.
 * @param h The (initial) timestep to be used for integration
 * @param tmax The total time to be simulated
 * @param ADAPTIVE_TIME_STEP Whether or not to use an adaptive timestep. This works for all integrators and doubles or halves the timestep depending on the error.
 * @param ADAPTIVE_RK45 When this is true RK45 with embedded adaptive timestep will be used. In order for this to run, `ADAPTIVE_TIME_STEP` must be set to false.
 * @return The files with the trajectories and energies are saved in the `traj` and `energy` folder respectively. The format of the file has been explained in the README file. The naming convetion is:
 * `SystemName_integrator_tmax_h(_adaptive).txt`
 */
void integrate (std::string in_cond, std::string integrator, double h, double tmax, bool ADAPTIVE_TIME_STEP, bool ADAPTIVE_RK45){
    // Each new integrator must be added to this map
    std::unordered_map<std::string, integ> functions = {
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
        {"Yoshida 4", Yoshida_4_friend} // 3 driver evaluations
    };

    double t = 0; // Starting time

    double Delta_max = pow(10, -10); // Parameter for when ADAPTIVE_TIME_STEP = true
    double Delta_min = pow(10, -11); // Parameter for when ADAPTIVE_TIME_STEP = true

    std::string SystemName = in_cond.substr(0, in_cond.size()-4); // System name based on initial conditions

    // Create folders for trajectories and energies of the system
    std::string outdir_traj = "traj/" + SystemName;
    fs::create_directories(outdir_traj);

    std::string outdir_energy = "energy/" + SystemName;
    fs::create_directories(outdir_energy);

    // start the execution
    int time_start = time(NULL);

    NSystem z = getvalues("Initial_conditions/" + in_cond); // Create an object of the class NSystem that correctly initialises the different bodies, based on the initial conditions.
    int N = z.n(); // N is the number of bodies

    std::string filename = SystemName + "_" + integrator + "_" + std::to_string(tmax) + "_" + std::to_string(h);
    if (ADAPTIVE_TIME_STEP){
        filename += "_adaptive"; // If ADAPTIVE_TIME_STEP = true, this is added to the file name.
    }

    std::ofstream outfile(outdir_traj + "/" + filename + ".txt");
    outfile << std::setprecision(12);
    outfile << t;
    for (int body_number = 0; body_number < N; body_number++){
        // add the positions to the trajectory file
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z() << ' ';
        }
    outfile << '\n';


    std::ofstream outfile_energy(outdir_energy + "/" + filename + ".txt");
    outfile_energy << std::setprecision(16);
    outfile_energy << t << ' ' << z.get_energy() << '\n';

    std::ofstream outfile_distances("temperatures/" + filename + ".txt"); // This is the file where a measure of temperature will be saved.

    integ integrator_function = functions[integrator]; // Select the integrator.

    int number_of_iterations = 0; // The number of iterations
    int driver_evaluations = get_driver_evaluations(integrator); // The number of driver evaluations, to be used in the analysis of the accuracy vs cost plot
    std::cout << "Number of driver evaluations per time step is " << driver_evaluations << std::endl;

    while(t < tmax){
        t+= h; // Progress the time
        number_of_iterations++;

        if (ADAPTIVE_TIME_STEP){ // Progress the number of iterations
            NSystem y = z; // Make a copy of the system, such that the system can be updated twice and independently. This is necessary to compute the error.

            integrator_function(z, h); // Update the system once with time step h
            integrator_function(y, h/2); // Update the system twice with time step h/2
            integrator_function(y, h/2);

            double error = compare_solutions(z, y); // Compute the error between the two solutions
            if (error > Delta_max) {
                h /= 2; // If the error is too large, the timestep is halved.
            } else if (error < Delta_min) {
                h *= 2; // If the error is too small, the timestep is doubled.
            } // If neither is true, the timestep remains fixed.
        } else if (ADAPTIVE_RK45) {
            RK45_step(z, h, 1e-5); // This updates both the system and the timestep using the RK45 method. 1e-5 is the error parameter.
        } else {
            integrator_function(z, h); // Using a fixed timestep, the system z is updated over a time h.
        }

        // The following lines write the positions of the different bodies to the correct file.
        outfile << t;
        for (int body_number = 0; body_number < N; body_number++){
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z() << ' ';
        }
        outfile << '\n';
        
        outfile_energy << t << ' ' << z.get_energy() << '\n'; // This writes the total energy of all bodies together to the correct file.

        // The following lines write the measure of temperature of the planet to the correct file.
        // This assumes that body 3 is a planet (of which the temperature is calculated), and that bodies 0,1, and 2 are suns.
        // Only relevent for certain initial conditions. SHould be ignored.
        double temperature = 0;
        for (int body_number = 0; body_number < N-1; body_number++){
            temperature += 1/(z.positions()[body_number] - z.positions()[3]).norm2();
        }
        outfile_distances << t << ' ' << pow(temperature,0.25) << '\n';

        // Every 1000 iterations, the number of iterations, together with the current value of t and h are displayed. This serves as an indication of how far the simulation has progressed.
        if (number_of_iterations % 1000 == 0){
            std::cout << "iterations = " << number_of_iterations << ", t = " << t << ", h = " << h << std::endl;
        }

    };
    // The following lines state that the simulation has finished and give some simple details on the total run time and the number of iterations and driver evaluations.
    int time_stop = time(NULL);
    std::cout << "Run Done" << std::endl;
    int time_diff = time_stop - time_start;
    std::cout << "Run Time = " << time_diff << " seconds" << std::endl;
    std::cout << "Total number of driver evaluations per time step is \n" << number_of_iterations << " iterations with " << driver_evaluations << " driver evaluations per time step gives " << driver_evaluations*number_of_iterations << " driver evaluations in total." << std::endl;
}

/**
 * This function takes some given N-body initial conditions and calculates their trajectories and total energy for a given maximum time using a given integrator and initial timestep. All integrators can be used with an adaptive timestep. This is essentially the same function as `integrate` with the only difference being that here the integrators are defined using Butcher tableau's.
 *
 * @param in_cond A string with the filename of the N-body initial condition to be read from the `initial_conditions` folder. The format of the file has been explained in the README file.
 * @param integrator A string with the name of the integrator to be used. This should be one of the `General_integrator` class defined in `classes.cpp`. The available integrators are listed in the README file.
 * @param h The (initial) timestep to be used for integration
 * @param tmax The total time to be simulated
 * @param ADAPTIVE_TIME_STEP Whether or not to use an adaptive timestep. This works for all integrators and doubles or halves the timestep depending on the error.
 * @return The files with the trajectories and energies are saved in the `traj` and `energy` folder respectively. The format of the file has been explained in the README file. The naming convetion is:
 * `SystemName_integrator-general_tmax_h(_adaptive).txt`
 */
void integrate_general (std::string in_cond, std::string integrator, double h, double tmax, bool ADAPTIVE_TIME_STEP){
    double t = 0; // Starting time

    double Delta_max = pow(10, -10); // Parameter for when ADAPTIVE_TIME_STEP = true
    double Delta_min = pow(10, -11); // Parameter for when ADAPTIVE_TIME_STEP = true

    std::string SystemName = in_cond.substr(0, in_cond.size()-4); // Systemname based on initial conditions

    // Create folders for trajectories and energies of the system
    std::string outdir_traj = "traj/" + SystemName;
    fs::create_directories(outdir_traj);

    std::string outdir_energy = "energy/" + SystemName;
    fs::create_directories(outdir_energy);

    // start the execution
    int time_start = time(NULL);

    NSystem z = getvalues("Initial_conditions/" + in_cond); // Create an object of the class NSystem that correctly initialises the different bodies, based on the initial conditions.
    int N = z.n(); // N is the number of bodies

    std::string filename = SystemName + "_" + integrator + "-general" + "_" + std::to_string(tmax) + "_" + std::to_string(h);
    if (ADAPTIVE_TIME_STEP){
        filename += "_adaptive"; // If ADAPTIVE_TIME_STEP = true, this is added to the file name.
    }

    std::ofstream outfile(outdir_traj +"/" + filename + ".txt");
    outfile << std::setprecision(12);
    outfile << t;
    for (int body_number = 0; body_number < N; body_number++){
            // add the positions to the trajectory file
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z() << ' ';
        }
    outfile << '\n';


    std::ofstream outfile_energy(outdir_energy +"/" + filename + ".txt");
    outfile_energy << std::setprecision(16);
    outfile_energy << t << ' ' << z.get_energy() << '\n';

    std::ofstream outfile_distances("temperatures/" + filename + ".txt"); // This is the file where a measure of temperature will be saved.

    General_integrator integrator_function = General_integrator(integrator); // Get the integrator function.

    int number_of_iterations = 0; // The number of iterations
    int driver_evaluations = get_driver_evaluations(integrator); // The number of driver evaluations, to be used in the analysis of the accuracy vs cost plot
    std::cout << "Number of driver evaluations per time step is " << driver_evaluations << std::endl;

    while(t < tmax){
        t+= h; // Progress the time
        number_of_iterations++; // Progress the number of iterations

        if (ADAPTIVE_TIME_STEP){
            NSystem y = z; // Make a copy of the system, such that the system can be updated twice and independently. This is necessary to compute the error.

            integrator_function(z, h); // Update the system once with time step h
            integrator_function(y, h/2); // Update the system twice with time step h/2
            integrator_function(y, h/2); // Compute the error between the two solutions

            double error = compare_solutions(z, y);
            if (error > Delta_max) {
                h /= 2; // If the error is too large, the timestep is halved.
            } else if (error < Delta_min) {
                h *= 2; // If the error is too small, the timestep is doubled.
            } // If neither is true, the timestep remains fixed.
        }else {
            integrator_function(z, h);  // Using a fixed timestep, the system z is updated over a time h.
        }

        // The following lines write the positions of the different bodies to the correct file.
        outfile << t;
        for (int body_number = 0; body_number < N; body_number++){
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z() << ' ';
        }
        outfile << '\n';
        
        outfile_energy << t << ' ' << z.get_energy() << '\n'; // This writes the total energy of all bodies together to the correct file.

        // The following lines write the measure of temperature of the planet to the correct file.
        // This assumes that body 3 is a planet (of which the temperature is calculated), and that bodies 0,1, and 2 are suns.
        // Only to be used for certain initial conditions. Should be ignored for others.
        double temperature = 0;
        for (int body_number = 0; body_number < N-1; body_number++){
            temperature += 1/(z.positions()[body_number] - z.positions()[3]).norm2();
        }
        outfile_distances << t << ' ' << pow(temperature,0.25) << '\n';

        // Every 1000 iterations, the number of iterations, together with the current value of t and h are displayed. This serves as an indication of how far the simulation has progressed.
        if (number_of_iterations % 1000 == 0){
            std::cout << "iterations = " << number_of_iterations << ", t = " << t << ", h = " << h << std::endl;
        }

    };

    // The following lines state that the simulation has finished and give some simple details on the total run time and the number of iterations and driver evaluations.
    int time_stop = time(NULL);
    std::cout << "Run Done" << std::endl;
    int time_diff = time_stop - time_start;
    std::cout << "Run Time = " << time_diff << " seconds" << std::endl;
    std::cout << "Total number of driver evaluations per time step is \n" << number_of_iterations << " iterations with " << driver_evaluations << " driver evaluations per time step gives " << driver_evaluations*number_of_iterations << " driver evaluations in total." << std::endl;
}

/**
 * This function takes some given N-body initial conditions and runs the `integrate` function for fixed timestep using different timesteps in the given range.
 *
 * @param in_cond A string with the filename of the N-body initial condition to be read from the `initial_conditions` folder. The format of the file has been explained in the README file.
 * @param integrator A string with the name of the integrator to be used. This should be one of the friend void integrators defined in `classes.cpp`. The available integrators are listed in the README file.
 * @param tmax The total time to be simulated
 * @param hmin The minimum timestep to be used
 * @param hmax The maximum timestep to be used
 * @param step The factor by which to loop from `hmin` to `hmax`. Should be greater than 1.
 * @return The files with the trajectories and energies are saved in the `traj` and `energy` folder respectively. The format of the file has been explained in the README file. The naming convetion is:
 * `SystemName_integrator_tmax_h(_adaptive).txt`
 */
void loop_h (std::string in_cond, std::string integrator, double tmax, double hmin, double hmax, int step){
    std::list<double> H;
    double t = hmin;
    while (t <= hmax)
    {
        integrate (in_cond, integrator, t, tmax, false, false);
        t *= step;
    }
}

/**
 * This function takes some given N-body initial conditions and runs the `integrate_general` function for fixed timestep using different timesteps in the given range.
 *
 * @param in_cond A string with the filename of the N-body initial condition to be read from the `initial_conditions` folder. The format of the file has been explained in the README file.
 * @param integrator A string with the name of the integrator to be used. This should be one of the friend void integrators defined in `classes.cpp`. The available integrators are listed in the README file.
 * @param tmax The total time to be simulated
 * @param hmin The minimum timestep to be used
 * @param hmax The maximum timestep to be used
 * @param step The factor by which to loop from `hmin` to `hmax`. Should be greater than 1.
 * @return The files with the trajectories and energies are saved in the `traj` and `energy` folder respectively. The format of the file has been explained in the README file. The naming convetion is:
 * `SystemName_integrator_tmax_h(_adaptive).txt`
 */
void loop_h_general (std::string in_cond, std::string integrator, double tmax, double hmin, double hmax, int step){
    std::list<double> H;
    double t = hmin;
    while (t <= hmax)
    {
        integrate_general (in_cond, integrator, t, tmax, false);
        t *= step;
    }
}

int main(){
    std::string in_cond = "Burrau.txt"; // File with the initial conditions to be read from the `Initial_conditions` folder
    std::string type_integ = "general"; // The type of integrator to be used. By default the `friend void` type integrators are used.
    std::string integrator = "RK4"; // he type of integrator to be used. For each type, available integrators are listed in the README file.
    double h = 0.001; // The (initial) timestep.
    double tmax = 70; // Total time considered in the simulation
    bool ADAPTIVE_TIME_STEP = true; // Whether or not to use an adaptive timestep
    bool ADAPTIVE_RK45 = false; // Whether or not to use RK45 embedded adaptive timestep. If you want to use this, ADAPTIVE_TIME_STEP should be set to false. Only implemented in the `integrate` function.

    if (type_integ == "general"){
        integrate_general(in_cond, integrator, h, tmax, ADAPTIVE_TIME_STEP);
    } else{
        integrate(in_cond, integrator, h, tmax, ADAPTIVE_TIME_STEP, ADAPTIVE_RK45);
    }

    // This can be commented out to test the `loop_h` function
    // loop_h(in_cond, integrator, tmax, h, 0.1, 10);
}