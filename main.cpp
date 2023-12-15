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

typedef void (*integ) (NSystem&, double);


int main(){

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

    std::string in_cond; // initial conditions
    std::string integrator; // integrator used for the timesteps
    double h; // dt
    // int iter; // Number of iterations. By default, tmax is used instead
    double tmax; // Total time considered in the simulation
    double t = 0; // Starting time
    bool ADAPTIVE_TIME_STEP; // Whether or not to use an adaptive timestep
    bool ADAPTIVE_RK45; // Whether or not to use RK45 embedded adaptive timestep. If you want to use this, ADAPTIVE_TIME_STEP should be set to false.

    h = 0.001;
    // iter = 50000;
    tmax = 1000;
    integrator = "PEFRL";
    in_cond = "Solar-System-New.txt";
    ADAPTIVE_TIME_STEP = false;
    ADAPTIVE_RK45 = false; // In order for this to be run, ADAPTIVE_TIME_STEP must be set to false.

    double Delta_max = pow(10, -10); // Parameter for when ADAPTIVE_TIME_STEP = true
    double Delta_min = pow(10, -11); // Parameter for when ADAPTIVE_TIME_STEP = true

    std::string SystemName = in_cond.substr(0, in_cond.size()-4); // Systemname based on initial conditions

    // start the execution
    int time_start = time(NULL);

    NSystem z = getvalues("Initial_conditions/" + in_cond); // Create an object of the class NSystem that correctly initialises the different bodies, based on the initial conditions.
    int N = z.n(); // N is the number of bodies

    std::string filename = SystemName + "_" + integrator + "_" + std::to_string(tmax) + "_" + std::to_string(h);
    if (ADAPTIVE_TIME_STEP){
        filename += "_adaptive"; // If ADAPTIVE_TIME_STEP = true, this is added to the file name.
    }

    std::ofstream outfile("traj/" + filename + ".txt");
    outfile << std::setprecision(12);
    outfile << t;
    for (int body_number = 0; body_number < N; body_number++){ // add the positions to the trajectory file
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z() << ' ';
        }
    outfile << '\n';


    std::ofstream outfile_energy("energy/" + filename + ".txt");
    outfile_energy << std::setprecision(16);
    outfile_energy << t << ' ' << z.get_energy() << '\n';

    std::ofstream outfile_temps("temperatures/" + filename + ".txt"); // This is the file where a measure of temperature will be saved.
    std::ofstream outfile_distance_planets("dist_mercury/" + filename + ".txt");
    outfile_distance_planets << std::setprecision(8); // This is the file where the distances of all bodies to the first body in function of time will be saved. This is included for the analysis of the solar system 

    // ***** Important: If the integrator is implemented using functions, the second line should be commented out.                   ***** //
    // ***** Important: If the integrator is implemented using the class General Integrator, the first line should be commented out. ***** //
    integ integrator_function = functions[integrator];
    // General_integrator integrator_function = General_integrator(integrator);

    int number_of_iterations = 0; // The number of iterations
    int driver_evaluations = get_driver_evaluations(integrator); // The number of driver evaluations, to be used in the analysis of the accuracy vs cost plot
    std::cout << "Number of driver evaluations per time step is " << driver_evaluations << std::endl;

    Vec sun_mercury_start = Vec(1, 0, 0); // For the analysis of the perihelion precession of mercury: 
    //sun_mercury_start is the vector with respect to which the vectors connecting mercury and the sun (in function of time) will be calculated. 
    // Using this vector, angles can be calculated, which will be necessary to calculate the angle for which the perihelion occurs.
    // This line initializes this variable, but will later be overwritten. 

    while(t < tmax){
        t+= h; // Progress the time
        number_of_iterations++; // Progress the number of iterations

        if (ADAPTIVE_TIME_STEP){
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
        double temperature = 0;
        for (int body_number = 0; body_number < N-1; body_number++){
            temperature += 1/(z.positions()[body_number] - z.positions()[3]).norm2();
        }
        outfile_temps << t << ' ' << pow(temperature,0.25) << '\n';

        if (number_of_iterations == 1) {
            sun_mercury_start = z.positions()[1]-z.positions()[0]; // Overwrite the initial value (Forthe analysis of the perihelion precession of mercury)
        }

        // The following lines write the positions to the sun and their angles wrt sun_mercury_start to the correct file.
        outfile_distance_planets << t << ' '; 
        for (int body_number = 1; body_number < N; body_number++){
            outfile_distance_planets << (z.positions()[body_number]-z.positions()[0]).norm() << ' ' << (z.positions()[body_number]-z.positions()[0]).angle(sun_mercury_start) << ' ';
        }
        outfile_distance_planets << std::endl;

        // Every 1000 iterations, the number of iterations, together with the current value of t and h are displayed. This serves as an indication of how far the simulation has progressed.
        if (number_of_iterations % 1000 == 0){
            std::cout << "iterations = " << number_of_iterations << ", t = " << t << ", h = " << h << std::endl;
        }
    };
    
    // ***** These lines can be commented out (and the lines above should be commented) when you do not want to work with a total simulation time,
    // but rather with a fixed total number of iterations. This is not recommended.
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
        outfile_temps << t << ' ' << pow(temperature,0.25) << '\n';

    }
    */

    // The following lines state that the simulation has finished and give some simple details on the total run time and the number of iterations and driver evaluations.
    int time_stop = time(NULL);
    std::cout << "Run Done" << std::endl;
    int time_diff = time_stop - time_start;
    std::cout << "Run Time = " << time_diff << " seconds" << std::endl;
    std::cout << "Total number of driver evaluations per time step is \n" << number_of_iterations << " iterations with " << driver_evaluations << " driver evaluations per time step gives " << driver_evaluations*number_of_iterations << " driver evaluations in total." << std::endl;

    return 1;
}