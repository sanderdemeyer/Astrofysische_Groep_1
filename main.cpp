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

#include "classes.cpp"

// --------------------TO DO--------------------
// switches (maak de main, plot, ... generiek): 
//      integrator -> wordt ook gebruikt in naam output (incl. banen, energies, runtime, ...)
// Nieuwe integratoren (VEEEEL)
// Niet alle timesteps plotten (overkill)
// Variable timestep

typedef void (*integ) (NSystem&, double);

int main(){

    // Each new integrator must be added to this map
    std::unordered_map<std::string, integ> functions ={
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
    iter = 40000;
    integrator = "RK4";
    in_cond = "initial_conditions.txt";

   
    // start the execution
    int time_start = time(NULL);

    NSystem z = getvalues(in_cond);
    int N = z.n();
    std::string filename= std::to_string(N)+ "_body_" + integrator + ".txt";

    std::ofstream outfile("traj/" + filename);
    outfile << std::setprecision(8);
    outfile << t;
    for (int body_number = 0; body_number < N; body_number++){
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z() << ' ';
        }
        outfile << '\n';


    std::ofstream outfile_energy("energy/" + filename);
    outfile_energy << std::setprecision(8);
    outfile_energy << t << ' ' << z.get_energy() << '\n';

    for (int i = 0; i <= iter; i++){
        t+= h;
        //z = RK4_step(z, h);
        //z = functions[integrator](z, h);
        functions[integrator](z, h);

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