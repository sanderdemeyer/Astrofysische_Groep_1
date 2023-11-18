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

typedef NSystem (*integ) (NSystem, double);

int main(){

    // Each new integrator must be added to this map
    std::map<std::string, integ> functions ={
        {"RK4", RK4_step},
        {"Forest Ruth", Forest_Ruth},
        {"PEFRL", PEFRL},
        {"Velocity Verlet", Velocity_Verlet},
        {"Position Verlet", Position_Verlet},
        {"Leapfrog", Leapfrog},
        {"Yoshida_4", Yoshida_4}
    };

    std::string in_cond;
    std::string integrator;
    double h;
    double t = 0;
    int iter;
    
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

    // start the execution
    int time_start = time(NULL);

    NSystem z = getvalues(in_cond);
    int N = z.n();
    std::string filename= std::to_string(N)+ "_body_" + integrator + ".txt";

    std::ofstream outfile("traj/" + filename);
    outfile << std::setprecision(8);
    outfile << t;
    for (int body_number = 0; body_number < N; body_number++){
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z();
        }
        outfile << '\n';


    std::ofstream outfile_energy("energy/" + filename);
    outfile_energy << std::setprecision(8);
    outfile_energy << t << ' ' << z.get_energy() << '\n';

    double h = 0.001;

    NSystem z = getvalues("initial_conditions.txt");

    for (int i = 0; i < 30000; i++){

        //z = RK4_step(z, h);
        z = Yoshida_4(z, h);

        outfile << i;
        for (int body_number = 0; body_number < 3; body_number++){
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z();
            //outfile << ' ' << x[body_number].x() << ' ' << x[body_number].y() << ' ' << x[body_number].z();
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