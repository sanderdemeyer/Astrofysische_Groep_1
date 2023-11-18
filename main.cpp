#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream> 
#include <ctime>
#include <vector>

#include "classes.cpp"

// --------------------TO DO--------------------
// switches (maak de main, plot, ... generiek): 
//      integrator -> wordt ook gebruikt in naam output (incl. banen, energies, runtime, ...)
// Nieuwe integratoren (VEEEEL)
// Niet alle timesteps plotten (overkill)
// Variable timestep

int main(){
    
    std::cout << "Run Start" << std::endl;
    int time_start = time(NULL);

    std::ofstream outfile("three_body_motion_RK4.txt");
    outfile << std::setprecision(8);
    std::ofstream outfile_energy("three_body_motion_energy_RK4.txt");
    outfile_energy << std::setprecision(8);

    double h = 0.001;

    NSystem z = getvalues("initial_conditions.txt");
    //std::vector<Vec> x = z.positions();
    //std::vector<Vec> v = z.velocities();
    //std::vector<double> masses = z.masses();

    for (int i = 0; i < 100000; i++){

        //z = RK4_step(z, h);
        //Yoshida_4_new(x, v, masses, h);
        z = Yoshida_4(z, h);
        //Yoshida_friend(z, h);

        outfile << i;
        for (int body_number = 0; body_number < 3; body_number++){
            outfile << ' ' << z.positions()[body_number].x() << ' ' << z.positions()[body_number].y() << ' ' << z.positions()[body_number].z();
            //outfile << ' ' << x[body_number].x() << ' ' << x[body_number].y() << ' ' << x[body_number].z();
        }
        outfile << '\n';
        
        outfile_energy << z.get_energy() << '\n';

    }

    int time_stop = time(NULL);
    std::cout << "Run Done" << std::endl;
    int time_diff = time_stop - time_start;
    std::cout << "Run Time = " << time_diff << " seconds" << std::endl;
    return 1;
}