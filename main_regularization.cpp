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

    std::ofstream outfile("two_body_motion_RK4.txt");
    outfile << std::setprecision(8);
    std::ofstream outfile_energy("two_body_motion_energy_RK4.txt");
    outfile_energy << std::setprecision(8);

    double h = 0.001;

    NSystem z_help = getvalues("initial_conditions_two_body.txt");
    std::vector<double> initial_masses;
    Regularized_coo regul_coo;
    NSystem_reg z = NSystem_reg(z_help, false, regul_coo, initial_masses);
    //std::vector<Vec> x = z.positions();
    //std::vector<Vec> v = z.velocities();
    //std::vector<double> masses = z.masses();

    int NoB = z.nsystem().positions().size(); // NoB = Number of Bodies

    bool regularized = false;
    double transform_distance = 0.9;
    double dtau = h/transform_distance;
    Vec u, r;

    for (int i = 0; i < 100000; i++){

        bool should_be_regularized = z.check_separation(transform_distance);
        if (should_be_regularized && (!regularized)){
            z = z.transform_forward(0, 1);
            regularized = !regularized;
        }
        else if ((!should_be_regularized) && (regularized)){
            z = z.transform_backward();
            regularized = !regularized;
        }
        //z = RK4_step(z, h);
        //Yoshida_4_new(x, v, masses, h);
        //PEFRL_friend(z, h);
        Leapfrog_reg(z, h);
        //RK4_step(z, h);
        //Yoshida_friend(z, h);

        outfile << i;
        for (int body_number = 0; body_number < NoB; body_number++){
            outfile << ' ' << z.nsystem().positions()[body_number].x() << ' ' << z.nsystem().positions()[body_number].y() << ' ' << z.nsystem().positions()[body_number].z();
            //outfile << ' ' << x[body_number].x() << ' ' << x[body_number].y() << ' ' << x[body_number].z();
        }
        outfile << '\n';
        
        outfile_energy << z.nsystem().get_energy() << '\n';

        if (i % 1000 == 0){
            std::cout << regularized << should_be_regularized << z.nsystem().positions().size() << std::endl;
        }
    }

    int time_stop = time(NULL);
    std::cout << "Run Done" << std::endl;
    int time_diff = time_stop - time_start;
    std::cout << "Run Time = " << time_diff << " seconds" << std::endl;
    return 1;
}