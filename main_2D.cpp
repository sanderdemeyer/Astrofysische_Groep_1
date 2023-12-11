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

#define ADAPTIVE_TIME_STEP false


typedef void (*integ) (NSystem&, double);

int main(){

    // Each new integrator must be added to this map
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
    iter = 1000000;
    integrator = "RK4";
    // in_cond = "perturbed_criss_cross.txt";
    in_cond = "Solar-System.txt";
    in_cond = "two-body-plane.txt";
    in_cond = "Burrau_scaled.txt";
    in_cond = "collision.txt";
    //in_cond = "100gauss.txt";

    double Delta_max = pow(10, -10);
    double Delta_min = pow(10, -15);

    bool regularized = false;
    double transform_distance = 0.5;
    double dtau = h/transform_distance;

    std::string SystemName = in_cond.substr(0, in_cond.size()-4);

    // start the execution
    int time_start = time(NULL);

    NSystem z_help = getvalues("Initial_conditions/" + in_cond);
    int N = z_help.n();
    std::vector<double> initial_masses;
    Regularized_coo_2D regul_coo;
    NSystem_reg_2D z = NSystem_reg_2D(z_help, false, regul_coo, initial_masses);

    int NoB = z.nsystem().positions().size(); // NoB = Number of Bodies
    std::cout << "size is " << NoB << std::endl;

    //std::string filename= std::to_string(N)+ "_body_" + integrator + ".txt";
    std::string filename = "Reg_2D_" + SystemName + "_" + integrator + "_" + std::to_string(iter) + "_" + std::to_string(h) + "_" + std::to_string(transform_distance);
    if (ADAPTIVE_TIME_STEP){
        filename += "_adaptive";
    }

    std::ofstream outfile("traj_reg/" + filename + ".txt");
    outfile << std::setprecision(8);
    outfile << t;
    for (int body_number = 0; body_number < N; body_number++){
            outfile << ' ' << z.nsystem().positions()[body_number].x() << ' ' << z.nsystem().positions()[body_number].y() << ' ' << z.nsystem().positions()[body_number].z() << ' ';
        }
        outfile << '\n';


    std::ofstream outfile_energy("energy_reg/" + filename + ".txt");
    outfile_energy << std::setprecision(8);
    outfile_energy << z.nsystem().get_energy() << '\n';

    Vec u, r;
    std::vector<int> should_be_regularized;

    for (int i = 0; i < iter; i++){ // 5540
        should_be_regularized = z.check_separation(transform_distance);
        //should_be_regularized = (i < 6003) && (i > 5997);
        //should_be_regularized = false;

        if (should_be_regularized[0] && (!regularized)){
            std::cout << "Forward for i = " << i << " and bodies " << should_be_regularized[1] << " and " << should_be_regularized[2] << std::endl;
            z = z.transform_forward(should_be_regularized[1], should_be_regularized[2]);
            regularized = !regularized;
        }
        else if ((!should_be_regularized[0]) && (regularized)){
            std::cout << "Backward for i = " << i << std::endl;
            z = z.transform_backward();
            regularized = !regularized;
        }
        //z = RK4_step(z, h);
        //Yoshida_4_new(x, v, masses, h);
        //PEFRL_friend(z, h);
        //RK4_step(z, h);
        //Yoshida_friend(z, h);
        if ((ADAPTIVE_TIME_STEP) && (!regularized)) {
            NSystem_reg_2D y = z;

            z.timestep(h, dtau);
            y.timestep(h/2, dtau);
            y.timestep(h/2, dtau);

            double error = compare_solutions(z, y);
            //std::cout << "For h = " << h << ", the error is " << error << std::endl;
            if (error > Delta_max) {
                h /= 2;
            } else if (error < Delta_min) {
                h *= 2;
            }
        } else {
            z.timestep(h, dtau);
        }


        if (!regularized){
            outfile << i;
            for (int body_number = 0; body_number < NoB; body_number++){
                outfile << ' ' << z.nsystem().positions()[body_number].x() << ' ' << z.nsystem().positions()[body_number].y() << ' ' << z.nsystem().positions()[body_number].z();
                //std::cout << ' ' << z.nsystem().positions()[body_number].x() << ' ' << z.nsystem().positions()[body_number].y() << ' ' << z.nsystem().positions()[body_number].z() << std::endl;
                //outfile << ' ' << x[body_number].x() << ' ' << x[body_number].y() << ' ' << x[body_number].z();
            }
            outfile << '\n';
            
            outfile_energy << z.nsystem().get_energy() << '\n';
        } 

        if (i < 8001 && i > 7999){
            std::cout << "Started for i = \n" << i;
            for (int body_number = 0; body_number < NoB; body_number++){
                std::cout << ' ' << z.nsystem().positions()[body_number].x() << ' ' << z.nsystem().positions()[body_number].y() << ' ' << z.nsystem().positions()[body_number].z() << std::endl;
                //std::cout << ' ' << z.nsystem().positions()[body_number].x() << ' ' << z.nsystem().positions()[body_number].y() << ' ' << z.nsystem().positions()[body_number].z() << std::endl;
                //outfile << ' ' << x[body_number].x() << ' ' << x[body_number].y() << ' ' << x[body_number].z();
            }
            z.nsystem().print_positions();

        }
    }

    int time_stop = time(NULL);
    std::cout << "Run Done" << std::endl;
    int time_diff = time_stop - time_start;
    std::cout << "Run Time = " << time_diff << " seconds" << std::endl;
    return 1;
}