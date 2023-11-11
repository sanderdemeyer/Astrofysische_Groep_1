#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <vector>

#include "classes.cpp"


int main(){
    
    std::cout << "Run Start" << std::endl;
    int time_start = time(NULL);

    std::ofstream outfile("two_body_motion_RK2.txt");
    outfile << std::setprecision(8);
    std::ofstream outfile_energy("two_body_motion_energy_RK2.txt");
    outfile_energy << std::setprecision(8);

    double h = 0.01;
    double m1 = 0.005;
    double m2 = 0.005;

    Vec pos1 = Vec(0.5, 0, 0);
    Vec pos2 = Vec(1, 0, 0);
    Vec vel1 = Vec(0.01, -0.06, 0);
    Vec vel2 = Vec(-0.01, 0.06, 0);

    //System y = System(pos1, pos2, vel1, vel2, m1, m2);

    std::vector<Vec> positions = {pos1, pos2};
    std::vector<Vec> velocities = {vel1, vel2};
    std::vector<double> masses = {m1, m2};
    NSystem z = NSystem(positions, velocities, masses);

    /*
    for (int i = 0; i < 20000; i++){

        y = RK4_step(y, h);
        outfile << i << ' ' << y.pos1().x() << ' ' << y.pos1().y()<< ' ' << y.pos2().x() << ' ' << y.pos2().y() << '\n';
        outfile_energy << y.get_energy() << '\n';

    }
    */
    
    
    for (int i = 0; i < 20000; i++){

        z = RK4_step(z, h);
        outfile << i << ' ' << z.positions()[0].x() << ' ' << z.positions()[0].y()<< ' ' << z.positions()[1].x() << ' ' << z.positions()[1].y() << '\n';
        outfile_energy << z.get_energy() << '\n';

    }

    int time_stop = time(NULL);
    std::cout << "Run Done" << std::endl;
    int time_diff = time_stop - time_start;
    std::cout << "Run Time = " << time_diff << " seconds" << std::endl;
    return 1;
}