#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "classes.cpp"


int main(){
    std::cout << "Started" << std::endl;

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

    Body b1 = Body(pos1, vel1, m1);
    Body b2 = Body(pos2, vel2, m2);

    b2.Print();
    std::cout << vel2.norm3() << std::endl;
    print(b2.drive_function(b1.pos()));

    System y = System(pos1, pos2, vel1, vel2, m1, m2);

    for (int i = 0; i < 20000; i++){
    y = RK4_step(y, h);
    // std::cout << "Body 1 is " << std::endl;
    // print(y.pos1()); 
    // std::cout << "Body 2 is " << std::endl;
    // print(y.pos2()); 
    // std::cout << y.get_energy() << std::endl;
    outfile << i << ' ' << y.pos1().x() << ' ' << y.pos1().y()<< ' ' << y.pos2().x() << ' ' << y.pos2().y() << '\n';
    outfile_energy << y.get_energy() << '\n';
    }

    std::cout << "Done" << std::endl;
    return 1;
}