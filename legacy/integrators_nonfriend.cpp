#define CONSTANT_G 1 //6.6743e-11
#define THETA 1.35120719195966
#define XI_PEFRL 0.1786178958448091
#define LAMBDA_PEFRL -0.2123418310626054
#define CHI_PEFRL -0.06626458266981849
#define YOSHIDA_W0 -1.702414384
#define YOSHIDA_W1 1.351207192

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream> 
#include <ctime>
#include <vector>
#include "classes.cpp"

NSystem RK4_step_old(NSystem y_n, double h){
    NSystem k1 = y_n.evaluate_g() * h;
    NSystem k2 = (y_n + k1*0.5).evaluate_g()*h;
    NSystem k3 = (y_n + k2*0.5).evaluate_g()*h;
    NSystem k4 = (y_n + k3).evaluate_g()*h;
    return y_n + k1/6 + k2/3 + k3/3 + k4/6;
    return y_n + k1/6 + k2/3 + k3/3 + k4/6;
}

NSystem Forest_Ruth_old(NSystem y_n, double h){
    std::vector<Vec> x = y_n.positions();
    std::vector<Vec> v = y_n.velocities();
    std::vector<double> masses = y_n.masses();

    x = x + (THETA*h/2)*v;
    v = v + THETA*h*evaluate_a(x, masses);
    x = x + ((1-THETA)*h/2)*v;
    v = v + ((1-2*THETA)*h)*evaluate_a(x, masses);
    x = x + ((1-THETA)*h/2)*v;
    v = v + THETA*h*evaluate_a(x, masses);
    x = x + (THETA*h/2)*v;

    return NSystem(x, v, masses);
}

NSystem PEFRL_old(NSystem y_n, double h){
    std::vector<Vec> x = y_n.positions();
    std::vector<Vec> v = y_n.velocities();
    std::vector<double> masses = y_n.masses();
    x = x + XI_PEFRL*h*v;
    v = v + ((1-2*LAMBDA_PEFRL)*h/2)*evaluate_a(x, masses);
    x = x + CHI_PEFRL*h*v;
    v = v + LAMBDA_PEFRL*h*evaluate_a(x, masses);
    x = x + (1-2*(CHI_PEFRL+XI_PEFRL))*h*v;
    v = v + LAMBDA_PEFRL*h*evaluate_a(x, masses);
    x = x + CHI_PEFRL*h*v;
    v = v + ((1-2*LAMBDA_PEFRL)*h/2)*evaluate_a(x, masses);
    x = x + XI_PEFRL*h*v;
    return NSystem(x, v, masses);
}

NSystem Velocity_Verlet(NSystem y_n, double h){
    std::vector<Vec> x = y_n.positions();
    std::vector<Vec> v = y_n.velocities();
    std::vector<double> masses = y_n.masses();
    v = v + (h/2)*evaluate_a(x, masses);
    x = x + h*v;
    v = v + (h/2)*evaluate_a(x, masses);
    return NSystem(x, v, masses);
}

NSystem Position_Verlet(NSystem y_n, double h){
    std::vector<Vec> x = y_n.positions();
    std::vector<Vec> v = y_n.velocities();
    std::vector<double> masses = y_n.masses();
    x = x + (h/2)*v;
    v = v + h*evaluate_a(x, masses);
    x = x + (h/2)*v;
    return NSystem(x, v, masses);
}

NSystem Leapfrog(NSystem y_n, double h){
    std::vector<Vec> x = y_n.positions();
    std::vector<Vec> v = y_n.velocities();
    std::vector<double> masses = y_n.masses();
    x = x + h*v;
    v = v + h*evaluate_a(x, masses);
    return NSystem(x, v, masses);
}

NSystem Yoshida_4(NSystem& y_n, double h){
    std::vector<Vec> x = y_n.positions();
    std::vector<Vec> v = y_n.velocities();
    std::vector<double> masses = y_n.masses();

    x = x + (YOSHIDA_W1/2)*h*v;
    v = v + YOSHIDA_W1*h*evaluate_a(x, masses);
    x = x + (YOSHIDA_W0+YOSHIDA_W1)*(h/2)*v;
    v = v +  YOSHIDA_W0*h*evaluate_a(x, masses);
    x = x + (YOSHIDA_W0+YOSHIDA_W1)*(h/2)*v;
    v = v + YOSHIDA_W1*h*evaluate_a(x, masses);
    x = x + (YOSHIDA_W1/2)*h*v;
    return NSystem(x, v, masses);
}

void Yoshida_4_new(std::vector<Vec>& x, std::vector<Vec>& v, std::vector<double> masses, double h){
    x = x + (YOSHIDA_W1/2)*h*v;
    v = v + YOSHIDA_W1*h*evaluate_a(x, masses);
    x = x + (YOSHIDA_W0+YOSHIDA_W1)*(h/2)*v;
    v = v +  YOSHIDA_W0*h*evaluate_a(x, masses);
    x = x + (YOSHIDA_W0+YOSHIDA_W1)*(h/2)*v;
    v = v + YOSHIDA_W1*h*evaluate_a(x, masses);
    x = x + (YOSHIDA_W1/2)*h*v;
}
