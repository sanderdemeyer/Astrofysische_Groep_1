#define _USE_MATH_DEFINES
#include <string.h>
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

// represents a 4-D vector
class Vec4 {
    double _x;
    double _y;
    double _z;
    double _a;

    public:
    Vec4() { _x=0; _y=0; _z=0; _a = 0;}
    Vec4(double x, double y, double z, double a) { _x = x; _y = y; _z = z; _a = a;}
    
    double x() const { return _x; }
    double y() const { return _y; }
    double z() const { return _z; }
    double a() const { return _a; }

    double norm() const {return sqrt(_x*_x + _y*_y + _z*_z + _a*_a);}
    double norm2() const {return _x*_x + _y*_y + _z*_z + _a*_a;}
    double norm3() const {return pow(_x*_x + _y*_y + _z*_z + _a*_a, 1.5);}

    friend Vec transform_vector_backward(Vec4 r);
    friend Vec4 transform_vector_forward(Vec r);

    Vec4& operator*=(double s) {
        _x *= s;
        _y *= s;
        _z *= s;
        _a *= s;
        return *this;
    }

    Vec4& operator/=(double s) {
        _x /= s;
        _y /= s;
        _z /= s;
        _a /= s;
        return *this;
    }

    Vec4& operator+=(Vec4 v) {
        _x += v._x;
        _y += v._y;
        _z += v._z;
        _a += v._a;
        return *this;
    }

    Vec4& operator-=(Vec4 v) {
        _x -= v._x;
        _y -= v._y;
        _z -= v._z;
        _a -= v._a;
        return *this;
    }

};

Vec4 operator*(Vec4 a, double s) { return a *= s; }
Vec4 operator/(Vec4 a, double s) { return a /= s; }
Vec4 operator*(double s, Vec4 a) { return a *= s; }
Vec4 operator+(Vec4 a, Vec4 b) { return a += b; }
Vec4 operator-(Vec4 a, Vec4 b) { return a -= b; }

// represents a 2-D vector
class Vec2 {
    double _x;
    double _y;

    public:
    Vec2() {_x=0; _y=0;}
    Vec2(double x, double y) { _x = x; _y = y;}
    
    double x() const { return _x; }
    double y() const { return _y; }

    double norm() const {return sqrt(_x*_x + _y*_y);}
    double norm2() const {return _x*_x + _y*_y;}
    double norm3() const {return pow(_x*_x + _y*_y, 1.5);}

    friend Vec2 transform_vector_backward(Vec2 r);
    friend Vec2 transform_vector_forward(Vec2 r);

    Vec2& operator*=(double s) {
        _x *= s;
        _y *= s;
        return *this;
    }

    Vec2& operator/=(double s) {
        _x /= s;
        _y /= s;
        return *this;
    }

    Vec2& operator+=(Vec2 v) {
        _x += v._x;
        _y += v._y;
        return *this;
    }

    Vec2& operator-=(Vec2 v) {
        _x -= v._x;
        _y -= v._y;
        return *this;
    }

};

Vec2 operator*(Vec2 a, double s) { return a *= s; }
Vec2 operator/(Vec2 a, double s) { return a /= s; }
Vec2 operator*(double s, Vec2 a) { return a *= s; }
Vec2 operator+(Vec2 a, Vec2 b) { return a += b; }
Vec2 operator-(Vec2 a, Vec2 b) { return a -= b; }


Vec2 transform_vector_forward_2D(Vec r){
    // Transforms from unregularized to regularized coordinates
    double expr;
    expr = r.x() + r.norm();
    
    double u1 = sqrt(expr/2);
    double u2 = r.y()/sqrt(2*expr);
    return Vec2(u1, u2);
}

Vec transform_vector_backward_2D(Vec2 r){
    // Transforms from regularized to unregularized coordinates
    double R1 = pow(r.x(), 2) - pow(r.y(), 2);
    double R2 = 2*r.x()*r.y();
    return Vec(R1, R2, 0);
}

class Regularized_coo {
// Defines the class dealing with the regularized coordinates
private:
    Vec4 _u;
    Vec4 _v;
    double _mu;

public: 
    // Constructors
    Regularized_coo(Vec4 u, Vec4 v, double mu){
    _u = u;
    _v = v;
    _mu = mu;
    }
    Regularized_coo(){
        _u = Vec4(0.0, 0.0, 0.0, 0.0);
        _v = Vec4(0.0, 0.0, 0.0, 0.0);
        _mu = 0;
    }

    // Friend declarations necessary for integrators
    friend double Leapfrog_reg(Regularized_coo& y_old, double h);
    friend double RK4_step_reg(Regularized_coo& y_n, double h);

    // Getter functions
    Vec4 u() const { return _u; }
    Vec4 v() const { return _v; }
    double mu() const { return _mu; }

    // Some basic arithmetic
    Regularized_coo& operator*=(double s) {
        _u *= s;
        _v *= s;
        return *this;
    }

    Regularized_coo& operator/=(double s) {
        _u /= s;
        _v /= s;
        return *this;
    }

    Regularized_coo& operator+=(Regularized_coo b) {
        _u += b.u();
        _v += b.v();
        return *this;
    }

    double u_squared(){
        return _u.norm2();
    }

    // Evaluate the driver function
    Regularized_coo evaluate_g_reg(){
        double E = (2*_v.norm2() - _mu)/(_u.norm2());
        Vec4 u_prime = _v;
        Vec4 v_prime = 0.5*E*_u;
        return Regularized_coo(u_prime, v_prime, _mu);
    }

};

Regularized_coo operator*(Regularized_coo a, double s) { return a *= s; }
Regularized_coo operator*(double s, Regularized_coo a) { return a *= s; }
Regularized_coo operator/(Regularized_coo a, double s) { return a /= s; }
Regularized_coo operator+(Regularized_coo a, Regularized_coo b) { return a += b; }

class Regularized_coo_2D {
// Does the same as the above class, but in 2D
private:
    Vec2 _u;
    Vec2 _v;
    double _mu;

public: 
    Regularized_coo_2D(Vec2 u, Vec2 v, double mu){
    _u = u;
    _v = v;
    _mu = mu;
    }
    Regularized_coo_2D(){
        _u = Vec2(0.0, 0.0);
        _v = Vec2(0.0, 0.0);
        _mu = 0;
    }

    friend double Leapfrog_reg_2D(Regularized_coo_2D& y_old, double h);
    friend double RK4_step_reg_2D(Regularized_coo_2D& y_n, double h);
    friend double Yoshida_4_step_reg_2D(Regularized_coo_2D& y_n, double h);

    Vec2 u() const { return _u; }
    Vec2 v() const { return _v; }
    double mu() const { return _mu; }

    Regularized_coo_2D& operator*=(double s) {
        _u *= s;
        _v *= s;
        return *this;
    }

    Regularized_coo_2D& operator/=(double s) {
        _u /= s;
        _v /= s;
        return *this;
    }

    Regularized_coo_2D& operator+=(Regularized_coo_2D b) {
        _u += b.u();
        _v += b.v();
        return *this;
    }

    double u_squared(){
        return _u.norm2();
    }

    Regularized_coo_2D evaluate_g_reg_2D(){
        double E = (2*_v.norm2() - _mu)/(_u.norm2());
        Vec2 u_prime = _v;
        Vec2 v_prime = 0.5*E*_u;
        return Regularized_coo_2D(u_prime, v_prime, _mu);
    }


};

Regularized_coo_2D operator*(Regularized_coo_2D a, double s) { return a *= s; }
Regularized_coo_2D operator*(double s, Regularized_coo_2D a) { return a *= s; }
Regularized_coo_2D operator/(Regularized_coo_2D a, double s) { return a /= s; }
Regularized_coo_2D operator+(Regularized_coo_2D a, Regularized_coo_2D b) { return a += b; }


Vec4 rotate_vec4(Vec4 u, double theta){
    // ***** LEGACY ***** //
    // Rotate a Vec4 in the plane of quaternion k.
    double comp1 = u.x()*cos(theta);
    double comp2 = u.y()*cos(theta);
    double comp3 = u.z()*cos(theta);
    double comp4 = u.a()*cos(theta);

    comp1 += -u.a()*sin(theta);
    comp2 += u.z()*sin(theta);
    comp3 += -u.y()*sin(theta);
    comp4 += u.x()*sin(theta);

    return Vec4(comp1, comp2, comp3, comp4);
}

Vec4 transform_vector_forward(Vec r){
    // ***** LEGACY ***** //
    double expr;
    expr = r.x() + r.norm();
    
    double u1 = -sqrt(expr/2);
    double u2 = r.y()/(2*u1);
    double u3 = r.z()/(2*u1);
    // double u2 = r.y()/sqrt(2*expr);
    // double u3 = r.z()/sqrt(2*expr);
    return Vec4(u1, u2, u3, 0);
}


Regularized_coo regularize_together(Vec u_old, Vec v_old, double Energy, double reduced_mass){
    // ***** LEGACY ***** //
    Vec4 u_new = transform_vector_forward(u_old);
    Vec4 v_new = transform_vector_forward(v_old);

    Vec4 u, v;
    double check = 10000;
    double theta1_best = 0;
    double theta2_best = 0;

    for (int i = 0; i < 1000; i++){
        for (int j = 0; j < 1000; j++){
            double theta1 = M_PI*i/1000;
            double theta2 = M_PI*j/1000;
            u = rotate_vec4(u_new, theta1);
            v = rotate_vec4(v_new, theta2);
            double check1 = v.x()*u.a() - v.a()*u.x() - v.y()*u.z() + v.z()*u.y();
            double check2 = Energy - 2*v.norm2()/u.norm2();
            if (pow(check1,2)+pow(check2,2) < check){
                theta1_best = theta1;
                theta2_best = theta2;
                check = pow(check1,2)+pow(check2,2);
            }
        }
    }
    return Regularized_coo(rotate_vec4(u_new, 0), rotate_vec4(v_new, 0), reduced_mass);
    return Regularized_coo(rotate_vec4(u_new, theta1_best), rotate_vec4(v_new, theta2_best), reduced_mass);
}


Regularized_coo_2D regularize_together_2D(Vec u_old, Vec v_old, double reduced_mass){
    // Transforms u and v to regularized coordinates
    Vec2 u_new = transform_vector_forward_2D(u_old);
    double v1_new = (v_old.x()*u_new.x() + v_old.y()*u_new.y())/(2*u_new.norm2());//*u_new.norm2();
    double v2_new = (v_old.y()*u_new.x() - v_old.x()*u_new.y())/(2*u_new.norm2());//*u_new.norm2();

    return Regularized_coo_2D(u_new, Vec2(v1_new, v2_new), reduced_mass);
}

Regularized_coo_2D regularize_together_2D_backward(Regularized_coo_2D reg_coo){
    // Performs the backward calculation as the function above
    Vec R_new = transform_vector_backward_2D(reg_coo.u());
    double V1 = (2*reg_coo.u().x()*reg_coo.v().x() - 2*reg_coo.u().y()*reg_coo.v().y());///reg_coo.u().norm2();
    double V2 = (2*reg_coo.v().x()*reg_coo.u().y() + 2*reg_coo.u().x()*reg_coo.v().y());///reg_coo.u().norm2();

    return Regularized_coo_2D(Vec2(R_new.x(), R_new.y()), Vec2(V1, V2), reg_coo.mu());
}

// Integrators for the regularized coordinates
double Leapfrog_reg(Regularized_coo& y_old, double dtau){
    double E = (2*y_old._v.norm2() - y_old._mu)/y_old._u.norm2();
    y_old._u = y_old._u + y_old._v*dtau;
    y_old._v = y_old._v + (0.5*E*y_old._u);
    return y_old._u.norm2();
}


double RK4_step_reg(Regularized_coo& y_n, double h){
    Regularized_coo k1 = y_n.evaluate_g_reg() * h;
    Regularized_coo k2 = (y_n + k1*0.5).evaluate_g_reg()*h;
    Regularized_coo k3 = (y_n + k2*0.5).evaluate_g_reg()*h;
    Regularized_coo k4 = (y_n + k3).evaluate_g_reg()*h;
    y_n = y_n + k1/6 + k2/3 + k3/3 + k4/6;
    return y_n._u.norm();
}

double RK4_step_reg_2D(Regularized_coo_2D& y_n, double h){
    Regularized_coo_2D k1 = y_n.evaluate_g_reg_2D() * h;
    Regularized_coo_2D k2 = (y_n + k1*0.5).evaluate_g_reg_2D()*h;
    Regularized_coo_2D k3 = (y_n + k2*0.5).evaluate_g_reg_2D()*h;
    Regularized_coo_2D k4 = (y_n + k3).evaluate_g_reg_2D()*h;
    y_n = y_n + k1/6 + k2/3 + k3/3 + k4/6;
    return y_n._u.norm2();
}

double Yoshida_4_step_reg_2D(Regularized_coo_2D& y_n, double h){
    y_n._u = y_n._u + (YOSHIDA_W1/2)*h*y_n._v;
    y_n._v = y_n._v + YOSHIDA_W1*h*y_n.evaluate_g_reg_2D()._v;
    y_n._u = y_n._u + (YOSHIDA_W0+YOSHIDA_W1)*(h/2)*y_n._v;
    y_n._v = y_n._v + YOSHIDA_W0*h*y_n.evaluate_g_reg_2D()._v;
    y_n._u = y_n._u + (YOSHIDA_W0+YOSHIDA_W1)*(h/2)*y_n._v;
    y_n._v = y_n._v + YOSHIDA_W1*h*y_n.evaluate_g_reg_2D()._v;
    y_n._u = y_n._u + (YOSHIDA_W1/2)*h*y_n._v;
    return y_n._u.norm2();
}


Vec transform_vector_backward(Vec4 u){
    // ***** LEGACY ***** //
    double r1 = pow(u._x, 2) - pow(u._y, 2) - pow(u._z, 2) + pow(u._a, 2);
    double r2 = 2*(u._x*u._y - u._z*u._a);
    double r3 = 2*(u._x*u._z + u._y*u._a);
    return Vec(r1, r2, r3);
}




class NSystem_reg {
// ***** LEGACY ***** //
// Defines a class with an NSystem for the bodies in unregularized coordiantes, together with (possibly) 2 regularized bodies.
private:
    NSystem _nsystem;
    bool _regularized;
    Regularized_coo _reg_coo;
    std::vector<double> _initial_masses;

public:
    NSystem_reg(NSystem nsystem, bool regularized, Regularized_coo reg_coo, std::vector<double> initial_masses)
    : _nsystem(nsystem), _regularized(regularized), _reg_coo(reg_coo), _initial_masses(initial_masses) {}

    NSystem nsystem() const { return _nsystem; }
    bool regularized() const { return _regularized; }
    Regularized_coo reg_coo() const { return _reg_coo; }
    std::vector<double> initial_masses() const { return _initial_masses; }

    //friend void Leapfrog_reg(NSystem_reg& y_n, double dtau);

    bool check_separation(double transform_distance){
        if (_regularized){
            return (_reg_coo.u_squared() < transform_distance);
        } else {
            return (((_nsystem.positions()[0] - _nsystem.positions()[1])).norm() < transform_distance);
        }
    }

    NSystem_reg transform_forward(int body1, int body2){
        if(body1 > body2){
            int temp = body2;
            body2 = body1;
            body1 = temp;
        }
        std::vector<Vec> pos_new = _nsystem.positions();
        std::vector<Vec> vel_new = _nsystem.velocities();
        std::vector<double> masses_new = _nsystem.masses();

        Vec relative_pos = Vec(pos_new[body2] - pos_new[body1]);
        Vec com_pos = Vec((masses_new[body2]*pos_new[body2] + masses_new[body1]*pos_new[body1])/(masses_new[body1]+masses_new[body2]));

        Vec relative_vel = Vec(vel_new[body2] - vel_new[body1]);
        Vec com_vel = Vec((masses_new[body2]*vel_new[body2] + masses_new[body1]*vel_new[body1])/(masses_new[body1]+masses_new[body2]));

        double joint_mass = masses_new[body2] + masses_new[body1];
        double reduced_mass = 1/(1/masses_new[body2] + 1/masses_new[body1]);

        std::vector<double> masses_old;
        masses_old.push_back(masses_new[body1]);
        masses_old.push_back(masses_new[body2]);

        pos_new.erase(pos_new.begin() + body2);
        vel_new.erase(vel_new.begin() + body2);
        masses_new.erase(masses_new.begin() + body2);

        pos_new.erase(pos_new.begin() + body1);
        vel_new.erase(vel_new.begin() + body1);
        masses_new.erase(masses_new.begin() + body1);

        pos_new.push_back(com_pos);
        vel_new.push_back(com_vel);
        masses_new.push_back(joint_mass);

        Regularized_coo reg_coo_rotated = regularize_together(relative_pos, relative_vel, _nsystem.get_energy(), reduced_mass);

        return NSystem_reg(NSystem(pos_new, vel_new, masses_new), true, reg_coo_rotated, masses_old);
    }
    
    
    NSystem_reg transform_backward(){
        std::vector<Vec> pos_new = _nsystem.positions();
        std::vector<Vec> vel_new = _nsystem.velocities();
        std::vector<double> masses_new = _nsystem.masses();

        Vec R_com = pos_new.back();
        Vec V_com = vel_new.back();
        Vec R_relative = transform_vector_backward(_reg_coo.u());
        Vec V_relative = transform_vector_backward(_reg_coo.v());

        double sum_masses = _initial_masses[0] + _initial_masses[1];
        Vec R1 = R_com - _initial_masses[1]*R_relative/sum_masses;
        Vec R2 = R1 + R_relative;
        Vec V1 = V_com - _initial_masses[1]*V_relative/sum_masses;
        Vec V2 = V1 + V_relative;

        pos_new.pop_back();
        vel_new.pop_back();
        masses_new.pop_back();
        pos_new.push_back(R1);
        pos_new.push_back(R2);
        vel_new.push_back(V1);
        vel_new.push_back(V2);
        masses_new.push_back(_initial_masses[0]);
        masses_new.push_back(_initial_masses[1]);

        return NSystem_reg(NSystem(pos_new, vel_new, masses_new), false, _reg_coo, _initial_masses);
    }

    void timestep(double dtau){
        if (_regularized){
            double u_squared = RK4_step_reg(_reg_coo, dtau);
            Yoshida_4_friend(_nsystem, dtau*u_squared);
        } else {
            Yoshida_4_friend(_nsystem, dtau);
        }
    }
};


class NSystem_reg_2D {
// Defines a class with an NSystem for the bodies in unregularized coordiantes, together with (possibly) 2 regularized bodies.
// This is implemented in 2D. 3D does not work.
private:
    NSystem _nsystem;
    bool _regularized;
    Regularized_coo_2D _reg_coo;
    std::vector<double> _initial_masses;

public:
    NSystem_reg_2D(NSystem nsystem, bool regularized, Regularized_coo_2D reg_coo, std::vector<double> initial_masses)
    : _nsystem(nsystem), _regularized(regularized), _reg_coo(reg_coo), _initial_masses(initial_masses) {}

    NSystem nsystem() const { return _nsystem; }
    bool regularized() const { return _regularized; }
    Regularized_coo_2D reg_coo() const { return _reg_coo; }
    std::vector<double> initial_masses() const { return _initial_masses; }

    //friend void Leapfrog_reg(NSystem_reg& y_n, double dtau);

    std::vector<size_t> check_separation(double transform_distance){
        // Check whether the system should be regularized based on the distances between the bodies.
        if (_regularized){
            if (_reg_coo.u_squared() < transform_distance) {
                std::vector<size_t> return_values = {1, 0, 0};
                return return_values;
            } else {
                std::vector<size_t> return_values = {0, 0, 0};
                return return_values;
            }
        } else {
            for (size_t body_1 = 0; body_1 < _nsystem.positions().size(); body_1++) {
                for (size_t body_2 = 0; body_2 < _nsystem.positions().size(); body_2++) {
                    if ((body_1 != body_2) && (((_nsystem.positions()[body_1] - _nsystem.positions()[body_2])).norm() < transform_distance)) {
                        std::vector<size_t> return_values = {1, body_1, body_2};
                        return return_values;
                    }
                }
            }
            std::vector<size_t> return_values = {0, 0, 0};
            return return_values;
        }
    }

    NSystem_reg_2D transform_forward(int body1, int body2){
        // Transforms the system to regularized coordinates for bodies 1 and 2
        if(body1 > body2){
            int temp = body2;
            body2 = body1;
            body1 = temp;
        }
        std::vector<Vec> pos_new = _nsystem.positions();
        std::vector<Vec> vel_new = _nsystem.velocities();
        std::vector<double> masses_new = _nsystem.masses();

        Vec relative_pos = Vec(pos_new[body2] - pos_new[body1]);
        Vec com_pos = Vec((masses_new[body2]*pos_new[body2] + masses_new[body1]*pos_new[body1])/(masses_new[body1]+masses_new[body2]));

        Vec relative_vel = Vec(vel_new[body2] - vel_new[body1]);
        Vec com_vel = Vec((masses_new[body2]*vel_new[body2] + masses_new[body1]*vel_new[body1])/(masses_new[body1]+masses_new[body2]));

        double joint_mass = masses_new[body2] + masses_new[body1];
        double reduced_mass = 1/(1/masses_new[body2] + 1/masses_new[body1]);

        std::vector<double> masses_old;
        masses_old.push_back(masses_new[body1]);
        masses_old.push_back(masses_new[body2]);

        pos_new.erase(pos_new.begin() + body2);
        vel_new.erase(vel_new.begin() + body2);
        masses_new.erase(masses_new.begin() + body2);

        pos_new.erase(pos_new.begin() + body1);
        vel_new.erase(vel_new.begin() + body1);
        masses_new.erase(masses_new.begin() + body1);

        pos_new.push_back(com_pos);
        vel_new.push_back(com_vel);
        masses_new.push_back(joint_mass);

        Regularized_coo_2D reg_coo = regularize_together_2D(relative_pos, relative_vel, reduced_mass);

        return NSystem_reg_2D(NSystem(pos_new, vel_new, masses_new), true, reg_coo, masses_old);
    }
    
    
    NSystem_reg_2D transform_backward(){
        // Performs the backward calculation as the above function.

        std::vector<Vec> pos_new = _nsystem.positions();
        std::vector<Vec> vel_new = _nsystem.velocities();
        std::vector<double> masses_new = _nsystem.masses();

        Vec R_com = pos_new.back();
        Vec V_com = vel_new.back();
        Regularized_coo_2D backwards_transformed = regularize_together_2D_backward(_reg_coo);
        Vec R_relative = Vec(backwards_transformed.u().x(), backwards_transformed.u().y(), 0);
        Vec V_relative = Vec(backwards_transformed.v().x(), backwards_transformed.v().y(), 0);

        double sum_masses = _initial_masses[0] + _initial_masses[1];
        Vec R1 = R_com - _initial_masses[1]*R_relative/sum_masses;
        Vec R2 = R1 + R_relative;
        Vec V1 = V_com - _initial_masses[1]*V_relative/sum_masses;
        Vec V2 = V1 + V_relative;

        pos_new.pop_back();
        vel_new.pop_back();
        masses_new.pop_back();
        pos_new.push_back(R1);
        pos_new.push_back(R2);
        vel_new.push_back(V1);
        vel_new.push_back(V2);
        masses_new.push_back(_initial_masses[0]);
        masses_new.push_back(_initial_masses[1]);

        return NSystem_reg_2D(NSystem(pos_new, vel_new, masses_new), false, _reg_coo, _initial_masses);
    }

    void timestep(double h, double dtau){
        // Performs a single timestep
        if (_regularized){
            double u_squared = RK4_step_reg_2D(_reg_coo, dtau);
            
            Yoshida_4_friend(_nsystem, dtau*u_squared);
        } else {
            Yoshida_4_friend(_nsystem, h);
        }
    }
};


double compare_solutions(NSystem_reg_2D a, NSystem_reg_2D b){
    // Compare two solutions a and b.
    double error = 0;
    for (size_t i = 0; i < a.nsystem().positions().size(); i++){
        error += (a.nsystem().positions()[i] - b.nsystem().positions()[i]).norm2();
        error += (a.nsystem().velocities()[i] - b.nsystem().velocities()[i]).norm2();
    }
    return sqrt(error);
}
