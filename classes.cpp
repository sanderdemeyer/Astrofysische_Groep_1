#define CONSTANT_G 39.473107
#define THETA 1.35120719195966
#define XI_PEFRL 0.1786178958448091
#define LAMBDA_PEFRL -0.2123418310626054
#define CHI_PEFRL -0.06626458266981849
#define YOSHIDA_W0 -1.702414384
#define YOSHIDA_W1 1.351207192
#define EPSILON 0.001

#include <string.h>
#include <assert.h>
#include <vector>
#include <algorithm>

// represents a 3-D vector
// This implements some simple Arithmetic like the norm, the angle, scalar multiplication and division, addition, and subtraction.
class Vec {
    double _x;
    double _y;
    double _z;

    public:
    Vec() { _x=0; _y=0; _z=0;}
    Vec(double x, double y, double z) { _x = x; _y = y; _z = z;}
    
    double x() const { return _x; }
    double y() const { return _y; }
    double z() const { return _z; }

    double norm() const {return sqrt(_x*_x + _y*_y + _z*_z);}
    double norm2() const {return _x*_x + _y*_y + _z*_z;}
    double norm3() const {return pow(_x*_x + _y*_y + _z*_z, 1.5);}

    double angle() const {return acos((_x)/norm());}
    double angle(const Vec other) const {return acos((_x*other.x() + _y*other.y() + _z*other.z())/(norm()*other.norm()));}

    Vec& operator*=(double s) {
        _x *= s;
        _y *= s;
        _z *= s;
        return *this;
    }

    Vec& operator/=(double s) {
        _x /= s;
        _y /= s;
        _z /= s;
        return *this;
    }

    Vec& operator+=(Vec v) {
        _x += v._x;
        _y += v._y;
        _z += v._z;
        return *this;
    }

    Vec& operator-=(Vec v) {
        _x -= v._x;
        _y -= v._y;
        _z -= v._z;
        return *this;
    }
};

Vec operator*(Vec a, double s) { return a *= s; }
Vec operator/(Vec a, double s) { return a /= s; }
Vec operator*(double s, Vec a) { return a *= s; }
Vec operator+(Vec a, Vec b) { return a += b; }
Vec operator-(Vec a, Vec b) { return a -= b; }


void print(Vec a){ 
    std::cout << "(" << a.x() << ", " << a.y() << ", " << a.z() << ")" << std::endl; 
}


std::vector<Vec> operator*(double s, std::vector<Vec> a){
    for (size_t i=0; i!=a.size(); ++i) {
        a[i] *= s;
    }
    return a;
}

std::vector<Vec> operator+(std::vector<Vec> a, std::vector<Vec> b){
    for (size_t i=0; i!=a.size(); ++i) {
        a[i] += b[i];
    }
    return a;
}

class NSystem {
// This is the class which represents the current positions and velocities of all the bodies present in the system, together with their masses.
private:
    std::vector<Vec> _positions; // Vector containing all the positions
    std::vector<Vec> _velocities; // Vector containing all the velocites
    std::vector<double> _masses; // Vector containing all the masses

public:
    NSystem(std::vector<Vec> positions, std::vector<Vec> velocities, std::vector<double> masses){
    _positions = positions;
    _velocities = velocities;
    _masses = masses;
    }

    NSystem() {
        std::vector<Vec> positions;
        std::vector<Vec> velocities;
        std::vector<double> masses;
        _positions = positions;
        _velocities = velocities;
        _masses = masses;
    }

    // Make the different integrators friend of this class, such that it can access the private properties.
    friend void Forward_Euler(NSystem& y_n, double h);
    friend void RK2_step(NSystem& y_n, double h);
    friend void Heun(NSystem& y_n, double h);
    friend void Heun3(NSystem& y_n, double h);
    friend void Ralston(NSystem& y_n, double h);
    friend void Ralston3(NSystem& y_n, double h);
    friend void RK3_step(NSystem& y_n, double h);
    friend void RK4_step(NSystem& y_n, double h);
    friend void Forest_Ruth_friend(NSystem& y_n, double h);
    friend void PEFRL_friend(NSystem& y_n, double h);
    friend void Velocity_Verlet_friend(NSystem& y_n, double h);
    friend void Position_Verlet_friend(NSystem& y_n, double h);
    friend void Leapfrog_friend(NSystem& y_n, double h);
    friend void Yoshida_4_friend(NSystem& y_n, double h);

    // Simple getter functions
    std::vector<Vec> positions() const { return _positions; }
    std::vector<Vec> velocities() const { return _velocities; }
    std::vector<double> masses() const { return _masses; }

    int n() const { return _masses.size(); } // Get the size of the system, namely the number of bodies

    // Underlying functions define some basic arithmetic.
    NSystem& operator*=(double s) {
        for (size_t i=0; i!=_positions.size(); ++i) {
            _positions[i] *= s;
            _velocities[i] *= s;
            }
        return *this;
    }

    NSystem& operator/=(double s) {
        for (size_t i=0; i!=_positions.size(); ++i) {
            _positions[i] /= s;
            _velocities[i] /= s;
            }
        return *this;
    }

    NSystem& operator+=(NSystem b) {
        for (size_t i=0; i!=_positions.size(); ++i) {
            _positions[i] += b.positions()[i];
            _velocities[i] += b.velocities()[i];
            }
        return *this;
    }

    NSystem& operator-=(NSystem b) {
        for (size_t i=0; i!=_positions.size(); ++i) {
            _positions[i] -= b.positions()[i];
            _velocities[i] -= b.velocities()[i];
            }
        return *this;
    }

    // This function evaluates the driver function. 
    NSystem evaluate_g(){
        std::vector<Vec> gs;
        for (size_t i=0; i!=_positions.size(); ++i) {
            Vec g = Vec(0, 0, 0);
            for (size_t j=0; j!=_positions.size(); ++j) {
                if (i != j) {
                    //g += (-_masses[j]*CONSTANT_G/(_positions[i]-_positions[j]).norm3()) * (_positions[i]-_positions[j]);
                    g += (-_masses[j]*CONSTANT_G/pow(pow(EPSILON,2)+(_positions[i]-_positions[j]).norm2(),1.5)) * (_positions[i]-_positions[j]);
                }
            }
            gs.push_back(g);
        }
        return NSystem(_velocities, gs, _masses);
    }

    double get_energy(){ 
        // This returns the total energy of the system.
        double E_kin = 0.0;
        double E_pot = 0.0;
        for (size_t i=0; i!=_positions.size(); ++i) {
            E_kin += (0.5*_masses[i])*_velocities[i].norm2();
            for (size_t j=0; j!=_positions.size(); ++j) {
                if (i != j) {
                E_pot += -0.5*CONSTANT_G*_masses[i]*_masses[j]/((_positions[i]-_positions[j]).norm());
                }
            }
        }
        return E_kin + E_pot;
    }

    void print_positions(){
        // Print all the positions of the different bodies.
        std::cout << "Positions are ";
        for (size_t i = 0; i < _positions.size(); i++){
            print(_positions[i]); 
        }
    }

};

// The underlying lines define some extra basic arithmetic.
NSystem operator*(NSystem a, double s) { return a *= s; }
NSystem operator*(double s, NSystem a) { return a *= s; }
NSystem operator/(NSystem a, double s) { return a /= s; }
NSystem operator+(NSystem a, NSystem b) { return a += b; }
NSystem operator-(NSystem a, NSystem b) { return a -= b; }

double compare_solutions(NSystem a, NSystem b){
    // Compute the 'distance' between two solutions. This is used for the adaptive timestep methods.
    double error = 0;
    for (size_t i = 0; i < a.positions().size(); i++){
        error += (a.positions()[i] - b.positions()[i]).norm2();
        error += (a.velocities()[i] - b.velocities()[i]).norm2();
    }
    return sqrt(error);
}

std::vector<Vec> evaluate_a(std::vector<Vec> positions, std::vector<double> masses){
    // This does the same as the function evaluate_g, but considers only the driver functions and not r' = v.
    // This returns a vector of Vec objects containing the driver functions, instead of a new object of the class NSystem.
    std::vector<Vec> gs;
    for (size_t i=0; i!=positions.size(); ++i) {
        Vec g = Vec(0, 0, 0);
        for (size_t j=0; j!=positions.size(); ++j) {
            if (i != j) {
                g += (-CONSTANT_G*masses[j]/(positions[i]-positions[j]).norm3()) * (positions[i]-positions[j]);
            }
        }
        gs.push_back(g);
    }
    return gs;
}

// The underlying functions are the integrator functions.
// These functions take an object of NSystem and a timestep h, and updates the NSystem over a time h.
// The return type is void, since the object is passed by reference and updated there.

void Forward_Euler(NSystem& y_n, double h){
    NSystem k1 = y_n.evaluate_g() * h;
    y_n = y_n + k1;
}

void RK2_step(NSystem& y_n, double h){
    NSystem k1 = y_n.evaluate_g() * h;
    NSystem k2 = (y_n + k1*0.5).evaluate_g()*h;
    y_n = y_n + k2;
}

void Heun(NSystem& y_n, double h){
    NSystem k1 = y_n.evaluate_g() * h;
    NSystem k2 = (y_n + k1*0.5).evaluate_g()*h;
    y_n = y_n + k1/2 + k2/2;
}

void Heun3(NSystem& y_n, double h){
    NSystem k1 = y_n.evaluate_g() * h;
    NSystem k2 = (y_n + k1*0.5).evaluate_g()*h;
    NSystem k3 = (y_n + k2*0.5).evaluate_g()*h;
    y_n = y_n + k1/4 + 3*k3/4;
}

void Ralston(NSystem& y_n, double h){
    NSystem k1 = y_n.evaluate_g() * h;
    NSystem k2 = (y_n + k1*0.5).evaluate_g()*h;
    y_n = y_n + k1/4 + 3*k2/4;
}

void Ralston3(NSystem& y_n, double h){
    NSystem k1 = y_n.evaluate_g() * h;
    NSystem k2 = (y_n + k1*0.5).evaluate_g()*h;
    NSystem k3 = (y_n + k2*0.5).evaluate_g()*h;
    y_n = y_n + 2*k1/9 + k2/3 + 4*k3/9;
}

void RK3_step(NSystem& y_n, double h){
    NSystem k1 = y_n.evaluate_g() * h;
    NSystem k2 = (y_n + k1*0.5).evaluate_g()*h;
    NSystem k3 = (y_n + k2*0.5).evaluate_g()*h;
    y_n = y_n + k1/6 + 2*k2/3 + k3/6;
}

void RK4_step(NSystem& y_n, double h){
    NSystem k1 = y_n.evaluate_g() * h;
    NSystem k2 = (y_n + k1*0.5).evaluate_g()*h;
    NSystem k3 = (y_n + k2*0.5).evaluate_g()*h;
    NSystem k4 = (y_n + k3).evaluate_g()*h;
    y_n = y_n + k1/6 + k2/3 + k3/3 + k4/6;
}


void Forest_Ruth_friend(NSystem& y_n, double h){
    y_n._positions = y_n._positions + (THETA*h/2)*y_n._velocities;
    y_n._velocities = y_n._velocities + THETA*h*evaluate_a(y_n._positions, y_n._masses);
    y_n._positions = y_n._positions + ((1-THETA)*h/2)*y_n._velocities;
    y_n._velocities = y_n._velocities + ((1-2*THETA)*h)*evaluate_a(y_n._positions, y_n._masses);
    y_n._positions = y_n._positions + ((1-THETA)*h/2)*y_n._velocities;
    y_n._velocities = y_n._velocities + THETA*h*evaluate_a(y_n._positions, y_n._masses);
    y_n._positions = y_n._positions + (THETA*h/2)*y_n._velocities;
}

void PEFRL_friend(NSystem& y_n, double h){
    y_n._positions = y_n._positions + XI_PEFRL*h*y_n._velocities;
    y_n._velocities = y_n._velocities + ((1-2*LAMBDA_PEFRL)*h/2)*evaluate_a(y_n._positions, y_n._masses);
    y_n._positions = y_n._positions + CHI_PEFRL*h*y_n._velocities;
    y_n._velocities = y_n._velocities + LAMBDA_PEFRL*h*evaluate_a(y_n._positions, y_n._masses);
    y_n._positions = y_n._positions + (1-2*(CHI_PEFRL+XI_PEFRL))*h*y_n._velocities;
    y_n._velocities = y_n._velocities + LAMBDA_PEFRL*h*evaluate_a(y_n._positions, y_n._masses);
    y_n._positions = y_n._positions + CHI_PEFRL*h*y_n._velocities;
    y_n._velocities = y_n._velocities + ((1-2*LAMBDA_PEFRL)*h/2)*evaluate_a(y_n._positions, y_n._masses);
    y_n._positions = y_n._positions + XI_PEFRL*h*y_n._velocities;
}

void Velocity_Verlet_friend(NSystem& y_n, double h){
    y_n._velocities = y_n._velocities + (h/2)*evaluate_a(y_n._positions, y_n._masses);
    y_n._positions = y_n._positions + h*y_n._velocities;
    y_n._velocities = y_n._velocities + (h/2)*evaluate_a(y_n._positions, y_n._masses);
}

void Position_Verlet_friend(NSystem& y_n, double h){
    y_n._positions = y_n._positions + (h/2)*y_n._velocities;
    y_n._velocities = y_n._velocities + h*evaluate_a(y_n._positions, y_n._masses);
    y_n._positions = y_n._positions + (h/2)*y_n._velocities;
}

void Leapfrog_friend(NSystem& y_n, double h){
    y_n._positions = y_n._positions + h*y_n._velocities;
    y_n._velocities = y_n._velocities + h*evaluate_a(y_n._positions, y_n._masses);
}

void Yoshida_4_friend(NSystem& y_n, double h){
    y_n._positions = y_n._positions + (YOSHIDA_W1/2)*h*y_n._velocities;
    y_n._velocities = y_n._velocities + YOSHIDA_W1*h*evaluate_a(y_n._positions, y_n._masses);
    y_n._positions = y_n._positions + (YOSHIDA_W0+YOSHIDA_W1)*(h/2)*y_n._velocities;
    y_n._velocities = y_n._velocities + YOSHIDA_W0*h*evaluate_a(y_n._positions, y_n._masses);
    y_n._positions = y_n._positions + (YOSHIDA_W0+YOSHIDA_W1)*(h/2)*y_n._velocities;
    y_n._velocities = y_n._velocities + YOSHIDA_W1*h*evaluate_a(y_n._positions, y_n._masses);
    y_n._positions = y_n._positions + (YOSHIDA_W1/2)*h*y_n._velocities;
}

void Yoshida_4_new(std::vector<Vec>& x, std::vector<Vec>& v, std::vector<double> masses, double h){
    // This is equivalent to the function Yoshida_4_friend. By default, Yoshida_4_friend is used.
    x = x + (YOSHIDA_W1/2)*h*v;
    v = v + YOSHIDA_W1*h*evaluate_a(x, masses);
    x = x + (YOSHIDA_W0+YOSHIDA_W1)*(h/2)*v;
    v = v +  YOSHIDA_W0*h*evaluate_a(x, masses);
    x = x + (YOSHIDA_W0+YOSHIDA_W1)*(h/2)*v;
    v = v + YOSHIDA_W1*h*evaluate_a(x, masses);
    x = x + (YOSHIDA_W1/2)*h*v;
}
 
void RK45_step(NSystem& y_n, double& h, double tolerance){
    // This is the function that defines the RK45 method, also called the Runge–Kutta–Fehlberg method
    // This function evaluates the error by comparing the solutions of a 4th order accurate and a 5th order accurate approximation.
    // This is more efficient than the other way that adaptive timestep control is implemented in this code, 
    // which updates the system twice and independently, and thus has more overhead than this method, which reuses some driver function evaluations.

    NSystem k1 = y_n.evaluate_g() * h;
    NSystem k2 = (y_n + 0.25*k1).evaluate_g()*h;
    NSystem k3 = (y_n + (3.0/32)*k1 + (9.0/32)*k2).evaluate_g()*h;
    NSystem k4 = (y_n + (1932.0/2197)*k1 + (-7200.0/2197)*k2 + (7296.0/2197)*k3).evaluate_g()*h;
    NSystem k5 = (y_n + (439.0/216)*k1 + (-8.0)*k2 + (3680.0/513)*k3 + (-845.0/4104)*k4).evaluate_g()*h;
    NSystem k6 = (y_n + (-8.0/27)*k1 + (2.0)*k2 + (-3544.0/2565)*k3 + (1859.0/4104)*k4 + (-11.0/40)*k5).evaluate_g()*h;
    
    NSystem z_new = y_n + (25.0/216)*k1 + (1408.0/2565)*k3 + (2197.0/4101)*k4 + (-1.0/5)*k5;
    y_n = y_n + (16.0/135)*k1 + (6656.0/12825)*k3 + (28561.0/56430)*k4 + (-9.0/50)*k5 + (2.0/55)*k6;
    std::vector<Vec> dif = (z_new - y_n).positions();
    double difference = 0;
    for (size_t i = 0; i < dif.size(); i++) {
        difference += dif[i].norm();
    }
    h *= pow(tolerance/(2*difference), 0.25);    
}



NSystem getvalues(std::string inputfile) {
    // This initializes an instance of the class NSystem based on a text file containing the initial conditions.
    std::ifstream MyreadFile(inputfile);

    std::vector<Vec> ini_positions;
    std::vector<Vec> ini_velocities;
    std::vector<double> masses;

    for (std::string line; getline(MyreadFile, line);) {
        std::istringstream iss(line);

        double mass;
        iss >> mass;
        masses.push_back(mass);

        double x, y, z;
        iss >> x >> y >> z;
        Vec r_i(x, y, z);
        ini_positions.push_back(r_i);

        double v_x, v_y, v_z;
        iss >> v_x >> v_y >> v_z;
        Vec v_i(v_x, v_y, v_z);
        ini_velocities.push_back(v_i);
    }

    MyreadFile.close();
    return NSystem(ini_positions, ini_velocities, masses);
}


class General_integrator{
    // This class defines a general integrator, based only on its Butcher Tableau.
    // Instead of defining the integrator to be a function, the integrator is an object of the class General Integrator,
    // for which the operator() is overloaded.
    // This is very useful for the higher-order integrators, for which implementing different functions is both cumbersome and error-prone.
    private: 
        int length; // Length of the Butcher Tableau
        std::vector<double> a_table; // Bulk of the Butcher tableau, which claculates the different values of k1, k2, etc...
        std::vector<double> b_table; // Last line of the Butcher tableau, which calculates the updated value of the system based on the values of k1, k2, etc...
    public:
        General_integrator(){} // Default constructor

        General_integrator(int _length, std::vector<double> _a_table, std::vector<double> _b_table) {
            // Constructor using a Butcher Tableau
            length = _length;
            a_table = _a_table;
            b_table = _b_table;
        }

        General_integrator(std::string integr) {
            // Constructor using a string, which is connected with a Butcher Tableau.
            if (integr.compare("RK4") == 0) {
                a_table = {0.0,0.0,0.0,0.0 , 0.5,0.0,0.0,0.0 , 0.0,0.5,0.0,0.0 , 0.0,0.0,1.0,0.0};
                b_table = {1.0/6, 1.0/3, 1.0/3, 1.0/6};
                length = 4;
            } else if (integr.compare("RK6") == 0) {
                a_table = {0.0,0.0,0.0,0.0,0.0,0.0,0.0 , 1.0/3,0.0,0.0,0.0,0.0,0.0,0.0 , 0.0,2.0/3,0.0,0.0,0.0,0.0,0.0,
                            1.0/12,1.0/3,-1.0/12,0.0,0.0,0.0,0.0 , -1.0/16,9.0/8,-3.0/16,-3.0/8,0.0,0.0,0.0 , 0.0,9.0/8,-3.0/8,-3.0/4,1.0/2,0.0,0.0,
                            9.0/44,-9.0/11,63.0/44,18.0/11,0.0,-16.0/11,0.0};
                b_table = {11.0/120, 0.0, 27.0/40, 27.0/40, -4.0/15, -4.0/15, 11.0/120};
                length = 7;
            } else if (integr.compare("Wray3") == 0) {
                a_table = {0.0,0.0,0.0 , 8.0/15,0.0,0.0 , 1.0/4,5.0/12,0.0};
                b_table = {1.0/4, 0.0, 3.0/4};
                length = 3;
            } else if (integr.compare("SSPRK3") == 0) {
                a_table = {0.0,0.0,0.0 , 1.0,0.0,0.0 , 1.0/4,1.0/4,0.0};
                b_table = {1.0/6, 1.0/6, 2.0/3};
                length = 3;
            } else if (integr.compare("3_over_8") == 0) {
                a_table = {0.0,0.0,0.0,0.0 , 1.0/3,0.0,0.0,0.0 , -1.0/3,1.0,0.0,0.0 , 1.0,-1.0,1.0,0.0};
                b_table = {1.0/8, 3.0/8, 3.0/8, 1.0/8};
                length = 4;
            } else if (integr.compare("Ralston4") == 0) {
                a_table = {0.0,0.0,0.0,0.0 , 0.4,0.0,0.0,0.0 , 0.29697761,0.15875964,0.0,0.0 , 0.21810040,-3.05096516,3.83286476,0.0};
                b_table = {0.17476028,-0.55148066,0.120553560,0.17118478};
                length = 4;
            } else if (integr.compare("RK8") == 0) {
                a_table = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 , 4.0/27,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 , 1.0/18,3.0/18,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                            1.0/12,0.0,3.0/12,0.0,0.0,0.0,0.0,0.0,0.0,0.0 , 1.0/8,0.0,0.0,3.0/8,0.0,0.0,0.0,0.0,0.0,0.0 , 13.0/54,0.0,-27.0/54,42.0/54,8.0/54,0.0,0.0,0.0,0.0,0.0,
                            389.0/4320,0.0,-54.0/4320,966.0/4320,-824.0/4320,243.0/4320,0.0,0.0,0.0,0.0 , -234.0/20,0.0,81.0/20,-1164.0/20,656.0/20,-122.0/20,800.0/20,0.0,0.0,0.0 , -127.0/288,0.0,18.0/288,-678.0/288,456.0/288,-9.0/288,576.0/288,4.0/288,0.0,0.0,
                            1481.0/820,0.0,-81.0/820,7104.0/820,-3376.0/820,72.0/820,-5040.0/820,-60.0/820,720.0/820,0.0};
                b_table = {41.0/840,0.0,0.0,27.0/840,272.0/840,27.0/840,216.0/840,0.0,216.0/840,41.0/840};
                length = 10;
            } else if (integr.compare("RK5_wrong") == 0) {
                a_table = {0.0,0.0,0.0,0.0,0.0,0.0,0.0 , 1.0/5,0.0,0.0,0.0,0.0,0.0,0.0 , 3.0/40,9.0/40,0.0,0.0,0.0,0.0,0.0,
                            44.0/45,-56.0/15,32.0/9,0.0,0.0,0.0,0.0 , 19372.0/6561,-25360.0/2187,64448.0/6561,-212.0/729,0.0,0.0,0.0 , -9017.0/3168,-355.0/33,46732.0/5247,49.0/176,-5103.0/18656,0.0,0.0,
                            35.0/384,0.0,500.0/1113,125.0/192,-2187.0/6784,11.0/84,0.0};
                b_table = {5179.0/57600, 0.0, 7571.0/16695, 393.0/640, -92097.0/339200, 187.0/2100, 1.0/40};
                length = 7;
            } else if (integr.compare("RK5") == 0) {
                a_table = {0.0,0.0,0.0,0.0,0.0,0.0 , 1.0/4,0.0,0.0,0.0,0.0,0.0 , 1.0/8,1.0/8,0.0,0.0,0.0,0.0,
                            0.0,-0.5,1.0,0.0,0.0,0.0 , 3.0/16,0.0,0.0,9.0/16,0.0,0.0 , -3.0/7,2.0/7,12.0/7,-12.0/7,8.0/7,0.0};
                b_table = {7.0/90,0.0,32.0/90,12.0/90,32.0/90,7.0/90};
                length = 6;
            } else if (integr.compare("IRK5_wrong") == 0) {
                // Implementation of IRK5 using https://www.ajbasweb.com/old/ajbas/2012/March/97-105.pdf
                // This implementation does not give correct results. Using this intergrator is highly discouraged.
                a_table = {0.0,0.0,0.0,0.0,0.0 , 0.25,0.0,0.0,0.0,0.0,
                            -0.7272,0.7322,0.0,0.0,0.0 , 0.5734,-2.2485,3.344,0.0,0.0,
                            0.1750,0.0121,0.0559,0.5517,0.0};
                b_table = {1.0222, -0.0961, 0.0295, -0.1, 0.6444};
                length = 5; // wrong: 
            } else if (integr.compare("IRK5_wrong2") == 0) {
                // Implementation of IRK5 using https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/99611d53-5ab6-45bc-8743-1e1080025f83/36a3fccd-1914-485e-bb6d-1d30fcb7e387/images/screenshot.png
                // This implementation does not give correct results. Using this intergrator is highly discouraged.
                a_table = {0.0,0.0,0.0,0.0,0.0 , 1.0/4,0.0,0.0,0.0,0.0,
                            -1.0/125,259.0/1000,0.0,0.0,0.0 , 0.386,-0.531,0.644,0.0,0.0,
                            0.206,-0.9,0.892,0.552,0.0};
                b_table = {46.0/45, 1.0/25, -0.107, -0.1, 29.0/45};
                length = 5;
            } else {
                assert((0 == 1) && ("This is not implemented with a Butcher Tableau"));
            }
        }

    void operator()(NSystem& y_n, double h) {
        // This overloads the operator()-function, such that the same syntax as the function-based integrators can be used.
        NSystem kis [length]; // This is an array of the values of k1, k2, etc...
        for (int i = 0; i < length; i++) {
            NSystem argument = y_n;
            for (int j = 0; j < i; j++) {
                argument += h*a_table[i*length+j]*kis[j];
            }
            kis[i] = argument.evaluate_g();
        }
        for (int i = 0; i < length; i++){
            y_n += h*b_table[i]*kis[i]; // The updated value the system is calculated based on the values of k1, k2, etc...
        }
    }
};


int get_driver_evaluations(std::string integrator) {
    // For both function-based integrators as general integrators, the number of driver evaulations per timestep is returned.
    std::vector<std::string> driver_evaluations_1 = {"Forward Euler", "Position Verlet", "Leapfrog"};
    std::vector<std::string> driver_evaluations_2 = {"RK2", "Heun", "Ralston", "Velocity Verlet"};
    std::vector<std::string> driver_evaluations_3 = {"Heun3", "Ralston3", "RK3", "Forest Ruth", "Yoshida_4", "Wray3", "SSPRK3"};
    std::vector<std::string> driver_evaluations_4 = {"RK4", "PEFRL", "3_over_8", "Ralston4"};
    std::vector<std::string> driver_evaluations_5 = {"IRK5"};
    std::vector<std::string> driver_evaluations_6 = {"RK5"};
    std::vector<std::string> driver_evaluations_7 = {"RK6", "RK5_wrong"};
    std::vector<std::string> driver_evaluations_10 = {"RK8"};

    if (std::find(driver_evaluations_1.begin(), driver_evaluations_1.end(), integrator) != driver_evaluations_1.end()) {
        return 1;
    } else if (std::find(driver_evaluations_2.begin(), driver_evaluations_2.end(), integrator) != driver_evaluations_2.end()) {
        return 2;
    } else if (std::find(driver_evaluations_3.begin(), driver_evaluations_3.end(), integrator) != driver_evaluations_3.end()) {
        return 3;
    } else if (std::find(driver_evaluations_4.begin(), driver_evaluations_4.end(), integrator) != driver_evaluations_4.end()) {
        return 4;
    } else if (std::find(driver_evaluations_5.begin(), driver_evaluations_5.end(), integrator) != driver_evaluations_5.end()) {
        return 5;
    } else if (std::find(driver_evaluations_6.begin(), driver_evaluations_6.end(), integrator) != driver_evaluations_6.end()) {
        return 6;
    } else if (std::find(driver_evaluations_7.begin(), driver_evaluations_7.end(), integrator) != driver_evaluations_7.end()) {
        return 7;
    } else if (std::find(driver_evaluations_10.begin(), driver_evaluations_10.end(), integrator) != driver_evaluations_10.end()) {
        return 7;
    } else {
        return -1;
    }
}

