#define CONSTANT_G 1
#define THETA 1.35120719195966
#define XI_PEFRL 0.1786178958448091
#define LAMBDA_PEFRL -0.2123418310626054
#define CHI_PEFRL -0.06626458266981849
#define YOSHIDA_W0 -1.702414384
#define YOSHIDA_W1 1.351207192

// represents a 3-D vector
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

    //friend Vec4 transform_vector_forward(Vec r);
    //friend Vec transform_vector_backward(Vec4 r);

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


// mag weg?
class Body{
    Vec _pos;
    Vec _vel;
    double _mass;

    public:
    Body(Vec init_pos, Vec init_vel, double mass) {_pos = init_pos; _vel = init_vel; _mass = mass;}
    Body(const Body& a){_pos = a.pos(); _vel = a.vel(); _mass = a.mass();}

    Vec pos() const { return _pos; }
    Vec vel() const { return _vel; }
    double mass() const { return _mass; }

    Body& operator+=(Body v) {
        _pos += v.pos();
        _vel += v.vel();
        return *this;
    }

    Vec calculate_a(){
        return _pos * (-_mass/pow(_pos.norm(),3));
    }

    double energy(double h){
        Vec posn = _pos + (h/2)*_vel;
        return 0.5*pow(_vel.norm(),2) - CONSTANT_G*_mass/posn.norm();
    }

    void initialize(double h){
        Vec a_0 = calculate_a();
        print(a_0);
        _pos = _pos + (1/2*h)*_vel + (pow(h,2)/8)*a_0;
        Vec a_12 = calculate_a();
        print(a_12);
        _vel += h*a_12;
    }

    Vec drive_function(Vec a){
        std::cout << (_mass/(a - _pos).norm3()) << std::endl;
        return (_mass/(a - _pos).norm3()) * (a - _pos);
    }

    void timestep(double h) {
        // starts with r_{n-1/2} and v_{n}
        // returns r_{n+1/2} and v_{n+1}
        // Vec a_pre = calculate_a();
        _pos += h*_vel;
        Vec a = calculate_a();
        _vel += h*a;
    }

    void Print(){
        std::cout << "Position: (" << _pos.x() << ", " << _pos.y() << ", " << _pos.z() << "). Velocity: (" << _vel.x() << ", " << _vel.y() << ", " << _vel.z() << ")." << std::endl;
    }
};

// Vec operator+(Body a, Body b) { return a += b; }
// mag weg?
class System{
    Vec _pos1;
    Vec _pos2;
    Vec _vel1;
    Vec _vel2;
    double _m1;
    double _m2;

    public:
    System(Vec pos1, Vec pos2, Vec vel1, Vec vel2, double m1, double m2){
    _pos1 = pos1;
    _pos2 = pos2;
    _vel1 = vel1;
    _vel2 = vel2;
    _m1 = m1;
    _m2 = m2;
    }

    Vec pos1() const { return _pos1; }
    Vec pos2() const { return _pos2; }
    Vec vel1() const { return _vel1; }
    Vec vel2() const { return _vel2; }

    System& operator*=(double s) {
        _pos1 *= s;
        _pos2 *= s;
        _vel1 *= s;
        _vel2 *= s;
        return *this;
    }

    System& operator/=(double s) {
        _pos1 /= s;
        _pos2 /= s;
        _vel1 /= s;
        _vel2 /= s;
        return *this;
    }

    System& operator+=(System b) {
        _pos1 += b.pos1();
        _pos2 += b.pos2();
        _vel1 += b.vel1();
        _vel2 += b.vel2();
        return *this;
    }

    System evaluate_g(){
        Vec g1 = (-_m2*CONSTANT_G/(_pos1-_pos2).norm3()) * (_pos1-_pos2);
        Vec g2 = (-_m1*CONSTANT_G/(_pos2-_pos1).norm3()) * (_pos2-_pos1);
        return System(_vel1, _vel2, g1, g2, _m1, _m2);
    }

    double get_energy(){
        double E_kin = (0.5*_m1)*_vel1.norm2() + (0.5*_m2)*_vel2.norm2();
        double E_pot = -0.5*CONSTANT_G*_m1*_m2/((_pos1-_pos2).norm());
        return E_kin + E_pot;
    }
};

System operator*(System a, double s) { return a *= s; }
System operator/(System a, double s) { return a /= s; }
System operator+(System a, System b) { return a += b; }

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

/*
std::vector<Vec>& operator+=(const std::vector<Vec>& a, const std::vector<Vec>& b){
    for (size_t i=0; i!=a.size(); ++i) {
        a[i] = a[i] + b[i];
    }
    return a;
}
*/

class NSystem {

private:
    std::vector<Vec> _positions;
    std::vector<Vec> _velocities;
    std::vector<double> _masses;

public:
    NSystem(std::vector<Vec> positions, std::vector<Vec> velocities, std::vector<double> masses){
    _positions = positions;
    _velocities = velocities;
    _masses = masses;
    }

    friend void Forest_Ruth_friend(NSystem& y_n, double h);
    friend void PEFRL_friend(NSystem& y_n, double h);
    friend void Velocity_Verlet_friend(NSystem& y_n, double h);
    friend void Position_Verlet_friend(NSystem& y_n, double h);
    friend void Leapfrog_friend(NSystem& y_n, double h);
    friend void Yoshida_4_friend(NSystem& y_n, double h);

    friend void Leapfrog_regularized_friend(NSystem& y_n, double dtau);

    std::vector<Vec> positions() const { return _positions; }
    std::vector<Vec> velocities() const { return _velocities; }
    std::vector<double> masses() const { return _masses; }
    int n() const { return _masses.size(); }

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

    NSystem evaluate_g(){
        std::vector<Vec> gs;
        for (size_t i=0; i!=_positions.size(); ++i) {
            Vec g = Vec(0, 0, 0);
            for (size_t j=0; j!=_positions.size(); ++j) {
                if (i != j) {
                    g += (-_masses[j]*CONSTANT_G/(_positions[i]-_positions[j]).norm3()) * (_positions[i]-_positions[j]);
                }
            }
            gs.push_back(g);
        }
        return NSystem(_velocities, gs, _masses);
    }

    double get_energy(){
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
    double get_energy_regularized(){
        return 0; // TO DO
    }

};

NSystem operator*(NSystem a, double s) { return a *= s; }
NSystem operator*(double s, NSystem a) { return a *= s; }
NSystem operator/(NSystem a, double s) { return a /= s; }
NSystem operator+(NSystem a, NSystem b) { return a += b; }

std::vector<Vec> evaluate_a(std::vector<Vec> positions, std::vector<double> masses){
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

void Leapfrog_regularized_friend(NSystem& y_n, double dtau){
    double E_first = y_n.get_energy_regularized();
    y_n._positions = y_n._positions + dtau*y_n._velocities;
    double E_second = y_n.get_energy_regularized();
    y_n._velocities = y_n._velocities + (dtau/4)*(E_first+E_second)*y_n._positions;
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

NSystem getvalues(std::string inputfile) {
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



class Regularized_coo {

private:
    Vec4 _u;
    Vec4 _v;
    double _mu;

public: 
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

    friend double Leapfrog_reg(Regularized_coo& y_old, double h);
    friend double RK4_step_reg(Regularized_coo& y_n, double h);

    Vec4 u() const { return _u; }
    Vec4 v() const { return _v; }
    double mu() const { return _mu; }

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

    Regularized_coo evaluate_g_reg(){
        double E = (2*_v.norm2() - _mu)/(_u.norm2());
        //std::cout << "energy is " << E << std::endl;
        Vec4 u_prime = _v;
        Vec4 v_prime = 0.5*E*_u;
        return Regularized_coo(u_prime, v_prime, _mu);
    }

};

Regularized_coo operator*(Regularized_coo a, double s) { return a *= s; }
Regularized_coo operator*(double s, Regularized_coo a) { return a *= s; }
Regularized_coo operator/(Regularized_coo a, double s) { return a /= s; }
Regularized_coo operator+(Regularized_coo a, Regularized_coo b) { return a += b; }


double Leapfrog_reg(Regularized_coo& y_old, double dtau){
    double E = (2*y_old._v.norm2() - y_old._mu)/y_old._u.norm2();
    std::cout << "Energy is " << E << std::endl;
    std::cout << "teller is " << (2*y_old._v.norm2() - y_old._mu) << std::endl;
    std::cout << "noemer is " << y_old._u.norm2() << std::endl;
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

Vec4 transform_vector_forward(Vec r){
    double expr;
    if (r.x() > 0) {
        expr = r.x() + r.norm();
    } else {
        expr = r.x() - r.norm();
    }
    double u1 = sqrt(expr/2);
    double u2 = r.y()/(2*u1);
    double u3 = r.z()/(2*u1);
    // double u2 = r.y()/sqrt(2*expr);
    // double u3 = r.z()/sqrt(2*expr);
    std::cout << "Teller of transform is " << expr << std::endl;
    std::cout << "u1 = " << u1 << std::endl;
    std::cout << "u2 = " << u2 << std::endl;
    std::cout << "u3 = " << u3 << std::endl;
    return Vec4(u1, u2, u3, 0);
}


Vec transform_vector_backward(Vec4 u){
    double r1 = pow(u._x, 2) - pow(u._y, 2) - pow(u._z, 2) + pow(u._a, 2);
    double r2 = 2*(u._x*u._y - u._z*u._a);
    double r3 = 2*(u._x*u._z + u._y*u._a);
    return Vec(r1, r2, r3);
}




class NSystem_reg {

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
        // only works for two bodies for the moment
        if (_regularized){
            return (_reg_coo.u_squared() < transform_distance);
        } else {
            //std::cout << "distance is " << ((_nsystem.positions()[0] - _nsystem.positions()[1])).norm() << std::endl;
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

        return NSystem_reg(NSystem(pos_new, vel_new, masses_new), true, Regularized_coo(transform_vector_forward(relative_pos), transform_vector_forward(relative_vel), reduced_mass), masses_old);
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
            double dt = RK4_step_reg(_reg_coo, dtau);
            Yoshida_4_friend(_nsystem, dt);
            //Yoshida_4_friend(_nsystem, dtau);
        } else {
            Yoshida_4_friend(_nsystem, dtau);
        }
    }
};

/*
void Leapfrog_reg(NSystem_reg& y_n, double dtau){
    if (y_n._regularized){
        Leapfrog_friend(y_n._nsystem, dtau);
    } else {
        Leapfrog_friend(y_n._nsystem, dtau);
    }
}
*/

NSystem RK4_step_old(NSystem y_n, double h){
    NSystem k1 = y_n.evaluate_g() * h;
    NSystem k2 = (y_n + k1*0.5).evaluate_g()*h;
    NSystem k3 = (y_n + k2*0.5).evaluate_g()*h;
    NSystem k4 = (y_n + k3).evaluate_g()*h;
    return y_n + k1/6 + k2/3 + k3/3 + k4/6;
    return y_n + k1/6 + k2/3 + k3/3 + k4/6;
}

NSystem Forest_Ruth(NSystem y_n, double h){
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

NSystem PEFRL(NSystem y_n, double h){
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

