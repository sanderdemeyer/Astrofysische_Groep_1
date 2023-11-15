#define CONSTANT_G 1 //6.6743e-11
#define THETA 1.35120719195966

// represents a 3-D vector
class Vec {
    double _x;
    double _y;
    double _z;

    public:
    Vec() { _x=0; _y=0; _z=0;}
    Vec(double x, double y, double z)
        { _x = x; _y = y; _z = z;}
    double x() const { return _x; }
    double y() const { return _y; }
    double z() const { return _z; }

    double norm() const {return sqrt(_x*_x + _y*_y + _z*_z);}
    double norm2() const {return _x*_x + _y*_y + _z*_z;}
    double norm3() const {return pow(_x*_x + _y*_y + _z*_z, 1.5);}

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

    std::vector<Vec> positions() const { return _positions; }
    std::vector<Vec> velocities() const { return _velocities; }
    std::vector<double> masses() const { return _masses; }

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

};

NSystem operator*(NSystem a, double s) { return a *= s; }
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

System RK2_step(System y_n, double h){
    System k1 = y_n.evaluate_g() * h;
    System k2 = (y_n + k1*0.5).evaluate_g()*h;
    return y_n + k2;
}

System RK4_step(System y_n, double h){
    System k1 = y_n.evaluate_g() * h;
    System k2 = (y_n + k1*0.5).evaluate_g()*h;
    System k3 = (y_n + k2*0.5).evaluate_g()*h;
    System k4 = (y_n + k3).evaluate_g()*h;
    return y_n + k1/6 + k2/3 + k3/3 + k4/6;
}


NSystem RK4_step(NSystem y_n, double h){
    NSystem k1 = y_n.evaluate_g() * h;
    NSystem k2 = (y_n + k1*0.5).evaluate_g()*h;
    NSystem k3 = (y_n + k2*0.5).evaluate_g()*h;
    NSystem k4 = (y_n + k3).evaluate_g()*h;
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

// Vectors vs arrays

// + en * operator overloaden op niveau van vectoren?