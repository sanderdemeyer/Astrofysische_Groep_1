#define CONSTANT_G 1 //6.6743e-11

// #include "classes.h"

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
        return *this;
    }

    Vec& operator/=(double s) {
        _x /= s;
        _y /= s;
        return *this;
    }

    Vec& operator+=(Vec v) {
        _x += v._x;
        _y += v._y;
        return *this;
    }

    Vec& operator-=(Vec v) {
        _x -= v._x;
        _y -= v._y;
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
        return 0.5*pow(_vel.norm(),2) - _mass/posn.norm();
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

class System_3{
    Vec _pos[3];
    Vec _vel[3];
    double _mass[3];

    public:
    System_3(Vec pos[], Vec vel[], double mass[]){
        for (int i = 0; i < 3; i++){
        _pos[i] = pos[i];
        _vel[i] = vel[i];
        _mass[i] = mass[i];
        }
    }

    Vec[] pos() const { return _pos; }
    Vec[] vel() const { return _vel; }

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