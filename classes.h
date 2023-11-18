#include <vector>

class Vec {
    double _x;
    double _y;
    double _z;

    public:
    Vec();
    Vec(double x, double y, double z);
    double x() const;
    double y() const;
    double z() const;

    double norm() const;
    double norm2() const;
    double norm3() const;
    
};

class System{
    Vec _pos1;
    Vec _pos2;
    Vec _vel1;
    Vec _vel2;
    double _m1;
    double _m2;

    public:
    System(Vec pos1, Vec pos2, Vec vel1, Vec vel2, double m1, double m2);

    Vec pos1() const ;
    Vec pos2() const;
    Vec vel1() const;
    Vec vel2() const;

    System evaluate_g();

    System RK2_step(double h);
};

class NSystem {

private:
    std::vector<Vec> _positions;
    std::vector<Vec> _velocities;
    std::vector<double> _masses;

public:
    NSystem(std::vector<Vec> positions, std::vector<Vec> velocities, std::vector<double> masses);

    friend void Forest_Ruth_friend(NSystem& y_n, double h);
    friend void PEFRL_friend(NSystem& y_n, double h);
    friend void Velocity_Verlet_friend(NSystem& y_n, double h);
    friend void Position_Verlet_friend(NSystem& y_n, double h);
    friend void Leapfrog_friend(NSystem& y_n, double h);
    friend void Yoshida_4_friend(NSystem& y_n, double h);
    friend void Leapfrog_regularized_friend(NSystem& y_n, double dtau);

    std::vector<Vec> positions() const;
    std::vector<Vec> velocities() const;
    std::vector<double> masses() const;

    NSystem& operator*=(double s); 
    NSystem& operator/=(double s);
    NSystem& operator+=(NSystem b);
    NSystem evaluate_g();
    bool check_separation(double transform_distance);
    double get_energy();
    double get_energy_regularized();

};
