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