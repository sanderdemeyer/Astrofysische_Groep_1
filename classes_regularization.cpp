class NSystem_reg {

private:
    NSystem _nsystem;
    bool _regularized;
    Vec _u;

public:
    NSystem_reg(NSystem nsystem, bool regularized, Vec u){
    _nsystem = nsystem;
    _regularized = regularized;
    _u = u;
    }
};
