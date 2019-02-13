#ifndef NEWTON_BUBBLE__
#define NEWTON_BUBBLE__

#include "levmar2.h"
#include "bubble_grid.h"
#include "lapack.h"
class newton_bubble : public levmar<Real, double> {
 public:
    valarray<Real> r0, r1, r2; // ... initial residual?
    //valarray<Real> rtmp_M, rtmp_N; // ... temporary residual
    valarray<Real> x0; //initial guess
    bubble_grid grid;
    int m1;
    int n1;
    int a;
    const Real alpha = 6.0; //surface tension
    Real F=2.0; // Froude number
    mp_mat<int> mat_vec;
    mp_mat<int> vec_mat;

    
    valarray<Real> ff, f_cut;
    valarray<Real> f1, f1_cut;
    valarray<Real> f11, f11_cut;
    valarray<Real> f111, f111_cut;
    valarray<Real> f2, f2_cut;
    valarray<Real> f12, f12_cut;
    valarray<Real> f112, f112_cut;
    valarray<Real> q;
    valarray<Real> qphi;
    valarray<Real> Kphi;

    newton_bubble(bubble_grid &_grid);
    void initial();
    void set_index();
    void cut(mp_mat<Real> &A, const mp_mat<Real> &B);
    void cut(valarray<Real> &u, const valarray<Real> &v);
    void fs(mp_mat<Real> &_f);
    void fcuts();
    void compute_q();
    void compute_qphi();
    void compute_Kphi();
    valarray<Real> compute_r(mp_mat<Real> &_f);
    virtual void compute_r();
    virtual void compute_J();
};

#endif
