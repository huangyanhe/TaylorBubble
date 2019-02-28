#ifndef NEWTON_BUBBLE__
#define NEWTON_BUBBLE__

#include "levmar2.h"
#include "bubble_grid.h"
#include "lapack.h"
class newton_bubble : public levmar<Real, double> {
 public:
    valarray<Real> r0, r1, r2; // ... initial residual?
    valarray<Real> x0; //initial guess
    bubble_grid grid;
    int m1;
    int m2;
    int n1;
    const Real alpha = 50.0; //surface tension
    Real F=0.49; // Froude number
    mp_mat<int> vec_mat, mat_vec;

    // for upper and lowr region fphi
    mp_mat<Real> fphi1;
    mp_mat<Real> fphiphi1;
    mp_mat<Real> fphi2;
    mp_mat<Real> fphiphi2;
    // all these are only for lowr-left boundary
    valarray<Real> ff;
    valarray<Real> f1;
    valarray<Real> f11;
    valarray<Real> f111;
    valarray<Real> f2;
    valarray<Real> f12;
    valarray<Real> f112;
    valarray<Real> q;
    valarray<Real> qphi;
    valarray<Real> Kphi;

    newton_bubble(bubble_grid &_grid);
    void initial();
    void set_index();
    void cut(mp_mat<Real> &A, const mp_mat<Real> &B);
    void cut(valarray<Real> &u, const valarray<Real> &v);
    void fs(mp_mat<Real> &_f);
    void compute_q();
    void compute_qphi();
    void compute_Kphi();
    void update_f();
    void shape();
    valarray<Real> compute_r(mp_mat<Real> &_f);
    virtual void compute_r();
    virtual void compute_J();
};

#endif
