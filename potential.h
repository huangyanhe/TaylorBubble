#ifndef POTENTIAL__
#define POTENTIAL__

#include "levmar2.h"
//#include "poly_grid.h"
#include "fftjw2.h"
#include "lapack.h"

class potential : public levmar<Real, double> {
 public:
  int N, M; // number of grid points in S1, S2
  int n; // number of quadrature points
  valarray<Real> theta1, theta2; // stores theta_n, theta_m
  valarray<Real> r0, r1, r2; // initial bubble shape guess r0, first and second derivatives r1, r2
  valarray<cx> c0;
  mp_mat<cx> c1, c2;
  Real V0; // bubble volume
  const Real alpha = 50; // surface tension
  Real F=0.49; //Froude number
  fft_class<Real> fft1, fft2; // fft on S1, S2
  valarray<Real> mu1, mu2;
  Real a, b, c, d, e, f, g, h, G; // stores G[a,b,c,d]
  valarray x1, x2, w1, w2, xx, ww; //classic&generalized guassian quadrature
  Real m;
  Real X; // for computing B, D when m>0.9
  valarray<int> m0, index_B, index_D;
  mp_mat<Real> B, D; // stores coeffs of B(m), D(m) Taylor expansion
  mp_mat<Real> A; // stores the lhs matrix when computing mu
  valarray<Real> rhs; // stores the rhs when computing mu
  potential(int N_, int M_);

  void set_theta();
  void set_c12();
  void set_r0();
  
  void compute_r12();
  void compute_V0();
  void compute_xw();
  void set_BD();
  Real compute_B(Real _m);
  Real compute_D(Real _m);
  void compute_G(Real _a, Real _b, Real _c, Real _d);
  void compute_mu();
  void update_shape();
  virtual void compute_r();
  virtual void compute_J();

};

#endif
  
