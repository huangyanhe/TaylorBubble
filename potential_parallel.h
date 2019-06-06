#ifndef POTENTIAL__
#define POTENTIAL__

#include "levmar2.h"
#include "real_def.h"
#include "fftjw2.h"
#include "lapack.h"
#include "str_double.h"

class potential : public levmar<Real, double> {
 public:
  int N, M; // number of grid points in S1, S2
  int n1; // number of quadrature points
  valarray<Real> theta1, theta2; // stores theta_n, theta_m
  //  valarray<Real> r0;//, r1, r2; // initial bubble shape guess r0, first and second derivatives r1, r2
  //  valarray<cmplx<Real> > c0;
  mp_mat<cmplx<Real> > c1, c2;
  Real V0; // bubble volume
  const Real U; //rising velocity
  const Real alpha; // surface tension
  const Real F; //Froude number
  fft_class<Real> fft1, fft2; // fft on S1, S2
  valarray<Real> mu1, mu2;
  // Real a, b, c, d, e, f, g, h, G; // stores G[a,b,c,d]
  valarray<Real> x1, x2, w1, w2; //classic&generalized guassian quadrature
  mp_mat<Real> xx1_1, xx1_2, xx1_3, ww1_1, ww1_2, ww1_3, xx2_1, xx2_2, xx2_3, ww2_1, ww2_2, ww2_3; // stores all the quadrature nodes&weights
  //Real X; // for computing B, D when m>0.9
  valarray<Real> m0; //stores the midpoint values of 10 division in calculating B, D
  valarray<int> index_B, index_D;
  mp_mat<Real> B, D; // stores coeffs of B(m), D(m) Taylor expansion
  //  mp_mat<Real> A; // stores the lhs matrix when computing mu
  //valarray<Real> rhs; // stores the rhs when computing mu
  //  valarray<Real> _r; // stores temporary residue.
  static const Real pi;
  
  potential(int N_, int M_);

  void set_theta();
  void set_c12();
  void set_x0();
  void compute_c0r12(valarray<Real> &_r0, valarray<cmplx<Real> > &_c0, valarray<Real> &_r1, valarray<Real> &_r2);
  void compute_V0();
  void compute_xw();
  void set_BD();
  void initial();
  
  Real compute_B(Real _m, Real _multi=0.0);
  Real compute_D(Real _m, Real _multi=0.0);
  Real compute_G(Real _a, Real _b, Real _c, Real _d);
  void compute_mu(valarray<Real> &r0, valarray<Real> &mu1, valarray<Real> &mu2);
  //void update_shape();
  //void update_r(vector<Real> &_r, valarray<Real> _rr);
  valarray<Real> compute_rr(valarray<Real> &_r0);
  
  virtual void compute_r();
  virtual void compute_J();
  void interpolate(int l);

};

#endif
  
