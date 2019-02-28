#ifndef BUBBLE_GRID__
#define BUBBLE_GRID__

#include "fft_template.h"
#include "mp_mat.h"
#include "real_def.h"
class bubble_grid {
 public:
  int m1, m2;
  int n;
  int N1, N2, N; // '' lower domain ''
  Real phi1;
  Real phi2; //phi domain is [-phi1, phi2]
  valarray<Real> x1, xx1, x2, xx2; // two grids for phi, xx:real grid, x: reference
  valarray<Real> y, yy; // grid for psi
  mp_mat<Real> D1_1, D1_2; // phi direction order 1&2&3 differentiation matrix
  mp_mat<Real> D2_1, D2_2, D2_3;
  mp_mat<Real> D3_1, D3_2; // psi direction...
  Real map1, map2, map3; 
  mp_mat<Real> f; // f evaluation
  mp_mat<Real> tmp11, tmp12, tmp13, tmp21, tmp22, tmp23;
  const Real pi = 4.0*atan(1.0);


  bubble_grid(int _m1, int _m2, int _n); //Chebyshev-Labatto grid
  void set_D();
  void set_xy();
  void set_f();
  //int delta(int _i, int _j);
};

#endif
