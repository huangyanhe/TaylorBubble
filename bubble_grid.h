#ifndef BUBBLE_GRID__
#define BUBBLE_GRID__

#include "fft_template.h"
#include "mp_mat.h"
#include "real_def.h"
class bubble_grid {
 public:
  int m; // might need to split m if we use domain decomposition
  int n;
  int A; //starting from here pde boundary condition
  int M; // number of all grid points (=(m+1)(n+1))
  int N; // number of all unkown grid values(solving N equations)
  Real phi1;
  Real phi2; //phi domain is [-phi1, phi2]
  
  valarray<Real> x, x1; // grid for phi
  valarray<Real> y, y1; // grid for psi
  mp_mat<Real> D1_1, D1_2, D1_3; // phi direction order 1&2 differentiation matrix
  mp_mat<Real> D2_1, D2_2; // psi direction...
  Real map1, map2; // change of variable of map to default domain [-1,1]*[-1,1],  4*psi-1, 2*phi/(phi2+phi1)-(phi2-phi1)/(phi2+phi1)
  mp_mat<Real> f; // update the known values at boundary
  mp_mat<Real> tmp, tmp1, tmp2, tmp3;
  const Real pi = 4.0*atan(1.0);

  bubble_grid(int _m, int _n);  // Gauss-Labatto mesh
  bubble_grid(int _m, int _n, int flag); //just to get small grid to create levmar.
  void set_D();
  void set_xy();
  void set_A();
  void set_N();
  void set_f();
  int delta(int _i, int _j);

};

#endif
