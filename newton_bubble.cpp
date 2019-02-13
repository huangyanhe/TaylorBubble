#include "newton_bubble.h"
newton_bubble::newton_bubble(bubble_grid &_grid) :
  levmar(_grid.N, _grid.N),
  r0(n),
  r1(n),
  r2(n),
  x0(n),
  grid(_grid),
  m1(_grid.m),
  n1(_grid.n),
  a(_grid.A),
  mat_vec(m1+1, n1+1, 0.0),
  vec_mat(n, 2, 0.0),
  ff(m1+1),
  f_cut(m1+1-a),
  f1(m1+1),
  f1_cut(m1+1-a),
  f11(m1+1),
  f11_cut(m1+1-a),
  f111(m1+1),
  f111_cut(m1+1-a),
  f2(m1+1),
  f2_cut(m1+1-a),
  f12(m1+1),
  f12_cut(m1+1-a),
  f112(m1+1),
  f112_cut(m1+1-a),
  q(m1+1-a),
  qphi(m1+1-a),
  Kphi(m1+1-a)
{
  // add more things here
  set_index();
  
}

// need to think about what f to choose. first update bubble_grid::set_f, then tranform to here.
void newton_bubble::initial() {
  for(int i=0; i<n; i++){
    
  }
}

void newton_bubble::set_index() {
  int k = 0;
  for(int j=1; j<n1; j++)
    for(int i=1; i<=m1; i++) {
      vec_mat(k,0) = i;
      vec_mat(k,1) = j;
      mat_vec(i,j) = k;
      k++;
    }
  for(int i=a; i<=m1; i++) {
    vec_mat(k,0) = i;
    vec_mat(k,1) = 0;
    mat_vec(i,0) = k;
    k++;
  }
  if (k!= n) throw gen_err("set_index dimension wrong");
}

void newton_bubble::cut(mp_mat<Real> &A, const mp_mat<Real> &B) {
  if(A.m != B.m-1 || A.n != B.n-2)
    throw gen_err("matrix dimension doesn't match in matrix cut");
  for(int i=1; i<B.m; i++)
    for(int j=1; j<B.n-1; j++)
      A(i-1, j-1) = B(i,j);
}

void newton_bubble::cut(valarray<Real> &u, const valarray<Real> &v) {
  if(v.size() != m1+1 || u.size() != m1+1-a)
     throw gen_err("vector dimension doesn't match in vector cut");
  for(int i=a; i<=m1; i++)
    u[i-a] = v[i];
}
   

// need to test
void newton_bubble::fs(mp_mat<Real> &_f) {
  mp_mat<Real> tmp_mat(m1+1, n1+1), tmp1_mat(m1+1, n1+1);
  ff = _f.extract_column(0);
  dgemv('N',m1+1,m1+1,1.0,grid.D1_1.p,m1+1,&ff[0],1,0.0,&f1[0],1);
  dgemv('N',m1+1,m1+1,1.0,grid.D1_2.p,m1+1,&ff[0],1,0.0,&f11[0],1);
  dgemv('N',m1+1,m1+1,1.0,grid.D1_3.p,m1+1,&ff[0],1,0.0,&f111[0],1);
  dgemm('N','N',m1+1,n1+1,n1+1,1.0,grid.f.p,m1+1,grid.D2_1.p,n1+1,0.0,tmp_mat.p,m1+1);
  f2 = tmp_mat.extract_column(0);
  dgemm('N','N',m1+1,n1+1,m1+1,1.0,grid.D1_1.p,m1+1,tmp_mat.p,m1+1,0.0,tmp1_mat.p,m1+1);
  f12 = tmp1_mat.extract_column(0);
  dgemm('N','N',m1+1,n1+1,m1+1,1.0,grid.D1_2.p,m1+1,tmp_mat.p,m1+1,0.0,tmp1_mat.p,m1+1);
  f112 = tmp1_mat.extract_column(0);
}

void newton_bubble::fcuts() {
  cut(f_cut, ff);
  cut(f1_cut, f1);
  cut(f11_cut, f11);
  cut(f111_cut, f111);
  cut(f2_cut, f2);
  cut(f12_cut, f12);
  cut(f112_cut, f112);
}

// need to test on real data
void newton_bubble::compute_q() {
  valarray<Real> temp0(m1+1-a, 0.0);
  temp0 = f1_cut*f1_cut/f_cut+f2_cut*f2_cut;
  q = 2.0/sqrt(temp0);
}

void newton_bubble::compute_qphi() {
  valarray<Real> temp0(m1+1-a, 0.0);
  qphi = -q/2.0;
  qphi = qphi*qphi*qphi;

  temp0 = f1_cut*(2.0*f11_cut*f_cut-f1_cut*f1_cut);
  temp0 /= f_cut*f_cut;
  temp0 += 2.0*f2_cut*f12_cut;
  qphi += temp0;			
}

void newton_bubble::compute_Kphi() {
  valarray<Real> temp0(m1+1-a, 0.0), temp1(m1+1-a, 0.0);
  Kphi = -(f12_cut*q+f2_cut*qphi)*sqrt(f_cut)+f2_cut*q/2.0/sqrt(f_cut);
  Kphi /= 2*f_cut;
  
  temp0 = 3.0*q*q*qphi*sqrt(f_cut);
  temp0 -= q*q*q/2.0/sqrt(f_cut);
  temp0 /= 4.0*f_cut;
  temp0 *= f2_cut*f11_cut-f1_cut*f12_cut-0.5*f2_cut*f1_cut*f1_cut/f_cut;
  Kphi = Kphi + temp0;
  
  temp0 = 0.25*q*q*q/sqrt(f_cut);
  temp1 = f12_cut*f11_cut + f2_cut*f111_cut-f11_cut*f12_cut-f1_cut*f112_cut;
  temp1 -= (f12_cut*f1_cut*f1_cut+2*f2_cut*f1_cut*f11_cut)*f_cut - f2_cut*f1_cut*f1_cut*f1_cut/2.0/f_cut/f_cut;
  Kphi = Kphi + temp0 * temp1;
}

valarray<Real> newton_bubble::compute_r(mp_mat<Real> &_f) {
  if (_f.m != m1+1 || _f.n != n1+1)
    throw gen_err("grid dimension mismatch");
  cut(grid.tmp1, _f*grid.D2_2);
  cut(grid.tmp2, grid.D1_2*_f);
  cut(grid.tmp3, _f);
  grid.tmp1 += grid.tmp2/grid.tmp3;
  cut(grid.tmp2, grid.D1_1*_f);
  grid.tmp2 = grid.tmp2/grid.tmp3;
  grid.tmp1 -= grid.tmp2<grid.tmp2;
  int m2 = m1*(n1-1);
  valarray<Real> _r(n);
  for(int i=0; i<m2; i++) {
    _r[i] = grid.tmp1(vec_mat(i,0)-1, vec_mat(i,1)-1);
  }
  
  fs(_f);
  fcuts();
  compute_q();
  compute_qphi();
  compute_Kphi();
  valarray<Real> temp0(m1+1-a);
  temp0 = q*qphi-1.0/(2.0*F*F)*f2_cut+1.0/alpha*Kphi;
  for(int i=m2; i<n; i++)
    _r[i] = temp0[i-m2];
  return _r;
}

void newton_bubble::compute_r() {
  r0 =  compute_r(grid.f);
  for(int i=0; i<n; i++)
    r[i] = r0[i];
}
void newton_bubble::compute_J() {
  mp_mat<Real> f1(m1, n1), f2(m1, n1);
  const Real eps = 0.00001; // check
  int i, j;
  for(int k=0; k<n; k++) {
    i = vec_mat(k, 0);
    j = vec_mat(k, 1);
    f1 = grid.f;
    f1(i, j) += eps;
    f2 = grid.f;
    f2(i, j) -= eps;
    r1 = compute_r(f1);
    r2 = compute_r(f2);
    r1 -= r2;
    for(int l=0; l<n; l++) 
      J(l, k) = r1[l]/(2*eps);
  }
  J.dump("J", 17, '%');
}


