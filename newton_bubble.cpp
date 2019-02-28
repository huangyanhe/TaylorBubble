#include "newton_bubble.h"
#include "str_double.h"
newton_bubble::newton_bubble(bubble_grid &_grid) :
  levmar(_grid.N, _grid.N),
  r0(n),
  r1(n),
  r2(n),
  grid(_grid),
  m1(_grid.m1),
  m2(_grid.m2),
  n1(_grid.n),
  vec_mat(n, 2, 0),
  mat_vec(m1+m2+1, n1, -1),
  fphi1(m1+1, n1+1, 0.0),
  fphiphi1(m1+1, n1+1, 0.0),
  fphi2(m2+1, n1+1, 0.0),
  fphiphi2(m2+1, n1+1, 0.0),
  ff(m2),
  f1(m2),
  f11(m2),
  f111(m2),
  f2(m2),
  f12(m2),
  f112(m2),
  q(m2),
  qphi(m2),
  Kphi(m2)
{
  // add more things here
  set_index();
  initial();
}

// need to think about what f to choose. first update bubble_grid::set_f, then tranform to here.

void newton_bubble::set_index() {
  int k=0;
  //upper inner
  for(int j=1; j<n1; j++) 
      for(int i=1; i<m1; i++) {
	vec_mat(k, 0) = i;
	vec_mat(k, 1) = j;
	mat_vec(i, j) = k;
	k++;
      }
  if (k!= grid.N1)
    throw gen_err("index set wrong, the upper domain");
  // lower inner(including the last row)
  for(int j=1; j<n1; j++)
    for(int i=1; i<=m2; i++) {
      vec_mat(k, 0) = m1+i;
      vec_mat(k, 1) = j;
      mat_vec(m1+i, j) = k;
      k++;
    }
  // interface
  for(int j=1; j<n1; j++) {
    vec_mat(k, 0) = m1;
    vec_mat(k, 1) = j;
    mat_vec(m1, j) = k;
    k++;
  }
  // left 
  for(int i=1; i<=m2; i++) {
    vec_mat(k, 0) = i+m1;
    vec_mat(k, 1) = 0;
    mat_vec(i+m1, 0) = k;
    k++;
  }
  if (k != n)
    throw gen_err("index set wrong, total");
}

void newton_bubble::initial() {
  for(int k=0; k<n; k++)
    x[k] = grid.f(vec_mat(k, 0), vec_mat(k, 1));
}

void newton_bubble::cut(mp_mat<Real> &A, const mp_mat<Real> &B) {
  int m3=A.m, n3=A.n;
  if(m3 == B.m-2 && n3 == B.n-2) {
    for(int i=0; i<m3; i++)
      for(int j=0; j<n3; j++)
	A(i, j) = B(i+1,j+1);
  }
  else if (m3 == B.m-1 && n3 == B.n-2) {
    for(int i=0; i<m3; i++)
      for(int j=0; j<n3; j++)
	A(i, j) = B(i+1, j+1);
  }
  else
    throw gen_err("matrix dimension doesn't satisfy 1 or 2");
}

void newton_bubble::cut(valarray<Real> &u, const valarray<Real> &v) {
  if (u.size() != v.size()-1)
    throw gen_err("cut vector size wrong");
  int m3= u.size();
  for(int i=0; i<m3; i++)
    u[i] = v[i+1];
}

//for lower left boundary
void newton_bubble::fs(mp_mat<Real> &_f) {
  mp_mat<Real> tmp_mat(m2+1, n1+1);
  valarray<Real> haha(m2+1), hehe(m2+1);
  haha = _f.extract_column(0);
  cut(ff, haha);
  dgemv('N',m2+1,m2+1,1.0,grid.D2_1.p,m2+1,&haha[0],1,0.0,&hehe[0],1);
  cut(f1, hehe);
  dgemv('N',m2+1,m2+1,1.0,grid.D2_2.p,m2+1,&haha[0],1,0.0,&hehe[0],1);
  cut(f11, hehe);
  dgemv('N',m2+1,m2+1,1.0,grid.D2_3.p,m2+1,&haha[0],1,0.0,&hehe[0],1);
  cut(f111, hehe);
  dgemm('N','N',m2+1,n1+1,n1+1,1.0,_f.p,m2+1,grid.D3_1.p,n1+1,0.0,tmp_mat.p,m2+1);
  hehe = tmp_mat.extract_column(0);
  cut(f2, hehe);
  //dgemm('N','N',m2+1,n1+1,m2+1,1.0,grid.D2_1.p,m2+1,tmp_mat.p,m2+1,0.0,tmp1_mat.p,m2+1);
  //f12 = tmp1_mat.extract_column(0);
  dgemv('N',m2+1,m2+1,1.0,grid.D2_1.p,m2+1,&hehe[0],1,0.0,&haha[0],1);
  cut(f12, haha);
  //dgemm('N','N',m2+1,n1+1,m2+1,1.0,grid.D2_2.p,m2+1,tmp_mat.p,m2+1,0.0,tmp1_mat.p,m2+1);
  //f112 = tmp1_mat.extract_column(0);
  dgemv('N',m2+1,m2+1,1.0,grid.D2_2.p,m2+1,&hehe[0],1,0.0,&haha[0],1);
  cut(f112, haha);
}


void newton_bubble::compute_q() {
  q = f1*f1/ff+f2*f2;
  q = 2.0/sqrt(q);
}


void newton_bubble::compute_qphi() {
  valarray<Real> temp0(m2, 0.0);
  qphi = -pow(q/2,3);
  temp0 = f1*(2.0*f11*ff-f1*f1);
  temp0 /= ff*ff;
  temp0 += 2.0*f2*f12;
  qphi *= temp0;
}

void newton_bubble::compute_Kphi() {
  valarray<Real> temp0(m2, 0.0), temp1(m2, 0.0);
  Kphi = -(f12*q+f2*qphi)*sqrt(ff)+f2*q*f1/2.0/sqrt(ff);
  Kphi /= 2*ff;
  
  temp0 = 3.0*qphi*sqrt(ff);
  temp0 -= q*f1/2.0/sqrt(ff);
  temp0 *= q*q;
  temp0 /= 4.0*ff;
  temp0 *= f2*f11-f1*f12-0.5*f2*f1*f1/ff;
  Kphi += temp0;
  
  temp0 = 0.25*pow(q,3)/sqrt(ff);
  temp1 = f2*f111-f1*f112;
  temp1 -= ((f12*f1*f1+2*f2*f1*f11)*ff - f2*pow(f1,3))/2.0/ff/ff;
  Kphi += temp0 * temp1;
}

valarray<Real> newton_bubble::compute_r(mp_mat<Real> &_f) {
  mp_mat<Real> _f1(m1+1, n1+1, 0.0), _f2(m2+1, n1+1, 0.0);
  for(int i=0; i<=m1; i++)
    for(int j=0; j<=n1; j++)
      _f1(i, j) = _f(i, j);
  for(int i=0; i<=m2; i++)
    for(int j=0; j<=n1; j++)
      _f2(i, j) = _f(i+m1, j);
  //upper inner variable
  fphi1 = grid.D1_1*_f1;
  fphiphi1 = grid.D1_2*_f1;
  cut(grid.tmp11, _f1*grid.D3_2);
  cut(grid.tmp12, fphiphi1);
  cut(grid.tmp13, _f1);
  grid.tmp11 += grid.tmp12/grid.tmp13;
  cut(grid.tmp12, fphi1);
  grid.tmp12 = grid.tmp12/grid.tmp13;
  grid.tmp11 -= grid.tmp12<grid.tmp12;
  valarray<Real> _r(n);
  
  int k=0;
  for(int j=0; j<n1-1; j++)
    for(int i=0; i<m1-1; i++){
      _r[k] = grid.tmp11(i, j);
      k++;
    }
  if (k!= grid.N1)
    throw gen_err("compute_r something wrong");

  // lower inner variable
  fphi2 = grid.D2_1*_f2;
  fphiphi2 = grid.D2_2*_f2;
  cut(grid.tmp21, _f2*grid.D3_2);
  cut(grid.tmp22, fphiphi2);
  cut(grid.tmp23, _f2);
  grid.tmp21 += grid.tmp22/grid.tmp23;
  cut(grid.tmp22, fphi2);
  grid.tmp22 = grid.tmp22/grid.tmp23;
  grid.tmp21 -= grid.tmp22<grid.tmp22;
 
  for(int j=0; j<n1-1; j++)
    for(int i=0; i<m2; i++){
      _r[k] = grid.tmp21(i, j);
      k++;
    }
  // interface variable
  for(int j=1; j<n1; j++){
    _r[k] = fphi1(m1, j) - fphi2(0, j);
    k++;
  }

  //lower left
  fs(_f2);
  compute_q();
  compute_qphi();
  compute_Kphi();
  valarray<Real> temp0(m1);
  temp0 = q*qphi-f2/(2.0*F*F)+Kphi/alpha;
  for(int i=0; i<m2; i++){
    _r[k] = temp0[i];
    k++;
  }
  if (k != n)
    throw gen_err("residual not all updated");
  return _r;
}

void newton_bubble::update_f() {
  for(int k=0; k<n; k++) 
    grid.f(vec_mat(k, 0), vec_mat(k, 1)) = x[k];
}

void newton_bubble::compute_r() {
  update_f();
  r0 =  compute_r(grid.f);
  for(int i=0; i<n; i++)
    r[i] = r0[i];
}
 
void newton_bubble::compute_J() {

  J.reset_values(0.0);
  // upper
  int k, k1;
  Real value;
  for(int j=1; j<n1; j++) {
    for(int i=1; i<m1; i++) {
      k = mat_vec(i, j);
      value = grid.f(i,j);
      for(int t=1; t<n1; t++) {
	k1 = mat_vec(i, t);
	if (t != j)
	  J(k, k1) = grid.D3_2(t, j);
	else 
	  J(k, k1) = -fphiphi1(i, j)/value/value + grid.D1_2(i, i)/value - 2*fphi1(i,j)/value*(grid.D1_1(i,i)/value - fphi1(i,j)/value/value) + grid.D3_2(j,j);
      }
      for(int s=1; s<=m1; s++) {
	k1 = mat_vec(s, j);
	if (s != i)
	  J(k, k1) = grid.D1_2(i,s)/value - 2*fphi1(i,j)*grid.D1_1(i,s)/value/value;
      }
    }
  }
  // lower
  for(int j=1; j<n1; j++) {
    for(int i=1; i<=m2; i++) {
      k = mat_vec(i+m1, j);
      value = grid.f(i+m1, j);
      for(int t=0; t<n1; t++) {
	k1 = mat_vec(i+m1, t);
	if (t != j)
	  J(k, k1) = grid.D3_2(t, j);
	else
	  J(k, k1) = -fphiphi2(i, j)/value/value + grid.D2_2(i, i)/value - 2*fphi2(i,j)/value*(grid.D2_1(i,i)/value - fphi2(i,j)/value/value) + grid.D3_2(j,j);
      }
      for(int s=0; s<=m2; s++) {
	k1 = mat_vec(s+m1, j);
	if (s != i)
	  J(k, k1) = grid.D2_2(i,s)/value - 2*fphi2(i,j)*grid.D2_1(i,s)/pow(value,2);
      }
    }
  }
  // interface
  for(int j=1; j<n1; j++) {
    k = mat_vec(m1, j);
    for(int s=1; s<=m1; s++) {
      k1 = mat_vec(s, j);
      J(k, k1) += grid.D1_1(m1, s);
    }
    for(int s=0; s<=m2; s++) {
      k1 = mat_vec(s+m1, j);
      J(k, k1) -= grid.D2_1(0, s);
    }
  }
  // lower left
  Real tempp1=0.0, tempp2 =0.0, q_st=0.0, qphi_st=0.0, Kphi_st=0.0;
  for(int i=1; i<=m2; i++) {
    k = mat_vec(i+m1, 0);
    for(int s=1; s<=m2; s++) {
      q_st = -pow(q[i-1]/2, 3);
      k1 = mat_vec(s+m1, 0);
      if (s == i) {
	q_st *= f1[s-1]*(2*grid.D2_1(i,i)*ff[s-1]-f1[s-1])/ff[s-1]/ff[s-1] + 2*f2[s-1]*grid.D3_1(0,0);
	//qphi_st = -q_st*q[s-1]*q[s-1]*1.5/4*(f1[s-1]*(2*f11[s-1]*ff[s-1]-f1[s-1]*f1[s-1])/ff[s-1]/ff[s-1] + 2*f2[s-1]*f12[s-1]);
	qphi_st = q_st*qphi[i-1]*3.0/q[i-1];
	tempp1 = (grid.D2_1(i,i)*(2*f11[i-1]*ff[i-1]-f1[i-1]*f1[i-1]) + f1[s-1]*(2*grid.D2_2(i,i)*ff[s-1]+2*f11[s-1]-2*f1[s-1]*grid.D2_1(i,i)))/ff[s-1]/ff[s-1];
	tempp1 -= 2*f1[s-1]*(2*f11[s-1]*ff[i-1]-f1[s-1]*f1[s-1])/pow(ff[s-1], 3);
	tempp1 += 2*grid.D3_1(0,0)*f12[s-1]+2*f2[s-1]*grid.D2_1(i,i)*grid.D3_1(0,0);
	tempp1 *= pow(q[i-1]/2, 3);
	qphi_st -= tempp1;

	Kphi_st = -(grid.D2_1(i,i)*grid.D3_1(0,0)*q[s-1]+f12[s-1]*q_st + grid.D3_1(0,0)*qphi[s-1]+f2[s-1]*qphi_st)*sqrt(ff[s-1]);
	Kphi_st -= (f12[s-1]*q[s-1]+f2[s-1]*qphi[s-1])/2/sqrt(ff[s-1]);
	Kphi_st += (grid.D2_1(i,i)*f2[s-1]*q[s-1]+f1[s-1]*grid.D3_1(0,0)*q[s-1]+f1[s-1]*f2[s-1]*q_st)/2/sqrt(ff[s-1]);
	Kphi_st -= f1[s-1]*f2[s-1]*q[s-1]/4/pow(ff[s-1],1.5);
	Kphi_st /= 2*ff[s-1];
	Kphi_st += (f12[s-1]*q[s-1]+f2[s-1]*qphi[s-1]*sqrt(ff[s-1]) - f2[s-1]*f1[s-1]*q[s-1]/2/sqrt(ff[s-1]))/2/ff[s-1]/ff[s-1];
	
	tempp1 = 6*q[s-1]*q_st*qphi[s-1]*sqrt(ff[s-1]) + 3*q[s-1]*q[s-1]*qphi_st*sqrt(ff[s-1]) + 3*q[s-1]*q[s-1]*qphi[s-1]/2/sqrt(ff[s-1]);
	tempp1 -= ((3*q[s-1]*q[s-1]*q_st*f1[s-1]+pow(q[s-1],3)*grid.D2_1(i,i))*sqrt(ff[s-1]) - pow(q[s-1],3)*f1[s-1]/2/sqrt(ff[s-1]))/2/ff[s-1];
	tempp1 *= f2[s-1]*f11[s-1]-f1[s-1]*f12[s-1]-f2[s-1]*f1[s-1]*f1[s-1]/2/ff[s-1];
	Kphi_st += tempp1/4/ff[i-1];
	tempp2 = (3*q[s-1]*q[s-1]*qphi[s-1]*sqrt(ff[s-1])-pow(q[s-1],3)*f1[s-1]/2/sqrt(ff[s-1]))/4/ff[s-1];
	tempp2 *= grid.D3_1(0,0)*f11[s-1]+f2[s-1]*grid.D2_2(i,i)-grid.D2_1(i,i)*f12[s-1]-f1[s-1]*grid.D2_1(i,i)*grid.D3_1(0,0)-((grid.D3_1(0,0)*f1[s-1]*f1[s-1]+2*f2[s-1]*f1[s-1]*grid.D2_1(i,i))*ff[s-1]-f2[s-1]*f1[s-1]*f1[s-1])/2/ff[s-1]/ff[s-1];
	Kphi_st += tempp2;
        tempp1 = (-3*q[s-1]*q[s-1]*qphi[s-1]*sqrt(ff[s-1])+pow(q[s-1],3)*f1[s-1]/2/sqrt(ff[s-1]))/4/ff[i-1]/ff[i-1];
	tempp1 *= f2[s-1]*f11[s-1]-f1[s-1]*f12[s-1]-f2[s-1]*f1[s-1]*f1[s-1]/2/ff[s-1];
	Kphi_st += tempp1;

	
	tempp1 = (3*q[s-1]*q[s-1]*q_st*sqrt(ff[s-1])-pow(q[s-1],3)/2/sqrt(ff[s-1]))/4/ff[s-1];
	tempp1 *= f2[s-1]*f111[s-1]-f1[s-1]*f112[s-1] - ((f12[s-1]*f1[s-1]*f1[s-1]+2*f2[s-1]*f1[s-1]*f11[s-1])*ff[s-1] - f2[s-1]*pow(f1[s-1],3))/2/ff[s-1]/ff[s-1];
	Kphi_st += tempp1;
	tempp2 = grid.D3_1(0,0)*f111[s-1]+f2[s-1]*grid.D2_3(i,i)-grid.D2_1(i,i)*f112[s-1]-f1[s-1]*grid.D2_2(i,i)*grid.D3_1(0,0);
	tempp1 = ff[s-1]*(grid.D2_1(i,i)*grid.D3_1(0,0)*f1[s-1]*f1[s-1]+2*f12[s-1]*f1[s-1]*grid.D2_1(i,i)+2*grid.D3_1(0,0)*f1[s-1]*f11[s-1]+2*f2[s-1]*grid.D2_1(i,i)*f11[s-1]+2*f1[s-1]*f2[s-1]*grid.D2_2(i,i));
	tempp1 -= grid.D3_1(0,0)*pow(f1[s-1],3)+3*f2[s-1]*f1[s-1]*f1[s-1]*grid.D2_1(i,i);
	tempp1 /= 2*ff[s-1]*ff[s-1];
	tempp2 -= tempp1;
	tempp2 += ((f12[s-1]*f1[s-1]*f1[s-1]+2*f2[s-1]*f1[s-1]*f11[s-1])*ff[s-1]-f2[s-1]*pow(f1[s-1],3))/pow(ff[s-1],3);
	tempp2 *= pow(q[s-1],3)/4/sqrt(ff[s-1]);
	Kphi_st += tempp2;
	J(k, k1) = q_st*qphi[s-1]+q[s-1]*qphi_st-grid.D3_1(0,0)/2/F/F+Kphi_st/alpha;
      }
      else {
	q_st *= f1[i-1]*2*grid.D2_1(i,s)/ff[i-1];
	
	qphi_st = q_st*qphi[i-1]*3.0/q[i-1];
	tempp1 = (grid.D2_1(i,s)*(2*f11[i-1]*ff[i-1]-f1[i-1]*f1[i-1]) + f1[i-1]*(2*grid.D2_2(i,s)*ff[i-1]-2*f1[i-1]*grid.D2_1(i,s)))/ff[i-1]/ff[i-1];
	tempp1 += 2*f2[i-1]*grid.D2_1(i,s)*grid.D3_1(0,0);
	tempp1 *= -pow(q[i-1]/2, 3);
	qphi_st += tempp1;
	
	Kphi_st = -(grid.D2_1(i,s)*grid.D3_1(0,0)*q[i-1]+f12[i-1]*q_st+f2[i-1]*qphi_st)*sqrt(ff[i-1]);
	Kphi_st += (grid.D2_1(i,s)*f2[i-1]*q[i-1]+f1[i-1]*f2[i-1]*q_st)/2/sqrt(ff[i-1]);
        Kphi_st /= 2*ff[i-1];
	
	tempp1 = 6*q[i-1]*q_st*qphi[i-1]*sqrt(ff[i-1]) + 3*q[i-1]*q[i-1]*qphi_st*sqrt(ff[i-1]);
	tempp1 -= (3*q[i-1]*q[i-1]*q_st*f1[i-1]+pow(q[i-1],3)*grid.D2_1(i,s))/sqrt(ff[i-1])/2;
	tempp1 *= (f2[i-1]*f11[i-1]-f1[i-1]*f12[i-1]-f2[i-1]*f1[i-1]*f1[i-1]/2/ff[i-1]);
	tempp1 /= 4*ff[i-1];
	tempp2 = (3*q[i-1]*q[i-1]*qphi[i-1]*sqrt(ff[i-1])-pow(q[i-1],3)*f1[i-1]/2/sqrt(ff[i-1]))/4/ff[i-1];
	tempp2 *= f2[i-1]*grid.D2_2(i,s)-grid.D2_1(i,s)*f12[i-1]-f1[i-1]*grid.D2_1(i,s)*grid.D3_1(0,0)-2*f2[i-1]*f1[i-1]*grid.D2_1(i,s)/2/ff[i-1];
	Kphi_st += tempp1+tempp2;
	tempp1 = 3*q[i-1]*q[i-1]*q_st/4/sqrt(ff[i-1]);
	tempp1 *= f2[i-1]*f111[i-1]-f1[i-1]*f112[i-1] - ((f12[i-1]*f1[i-1]*f1[i-1]+2*f2[i-1]*f1[i-1]*f11[i-1])*ff[i-1] - f2[i-1]*pow(f1[i-1],3))/2/ff[i-1]/ff[i-1];
	Kphi_st += tempp1;
	tempp2 = f2[i-1]*grid.D2_3(i,s)-grid.D2_1(i,s)*f112[i-1]-f1[i-1]*grid.D2_2(i,s)*grid.D3_1(0,0);	  
	tempp1 = ff[i-1]*(grid.D2_1(i,s)*grid.D3_1(0,0)*f1[i-1]*f1[i-1]+2*f12[i-1]*f1[i-1]*grid.D2_1(i,s)+2*f2[i-1]*grid.D2_1(i,s)*f11[i-1]+2*f1[i-1]*f2[i-1]*grid.D2_2(i,s));
	tempp1 -= 3*f2[i-1]*f1[i-1]*f1[i-1]*grid.D2_1(i,s);
	tempp1 /= 2*ff[i-1]*ff[i-1];
	tempp2 -= tempp1;
	tempp2 *= pow(q[i-1],3)/4/sqrt(ff[i-1]);
	Kphi_st += tempp2;
	J(k, k1) = q_st*qphi[i-1]+q[i-1]*qphi_st+Kphi_st/alpha;
      }
    }
    for(int t=1; t<n1; t++) {
      k1 = mat_vec(i+m1, t);
      q_st = -pow(q[i-1]/2, 3)*2*f2[i-1]*grid.D3_1(t,0);
      
      qphi_st = qphi[i-1]*q_st*3.0/q[i-1];
      tempp1 = 2*grid.D3_1(t,0)*f12[i-1]+2*f2[i-1]*grid.D2_1(i,i)*grid.D3_1(t,0);
      tempp1 *= -pow(q[i-1]/2, 3);
      qphi_st += tempp1;
      
      Kphi_st = -(grid.D2_1(i,i)*grid.D3_1(t,0)*q[i-1]+f12[i-1]*q_st + grid.D3_1(t,0)*qphi[i-1]+f2[i-1]*qphi_st)*sqrt(ff[i-1]);
      Kphi_st += (f1[i-1]*grid.D3_1(t,0)*q[i-1]+f1[i-1]*f2[i-1]*q_st)/2/sqrt(ff[i-1]);
      Kphi_st /= 2*ff[i-1];
      
      tempp1 = 6*q[i-1]*q_st*qphi[i-1]*sqrt(ff[i-1]) + 3*q[i-1]*q[i-1]*qphi_st*sqrt(ff[i-1]);
      tempp1 -= 3*q[i-1]*q[i-1]*q_st*f1[i-1]/2/sqrt(ff[i-1]);
      tempp1 *= (f2[i-1]*f11[i-1]-f1[i-1]*f12[i-1]-f2[i-1]*f1[i-1]*f1[i-1]/2/ff[i-1]);
      tempp2 = (3*q[i-1]*q[i-1]*qphi[i-1]*sqrt(ff[i-1])-pow(q[i-1],3)*f1[i-1]/2/sqrt(ff[i-1]))/4/ff[i-1];
      tempp2 *= grid.D3_1(t,0)*f11[i-1]-f1[i-1]*grid.D2_1(i,i)*grid.D3_1(t,0)-grid.D3_1(t,0)*f1[i-1]*f1[i-1]/2/ff[i-1];
      Kphi_st += tempp1/4/ff[i-1] + tempp2;

      
      tempp1 = 3*q[i-1]*q[i-1]*q_st/4/sqrt(ff[i-1]);
      tempp1 *= f12[i-1]*f11[i-1]+f2[i-1]*f111[i-1]-f11[i-1]*f12[i-1]-f1[i-1]*f112[i-1] - ((f12[i-1]*f1[i-1]*f1[i-1]+2*f2[i-1]*f1[i-1]*f11[i-1])*ff[i-1] - f2[i-1]*pow(f1[i-1],3))/2/ff[i-1]/ff[i-1];
      Kphi_st += tempp1;
      tempp2 = grid.D3_1(t,0)*f111[i-1]-f1[i-1]*grid.D2_2(i,i)*grid.D3_1(t,0);
      tempp1 = ff[i-1]*(grid.D2_1(i,i)*grid.D3_1(t,0)*f1[i-1]*f1[i-1]+2*grid.D3_1(t,0)*f1[i-1]*f11[i-1]);
      tempp1 -= grid.D3_1(t,0)*pow(f1[i-1],3);
      tempp1 /= 2*ff[i-1]*ff[i-1];
      tempp2 -= tempp1;
      tempp2 *= pow(q[i-1],3)/4/sqrt(ff[i-1]);
      Kphi_st += tempp2;
      J(k, k1) = q_st*qphi[i-1]+q[i-1]*qphi_st-grid.D3_1(t,0)/2/F/F+Kphi_st/alpha;
    }
  }

  for(int i=1; i<=m2; i++) {
    k = mat_vec(i+m1, 0);
    q_st = 0.0;
    for(int t=1; t<n1; t++) {
      for(int s=0; s<=m2; s++) {
	if ( s!=i ) {
	  k1 = mat_vec(s+m1, t);
	  qphi_st = -pow(q[i-1]/2, 3);
	  qphi_st *= 2*f2[i-1]*grid.D2_1(i,s)*grid.D3_1(t,0);
	  
	  Kphi_st = -(grid.D2_1(i,s)*grid.D3_1(t,0)*q[i-1]+f2[i-1]*qphi_st)/2/sqrt(ff[i-1]);
	  
	  tempp1 = 3*q[i-1]*q[i-1]*qphi_st/4/sqrt(ff[i-1]); 
	  tempp1 *= f2[i-1]*f11[i-1]-f1[i-1]*f12[i-1]-f2[i-1]*f1[i-1]*f1[i-1]/2/ff[i-1];
	  tempp2 = (3*q[i-1]*q[i-1]*qphi[i-1]*sqrt(ff[i-1])-pow(q[i-1],3)*f1[i-1]/2/sqrt(ff[i-1]))/4/ff[i-1];
	  tempp2 *= -f1[i-1]*grid.D2_1(i,s)*grid.D3_1(t,0);
	  Kphi_st += tempp1+ tempp2;

	  tempp2 = -f1[i-1]*grid.D2_2(i,s)*grid.D3_1(t,0);
	  tempp1 = grid.D2_1(i,s)*grid.D3_1(t,0)*f1[i-1]*f1[i-1]/2/ff[i-1];
	  tempp2 -= tempp1;
	  tempp2 *= pow(q[i-1],3)/4/sqrt(ff[i-1]);
	  Kphi_st += tempp2;
	  J(k, k1) = q_st*qphi[i-1]+q[i-1]*qphi_st+Kphi_st/alpha;
	}
      }
    }
  }
  /*
  J.dump("J", 17, '%');
  mp_mat<Real> g1(m1+m2+1, n1), g2(m1+m2+1, n1);
  const Real eps = 1.0e-7;
  int i, j;
  for(int k=0; k<n; k++) {
    i = vec_mat(k, 0);
    j = vec_mat(k, 1);
    g1 = grid.f;
    g1(i, j) += eps;
    g2 = grid.f;
    g2(i, j) -= eps;
    r1 = compute_r(g1);
    r2 = compute_r(g2);
    r1 -= r2;
    for(int l=0; l<n; l++) 
      J(l, k) = r1[l]/(2*eps);
  }
  J.dump("P", 17, '%');
  */
}


void newton_bubble::shape() {
  FILE *fp1 = fopen("shape","w");
    for(int k=0; k<m2; k++) {
      fprintf(fp1, "%23s\n", str(grid.xx2[k+1],0));
      fprintf(fp1, "%23s\n", str(sqrt(x[n-m2+k]),0));
    }
}
