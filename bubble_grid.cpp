#include "bubble_grid.h"
bubble_grid::bubble_grid(int _m1, int _m2, int _n) :
  m1(_m1),
  m2(_m2),
  n(_n),
  N1((m1-1)*(n-1)),
  N2((m2+1)*n-1),
  N(N1+N2),
  phi1(4.0),
  phi2(20.0),
  x1(0.0, m1+1),
  xx1(0.0, m1+1),
  x2(0.0, m2+1),
  xx2(0.0, m2+1),
  y(0.0, n+1),  //x for phi, y for psi
  yy(0.0, n+1), //x, y are for reference. x1 y1 are the real grid 
  D1_1(m1+1, m1+1, 0.0),
  D1_2(m1+1, m1+1, 0.0),
  D2_1(m2+1, m2+1, 0.0),
  D2_2(m2+1, m2+1, 0.0),
  D2_3(m2+1, m2+1, 0.0),
  D3_1(n+1, n+1, 0.0),
  D3_2(n+1, n+1, 0.0),
  map1(-2.0/phi1),
  map2(-2.0/phi2),
  map3(-4.0),
  f(m1+m2+1, n+1, 0.0), 
  tmp11(m1-1, n-1, 0.0),
  tmp12(m1-1, n-1, 0.0),
  tmp13(m1-1, n-1, 0.0),
  tmp21(m2, n-1, 0.0),
  tmp22(m2, n-1, 0.0),
  tmp23(m2, n-1, 0.0)
{ 
  set_xy();
  set_D();
  set_f();
}

// y = -(phi1+phi2)/2*(x-1)-phi1
void bubble_grid::set_xy() {
  for(int i=0; i<=m1; i++) {
    x1[i] = cos(pi*i/m1);
    xx1[i] = 1/map1*(x1[i]+1);
  }
  for(int i=0; i<=m2; i++) {
    x2[i] = cos(pi*i/m2);
    xx2[i] = 1/map2*(x2[i]-1);
  }
  for(int i=0; i<=n; i++) {
    y[i] = cos(pi*i/n);
    yy[i] = 1/map3*(y[i]-1);
  }
}

void bubble_grid::set_D() {
  valarray<Real> c1(1.0, m1+1), c2(1.0, m2+1), c3(1.0, n+1);
    c1[0] = c1[m1] = 2.0;
    c2[0] = c2[m2] = 2.0;
    c3[0] = c3[n] = 2.0;
    // D1_1
    for(int i=0; i<=m1; i++) {
      for(int j=0; j<=m1; j++) {
	if(i == j) {
	  if(i == 0) D1_1(0,0) = (2*m1*m1+1)/6.0;
	  else {
	    if (i == m1) D1_1(i,i) = -(2*m1*m1+1)/6.0;
	    else D1_1(i,j) = -x1[i]/2/(1-x1[i]*x1[i]);
	  }
	}
	else D1_1(i,j) = c1[i]*pow(-1, (i+j)%2)/c1[j]/(x1[i]-x1[j]);
      }
    }
    D1_1 *= map1;
    D1_2 = D1_1*D1_1;

    // D2_1
    for(int i=0; i<=m2; i++) {
      for(int j=0; j<=m2; j++) {
	if(i == j) {
	  if(i == 0) D2_1(0,0) = (2*m2*m2+1)/6.0;
	  else {
	    if (i == m2) D2_1(i,i) = -(2*m2*m2+1)/6.0;
	    else D2_1(i,j) = -x2[i]/2/(1-x2[i]*x2[i]);
	  }
	}
	else D2_1(i,j) = c2[i]*pow(-1, (i+j)%2)/c2[j]/(x2[i]-x2[j]);
      }
    }
    D2_1 *= map2;
    D2_2 = D2_1*D2_1;
    D2_3 = D2_1*D2_2;
    

    // D3_1^T
    for(int j=0; j<=n; j++) {
      for(int i=0; i<=n; i++) {
	if(j == i) {
	  if(j == 0) D3_1(0,0) = (2*n*n+1)/6.0;
	  else {
	    if (j == n) D3_1(j,j) = -(2*n*n+1)/6.0;
	    else D3_1(i,j) = -y[j]/2/(1-y[j]*y[j]);
	  }
	}
	else D3_1(i,j) = c3[j]*pow(-1, (i+j)%2)/c3[i]/(y[j]-y[i]);
      }
    }
    D3_1 *= map3;
    D3_2 = D3_1*D3_1;
}



// need to initialize all of them and then use mat_vec to initialize x.
void bubble_grid::set_f() {
  for (int i=0; i<=m1; i++) {
    f(i, n) = 1.0;
    f(i, 0) = 0.0;
  }
  
  Real temp;
  for (int j=1; j<n; j++) {
    temp = 2*yy[j];
    for(int i=0; i<m1; i++)
      f(i, j) = temp;
  }
  for(int i=1; i<=m2; i++)
    f(i+m1, n) = 1.0;
  
  for(int j=0; j<n; j++)
    for(int i=0; i<=m2; i++)
      //f(i+m1, j) = sqrt(xx2[i]/20)*(1-2*yy[j])+2*yy[j];
      //f(i+m1, j) = sqrt(xx2[i]/22.0)*(1-2*yy[j])+2*yy[j];
      f(i+m1, j) = 2.0/pi*atan(1.6*pow(xx2[i], 2.0/3))*(1-2*yy[j])+2*yy[j];
}


/*
int bubble_grid::delta(int _i, int _j) {
  if(_i == _j) return 1;
  return 0;
}
*/
