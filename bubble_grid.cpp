#include "bubble_grid.h"
bubble_grid::bubble_grid(int _m, int _n) :
  m(_m),
  n(_n),
  M((m+1)*(n+1)),
  phi1(4.0),
  phi2(20.0),
  x(0.0, m+1),
  x1(0.0, m+1),
  y(0.0, n+1),  //x for phi, y for psi
  y1(0.0, n+1), //x, y are for reference. x1 y1 are the real grid 
  D1_1(m+1, m+1, 0.0),
  D1_2(m+1, m+1, 0.0),
  D1_3(m+1, m+1, 0.0),
  D2_1(n+1, n+1, 0.0),
  D2_2(n+1, n+1, 0.0),
  map1(-4.0), // [1,-1] -> [-psi1,psi2]=[0,1/2]
  map2(-2/(phi1+phi2)), // [1,-1] -> [-phi1, phi2]:y = -(phi1+phi2)/2*(x-1)-phi1
  f(m+1, n+1, 10),// need to change here?
  tmp(m+1, n+1, 0.0),
  tmp1(m, n-1, 0.0),
  tmp2(m, n-1, 0.0),
  tmp3(m, n-1, 0.0)
{
  set_xy();
  set_A();
  set_N();
  set_D();
  set_f();
}

// this is only for getting N, M to build up levmar
bubble_grid::bubble_grid(int _m, int _n, int flag):
  m(_m),
  n(_n),
  M((m+1)*(n+1)),
  phi1(4.0),
  phi2(20.0),
  x(0.0, m+1),
  x1(0.0, m+1),
  y(0.0, n+1),
  y1(0.0, n+1),
  map1(-4.0),
  map2(-2/(phi1+phi2))
{
  set_xy();
  set_A();
  set_N();
}

// y = -(phi1+phi2)/2*(x-1)-phi1
void bubble_grid::set_xy() {
  for(int i=0; i<=m; i++) {
    x[i] = cos(pi*i/m);
    x1[i] = 1/map2*(x[i]-1)-phi1;
  }
  for(int i=0; i<=n; i++) {
    y[i] = cos(pi*i/n);
    y1[i] = 1/map1*(y[i]-1);
  }
}

void bubble_grid::set_A() {
  int i = 0;
  while (i<n+1) {
    if (x1[i] >= 0) break; // need to change here!!!
    i++;
  }
  A = i;
}

void bubble_grid::set_N() {
  N = m*(n-1)+m+1-A;
}

void bubble_grid::set_D() {
   valarray<Real> c1(1.0, m+1), c2(1.0, n+1);
    c1[0] = c1[m] = 2.0;
    c2[0] = c2[n] = 2.0;
    for(int i=0; i<=m; i++) {
      for(int j=0; j<=m; j++) {
	if(i == j) {
	  if(i == 0) D1_1(0,0) = (2*m*m+1)/6.0;
	  else {
	    if (i == m) D1_1(i,i) = -(2*m*m+1)/6.0;
	    else D1_1(i,j) = -x[i]/2/(1-x[i]*x[i]);
	  }
	}
	else D1_1(i,j) = c1[i]*pow(-1, (i+j)%2)/c1[j]/(x[i]-x[j]);
      }
    }
    D1_1 *= map1;
    
    for(int i=1; i<m; i++) {
      for(int j=0; j<=m; j++) {
        if (i != j) D1_2(i,j) = pow(-1.0, (i+j)%2)*(x[i]*x[i]+x[i]*x[j]-2)/c1[j]/(1-x[i]*x[i])/(x[i]-x[j])/(x[i]-x[j]);
	else D1_2(i,i) = -((m*m-1)*(1-x[i]*x[i])+3)/3/(1-x[i]*x[i])/(1-x[i]*x[i]);
      }
    }
    for(int j=0; j<m; j++) {
      D1_2(0,j+1) = 2*pow(-1, (j+1)%2)*((2*m*m+1)*(1-x[j+1])-6)/3/c1[j+1]/(1-x[j+1])/(1-x[j+1]);
      D1_2(m, j) = 2*pow(-1, (m+j)%2)*((2*m*m+1)*(1+x[j])-6)/3/c1[j]/(1+x[j])/(1+x[j]);
    }
    D1_2(0,0) = D1_2(m,m) = (pow(m, 4)-1)/5.0;
    D1_2 *= map1*map1;

    //D1_3 to be implemented here
    
    for(int i=1; i<n; i++) {
      for(int j=0; j<=n; j++) {
        if (i != j) D2_2(i,j) = pow(-1.0, (i+j)%2)*(x[i]*x[i]+x[i]*x[j]-2)/c2[j]/(1-x[i]*x[i])/(x[i]-x[j])/(x[i]-x[j]);
	else D2_2(i,i) = -((n*n-1)*(1-x[i]*x[i])+3)/3/(1-x[i]*x[i])/(1-x[i]*x[i]);
      }
    }
    for(int j=0; j<n; j++) {
      D2_2(0,j+1) = 2*pow(-1, (j+1)%2)*((2*n*n+1)*(1-x[j+1])-6)/3/c2[j+1]/(1-x[j+1])/(1-x[j+1]);
      D2_2(n, j) = 2*pow(-1, (n+j)%2)*((2*n*n+1)*(1+x[j])-6)/3/c2[j]/(1+x[j])/(1+x[j]);
    }
    D2_2(0,0) = D1_2(n,n) = (pow(n, 4)-1)/5.0;
    D2_2 *= map2*map2;
}



// need to initialize all of them and then use mat_vec to initialize f.
void bubble_grid::set_f() {
  for (int i=0; i<=m; i++) {
    f(i, n) = 1;
    if (i < A) f(i, 0) = 0;
  }
  for (int j=0; j<=n; j++) {
    f(0,j) = 2*y[j];
  }
}


int bubble_grid::delta(int _i, int _j) {
  if(_i == _j) return 1;
  return 0;
}
