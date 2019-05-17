#include "potential.h"
#include "str_double.h"

int main(int argc, char*argv[]) {
  potential bubble(20, 60);
  bubble.set_delta(1.0e-10);
  //Real bb = bubble.compute_B(1000000.0/1000001,1.0/1000001);
  //Real dd = bubble.compute_D(1000000.0/1000001,1.0/1000001);
  //printf("%23s, %23s\n", str( (bb+dd/1000001)*sqrt(1000001), 0), str((bb+dd)/sqrt(1000001),0 ));
  //bubble.compute_G(1,2,3,2.999999);
  //printf("%23s\n",str(bubble.G, 0) );
  bubble.compute_r();
  bubble.compute_J();
  //bubble.chkder();
  bubble.solve();
  FILE *fp = fopen("initial_new", "w");
  for(int i=0; i<=bubble.N; i++) 
    fprintf(fp, "%23s\n", str(bubble.x[i],0));
  //bubble.interpolate(50);

}

/*
#include "newton_bubble.h"
#include "str_double.h"

int main(int argc, char*argv[]) {
  
  if (argc != 3)
    throw gen_err("uasage: ./test_levmar F, alpha");
  Real F, alpha;
  str_to_real(argv[1], F);
  str_to_real(argv[2], alpha);
  
  bubble_grid grid(20, 80, 20);
  
  newton_bubble bubble(grid, 1);
  bubble.set_delta(1.0e-8);
  bubble.compute_r();
  bubble.compute_J();
  
  bubble.solve();
  bubble.compute_g();
  bubble.dump();
  bubble.shape();
  
  return 0;
}
*/
