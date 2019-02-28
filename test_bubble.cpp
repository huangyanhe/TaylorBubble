#include "newton_bubble.h"
#include "str_double.h"

int main(int argc, char*argv[]) {
  /*
  if (argc != 3)
    throw gen_err("uasage: ./test_levmar F, alpha");
  Real F, alpha;
  str_to_real(argv[1], F);
  str_to_real(argv[2], alpha);
  */
  
  bubble_grid grid(10, 40, 10);
  newton_bubble bubble(grid);
  bubble.set_delta(.6);
  bubble.compute_r();
  bubble.compute_J();
  
  bubble.solve();
  bubble.compute_g();
  bubble.dump();
  bubble.shape();
  return 0;
}
