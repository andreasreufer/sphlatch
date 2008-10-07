#include <iostream>
#include "quick.cc"

int main(void) {
  int  f[11] = {6, 3, 7, 2, 4, 9, 5, 6, 4, 3, 0};
  int  d[11], r[11];

  quickStart(f, d, 11);
  
  for (int i = 0; i < 11; i++) { 
    std::cout << f[i] << " ";
    r[d[i]] = i;
  }
  std::cout << std::endl << std::endl;
  for (int i = 0; i < 11; i++) std::cout << f[d[i]] << " "; std::cout << std::endl;
  for (int i = 0; i < 11; i++) std::cout << d[i] << " "; std::cout << std::endl;
  for (int i = 0; i < 11; i++) std::cout << r[i] << " "; std::cout << std::endl;
}
