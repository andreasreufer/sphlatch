/******************************************************************************
 * ParaSPH -- Version 24.11.2003                                              *
 *----------------------------------------------------------------------------*
 * File:      CommStrat.cc                                                    *
 * Purpose:   This class implements a communication strategy found in my      *
 *            thesis (chapter 4.1.2, pp. 38-41).                              *
 *            That was quite tough, so do not hope to understand at the first *
 *            glimpse. If you want to test, just uncomment the comments '//!',*
 *            compile as a normal program (input is the number of processes). *
 *****************************************************************************/
#ifndef COMMSTRAT_CC
#define COMMSTRAT_CC

//!#include <iostream>

class CommStrat {
private:
  int rank, npro;
  typedef signed char cstype;
  cstype **comm, rounds;
  
  void set(const int &round, const int &one, const int &two) {
    comm[round][one] = two; comm[round][two] = one;
  }

  // 
  void func(int from, int to) {
    int i, j, num = to - from, half = num/2, mid = from + half;

    if (num % 2 == 0) { 
      func(from, mid); func(mid, to);
      for (i = 0; i < half; i++) for (j = from; j < mid; j++) {
	if (comm[i][j] == j) set(i, j, j + half); 
	if (i >= half % 2) set(half + i - 1, j, (j + i) % half + mid);
      }
    } else 
      for (i = 0; i < num; i++) for (j = from; j < to; j++)
	set(i, j, from + (i + 1 - j + npro) % num);
  }
public:
  // Without initialization the program will die when calling get ...
  void init(const int &_rank, const int &_npro) {
    rank = _rank; npro = _npro;
    rounds = (npro - 1) / 2 * 2 + 1;
    comm  = new cstype*[rounds];
    for (int i = 0; i < rounds; i++) comm[i] = new cstype[npro];    
    func(0, npro);
//!    show();
  }

  // Get the process id to communicate with at a certain round
  cstype get(const int &round) { return comm[round][rank]; }

  // Get the number of communication rounds needed
  cstype getRounds() const { return rounds; }

//!  void show() {
//!    for (int i = 0; i < rounds; i++, std::cout << std::endl) 
//!      for (int j = 0; j < npro; j++)
//!     	  if (comm[i][j] == j) std::cout << "* ";
//!  	  else std::cout << (int)comm[i][j] << " ";
//!  }
};

#endif

//!int main(int argc, char *argv[]) { CommStrat c; c.init(0, atol(argv[1])); }
