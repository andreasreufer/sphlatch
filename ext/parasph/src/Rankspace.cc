// Input:  Particle.pos   (physical coordinates)
// Output: Particle.coord (rankcell coordinates)

#ifndef RANKSPACE_CC
#define RANKSPACE_CC

#include "Array.cc"
#include "quick.cc"
#include "Particle.cc"
#include "TList.cc"

class Rankspace {
private:
  Array<ftype>     input;
  int             *disp, *index, *point, size, *tot;
  ftype            slicem1, *sort;
  TList<Particle> &pList;
  
  void pSizes() {
    size = pList.getSize();
    input.setSize(size);

    MCW.Allgather(&size, 1, MPI::INT, tot, 1, MPI::INT);
    disp[0] = point[0] = 0;
    for (int i = 1; i < global::npro; i++) 
      disp[i] = point[i] = disp[i-1] + tot[i-1];
  }

public:
  Rankspace(TList<Particle> &_pList) : pList(_pList) {}

  void init() {
    slicem1 = 1. / (ftype)global::slice;

    index  = new int[global::totNumPart];
    sort   = new ftype[global::totNumPart];
    disp   = new int[global::npro];
    point  = new int[global::npro];
    tot    = new int[global::npro];
  }

  void pRanking(const int &co) {
    int   i, nr, p, boss = co % global::npro, win;
    ftype min;

    pSizes();
    for (p = 0; p < size; p++) input[p] = pList[p].getPos()[co];
    quickStart(&input[0], index, size);
    for (p = 0; p < size; p++) 
      input[p] = pList[index[p]].getPos()[co];
    MCW.Gatherv(&input[0], size, MPI_ftype, sort, tot, disp, MPI_ftype, boss);

    if (global::rank == boss)
      for (nr = 0; nr < global::totNumPart; nr++) {	
	for (i = 0, min = 1.e30, win = -1; i < global::npro; i++)
	  if (point[i] < (tot[i]+disp[i]) && sort[point[i]] < min) { 
	    min = sort[point[i]]; win = i;
	  }
	sort[point[win]++] = nr;
      } 
   
    MCW.Scatterv(sort, tot, disp, MPI_ftype, &input[0], size, MPI_ftype, boss);
    for (p = 0; p < size; p++)
      pList[index[p]].setCoord(co, (int)(input[p]*slicem1));
  }
};

/* To test if particles are distributed equally among cells
  int size = localList.getSize();
  int *cell, *All, i, p, gs = global::dim3;
  cell = new int[global::dim3]; All = new int[global::dim3];
  for (i = 0; i < global::dim3; i++) cell[i] = 0;

  Vector<int> c;
  for (p = 0; p < size; p++) {
    c = localList(p)->getCoord();
    cell[c[0]*global::dim2 + c[1]*global::dim + c[2]]++;
  }

  MCW.Allreduce(cell, All, global::dim3, MPI::INT, MPI::SUM);

  if (global::rank == 0) {
    int h = 100, *count = new int[h];
    for (i = 0; i < h; i++) count[i] = 0;
    for (i = 0; i < global::dim3; i++) count[All[i]]++;
    for (i = 0; i < h; i++) if (count[i] > 0) 
      std::cout << "count[" << i << "] = " << count[i] << std::endl;
  }
*/

#endif
