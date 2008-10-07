// Input:  Particle.coord (rankcell coordinates)
// Output: Particle.near  (range of neighbouring cells)

#ifndef CELLMAP_CC
#define CELLMAP_CC

#include "Def.cc"
#include "Particle.cc"
#include "TList.cc"
#include "Vector.cc"

class CellMap {
private:
  TList<Particle> &pList;
  Vector<ftype>   *glob, *mapmax, *mapmin;

public:
  CellMap(TList<Particle> &_pList) : pList(_pList) {}

  void init() {
    glob   = new Vector<ftype>[global::dim];
    mapmax = new Vector<ftype>[global::dim];
    mapmin = new Vector<ftype>[global::dim];
  }

  void pCreate() {
    int           c, co, p;
    Vector<ftype> pos;
    Vector<int>   coord;

    for (c = 0; c < global::dim; c++) { 
      mapmin[c] = maxVec; mapmax[c] = minVec; 
    }
    mapmin[0] = minVec; mapmax[global::dim-1] = maxVec;

    for (p = 0; p < pList.getSize(); p++) {
      coord = pList[p].getCoord();
      pos   = pList[p].getPos();
      for (co = 0; co < 3; co++) {
	mapmin[coord[co]][co] = Min(pos[co], mapmin[coord[co]][co]);
	mapmax[coord[co]][co] = Max(pos[co], mapmax[coord[co]][co]);
      }
    }

    MCW.Allreduce(mapmin, glob, 3*global::dim, MPI_ftype, MPI::MIN);
    for (c = 0; c < global::dim; c++) mapmin[c] = glob[c];
    MCW.Allreduce(mapmax, glob, 3*global::dim, MPI_ftype, MPI::MAX);
    for (c = 0; c < global::dim; c++) mapmax[c] = glob[c];

    /*
    // testing
    for (p = 0; p < pList.getSize(); p++) {
      Vector<ftype> pos   = pList[p].getPos();
      Vector<int>   coord = pList[p].getCoord();

      bool ok1 = true;
      for (int co = 0; co < 3; co++) {
	if (pos[co] < mapmin[coord[co]][co]) ok1 = false;
	if (pos[co] > mapmax[coord[co]][co]) ok1 = false;
      }
      if (!ok1) std::cout << "Particle " << pos << " not within bounds "
			 << coord << std::endl;
    }

    for (int co = 0; co < 3; co++) {
      for (int i = 0; i < global::dim; i++) {
	if (mapmin[i][co] > mapmax[i][co]) 
	  std::cout << "mapmin > mapmax" << std::endl;
	if (i < global::dim-1) if (mapmin[i][co] > mapmin[i+1][co]) 
	  std::cout << "mapmin i > mapmin i+1" << std::endl;
	if (i < global::dim-1) if (mapmax[i][co] > mapmax[i+1][co]) 
	  std::cout << "mapmax i > mapmax i+1" << std::endl;
      }
    }
    */
  }

  void findNear() {
    ftype         hmax;
    int           co, p;
    Particle     *part;
    Vector<ftype> pos;
    Vector<int>   tmp;
    Range<int>    near;

    for (p = 0; p < pList.getSize(); p++) {
      part = &pList[p];
      hmax = part->getH() * 2.;
      pos  = part->getPos();
      near.set(part->getCoord(), part->getCoord());	

      for (co = 0; co < 3; co++) {
	for (tmp = near.min; tmp[co] >= 0; tmp[co]--)
	  if (mapmax[tmp[co]][co] >= pos[co] - hmax)
	    near.min[co] = tmp[co]; else tmp[co] = -1;
	for (tmp = near.max; tmp[co] < global::dim; tmp[co]++)
	  if (mapmin[tmp[co]][co] <= pos[co] + hmax)
	    near.max[co] = tmp[co]; else tmp[co] = global::dim;
      }
      near.max = near.max + 1;
      part->setNear(near);
    }
  }
};

#endif
