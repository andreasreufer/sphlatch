// Input:  Particle.coord, Particle.near
// Output: Particle.next

#ifndef CELLLIST_CC
#define CELLLIST_CC

#include "Array.cc"
#include "Cell.cc"
#include "Debug.cc"
#include "Def.cc"
#include "Ghost.cc"
#include "Particle.cc"
#include "TList.cc"
#include "Tree.cc"

class CellList {
private:
  Array<int>         newtonCell;
  Cell              *cellList;
  int               *lRequest, *gRequest;
  TList<GhostRead>  *readGhostList,  tmpGhostRList;
  TList<GhostWrite> *writeGhostList, tmpGhostWList;
  TList<Particle>   &pList;
  Tree               tree;
  
  int calcIndex(const Vector<int> &c) {
    return c[0]*global::dim2 + c[1]*global::dim + c[2];
  }

  void calcRequest(Particle *part) {
    int        i, j, k, proc = 1 << global::rank;
    Range<int> near = part->getNear();

    near.min[0] *= global::dim2; near.max[0] *= global::dim2;   
    near.min[1] *= global::dim;  near.max[1] *= global::dim;   

    for (i = near.min[0]; i < near.max[0]; i += global::dim2)
      for (j = near.min[1]; j < near.max[1]; j += global::dim)
	for (k = near.min[2]; k < near.max[2]; k++)
	  lRequest[i + j + k] = proc;
#ifdef GRAV
    tree.findGNear(part, lRequest, Vector<int>(0, 0, 0), 0, proc);
#endif
  }

  template <class T>
  int pSwapTLists(TList<T> &from, TList<T> &to, const int &with) {
    int recv, size = from.getSize();
    MCW.Sendrecv(&size, 1, MPI::INT, with, 2*with+1,
                 &recv, 1, MPI::INT, with, 2*global::rank+1);
    to.setSize(recv);
    MCW.Sendrecv(&from[0], size*sizeof(T), MPI::BYTE, with, 2*with+1,
                 &to[0],   recv*sizeof(T), MPI::BYTE, with, 2*global::rank+1);
    return recv;
  }

public:
  CellList(TList<Particle> &_pList) : newtonCell(5000), pList(_pList) {}

  void init() {
    cellList       = new Cell[global::dim3];
    lRequest       = new int[global::dim3]; 
    gRequest       = new int[global::dim3];
    readGhostList  = new TList<GhostRead>[global::npro];
    writeGhostList = new TList<GhostWrite>[global::npro];
#ifdef GRAV
    tree.init();
#endif
  }

  void pFindNeighbourCells() {
    DEBUG("pFindNeighbourCells()", "");
    int         c, cellNr, p;
    Particle   *part;
#ifdef GRAV
    tree.clear();
#endif
    for (c = 0; c < global::dim3; c++) { 
      cellList[c].clear();
      lRequest[c] = gRequest[c] = 0; 
    }
    
    for (p = 0; p < pList.getSize(); p++) {
      part   = &pList[p];
      cellNr = calcIndex(part->getCoord());
      cellList[cellNr].addParticle(part, p, global::rank);
#ifdef GRAV      
      tree.addParticle(part);
#endif
    }
#ifdef GRAV
    tree.pBuild();
#endif

    for (p = 0; p < pList.getSize(); p++) calcRequest(&pList[p]);
    MCW.Allreduce(lRequest, gRequest, global::dim3, MPI::INT, MPI::BOR);
  }

  void pSpreadGhosts() {
    DEBUG("pSpreadGhosts()", "");
    if (global::npro == 1) return;

    int       bit, c, n, p, with;
    Link      l;
    GhostRead ghost, *tmp;

    for (n = 0; n < global::comm.getRounds(); n++) {
      with = global::comm.get(n); 
      if (with != global::rank) {
	tmpGhostRList.setSize(0);
	for (bit = 1 << with, c = 0; c < global::dim3; c++) 
	  if (gRequest[c] & bit)
	    for (l = cellList[c].first; l.nr > -1; )
	      if (l.proc == global::rank) {
		ghost = pList[l.nr];
		ghost.setNr(l.nr);
		ghost.setCell(c);
		tmpGhostRList.safeAppend(ghost);
		l = ghost.getNext();
	      } else l = readGhostList[l.proc][l.nr].getNext();
	
	pSwapTLists(tmpGhostRList, readGhostList[with], with);
        writeGhostList[with].setSize(readGhostList[with].getSize());

        for (p = 0; p < readGhostList[with].getSize(); p++) {
          tmp = &readGhostList[with][p];
          writeGhostList[with][p].setZero();
          writeGhostList[with][p].setNr(tmp->getNr());
	  cellList[tmp->getCell()].addParticle(tmp, p, with);
        }
      }
    }
  }

  void pCollectGhosts() {
    DEBUG("pCollectGhosts()", "");
    GhostWrite *ghost;
    int         c, n, with;
    Particle   *part;

    for (n = 0; n < global::comm.getRounds(); n++) {
      with = global::comm.get(n);
      if (with != global::rank) {
        pSwapTLists(writeGhostList[with], tmpGhostWList, with);

        for (c = 0; c < tmpGhostWList.getSize(); c++) {
          ghost = &tmpGhostWList[c];
          part  = &pList[ghost->getNr()];
          ghost->update(part);
        }
      }
    }
  }

#ifdef SPH
  void forcesSPH(Particle *part) {
    int        i, j, k;
    Link       l;
    Range<int> near = part->getNear();

    near.min[0] *= global::dim2; near.max[0] *= global::dim2;   
    near.min[1] *= global::dim;  near.max[1] *= global::dim;   

    for (i = near.min[0]; i < near.max[0]; i += global::dim2)
      for (j = near.min[1]; j < near.max[1]; j += global::dim)
	for (k = near.min[2]; k < near.max[2]; k++)
	  for (l = cellList[i + j + k].first; l.nr > -1;)
	    if (l.proc == global::rank) {
	      part->forcesSPH(&pList[l.nr], &pList[l.nr]);
	      l = pList[l.nr].getNext();
	    } else {
	      part->forcesSPH(&readGhostList[l.proc][l.nr], 
			      &writeGhostList[l.proc][l.nr]);
	      l = readGhostList[l.proc][l.nr].getNext();
	    }
  }
#endif
#ifdef CORR
  void corrtensor(Particle *part) {
    int        i, j, k;
    Link       l;
    Range<int> near = part->getNear();

    near.min[0] *= global::dim2; near.max[0] *= global::dim2;
    near.min[1] *= global::dim;  near.max[1] *= global::dim;

    for (i = near.min[0]; i < near.max[0]; i += global::dim2)
      for (j = near.min[1]; j < near.max[1]; j += global::dim)
        for (k = near.min[2]; k < near.max[2]; k++)
          for (l = cellList[i + j + k].first; l.nr > -1;)
            if (l.proc == global::rank) {
              part->corrtensor(&pList[l.nr], &pList[l.nr]);
              l = pList[l.nr].getNext();
            } else {
              part->corrtensor(&readGhostList[l.proc][l.nr],
                              &writeGhostList[l.proc][l.nr]);
              l = readGhostList[l.proc][l.nr].getNext();
            }
  }
#endif
#ifdef GRAV
  void forcesGrav(Particle *part) {
    int  i, sum = tree.forcesGrav(part, newtonCell);
    Link l;

    for (i = 0; i < newtonCell.getSize(); i++)
      for (l = cellList[newtonCell[i]].first; l.nr > -1; ++sum)
	if (l.proc == global::rank) {
	  part->forcesNewton(&pList[l.nr]);
	  l = pList[l.nr].getNext();
	} else {
	  part->forcesNewton(&readGhostList[l.proc][l.nr]); 
	  l = readGhostList[l.proc][l.nr].getNext();
	}
    //std::cout << global::rank << ": " << sum << std::endl;
  }
#endif

/*
  int findNeigh(Particle *part, Array<Link> *list, 
		 Array<GhostRead> *in, Array<GhostWrite> *out) {
    int        c = 0, i, j, k;
    Link       l;
    Range<int> near = part->getNear();

    near.min[0] *= global::dim2; near.max[0] *= global::dim2;   
    near.min[1] *= global::dim;  near.max[1] *= global::dim;   

    for (i = near.min[0]; i < near.max[0]; i += global::dim2)
      for (j = near.min[1]; j < near.max[1]; j += global::dim)
	for (k = near.min[2]; k < near.max[2]; k++) {
	  for (l = cellList[i + j + k].first; l.nr > -1;) {
	    if (l.proc == global::rank) {
	      if (part->isNeighbour(&pList[l.nr]))
		{ (*list)[c] = l; (*in)[c++] = pList[l.nr]; }
	      l = pList[l.nr].getNext();
	    } else {
	      if (part->isNeighbour(&readGhostList[l.proc][l.nr]))
		{ (*list)[c] = l; (*in)[c++] = readGhostList[l.proc][l.nr]; }
	      l = readGhostList[l.proc][l.nr].getNext();
	    }
	  }
	}
    return c;
  }

  void updateNeigh(Particle *part, Array<Link> *list, Array<GhostRead> *in, 
		  Array<GhostWrite> *out, const int &neigh) {
    int  c;
    Link l;
    for (c = 0, l = (*list)[c]; c < neigh; l = (*list)[c++])
      if (l.proc == global::rank) (*out)[c].update(&pList[l.nr]);
      else (*out)[c].update(&writeGhostList[l.proc][l.nr]);
  }
*/
};

#endif
