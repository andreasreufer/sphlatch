#ifndef TREE_CC
#define TREE_CC

#include "Def.cc"
#include "Cell.cc"
#include "Particle.cc"

class Tree {
private:
  Array<GravCell> gravList;
  int             depth, freeCell, *level, *size, *size2, totSize, *tree;

  template <class T>
  int pSwapArrays(T *from, const int size, T *to, const int &with) {
    int recv;
    MCW.Sendrecv(&size, 1, MPI::INT, with, 2*with+1,
                 &recv, 1, MPI::INT, with, 2*global::rank+1);
    MCW.Sendrecv(from, size*sizeof(T), MPI::BYTE, with, 2*with+1,
                 to,   recv*sizeof(T), MPI::BYTE, with, 2*global::rank+1);
    return recv;
  }

public:
  void init() {
    int d;
    depth = (int)ceil(log(global::dim)/log(2)) + 1;

    size = new int[depth];
    for (d = depth-2, size[depth-1] = (global::dim+1)>>1; d >= 0; --d) {
      size[d]     = (size[d+1] + 1) >> 1; 
      size[d+1] <<= 1;
    }

    size2 = new int[depth]; 
    for (d = 0, totSize = 0; d < depth; ++d) {
      size2[d] = size[d]*size[d];
      totSize += size[d]*size2[d];
    }

    tree = new int[totSize];
    level = new int[depth];
    for (d = 1, level[0] = 0; d < depth; ++d) 
      level[d] = level[d-1] + size[d-1]*size2[d-1];
    
    if (global::rank == 0) {
      std::cout << "Tree::depth    = " << depth << std::endl;
      std::cout << "Tree::totSize  = " << totSize << std::endl;
      for (d = 0; d < depth; ++d) 
	std::cout << "Tree::size[" << d << "]  = " << size[d] << "   size2[" 
		  << d << "] = " << size2[d] << std::endl;
      for (d = 0; d < depth; ++d) 
	std::cout << "Tree::level[" << d << "]  = " << level[d] << std::endl;
      std::cout << "global::dim = " << global::dim << std::endl;
    }
  }

  void clear() { 
    freeCell = 0;
    for (int i = 0; i < totSize; ++i) tree[i] = -1;
    gravList.setSize(global::dim3);
  }

  int findTreePos(const Vector<int> &c, const int &d) {
    return level[d] + c[0]*size2[d] + c[1]*size[d] + c[2];
  }

  GravCell &getCell(const Vector<int> &c, const int &d) {
    int cellNr = findTreePos(c, d);
    if (tree[cellNr] == -1) {
      tree[cellNr] = freeCell++; assert(freeCell <= global::dim3);
      gravList[tree[cellNr]].clear();
      gravList[tree[cellNr]].cellNr = cellNr;
      gravList[tree[cellNr]].coord  = c;
    }
    return gravList[tree[cellNr]];
  }

  void addParticle(Particle *part) {
    getCell(part->getCoord(), depth-1).addParticle(part);
  }

  void addCellToFather(GravCell *cell, const int &d) {
    getCell(cell->coord >> 1, d-1).addCell(cell);
  }

  void pBuild() {
    DEBUG("Tree::pBuild()", "si", "gravList.getSize() = ", gravList.getSize());
    for (int cell = 0, d = depth-1; d > 0; --d) 
      for (int end = freeCell; cell < end; cell++)
	addCellToFather(&gravList[cell], d);

    pExchange();
  }

  void pExchange() {
    DEBUG("Tree::pExchange()", "");
    int cell = freeCell, cellNr, dest = freeCell, n, tot, with;

    MCW.Allreduce(&freeCell, &tot, 1, MPI::INT, MPI::SUM);
    gravList.growTo(tot); gravList.setSize(tot);

    for (n = 0; n < global::comm.getRounds(); n++) {
      with = global::comm.get(n); 
      if (with != global::rank)
	dest += pSwapArrays(&gravList[0], freeCell, &gravList[dest], with);
    }    
    for (; cell < dest; cell++) {
      cellNr = gravList[cell].cellNr;
      if (tree[cellNr] == -1) tree[cellNr] = cell;
      else gravList[tree[cellNr]].addCell(&gravList[cell]);
    }    
  }

  void findGNear(Particle *part, int *lRequest, Vector<int> c, int d, 
		 const int &proc) {
    int cellNr = findTreePos(c, d);
    if (tree[cellNr] == -1) return;
    if (gravList[tree[cellNr]].tooCloseGrav(part->getPos()))
      if (d == depth-1) 
	lRequest[c[0]*global::dim2 + c[1]*global::dim + c[2]] = proc;
      else {
	c = c << 1;
	++d;    findGNear(part, lRequest, c, d, proc); 
	c[0]++; findGNear(part, lRequest, c, d, proc);
	c[1]++; findGNear(part, lRequest, c, d, proc);
	c[0]--; findGNear(part, lRequest, c, d, proc);
	c[2]++; findGNear(part, lRequest, c, d, proc);
	c[0]++; findGNear(part, lRequest, c, d, proc);
	c[1]--; findGNear(part, lRequest, c, d, proc);
	c[0]--; findGNear(part, lRequest, c, d, proc);
      }
  }

/*
  int forces(Particle *part, Array<int> &newtonCell, Vector<int> c, int d) {
    int cellNr = findTreePos(c, d);
    if (tree[cellNr] == -1) return 0;
    if (gravList[tree[cellNr]].tooCloseGrav(part->getPos())) {
      if (d == depth-1) {
	newtonCell.append(c[0]*global::dim2 + c[1]*global::dim + c[2]);
	return 0;
      } else {
	Vector<int> fa = c << 1; int sum;
	sum = forces(part, newtonCell, fa, d+1); 
	fa[0]++; sum += forces(part, newtonCell, fa, d+1);
	fa[1]++; sum += forces(part, newtonCell, fa, d+1);
	fa[0]--; sum += forces(part, newtonCell, fa, d+1);
	fa[2]++; sum += forces(part, newtonCell, fa, d+1);
	fa[0]++; sum += forces(part, newtonCell, fa, d+1);
	fa[1]--; sum += forces(part, newtonCell, fa, d+1);
	fa[0]--; sum += forces(part, newtonCell, fa, d+1);
	return sum;
      }
    } else {
      part->forcesQuad(&gravList[tree[cellNr]]);
      return gravList[tree[cellNr]].size;
    }
  }

  int forcesGrav(Particle *part, Array<int> &newtonCell) {
    newtonCell.setSize(0);
    return forces(part, newtonCell, Vector<int>(0, 0, 0), 0); 
  }
*/

  void forces(Particle *part, Array<int> &newtonCell, Vector<int> c, int d,
	      int &sum) {
    int cellNr = findTreePos(c, d);
    if (tree[cellNr] == -1) return;
    if (gravList[tree[cellNr]].tooCloseGrav(part->getPos())) {
      if (d == depth-1)
	newtonCell.append(c[0]*global::dim2 + c[1]*global::dim + c[2]);
      else {
	c = c << 1;
	++d;    forces(part, newtonCell, c, d, sum); 
	c[0]++; forces(part, newtonCell, c, d, sum);
	c[1]++; forces(part, newtonCell, c, d, sum);
	c[0]--; forces(part, newtonCell, c, d, sum);
	c[2]++; forces(part, newtonCell, c, d, sum);
	c[0]++; forces(part, newtonCell, c, d, sum);
	c[1]--; forces(part, newtonCell, c, d, sum);
	c[0]--; forces(part, newtonCell, c, d, sum);
      }
    } else {
      part->forcesQuad(&gravList[tree[cellNr]]);
      sum += gravList[tree[cellNr]].size;
    }
  }

  int forcesGrav(Particle *part, Array<int> &newtonCell) {
    int sum = 0;
    newtonCell.setSize(0);
    forces(part, newtonCell, Vector<int>(0, 0, 0), 0, sum); 
    return sum;
  }

/*  int forcesGrav2(Particle *part, Array<int> &newtonCell) {
    int         cellNr, d, sum = 0;
    Vector<int> c, fa;
    
    newtonCell.setSize(0);

    for (stack.push(Vector<int>(0, 0, 0), 0); !stack.empty; ) {
      stack.top(c, d); stack.pop();

      cellNr = findTreePos(c, d);
      if (tree[cellNr] > -1) 
	if (gravList[tree[cellNr]].tooCloseGrav(part->getPos())) {
	  if (d == depth-1)
	    newtonCell.append(c[0]*global::dim2 + c[1]*global::dim + c[2]);
	  else {
	    ++d; fa = c << 1; 
	    stack.push(fa, d);	        fa[0]++; stack.push(fa, d);
	    fa[1]++; stack.push(fa, d); fa[0]--; stack.push(fa, d);
	    fa[2]++; stack.push(fa, d); fa[0]++; stack.push(fa, d);
	    fa[1]--; stack.push(fa, d); fa[0]--; stack.push(fa, d);
	  }
	} else {
	  part->forcesQuad(&gravList[tree[cellNr]]);
	  sum += gravList[tree[cellNr]].size;
	}
    }
    return sum;
  }*/
};

#endif
