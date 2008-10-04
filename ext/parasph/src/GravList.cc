#ifndef GRAVLIST_CC
#define GRAVLIST_CC

#include <iomanip>

#include "Def.cc"
#include "Array.cc"
#include "Particle.cc"
#include "TList.cc"
#include "Vector.cc"

/*
class GravCell {
  friend class GravList;
private:
  ftype         Qxx, Qxy, Qxz, Qyy, Qyz, Qzz, mass, R;
  int           cellNr, size;
  Vector<ftype> pos;

public:
  void clear() { 
    Qxx = Qxy = Qxz = Qyy = Qyz = Qzz = mass = R = 0.; 
    pos.set(0., 0., 0.); 
    size = 0;
  }

  void addQ(Particle *part) {
    ftype         Mass = part->getMass();
    Vector<ftype> Pos  = part->getPos();
    if (mass == 0.) { pos = Pos; mass = Mass; size = 1; }
    else {
      Vector<ftype> d   = pos  - Pos;
      ftype         sum = mass + Mass, summ1 = 1./sum,
	            mu  = mass * Mass * summ1, zij = d.len() * summ1;
      R     = Max(mass * zij, Mass * zij + R);
      Qxx  += mu * d[0] * d[0]; Qxy  += mu * d[0] * d[1];
      Qxz  += mu * d[0] * d[2]; Qyy  += mu * d[1] * d[1];
      Qyz  += mu * d[1] * d[2]; Qzz  += mu * d[2] * d[2];
      pos   = (pos * mass + Pos * Mass) * summ1;
      mass  = sum;
      size++;
    }
  }

  void addQ(GravCell *other) {
    if (mass == 0.) {
      mass = other->mass; pos = other->pos; R   = other->R;
      Qxx  = other->Qxx;  Qxy = other->Qxy; Qxz = other->Qxz;
      Qyy  = other->Qyy;  Qyz = other->Qyz; Qzz = other->Qzz;
      size = other->size;
    } else {
      if (other->mass > 0.) {
        Vector<ftype> d   = pos - other->pos;
        ftype         M2  = other->mass, M3 = mass + M2, M3m1 = 1./M3,
                      mu  = mass * M2 * M3m1, zij = d.len() * M3m1;
        R     = Max(mass * zij + other->R, M2 * zij + R);
        Qxx  += other->Qxx + mu * d[0] * d[0];
        Qxy  += other->Qxy + mu * d[0] * d[1];
        Qxz  += other->Qxz + mu * d[0] * d[2];
        Qyy  += other->Qyy + mu * d[1] * d[1];
        Qyz  += other->Qyz + mu * d[1] * d[2];
        Qzz  += other->Qzz + mu * d[2] * d[2];
        pos   = (pos * mass + other->pos * M2) * M3m1;
        mass  = M3;
	size += other->size;
      }
    }
  }

  const ftype         &getMass() const { return mass; }
  const Vector<ftype> &getPos()  const { return pos; }
  const ftype         &getR()    const { return R; }
 
  Vector<ftype> Qtimesr(const Vector<ftype> &r) const {
    return Vector<ftype>(Qxx * r[0] + Qxy * r[1] + Qxz * r[2],
                         Qxy * r[0] + Qyy * r[1] + Qyz * r[2],
                         Qxz * r[0] + Qyz * r[1] + Qzz * r[2]);
  }
 
  ftype traceQ() const { return Qxx + Qyy + Qzz; }
 
  bool tooCloseGrav(const Vector<ftype> &partPos) const {
    return (R / (pos - partPos).len()) >= global::theta;
  }
};
*/
/*
class GravList {
private:
  Array<GravCell>  gravCell;
  int              depth, freeCell, *pre, *size, *size2, *tree;
  TList<Particle> *pList;

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
  void init(TList<Particle> *_pList) {
    pList = _pList;
    int d;
    depth = (int)ceil(log(global::dim)/log(2)) + 1;

    size = new int[depth]; size[depth-1]  = global::dim;
    for (d = depth-2; d >= 0; --d) size[d] = (size[d+1] + 1) >> 1; 
    size2 = new int[depth]; size2[depth-1] = global::dim * global::dim;
    for (d = 0; d < depth; ++d) size2[d] = size[d]*size[d];
    pre = new int[depth+1]; pre[0] = 0;
    for (d = 1; d <= depth; ++d) pre[d] = pre[d-1] + size[d-1]*size2[d-1];
    
    tree = new int[pre[depth]];
    // Warning: the following size could be too small -> investigate :-)
    gravCell.setSize(global::dim3);
  }

  void clear() { 
    freeCell = 0; 
    for (int i = 0; i < pre[depth]; ++i) tree[i] = -1;
    for (int i = 0; i < global::dim3; ++i) gravCell[i].clear();
  }

  void addParticle(const int &cellNr, Particle *part) {
    if (tree[cellNr] == -1) {
      tree[cellNr] = freeCell++;
      gravCell[tree[cellNr]].cellNr = cellNr;
      gravCell[tree[cellNr]].coord  = part->getCoord();
    }
    gravCell[tree[cellNr]].addParticle(part); 
  }

  void buildLocalTree() {
    for (int c = 0; c < freeCell; c++) {
      if (gravCell[i]
    }
  }

  void addCell(int x, int y, int z, GravCell *gCell) {
    int cellNr, d;
    for (d = depth-1; d >= 0; --d) {
      cellNr = pre[d] + x*size2[d] + y*size[d] + z;
      x >>= 1; y >>= 1; z >>= 1;
      if (tree[cellNr] == -1) {
	tree[cellNr] = freeCell++;
	gravCell[tree[cellNr]].cellNr = cellNr;
      }
      gravCell[tree[cellNr]].addQ(gCell); 
    }
  }

  void pExchange() {
    int cellNr, dest = freeCell, n, with;

    for (n = 0; n < global::comm.getRounds(); n++) {
      with = global::comm.get(n); 
      if (with != global::rank)
	dest += pSwapArrays(&gravCell[0], freeCell, &gravCell[dest], with);
    }
    
    for (; freeCell < dest; freeCell++) {
      cellNr = gravCell[freeCell].cellNr;
      if (tree[cellNr] == -1) tree[cellNr] = freeCell;
      else gravCell[tree[cellNr]].addQ(&gravCell[freeCell]);
    }    
  }

  void forces(Particle *part, int cellNr) {
    std::cout << "cellNr = " << cellNr << " | " << "tree[] = " << tree[cellNr]
	      << "   " << pre[depth-1] << std::endl;
    
    if (cellNr >= pre[depth]) {
      std::cout << "cellNr out of range " << cellNr << " " << pre[depth]
		<< std::endl;
      return;
    }
    if (tree[cellNr] == -1) return;
    if (gravCell[tree[cellNr]].tooCloseGrav(part->getPos())) {
      if (tree[cellNr] >= pre[depth-1]) {
	std::cout << "particle information needed" << std::endl;
	return;
      }
      cellNr <<= 3;
      forces(part, ++cellNr); 
      //forces(part, cellNr++);
      //forces(part, cellNr++); forces(part, cellNr++);
      //forces(part, cellNr++); forces(part, cellNr++);
      //forces(part, cellNr++); forces(part, cellNr);
    } else part->forcesQuad(&gravCell[tree[cellNr]]);
  }
};
*/
#endif
