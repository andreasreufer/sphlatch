#ifndef CELL_CC
#define CELL_CC

#include "Def.cc"
#include "Vector.cc"

class Cell {
  friend class CellList;
private:
  Link first;

public:
  Cell() { clear(); }

  void clear() { first.clear(); }

  template <class T>
  void addParticle(T *part, const int &nr, const int &proc) { 
    first.addParticle(part, nr, proc);
  }
};

class GravCell {
  friend class Tree;
private:
  ftype         Qxx, Qxy, Qxz, Qyy, Qyz, Qzz, mass, R;
  int           cellNr, size;
  Vector<ftype> pos;
  Vector<int>   coord;

public:
  GravCell() { clear(); }

  void clear() { 
    Qxx = Qxy = Qxz = Qyy = Qyz = Qzz = mass = R = 0.; 
    pos.set(0., 0., 0.); coord.set(-1, -1, -1); 
    size = 0;
  }

  const ftype         &getMass() const { return mass; }
  const Vector<ftype> &getPos()  const { return pos; }
  const ftype         &getR()    const { return R; }

  template <class T>
  void addParticle(T *part) { 
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

  void addCell(GravCell *other) {
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

  Vector<ftype> Qtimesr(const Vector<ftype> &r) const {
    return Vector<ftype>(Qxx * r[0] + Qxy * r[1] + Qxz * r[2],
                         Qxy * r[0] + Qyy * r[1] + Qyz * r[2],
                         Qxz * r[0] + Qyz * r[1] + Qzz * r[2]);
  }

  ftype traceQ() const { return Qxx + Qyy + Qzz; }
 
  bool tooCloseGrav(const Vector<ftype> &partPos) const {
    ftype dist = (pos - partPos).len();
    return (dist == 0. || R / dist >= global::theta);
  }
};

#endif
