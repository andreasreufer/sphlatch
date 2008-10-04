/******************************************************************************
 * ParaSPH -- Version 19.01.2004                                              *
 *----------------------------------------------------------------------------*
 * File:      Ghost.cc                                                        *
 * Purpose:   Sending ghosts between processes avoids sending a lot of data   *
 *            which is not really needed, but it makes your life somewhat     *
 *            more complicated. If you desire to change something, then its   *
 *            most probably because you would like to add a variable. Here is *
 *            what you should do:                                             *
 * 1) Do your changes in class 'Particle'.                                    *
 * 2) Make changes in class 'GhostRead' if and only if you need to read from  *
 *    your variable during 'forcesSPH'.                                       *
 * 2a) Add the variable to the variable block (between appropriate compiler   *
 *     switches).                                                             *
 * 2b) Add an assign operation in 'operation=' (or your ghost variable will   *
 *     never get initialized).   
 * 3) Make changes in class 'GhostWrite' if and only if you write to the new  *
 *    variable during 'forcesSPH'.                                            *
 * 3a) Add the variable to the variable block (between appropriate compiler   *
 *     switches).                                                             *
 * 3b) Add an addition operation (+=) in 'update' (or the source particle     *
 *     will never know about changes done by the alien process).              *
 *****************************************************************************/
#ifndef GHOST_CC
#define GHOST_CC

#include "Particle.cc"
#include "Vector.cc"

class GhostRead {
public:
  Vector<ftype> pos, v;
  ftype         h, id, mass;
  int           cell, nr;
  Link          next;
#ifdef SPH
  double        rankH;
  ftype         mat, p, rho, rhom1, u, vsound;
#ifdef SOLID
  ftype         reduce, sxx, sxy, sxz, syy, syz;
  bool          str;
#endif
#endif

public:
  GhostRead() { cell = -1; id = -1.; next.clear(); nr = -1; }

  const int           &getCell() { return cell; }
  const ftype         &getId()   { return id; }
  const ftype         &getMass()   const { return mass; }
  const Link          &getNext() { return next; }
  const int           &getNr()   { return nr; }
  const Vector<ftype> &getPos()    const { return pos; }

  void setCell(const int  &_cell) { cell = _cell; }
  void setNext(const Link &_next) { next = _next; }
  void setNr  (const int  &_nr)   { nr   = _nr; }

#ifdef SPH
  int  matnr()  const { return ((int)mat & 31); }
  int  bodyNr() const { return ((int)mat & 768) >> 8; }
  void divrho()   { rhom1 = 1. / rho; }
  void setRankH() { rankH = h*Particle::hMinm1 + id*Particle::totNumPartm1; }
#ifdef SOLID
  bool matstr() const { return str; }
#endif
#endif

  void operator=(const Particle &o) {
    pos    = o.pos;    v   = o.v;   h   = o.h;   id   = o.id;  mass = o.mass;
    next   = o.next; 
#ifdef SPH
    mat    = o.mat;    p   = o.p;   rho = o.rho; u    = o.u; 
    vsound = o.vsound;
    divrho(); setRankH();
#ifdef SOLID
    reduce = o.reduce; str = o.str; 
    sxx    = o.sxx;    sxy = o.sxy; sxz = o.sxz; syy  = o.syy; syz  = o.syz; 
#endif
#endif
  }
};


class GhostWrite {
public:
  Vector<ftype> f;
  ftype         work;
  int           nr;
#ifdef SPH
  Vector<ftype> deltav, tif;
  ftype         divv, du, hc, xmumax;
  int           neib;
#ifdef SOLID
  ftype         epsxx, epsxy, epsxz, epsyy, epsyz, epszz, rxy, rxz, ryz;
#endif
#endif

public:
  void setZero() {
    f      = Vector<ftype>(0., 0., 0.); work = 0.;
#ifdef SPH
    deltav = tif = Vector<ftype>(0., 0., 0.);
    divv   = du = hc = xmumax = 0.; neib  = 0;
#ifdef SOLID
    epsxx  = epsxy = epsxz = epsyy = epsyz = epszz = rxy = rxz = ryz = 0.;
#endif
#endif
  }

  const int &getNr() { return nr; }
  void setNr(const int &_nr) { nr = _nr; }

  template <class T>
  void update(T *o) {
    o->f      += f;      o->work += work;
#ifdef SPH    
    o->deltav += deltav; o->divv += divv; o->du += du; o->hc += hc; 
    o->neib   += neib;
    o->xmumax  = Max(o->xmumax, xmumax);
#ifdef SOLID
    o->epsxx  += epsxx;  o->epsxy += epsxy; o->epsxz += epsxz; 
    o->epsyy  += epsyy;  o->epsyz += epsyz; o->epszz += epszz; 
    o->rxy    += rxy;    o->rxz   += rxz;   o->ryz   += ryz;
#endif
#endif
#ifdef TENSILE
    o->tif += tif;
#endif
  }
};

#endif
