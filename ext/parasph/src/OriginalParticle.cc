#ifndef PARTICLE_CC
#define PARTICLE_CC

#include "Eos.cc"
#include "Kernel.cc"
#include "Manager.cc"
#include "Vector.cc"
#include "Xdr.cc"

class Particle {
  friend class Processor;
  friend class Xdr;

  friend class Cracks;
  friend class GhostRead;
  friend class GhostWrite;
  friend class Init;

private:
  // Variables to store
  Vector<ftype>   pos, v;
  ftype           h, id, mass;
#ifdef SPH
  ftype           mat, p, T, rho, u;
#ifdef SOLID
  ftype           dm, sLim, sxx, sxy, sxz, syy, syz;
#endif
#endif

  // Temporarily used variables
  float           work;
  int             proc;
  Link            next;         
  Vector<int>     coord;
  Range<int>      near;  

  // temp. SPH vars
  Vector<ftype>   deltav, f, sf;

#ifdef SPH
  double          rankH;
  ftype           divv, hc, rhom1, vsound, xmumax;
  ftype           dh, drho, du, sdh, sdrho, sdu;
  int             neib;
  Vector<ftype>   tif;

  static Eos      eos;
  static ftype    avAlpha, avBeta, dampFactor, dtFactor, xsphFactor;
  static ftype    heatconG1, heatconG2;
  static ftype    hMax, hMin, hMinm1, rhoMin, uMin;  

#ifdef SOLID
  bool            str;
  ftype           epsxx, epsyy, epszz, epsxy, epsxz, epsyz;
  ftype           reduce, rxy, rxz, ryz, vonmises;
  ftype           ddm,  dsxx,  dsxy,  dsxz,  dsyy,  dsyz;
  ftype           sddm, sdsxx, sdsyy, sdsxy, sdsxz, sdsyz;
  ftype           sig1, sig2, sig3;

  ftype           acoef, epsmin, grav, m, young;
  int             flaws;

  static Material<Mat> matTab;
  static ftype         dmMin, sMin, tiny;
#endif
#endif

  static ftype        G, gravHb1, gravHb12, gravHb1m1, gravHb1m2;
  static ftype           gravHb2, gravHb22, gravHb2m1, gravHb2m2;

  static double       totNumPartm1;
  static const ftype  third, sqthird, A, B, C;
  static const Kernel kern;

public:
  static bool         heatCon, xsph;
  static int          numVar;
  static std::string  varName;

#ifdef SPH
  class ErrorHTooBig {
  public:
    double h, hMax;
    ErrorHTooBig(const double &_h) { h = _h; hMax = Particle::hMax; }
  };
#endif

  Particle() { work = 1.; next.clear(); }

  static void init(Manager &man) {
    numVar   = 9;
    varName  = "x\ny\nz\nvx\nvy\nvz\nh\nid\nmass\n";
#ifdef SPH
    numVar  += 5;
    varName += "mat\np\nT\nrho\nu\n";

    int         eosNr;
    std::string eosDataFile, path;
    man.getValue("simulation.eos",    eosNr);
    man.getValue("simulation.path",   path);
    man.getValue("input.eosDataFile", eosDataFile);
    eosDataFile.insert(0, path + "/");
    eos.init(eosNr, global::noisy, eosDataFile);
    
    man.getValue("thres.hMin",   hMin);
    man.getValue("thres.rhoMin", rhoMin);
    man.getValue("thres.uMin",   uMin);

    man.getValue("SPH.avAlpha",    avAlpha);
    man.getValue("SPH.avBeta",     avBeta);
    man.getValue("SPH.dtFactor",   dtFactor);
    man.getValue("SPH.damping",    dampFactor);

    hMinm1 = 1. / hMin;
#ifdef SOLID
    numVar  += 7;
    varName += "dm\nsLim\nsxx\nsxy\nsxz\nsyy\nsyz\n";

    std::string matDataFile;
    man.getValue("input.matDataFile", matDataFile);
    matDataFile.insert(0, path + "/");
    matTab.read(global::noisy, matDataFile);

    man.getValue("thres.dmMin",  dmMin);
    man.getValue("thres.sMin",   sMin);
    man.getValue("thres.tiny",   tiny);
#endif
#endif
    man.getValue("physics.G", G);
    man.getValue("simulation.gravHb1", gravHb1);
    man.getValue("simulation.gravHb2", gravHb2);
 
#ifdef HEATCON
    man.getValue("simulation.heatCon", heatCon);
    man.getValue("SPH.heatconG1", heatconG1);
    man.getValue("SPH.heatconG2", heatconG2);
#endif
#ifdef XSPH
    man.getValue("simulation.XSPH",    xsph);
    man.getValue("SPH.xsphFactor", xsphFactor);
#endif
  }

  static void setGravH(const ftype &hMinGlob) {
    //gravH   *= hMinGlob;
    gravHb12   = gravHb1*gravHb1;
    gravHb1m1  = 1. / gravHb1;
    gravHb1m2  = gravHb1m1 * gravHb1m1;
    gravHb22   = gravHb2*gravHb2;
    gravHb2m1  = 1. / gravHb2;
    gravHb2m2  = gravHb2m1 * gravHb2m1;
  }

  static void calcHMin() {
    totNumPartm1 = 1. / (double)global::totNumPart;

#ifdef SPH
    // find the maximal h (for which the rankH of two particles still will
    // be different)
    bool   prec = true;
    double d1, d2;
    for (hMax = 1.; prec; hMax *= 10.) {
      d1   = hMax * hMinm1; d2 = hMax * hMinm1 + totNumPartm1;
      prec = (d2 > d1);
    }
#endif
  }

  // Get variables
  const Vector<int>   &getCoord()  const { return coord; }  
  const ftype         &getH()      const { return h; }
  const ftype         &getId()     const { return id; }
  const ftype         &getMass()   const { return mass; }
  const Range<int>    &getNear()   const { return near; }
  const Link          &getNext()   const { return next; }
  const Vector<ftype> &getPos()    const { return pos; }
  const int           &getProc()   const { return proc; }
  const float         &getWork()   const { return work; }
#ifdef SOLID
  const Mat           &getMatVar() const { return matTab.get(matNr()); }
#endif

  // Set variables
  void setCoord (const int &c, const int &x)    { coord[c] = x; }  
  void setNear  (const Range<int>    &_near)    { near     = _near; }
  void setNext  (const Link          &_next)    { next     = _next; }
  void setProc  (const int           &_proc)    { proc     = _proc; }
  void addToProc(const int           &_proc)    { proc    += _proc; }
  // The work is set to _work, only if it is 1. (default value),
  // otherwise it is changed slowly
  void setWork  (const float         &_work)    { 
    //if (work == 1.) work = _work; 
    //else {
     ftype ratio = _work/work; 
     ratio = Min(ratio,(float) 2.); ratio = Max(ratio, float (0.5));
     work = ratio * work;
    //}
  }

  // Read/write variables from/to a xdr-file
  // Changes: the first numVar variables of *this
  void read (Xdr &file) { file.read((ftype *)this); }
  void write(Xdr &file) { file.write((ftype *)this); }

  /***************************************************************************
   * Preforces                                                               *
   *-------------------------------------------------------------------------*
   * Insert here everything that has to be calculated before ghost exchange. *
   * Do not forget to call your new method within preForces().               *
   **************************************************************************/
#ifdef SPH
  int  matNr()  const { return ((int)mat & 31); }

  int  bodyNr() const { return ((int)mat & 3840) >> 8; }

  void divrho() { rhom1 = 1. / rho; }

  // This is a little trick: forces between two particles should only be
  // be computed once. So, just compute them if getRankH_i is larger than
  // getRankH_j, ie. if (h_i > h_j), or in case that (h_i == h_j) just
  // if (id_i > id_j).
  void setRankH() { rankH = h * hMinm1 + id * totNumPartm1; }

// equation of state
  void calleos() {
    eos.eos(rho, rhom1, u, matNr(), p, vsound, T);
  }
#ifdef SOLID
  // Von Mises criterion to set plastic yielding
  void plastic() {
    Mat   mat   = matTab.get(matNr());
    //ftype yst = mat.yield;
    ftype ratio = u / mat.umelt, yst;
    if (ratio < tiny) yst = mat.yield;
    else              yst = mat.yield * Max((1. - ratio) , 0.);
    ftype szz = -sxx - syy;
    ftype sqJ2 = sqrt(.5*(sxx*sxx+syy*syy+szz*szz)+sxy*sxy+sxz*sxz+syz*syz) + tiny;
    sLim = Min(sqthird*yst/sqJ2 , (ftype) 1.);
  }

  // Limit stresses
  void stressLimiters() {
    // Von Mises Yield criterion
    plastic();
  }

  // Reduce stresses due to fracture
  void damage() {
    reduce = sLim*(1. - dm*dm*dm);
    // limit negative pressure
    if (p < 0.) {p *= (1.0-dm*dm*dm);}
  }

  // particle flags
  void flagp() {str = reduce != 0.0;}
  bool matstr() const { return str; }
  bool strength(Particle *other) const {
    return matstr() && other->matstr() && (bodyNr() == other->bodyNr());
  }

#endif
#endif
  void setZero() { 
    f      = Vector<ftype>(0., 0., 0.);
#ifdef SPH
    deltav = tif = Vector<ftype>(0., 0., 0.);
    divv   = drho = du = hc = xmumax = 0.; neib  = 0;
#ifdef SOLID
    dsxx   = dsxy = dsxz = dsyy = dsyz = 0.;
    epsxx  = epsxy = epsxz = epsyy = epsyz = epszz = rxy = rxz = ryz = 0.;
#endif
#endif
  }

  void preForces() {
#ifdef SPH
    divrho();
    calleos();
    setRankH();
#ifdef SOLID
    stressLimiters();
    damage();
    flagp();
#endif
#endif
    setZero();
  }  

  /***************************************************************************
   * Forces                                                                  *
   *-------------------------------------------------------------------------*
   * Remember: everything you add to the following methods will be executed  *
   *           approx. 60 * N times (with N being the number of particles).  *
   *           Especially try to avoid 'sqrt' and division operations :-)    *
   **************************************************************************/
  template <class T>
  bool isNeighbourOld(T *other) {
    if (rankH < other->rankH) return false;

    ftype hmean   = .5*(h + other->h);
    ftype hmeanm1 = 1. / hmean, hmeanm2 = hmeanm1 * hmeanm1;
       
    Vector<ftype> d    = pos - other->pos;
    ftype         rij2 = d.len2();
    ftype         v2   = rij2 * hmeanm2;
 
    if (v2 < kern.v2max) return true; else return false;
  }

  template <class T>
  void forcesSPHOld(T *other) {
    if (rankH < other->rankH) return;

    ftype hmean   = .5*(h + other->h);
    ftype hmeanm1 = 1. / hmean, hmeanm2 = hmeanm1 * hmeanm1;
       
    Vector<ftype> d    = pos - other->pos;
    ftype         rij2 = d.len2();
    ftype         v2   = rij2 * hmeanm2;
 
    if (v2 < kern.v2max) {
      mass += 1.; // do some calculations
    }
  }

#ifdef SPH
  template <class T>
  bool isNeighbour(T *other) {
    if (rankH < other->rankH) return false;
    ftype hmean = (h + other->h);
    if ((pos - other->pos).len2() < hmean*hmean) return true;
    return false;
  }

  template <class IN, class OUT>
  void forcesSPH(IN *read, OUT *write) {
    if (!isNeighbour(read)) return;

    neib++; write->neib++;

    ftype hmean   = .5*(h + read->h);
    ftype hmean2  = hmean*hmean;
    ftype hmeanm1 = 1. / hmean, hmeanm2 = hmeanm1 * hmeanm1;
    ftype hmeanm3 = hmeanm2 * hmeanm1, hmeanm5 = hmeanm3 * hmeanm2;
      
    Vector<ftype> d    = pos - read->pos;
    ftype         rij2 = d.len2();
    ftype         v2   = rij2 * hmeanm2;
    //ftype         v2 = d.len2() / (hmean*hmean);

    // One of the following two lines has to be a comment
    // First line: Pi/rhoi^2 + Pj/rhoj^2
    // 2nd line:   (Pi + Pj) / (rhoi * rhoj)
    ftype rhoim2 = rhom1 * rhom1, rhojm2 = read->rhom1 * read->rhom1;
    //**//ftype rhoim2 = rhom1 * other->rhom1, rhojm2 = rhoim2;
      
    ftype gW         = kern.getgW(v2) * hmeanm5;
    // direct kernel computation (requires a sqrt)
    // ftype gW      = kern.fgW(sqrt(v2)) * hmeanm5;
      
    ftype gWmi       = mass*gW;
    ftype gWmj       = read->mass*gW;
    ftype gWmri      = gWmi*rhom1;
    ftype gWmrj      = gWmj*read->rhom1;
    
    Vector<ftype> dv = v - read->v;

    ftype exx        = dv[0]*d[0];
    ftype eyy        = dv[1]*d[1];
    ftype ezz        = dv[2]*d[2];
    ftype projv      = exx + eyy + ezz;
    
    ftype rhobar     = .5*(rho + read->rho);
    
    // Artificial viscosity
    ftype Piij, Q2 = 0.;
    if (projv < 0) {
      ftype vsbar   = .5*(vsound + read->vsound);
      ftype muijn   = rij2 + 0.05*hmean2;
      Q2            = projv / muijn; 
      ftype ff      = Q2 * hmean;
      xmumax        = Max(Fabs(ff), xmumax);
      write->xmumax = Max(Fabs(ff), write->xmumax);
      Piij          = (avBeta*ff - avAlpha*vsbar)*ff/rhobar;
    } else Piij = 0.0;
    
    // Pressure
    ftype pij  = p * rhoim2 + read->p * rhojm2;
    ftype diag = pij + Piij;
    
    f         -= d     * diag * gWmj;
    write->f  += d     * diag * gWmi;
    du        += projv * diag * gWmj;
    write->du += projv * diag * gWmi;
    
    // Divergence velocity
    divv        -= projv*gWmrj;
    write->divv -= projv*gWmri;
    
#ifdef SOLID
    // Material strength (if particles belong to same object)
    ftype sigxx, sigxy, sigxz, sigyy, sigyz, sigzz;

    if (matstr() && read->matstr() && (bodyNr() == read->bodyNr())) {
      ftype ror2i = reduce       * rhoim2;
      ftype ror2j = read->reduce * rhojm2;
      
      sigxx = ror2i*sxx        + ror2j*read->sxx;
      sigyy = ror2i*syy        + ror2j*read->syy;
      sigzz = ror2i*(-sxx-syy) + ror2j*(-read->sxx-read->syy);
      sigxy = ror2i*sxy        + ror2j*read->sxy;
      sigxz = ror2i*sxz        + ror2j*read->sxz;
      sigyz = ror2i*syz        + ror2j*read->syz;
      
      Vector<ftype> t;
      t[0] = sigxx*d[0] + sigxy*d[1] + sigxz*d[2];
      t[1] = sigxy*d[0] + sigyy*d[1] + sigyz*d[2];
      t[2] = sigxz*d[0] + sigyz*d[1] + sigzz*d[2];
      f         += t   * gWmj;
      write->f  -= t   * gWmi;
      ftype tdv  = t   * dv;
      du        -= tdv * gWmj;
      write->du -= tdv * gWmi;
      
      // Strain rate tensor
      epsxx        -= gWmrj*exx;
      epsyy        -= gWmrj*eyy;
      epszz        -= gWmrj*ezz;
      write->epsxx -= gWmri*exx;
      write->epsyy -= gWmri*eyy;
      write->epszz -= gWmri*ezz;
      ftype dvxy    = dv[0]*d[1], dvyx = dv[1]*d[0];
      ftype dvxz    = dv[0]*d[2], dvzx = dv[2]*d[0];
      ftype dvyz    = dv[1]*d[2], dvzy = dv[2]*d[1];
      ftype exy     = .5*(dvxy + dvyx);
      ftype exz     = .5*(dvxz + dvzx);
      ftype eyz     = .5*(dvyz + dvzy);
      epsxy        -= gWmrj*exy;
      epsxz        -= gWmrj*exz;
      epsyz        -= gWmrj*eyz;
      write->epsxy -= gWmri*exy;
      write->epsxz -= gWmri*exz;
      write->epsyz -= gWmri*eyz;
      
      // Rotation rate tensor
      ftype rrxy     = .5*(dvxy - dvyx);
      ftype rrxz     = .5*(dvxz - dvzx);
      ftype rryz     = .5*(dvyz - dvzy);
      rxy        -= gWmrj*rrxy;
      rxz        -= gWmrj*rrxz;
      ryz        -= gWmrj*rryz;
      write->rxy -= gWmri*rrxy;
      write->rxz -= gWmri*rrxz;
      write->ryz -= gWmri*rryz;
    } else sigxx = sigxy = sigxz = sigyy = sigyz = sigzz = 0.0;
#endif	

    // The following three code blocks are extensions to standard SPH
    // You can omit them by just not compiling them (ie. without compiler
    // switches). The if tests in XSPH and HEATCON are needed only if you
    // want to switch it off during one simulation. If not just comment
    // the if tests. For further information read chapter 2.3 of my thesis. 
    
#ifdef XSPH
    // XSPH
    if (xsph) {
      ftype W        = kern.getW(v2) * hmeanm3;
      deltav        -= dv * (read->mass * W / rhobar);
      write->deltav += dv * (mass       * W / rhobar);
    }
#endif
    
#ifdef HEATCON
    // Heat conduction
    // I did not need the second order terms. They are expensive :-)
    if (heatCon) {
      ftype help = hmean * rhobar * rhom1 * read->rhom1 * (u - read->u);
      ftype Qij  = heatconG1 * help  * (vsound + read->vsound);
      ftype Qij2 = heatconG2 * help  * hmean * Q2 / sqrt(rij2 + 0.01*hmean2);
      hc        += Qij   * gWmj;
      write->hc -= Qij   * gWmi;
      hc        -= Qij2  * gWmj;
      write->hc += Qij2  * gWmi;
    }
#endif
    
#ifdef TENSILE
    // Tensile instability
    ftype ti_eps = pow(kern.getW(v2)/kern.fW(1.), 4.) * 1.;
    if (ti_eps > 0.) {  
      if (sigxx < pij) sigxx = 0.; else sigxx -= pij; 
      if (sigyy < pij) sigyy = 0.; else sigyy -= pij;
      if (sigzz < pij) sigzz = 0.; else sigzz -= pij;
      if (sigxy < 0.)  sigxy = 0.;
      if (sigxz < 0.)  sigxz = 0.;
      if (sigyz < 0.)  sigyz = 0.;
      
      Vector<ftype> R;
      R[0] = sigxx*d[0] + sigxy*d[1] + sigxz*d[2];
      R[1] = sigxy*d[0] + sigyy*d[1] + sigyz*d[2];
      R[2] = sigxz*d[0] + sigyz*d[1] + sigzz*d[2];
      tif        -= R * ti_eps * gWmj;
      write->tif += R * ti_eps * gWmi;
    }
#endif
  }
#endif

  template <class T>
  void forcesNewton(T *other) {
    ftype gravH2, gravHm1, gravHm2;
    if(bodyNr() == 0) gravH2=gravHb12, gravHm1=gravHb1m1, gravHm2=gravHb1m2;
    else              gravH2=gravHb22, gravHm1=gravHb2m1, gravHm2=gravHb2m2;
    Vector<ftype> r       = pos - other->pos;
    ftype         rij2    = r.len2() + 0.01*gravH2; 
    ftype         v2      = rij2 * gravHm2;
    ftype         gravij  = G * other->mass / rij2;
 
    if (v2 < kern.v2max) f -= r * gravij * kern.getFm(v2) * gravHm1;
    else                 f -= r * (gravij / sqrt(rij2));
  }
 
  // Gravitational forces (using quadrupole moments)
  // T has to be a cell
  template <class T>
  void forcesQuad(T *other) {
//    std::cout << "forcesQuad called" << std::endl;
    Vector<ftype> r = pos - other->getPos();
    Vector<ftype> t = other->Qtimesr(r);

    // This should not happen ...
    //if (r.len() == 0.) return;

    double rijm1 = 1. / (r.len() + 0.01*gravHb1);    //   1 /  r
    double rijm2 = rijm1 * rijm1;                   //   1 /  r^2
    double C1    = rijm2 * rijm1 * G;               //   G /  r^3
    double C2    = rijm2 * C1 * -1.5;               //  3G / 2r^5
    double C3    = rijm2 * C2 * -5.;                // 15G / 2r^7
    double fac   = C1 * other->getMass() + C2 * other->traceQ() + t * r * C3;
    C2 *= 2.;

    /*double r0 = r[0] * fac, r1 = r[1] * fac, r2 = r[2] * fac;
    double t0 = t[0] * C2,  t1 = t[1] * C2,  t2 = t[2] * C2;
    
    if (isnan(r0) || isnan(r1) || isnan(r2)) {
      std::cout << "nan error: " << r << " " 
		<< C1 << " " << other->getMass() << " "
		<< C2 << " " << other->traceQ() << " " << t
		//<< C3 << " " << t * r << " " 
		//<< fac 
		<< std::endl;
      quit(0);
    }

    f[0] -= r0 + t0 ; f[1] -= r1 + t1; f[2] -= r2 + t2;
    */
    f -= r * fac + t * C2;
  }

  /***************************************************************************
   * Postforces                                                              *
   *-------------------------------------------------------------------------*
   * Insert here everything that has to be calculated after ghost collect.   *
   **************************************************************************/
#ifdef SPH
  void density() { drho = -rho * divv; }
  
  void energy()  { du   = .5 * (du + hc); }
  
  void damping() { f   -= v * dampFactor; }
  
  void addTif()  { f += tif; }

  ftype getEnergy() { return mass * (u + .5*v.len2()); }
  
#ifdef SOLID
  void deviator() {
    if (bodyNr() == 0 && dm < 1.) {
      ftype mu   = matTab.get(matNr()).mu;
      ftype szz  = - sxx - syy;
      ftype div3 = (epsxx + epsyy + epszz) * third;
      dsxx = 2.*(mu*(epsxx - div3) + sxy*rxy + sxz*rxz);
      dsyy = 2.*(mu*(epsyy - div3) - sxy*rxy + syz*ryz);
      dsxy = 2.*mu*epsxy + (syy-sxx)*rxy + sxz*ryz + syz*rxz;
      dsxz = 2.*mu*epsxz + (szz-sxx)*rxz - sxy*ryz + syz*rxy;
      dsyz = 2.*mu*epsyz + (szz-syy)*ryz - sxy*rxz - sxz*rxy;
    } else {
      dsxx = 0.0; dsxy = 0.0; dsxz = 0.0; dsyy = 0.0; dsyz = 0.0;
    }
  }

  // compute principal stresses
  void pstresses() {
    ftype Sxx = reduce*sxx-p-grav; ftype Syy = reduce*syy-p-grav;
    ftype Sxy = reduce*sxy; ftype Sxz = reduce*sxz; ftype Syz = reduce*syz;
    ftype Szz = reduce*(-sxx-syy)-p-grav;
    ftype Smax = Max(Fabs(Sxy), Fabs(Syz), Fabs(Syz));

    sig1 = sig2 = sig3 = 0.;
    if (Smax == 0.) {
       sig1 = sig2 = sig3 = -p-grav;
#ifdef TENSILE
       pa1  = Vector<ftype>(1.,0.,0.);
       pa2  = Vector<ftype>(0.,1.,0.);
       pa3  = Vector<ftype>(0.,0.,1.);
#endif
       }
    else {
      ftype Smax = Max(Smax, Fabs(Sxx), Fabs(Syy), Fabs(Szz));
      Sxx /= Smax; Syy /= Smax; Szz /= Smax;
      Sxy /= Smax; Sxz /= Smax; Syz /= Smax;
      ftype pp = - Sxx - Syy - Szz,
      q  = Sxx*Syy + Sxx*Szz + Syy*Szz - Sxy*Sxy - Sxz*Sxz - Syz*Syz,
      r  = Sxx*Syz*Syz + Syy*Sxz*Sxz + Szz*Sxy*Sxy - Sxx*Syy*Szz
           - 2*Sxy*Sxz*Syz,
      a  = third * (3.*q - pp*pp),
      b  = (2*pp*pp*pp - 9.*pp*q + 27.*r) / 27.,
      a1 = a*a*a / 27.;
      if (0.25 * b * b + a1 < 0.) {
         ftype t1 = 2.*sqrt(-third * a),
         p3   = third * pp,
         phi  = acos(-0.5 * b / sqrt(-a1)),
         phi1 = third * phi,
         phi2 = third * (phi +  6.28318530718),
         phi3 = third * (phi + 12.56637061436);
         sig1 = (t1*cos(phi1) - p3);
         sig2 = (t1*cos(phi2) - p3);
         sig3 = (t1*cos(phi3) - p3);
         sig1 *= Smax, sig2 *= Smax, sig3 *= Smax;
      }
    }
  }

  // compute fracture growth - Damage limited to 1/Nflaws
  void fracture() {
    ddm = 0.;
    if (dm < 1.) {
      pstresses();
      ftype youngi = young*(1.-dm*dm*dm) + tiny ;
      ftype rft    = Max(sig1, sig2, sig3) / (epsmin * youngi);
      if (rft > 1. && epsmin < 1.e10){
         ftype aflaws = pow(rft,m);
         if (aflaws > flaws) {
            ddm = acoef * pow(flaws,third);}
         else{
            ftype dmmax= aflaws/flaws;
            if (dm*dm*dm < dmmax){
               ddm = acoef * pow(aflaws,third);}
            else{
               ddm = 0.;}
         }
      }
    }
  }

  // compute fracture growth - Damage limited to 1/Nflaws
  /*void fracture() {
    ddm = 0.;
    if (dm < 1.) {
      pstresses();
      ftype youngi = young*(1.-dm*dm*dm) + tiny ;
      ftype rft    = Max(sig1, sig2, sig3) / (epsmin * youngi);
      if (rft > 1. && epsmin < 1.e10){
         ftype xflaws = Min((double)pow(rft,m),(double) flaws);
         ddm = acoef * pow(xflaws,third);
      }
    }
  }*/

#endif
#endif

  void postForces() {
#ifdef SPH
    density();
    energy();
    damping();
    addTif();
#ifdef SOLID
    deviator();
    fracture();
#endif 
#endif
  }

#ifdef SPH
  void maxdivv(ftype *divmax) {
    divmax[bodyNr()] = Max(divmax[bodyNr()], Fabs(divv));
  }

  void hdot(ftype *divmax) {
    dh = h * divv * third;
    
    ftype dhs   = -divmax[bodyNr()] * h * third;
    int   dnsup = Max(120 - neib, -50);
    int   dninf = Max(neib -  10, -50);
    if (dnsup < 10) {
      ftype wsex1 = exp(dnsup * 0.2);
      ftype wsex2 = 1. / wsex1;
      dh = (wsex1 * dh + wsex2 * dhs)/(wsex1 + wsex2);
    } else if (dninf < 20) {
      ftype wiex1 = exp(dninf * 0.2);
      ftype wiex2 = 1. / wiex1;
      dh = (wiex1 * dh - wiex2 * dhs)/(wiex1 + wiex2);
    }
  }
#endif

  /***************************************************************************
   * Timestep conditions                                                     *
   *-------------------------------------------------------------------------*
   **************************************************************************/
  void tsconds(dtcond &tc) {
    double help;

    // acceleration cond
    help = pow((double)(h*h / f.len2()), (double)0.25);
    if (help < tc.dt) tc.set(help, 1, (int)getId(), -1, -1);

#ifdef SPH
    // courant cond
    help = h / (vsound + 1.2 * avAlpha*(vsound + 2*xmumax));
    if (help < tc.dt) tc.set(help, 2, (int)getId(), -1, -1);

    // energy cond
    if (u > uMin) if (du != 0.0) {
      help = dtFactor*(u + uMin)/Fabs(du);
      if (help < tc.dt)	tc.set(help, 3, (int)getId(), u, du);
    }

    // smoothing cond
    if (h > hMin) if (dh != 0.0) {
      help = dtFactor*h/Fabs(dh);
      if (help < tc.dt) tc.set(help, 4, (int)getId(), h, dh);
    }

    // density cond
    if (rho > rhoMin) if (drho != 0.0) {
      help = dtFactor*(rho + rhoMin)/Fabs(drho);
      if (help < tc.dt) tc.set(help, 5, (int)getId(), rho, drho);
    }
#ifdef SOLID
    // damage cond
    if (dm > 2. * dmMin) if (ddm != 0.0) {
      help = dtFactor*(dm + dmMin)/Fabs(ddm);
      if (help < tc.dt) tc.set(help, 6, (int)getId(), dm, ddm);
    }
    
    // sxx cond
    if (Fabs(sxx) > 2.*sMin) if (dsxx != 0.0) {
      help = dtFactor*(Fabs(sxx) + sMin)/Fabs(dsxx);
      if (help < tc.dt) tc.set(help, 7, (int)getId(), sxx, dsxx);
    }
    
    // sxy cond
    if (Fabs(sxy) > 2.*sMin) if (dsxy != 0.0) {
      help = dtFactor*(Fabs(sxy) + sMin)/Fabs(dsxy);
      if (help < tc.dt) tc.set(help, 8, (int)getId(), sxy, dsxy);
    }
    
    // sxz cond
    if (Fabs(sxz) > 2.*sMin) if (dsxz != 0.0) {
      help = dtFactor*(Fabs(sxz) + sMin)/Fabs(dsxz);
      if (help < tc.dt) tc.set(help, 9, (int)getId(), sxz, dsxz);
    }
    
    // syy cond
    if (Fabs(syy) > 2.*sMin) if (dsyy != 0.0) {
      help = dtFactor*(Fabs(syy) + sMin)/Fabs(dsyy);
      if (help < tc.dt) tc.set(help, 10, (int)getId(), syy, dsyy);
    }
    
    // syz cond
    if (Fabs(syz) > 2.*sMin) if (dsyz != 0.0) {
      help = dtFactor*(Fabs(syz) + sMin)/Fabs(dsyz);
      if (help < tc.dt) tc.set(help, 11, (int)getId(), sxz, dsyz);
    }
#endif
#endif
  }

  /***************************************************************************
   * Integrator (predictor/corrector)                                        *
   *-------------------------------------------------------------------------*
   **************************************************************************/
  void predictor(const double &dt) {
    pos += (v + deltav*xsphFactor)*dt + f*dt*dt*0.5;
    v   += f*dt;
    sf   = f;
#ifdef SPH
    h    += dt*dh;
    rho  += dt*drho; 
    u    += dt*du;
    if (h > Particle::hMax) throw ErrorHTooBig(h);
    //h     = Max(h,   hMin,   "Particle::Warning: hMin");
    //rho   = Max(rho, rhoMin, "Particle::Warning: rhoMin");
    //u     = Max(u,   uMin,   "Particle::Warning: uMin");
    h    = Max(h,   hMin);
    rho   = Max(rho, rhoMin);
    u     = Max(u,   uMin);
    sdh   = dh; sdrho = drho; sdu = du;
#ifdef SOLID
    dm   += dt*ddm;
    sxx  += dt*dsxx;
    syy  += dt*dsyy;
    sxy  += dt*dsxy;
    sxz  += dt*dsxz;
    syz  += dt*dsyz;
    dm    = Min(dm, (ftype)1.0);
    sddm  = ddm;  sdsxx = dsxx; sdsyy = dsyy;
    sdsxy = dsxy; sdsxz = dsxz; sdsyz = dsyz;
#endif
#endif
  }

  void corrector(const double &dt) {
    Vector<ftype> df  = f - sf;

    //pos += df*A*.5*dt*dt;
    v   += df*B*dt;
#ifdef SPH
    ftype  Cdt = C * dt;

    h   += Cdt*(dh   - sdh);
    rho += Cdt*(drho - sdrho); 
    u   += Cdt*(du   - sdu);
    if (h > Particle::hMax) throw ErrorHTooBig(h);
    //h    = Max(h,   hMin,   "Particle::Warning: hMin");
    //rho  = Max(rho, rhoMin, "Particle::Warning: rhoMin");
    //u    = Max(u,   uMin,   "Particle::Warning: uMin");
    h    = Max(h,   hMin);
    rho  = Max(rho, rhoMin);
    u    = Max(u,   uMin);
#ifdef SOLID
    sxx += Cdt*(dsxx - sdsxx);
    syy += Cdt*(dsyy - sdsyy);
    sxy += Cdt*(dsxy - sdsxy);
    sxz += Cdt*(dsxz - sdsxz);
    syz += Cdt*(dsyz - sdsyz);
    ftype olddm = dm;
    if (dm <  1.) dm += Cdt*(ddm  - sddm);
    if (dm >= 1.) dm  = olddm;
    if (dm >= 1.) { 
      dm = 1.; 
      sxx  = sxy  = sxz  = syy  = syz = 0; 
      dsxx = dsxy = dsxz = dsyy = dsyz = 0.;
    }
#endif
#endif
  }
};

#ifdef SPH
Eos           Particle::eos;
ftype         Particle::avAlpha,  Particle::avBeta, Particle::dampFactor;
ftype         Particle::dtFactor, Particle::hMax,   Particle::hMin;
ftype         Particle::hMinm1,   Particle::rhoMin, Particle::uMin;
ftype         Particle::xsphFactor, Particle::heatconG1, Particle::heatconG2;
#ifdef SOLID
Material<Mat> Particle::matTab;
ftype         Particle::dmMin, Particle::sMin, Particle::tiny;
#endif
#endif

ftype         Particle::G;
ftype         Particle::gravHb1, Particle::gravHb12, Particle::gravHb1m1, Particle::gravHb1m2;
ftype         Particle::gravHb2, Particle::gravHb22, Particle::gravHb2m1, Particle::gravHb2m2;

const ftype   Particle::third     = 1./3.;
const ftype   Particle::sqthird   = sqrt(third);
const ftype   Particle::A         = third;       // corrector
const ftype   Particle::B         = .5;          // corrector
const ftype   Particle::C         = .5;          // corrector

const Kernel  Particle::kern;
double        Particle::totNumPartm1;
bool          Particle::heatCon, Particle::xsph;
int           Particle::numVar;
std::string   Particle::varName;

#endif
