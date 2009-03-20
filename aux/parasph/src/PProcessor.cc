#ifndef PROCESSOR_CC
#define PROCESSOR_CC

#include "CellList.cc"
#include "Debug.cc"
#include "Def.cc"
#include "DomainDecomp.cc"
#include "Manager.cc"
#include "CellMap.cc"
#include "ParaIO.cc"
#include "Particle.cc"
#include "Rankspace.cc"
#include "Timer.cc"
#include "TList.cc"

#include <fstream>

class Processor {
private:
  CellList        cellList;
  DomainDecomp    domainDecomp;
  double          time, dt;
  int             numBody;
  Manager         man;
  CellMap         cellMap;
  ParaIO          paraIO;
  Rankspace       rankspace;
  TList<Particle> pList;

  static double   courant, dtMax;

public:
  Processor(const std::string &pathName) : cellList(pList), 
					   domainDecomp(pList), 
					   man(pathName, true),
					   cellMap(pList),
					   paraIO(pList, pathName),
					   rankspace(pList) {
    global::init();
    Particle::init(man);
    time = paraIO.getInput(man);
#ifdef SPH
    numBody = paraIO.getNumBody();
#endif

    rankspace.init();
    cellMap.init();
    cellList.init();

    man.getValue("SPH.courant", courant);
    man.getValue("SPH.dtMax",   dtMax);
    man.getValue("simulation.theta", global::theta);
  }
  
  void bookKeeping() {
    DEBUG("bookKeeping()", "");
    rankspace.pRanking(0); rankspace.pRanking(1); rankspace.pRanking(2);
    DMORE("s", "rankspace done");
    cellMap.pCreate(); cellMap.findNear();
    DMORE("s", "cellMap done");
    domainDecomp.pDomainDecomp();
    DMORE("s", "domainDecomp done");
  }

  void deriv() {
    DEBUG("deriv()", "");
    int   p, size = pList.getSize();
    Timer t;

    for (p = 0; p < size; p++) pList[p].preForces();
#ifdef PHASESTAT
    const size_t noPhases = 8;
    ftype phase_histogram[noPhases];
    for (p = 0; p < size; p++) 
    {
      // quick hack: only do it for ice
      if ( lrint(pList[p].getMatId()) == 17 )
      {
        const size_t phase = pList[p].getPhase();
        phase_histogram[phase] += pList[p].getMass();
      }
    };  
      
    ftype phase_histogram_tot[noPhases];
    MCW.Allreduce(&phase_histogram,
                  &phase_histogram_tot,
                  noPhases, MPI_ftype, MPI_SUM);
    
    if (global::noisy)
    {
      std::fstream fout;
      fout.open("ice_phases_histogram.txt", std::ios::out | std::ios::app);
      
      std::cout << "  >>> ice phases: [";
      fout << time << "\t";

      for (size_t i = 0; i < noPhases; i++)
      {
        fout << phase_histogram_tot[i] << "\t";
        std::cout << phase_histogram_tot[i] << ",";
      }
      
      fout << "\n";
      fout.close();
      std::cout << "]\n";
    }
#endif

    cellList.pFindNeighbourCells();
    cellList.pSpreadGhosts();
#ifdef CORR
    for (p = 0; p < size; p++) cellList.corrtensor(&pList[p]);
    cellList.pCollectGhosts();

    for (p = 0; p < size; p++) pList[p].inverse();

    cellList.pFindNeighbourCells();
    cellList.pSpreadGhosts();
#endif
    DMORE("s", "calculating Forces ...");
    for (p = 0; p < size; p++) {
      t.reset(); t.start();
#ifdef SPH
      cellList.forcesSPH(&pList[p]);
#endif
#ifdef GRAV
      cellList.forcesGrav(&pList[p]);
#endif
      t.stop();
      pList[p].setWork(1. +  t.result());
    }

    cellList.pCollectGhosts();

    for (p = 0; p < size; p++) pList[p].postForces();

#ifdef SPH
    ftype *buf = new ftype[numBody], *divMax = new ftype[numBody];
    for (int i = 0; i < numBody; i++) buf[i] = divMax[i] = 0.;
    for (p = 0; p < size; p++) pList[p].maxdivv(buf);
    MCW.Allreduce(buf, divMax, numBody, MPI_ftype, MPI_MAX);
    for (p = 0; p < size; p++) pList[p].hdot(divMax);

    // the energy calculation is not correct if GRAVITY is specified. Sorry.
    ftype pE = 0., E;
    for (p = 0; p < size; p++) pE += pList[p].getEnergy();
    MCW.Allreduce(&pE, &E, 1, MPI_ftype, MPI_SUM);
    if (global::noisy) 
      std::cout << "  >>> Total energy E = " << E << " <<<" << std::endl;
#endif
  }

  /***************************************************************************
   * Timestep & integrator                                                   *
   *-------------------------------------------------------------------------*
   * integrator: reads first the input file, runs the firstStep method to    *
   *             get a first estimate of the derivations, and then starts    *
   *             the integration. The program will run until 'done' is set   *
   *             to 'true' which happens as soon as the last output file has *
   *             been written.                                               *
   * step:       advances the time, calls the predictor and corrector method,*
   *             starts the DD and tree building and computes the calls      *
   *             deriv. Time and timestep are printed to standard output.    *
   * firstStep:  the same as step, but without predictor and corrector.      *
   * pTimeStep:  computes the minimal timestep of all processes and prints   *
   *             that information to standard out.                           *
   **************************************************************************/
  double pTimestep() {
    DEBUG("pTimestep()", "");
    dtcond *tc = new dtcond[global::npro];
    struct { double dt; int rank; } loc, glob;

    for (int p = 0; p < pList.getSize(); p++) 
      pList[p].tsconds(tc[global::rank]);

    tc[global::rank].dt *= courant;
    if (dtMax < tc[global::rank].dt) tc[global::rank].set(dtMax, 0, -1, -1, -1);
    loc.dt = tc[global::rank].dt; loc.rank = global::rank;
    MCW.Allreduce(&loc, &glob, 1, MPI::DOUBLE_INT, MPI::MINLOC);

    if (global::rank == glob.rank) tc[global::rank].print(global::rank);
    return glob.dt;
  }


  void firstStep() {
    DEBUG("firstStep()", "");

    bookKeeping();
    deriv();

    dt = pTimestep();
    dt = Min(dt , 0.01*dtMax);
    if (global::noisy) 
      std::cout << "  step 0: t = " << time << ", dt = " << dt << std::endl;
  }

  void step(const int &intStep) {
    DEBUG("step(intStep):", "si", "intStep = ", intStep);
    int   p;
    Timer t; t.reset(); t.start();

    time += dt;

    for (p = 0; p < pList.getSize(); p++) pList[p].predictor(dt);

    bookKeeping();
    deriv();

    for (p = 0; p < pList.getSize(); p++) pList[p].corrector(dt);

    dt = pTimestep();
    if(intStep < 30){
       dt = Min(dt , 0.01*dtMax);
       if (global::noisy) std::cout << "dt limited to 0.01 dtmax for first 30 steps. Currently:  "<< intStep << std::endl;
    }


    t.stop();
    
    float *f = new float[global::npro], fl = t.result();
    MCW.Allgather(&fl, 1, MPI::FLOAT, f, 1, MPI::FLOAT);
    if (global::noisy) {
      std::cout << "  step " << intStep << ": t = " << time << ", dt = " << dt
		<< ", cpu(s) = ";
      for (int i = 0; i < global::npro; i++) std::cout << f[i] << " ";
      std::cout << std::endl;
    }
  }

  void pCalcGravH() {
    ftype hMin = 1.e30, hMinGlob;
    
    for (int p = 0; p < pList.getSize(); p++)
      hMin = Min(pList[p].getH(), hMin);
    MCW.Allreduce(&hMin, &hMinGlob, 1, MPI_ftype, MPI::MIN);
    
    Particle::setGravH(hMinGlob);
  }

/*
  int ox[4], oy[4];
  std::ofstream ps;
  
  void psHeader() {
    ps.open("test.ps");
    ps << "%!PS-Adobe-2.0\n%%Orientation: Portrait\n%%DocumentMedia: A4 596 "
       << "842\n%%Pages: (atend)\n%%EndComments\n\n";
    ps << "/point{\n  0.1 0 360 arc closepath\n  stroke\n} def\n\n";
  
    for (int i = 0; i < pList.getSize(); ++i)
      psDraw(pList[i].pos[0], pList[i].pos[1], i, false);
  }

  void psDraw(const ftype &x, const ftype &y, const int &n, bool ok) {
    int nx = (int)(x*2+300), ny = (int)(y*2+300);

    if (ok) ps << ox[n] << " " << oy[n] << " moveto\n"
	       << nx << " " << ny << " lineto\n\nstroke\n";

    ox[n] = nx; oy[n] = ny;
  }

  void psClose() { ps.close(); }
*/

  void integrator() {
    DEBUG("integrator()", "");
    bool   done;
    double lastSave = 0.;
    int    intStep;

    pCalcGravH();
    firstStep();
    

    //psHeader();
  
    for (done = false, intStep = 1; !done; intStep++) {
      step(intStep);

      //for (int i = 0; i < pList.getSize(); ++i)
	//psDraw(pList[i].pos[0], pList[i].pos[1], i, true);

      done = paraIO.produceOutput(man, time, lastSave);
//      if (intStep >= 1) done = true;
    }
    //psClose();
  }
};

double Processor::courant, Processor::dtMax;

#endif
