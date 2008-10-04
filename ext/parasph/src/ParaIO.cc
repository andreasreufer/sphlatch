#ifndef PARAIO_CC
#define PARAIO_CC

#include <string>

#include "Cracks.cc"
#include "Manager.cc"
#include "Particle.cc"
#include "TList.cc"
#include "Xdr.cc"

class ParaIO {
private:
  TList<Particle> &pList;
  std::string      xdrFile, title;

public:
  ParaIO(TList<Particle> &_pList, const std::string &_xdrFile) : 
    pList(_pList), xdrFile(_xdrFile) {}

  void readData(const std::string &filename, double &time) {
    int  numPart;
    long p, pos, size, start;
    Xdr  input(Particle::numVar, global::noisy);

    input.open(filename, 'r');
    pos   = input.readheader(numPart, time, title);
    global::setTotNumPart(numPart);
    Particle::calcHMin();
    start =  global::rank    * global::totNumPart/global::npro;
    size  = (global::rank+1) * global::totNumPart/global::npro - start;
    input.setpos(start * Particle::numVar * Xdr::varSize + pos);
    input.createXdr();
    pList.setSize(size);
    for (p = 0; p < size; p++) pList[p].read(input);
    input.close();
  }

#ifdef SPH
  int getNumBody() {
    // Getting the maximal number of different bodies
    int bodies = 0, numBody, p;
    for (p = 0; p < pList.getSize(); p++) 
      bodies = Max(bodies, pList[p].bodyNr());
    bodies++;
    MCW.Allreduce(&bodies, &numBody, 1, MPI_INT, MPI_MAX);
    if (global::noisy) std::cout << "Different bodies: " << numBody 
				 << std::endl;
    return numBody;
  }
#endif

#ifdef SOLID
  void readCrackFile(const std::string &crackFile) {
    Cracks cracks(crackFile, false, global::noisy);
    long   size, start;

    start =  global::rank    * global::totNumPart/global::npro;
    size  = (global::rank+1) * global::totNumPart/global::npro - start;
    cracks.readFlaws(&pList[0], start, size);
  }
#endif

  double getInput(Manager &man) { 
    double      time;
    std::string crackFile, path, xdrFile;
 
    man.getValue("simulation.path", path);
    man.getValue("input.xdrFile",   xdrFile);
    xdrFile.insert(0, path + "/");
    readData(xdrFile, time);
#ifdef SOLID
    man.getValue("input.crackFile", crackFile);
    crackFile.insert(0, path + "/");
    readCrackFile(crackFile);
#endif
    return time;
  }

  void pWriteData(const std::string &filename, const double &time) {
    long c, p, pos, *psize = new long[global::npro], size = pList.getSize();
    Xdr  output(Particle::numVar, global::noisy);

    output.open(filename, 'w');
    pos = output.writeheader(Particle::varName, global::totNumPart, time, 
			     title);
    
    MCW.Allgather(&size, sizeof(long), MPI::BYTE, 
		  psize, sizeof(long), MPI::BYTE);

    for (c = 0; c < global::rank; c++) 
      pos += psize[c] * Particle::numVar * Xdr::varSize;
    output.setpos(pos);
    output.createXdr();
    for (p = 0; p < size; p++) pList[p].write(output);
    output.close();
  }

  bool produceOutput(Manager &man, const double &time, double &lastSave) {
    char        num[3];
    int         nr;
    double      next;
    std::string out;
 
    if (man.getMinValueGreater("output.saveTime", lastSave, next)) {
      if (time >= next) {
        try { man.getValue("output.lastSave", nr); }
        catch (Manager::ErrorVarNotFound &e) { nr = -1; }
        sprintf(num, "%.3i", ++nr);
        man.getValue("simulation.path", out);
        out += "/out###.xdr";
        out.replace(out.find("###"), 3, num);
        pWriteData(out, time);
        lastSave = time;
        man.setValue("output.lastName", out);
        man.setValue("output.lastSave", num);
      }
    } else return true;
    return false;
  }
};

#endif
