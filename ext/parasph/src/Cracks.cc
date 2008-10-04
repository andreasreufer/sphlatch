/******************************************************************************
 * ParaSPH -- Version 24.11.2003                                              *
 *----------------------------------------------------------------------------*
 * File:      Cracks.cc                                                       *
 * Purpose:   
 *****************************************************************************/
#ifndef CRACKS_CC
#define CRACKS_CC

#include <fstream>
#include <stdlib.h>

#include "Def.cc"
#include "Material.cc"
#include "Particle.cc"

class Cracks {
private:
  std::fstream file;

public:
  class ErrorFileNotOpen {
  public:
    std::string file;
    ErrorFileNotOpen(const std::string &_file) { file = _file; }
  };

  Cracks(const std::string &crackPath, const bool &write, const bool &noisy) {
    if (write) {
      file.open(crackPath.c_str(), std::fstream::binary | 
		std::fstream::out | std::fstream::trunc);
      if (!file) throw ErrorFileNotOpen(crackPath); 
      if (noisy) std::cout << "Writing file '" << crackPath << "'" 
			   << std::endl;
    } else {
      file.open(crackPath.c_str(), std::fstream::binary | std::fstream::in);
      if (!file) throw ErrorFileNotOpen(crackPath);   
      if (noisy) std::cout << "Reading crack configuration from file '"
			   << crackPath << "'" << std::endl;
    }
  }

  ~Cracks() { file.close(); }

#ifdef SOLID
  void setFlaws(Particle *part, const int &from, const int &to,
		const int &iallow, const int &ishield, const ftype &radius,
		const ftype &bulk, const ftype &rho0) {
    double         cw, pw, vol = 0.;
    ftype          eps; 
    long           flawsmax, i, nr;
    Mat            m;

    // Assuming that each body consist of one material
    m = part[from].getMatVar();

    for (i = from; i < to; i++) {
      // Is the whole body of the same material
      if (part[i].mat != part[from].mat)
	std::cerr << "setflaws: different materials in same body!"
		  << std::endl;

      part[i].acoef = .4*sqrt((bulk+4.*m.mu/3.)/rho0) / (2.*part[i].h);
      part[i].young = 9.*bulk * m.mu  / (3.*bulk + m.mu);
      part[i].flaws = 0;

      vol += part[i].mass / part[i].rho;
    }
    
    // What is the volume of your body? If you are simulating the whole
    // body then the method right above may be correct. Otherwise you will
    // need to specify something like this:
    // vol = 4.18879*pow((double)radius, 3.); // 4*Pi/3*r^3

    if (ishield == 0) for (i = from; i < to; i++) part[i].grav = 0.;

    if (iallow == 0) {
      for (i = from; i < to; i++) {
	part[i].epsmin = 1.e20; part[i].m     = 1.; 
	part[i].flaws  = 0;     part[i].acoef = 0.;
      }
    } else { 
      cw = 1. /(m.cweib * vol);
      pw = 1. / m.pweib;
      // The number of flaws should be at least a factor of 10 higher than
      // the number of particles. The first line is the same as within the
      // old fortran code, the second is used, if you use millions of p.
      flawsmax = (int)1e7;
      //flawsmax = 10 * (to - from);

      for (nr = 0; nr < flawsmax; nr++) {
	i = from + (int) ((double)(to - from) * rand() / (RAND_MAX + 1.0) );
	eps = pow(cw*((double)nr+1.), pw);
	if (part[i].flaws == 0) part[i].epsmin = eps;
	part[i].flaws++;
	part[i].m = eps;
      }

      for (i = from; i < to; i++) {
	if (part[i].flaws == 0) {
	  part[i].flaws = 1;
	  part[i].epsmin = pow(cw*((double)flawsmax+1.), pw);
	}
	if (part[i].flaws == 1) part[i].m = 1.;
	else part[i].m = log((ftype)part[i].flaws) / log(part[i].m / part[i].epsmin);
      }
    }

    for (i = from; i < to; i++) {
      file.write((char *) &part[i].epsmin, sizeof(part[i].epsmin));
      file.write((char *) &part[i].m,      sizeof(part[i].m));
      file.write((char *) &part[i].flaws,  sizeof(part[i].flaws));
      file.write((char *) &part[i].acoef,  sizeof(part[i].acoef));
      file.write((char *) &part[i].young,  sizeof(part[i].young));
      file.write((char *) &part[i].grav,   sizeof(part[i].grav));
    }
  }

  // This mehtod has to work in parallel
  void readFlaws(Particle *part, const int &start, const int &numPart) {
    ftype fdum;
    int   idum;

    file.seekg(start*(5*sizeof(fdum)+sizeof(idum)));
    for (int i = 0; i < numPart; i++) {
      file.read((char *) &fdum, sizeof(fdum)); part[i].epsmin = fdum;
      file.read((char *) &fdum, sizeof(fdum)); part[i].m      = fdum;
      file.read((char *) &idum, sizeof(idum)); part[i].flaws  = idum;
      file.read((char *) &fdum, sizeof(fdum)); part[i].acoef  = fdum;
      file.read((char *) &fdum, sizeof(fdum)); part[i].young  = fdum;
      file.read((char *) &fdum, sizeof(fdum)); part[i].grav   = fdum;
    }
    //std::cout << "epsmin "<< part[1].epsmin << std::endl;
  }
#endif
};

#endif
