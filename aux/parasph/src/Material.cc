/******************************************************************************
 * ParaSPH -- Version 13.01.2004                                              *
 *----------------------------------------------------------------------------*
 * File:      Material.cc                                                     *
 * Purpose:   This class takes care of the various material constants like    *
 *            bulk modulus, initial density etc.                              *
 *            The file that has to be read has a very special format and is   *
 *            probably named like 'tillotson.input' and/or 'matdata.input'.   *
 *            You better leave it alone :-)                                   *
 *            'aneos.input' is not read with this class; ANEOS does this.     *
 *            The numbering in the file is fortran standard. Remember that    *
 *            the first material type in C has the number 0.                  *
 * Note:      The variable C_FORTRAN_MATERIAL_SHIFT can be 1 or 0. In case of *
 *            1, the materials are numbered in fortran style, ie. from 1 to   *
 *            nmat; otherwise (0) from 0 to nmat - 1.                         *
 *            Maybe this could go into '.config'. What do you think about it? *
 *****************************************************************************/
#ifndef MATERIAL_CC
#define MATERIAL_CC

// only set the shifting variable if no value has been given to the compiler
// 0: C style (0..nmat-1); 1: Fortran style (1..nmat); 
#ifndef C_FORTRAN_MATERIAL_SHIFT
#define C_FORTRAN_MATERIAL_SHIFT 1
#endif

#include <fstream>
#include <stdlib.h>
#include <string>

#include "Def.cc"

typedef struct {
  ftype rho0, A, B, a, b, alpha, beta, u0, Eiv, Ecv;
} eosMat;

typedef struct {
  ftype mu, umelt, yield, pweib, cweib;
} Mat;

template<class T>
class Material {
private:
  int  nmat;
  T   *konst;

public:
  class ErrorFileNotFound {
  public:
    std::string file;
    ErrorFileNotFound(const std::string &_file) { file = _file; }
  };
  class ErrorFileCorrupt  {
  public:
    std::string file;
    ErrorFileCorrupt(const std::string &_file) { file = _file; }
  };

  void read(const bool &noisy, const std::string &matDataFile) {
    int         var = sizeof(T)/sizeof(ftype);
    ftype      *fp;
    std::string line, text;

    std::ifstream in(matDataFile.c_str());
    if (in.fail()) throw ErrorFileNotFound(matDataFile);

    // The first lines of the material files contain a comment (usually 
    // 12 lines). You can use following line to jump over these lines or the
    // two lines after (default):
    // for (int i = 0; i < 12; i++) in.ignore(80, '\n'); 
    do { getline(in, line);
    } while (line.find("number of materials") == std::string::npos);

    in >> nmat;
    if (in.fail()) throw ErrorFileCorrupt(matDataFile);
    konst = new T[nmat];
    if (noisy) std::cout << "Reading material file '" << matDataFile << "': "
			 << C_FORTRAN_MATERIAL_SHIFT << ".." 
			 << nmat-1 + C_FORTRAN_MATERIAL_SHIFT << " different "
			 << "materials, " << var << " constants each."
			 << std::endl;

    for (int c = 0; c < var; c++) {
      for (int i = 0; i < 4; i++) in.ignore(80, '\n');
      for (int i = 0; i < nmat; i++) {
	in >> text;
	if (in.fail()) throw ErrorFileCorrupt(matDataFile);
	fp  = (ftype*)&konst[i]; fp += c; 
	*fp = atof(text.c_str());
      } 
    }
  }

  // Get a constant list for a material 'mat' (see material file)
  const T &get(const int &mat) const { 
    int matNr = mat - C_FORTRAN_MATERIAL_SHIFT;
    if (matNr < 0 || matNr >= nmat) std::cout << "Material error " << mat << " " << matNr << " " << nmat << std::endl;
    assert(0 <= matNr && matNr < nmat);
    return konst[matNr];
  }
};

#endif
