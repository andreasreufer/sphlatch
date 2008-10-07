/******************************************************************************
 * ParaSPH -- Version 24.11.2003                                              *
 *----------------------------------------------------------------------------*
 * File:      Xdr.cc                                                          *
 * Purpose:   This class reads and writes the XDR type files which are read   *
 *            by 'xplot'. The moment 'xplot' is not needed anymore this class *
 *            can be rewritten easily.                                        *
 *            To understand what is done here I recommend to have a look at a *
 *            xdr-file.                                                       *
 *****************************************************************************/
#ifndef XDR_CC
#define XDR_CC

#include <iostream>
#include <string>
// These two libraries are standard in Linux (meanwhile) and most other
// real operating systems :-)
#include <rpc/types.h>
#include <rpc/xdr.h>

class Xdr {
private:
  bool        noisy;
  double      dvar;
  FILE       *file;
  float       fvar;
  int         numvar;
  std::string name, tstr;
  void       (Xdr::*xRead)(XDR *, ftype *);
  void       (Xdr::*xWrite)(XDR *, ftype *);
  XDR         xdrs[1];
  xdr_op      xop;

  void doubleRead(XDR *xdrs, ftype *var) { 
    xdr_double(xdrs, &dvar); *var = (ftype)dvar;
  }
  void doubleWrite(XDR *xdrs, ftype *var) {
    dvar = (double)*var; xdr_double(xdrs, &dvar);
  }
  void singleRead(XDR *xdrs, ftype *var) { 
    xdr_float(xdrs, &fvar); *var = (ftype)fvar;
  }
  void singleWrite(XDR *xdrs, ftype *var) {
    fvar = (float)*var; xdr_float(xdrs, &fvar);
  }

public:
  static int varSize;

  class ErrorFileNotOpen {
  public:
    std::string file;
    ErrorFileNotOpen(const std::string &_file) { file = _file; }
  };
  class ErrorWrongVariableNumber {
    public:
    int         num;
    std::string file;
    ErrorWrongVariableNumber(const int &_num, const std::string &_file) {
      num = _num; file = _file;
    }
  };

  Xdr(const int &_numvar, const bool &_noisy) {
    numvar = _numvar; noisy = _noisy;
  }

  void open(const std::string &_name, const char &mode) {
    int count = 5;
    name = _name; 
    do { file = fopen(name.c_str(), &mode); } while (--count > 0 && !file);
    if (!file) throw ErrorFileNotOpen(name);
    if (mode == 'w') xop = XDR_ENCODE; else xop = XDR_DECODE;
  } 

  long writeheader(const std::string &varname, const int &numpart, 
		   const double &time, const std::string &title) {
    if (noisy) std::cout << "Writing '"<< name << "', " << numpart
			 << " parts., " << numvar << " vars., time = "
			 << time << ", title = " << title << std::endl;
    if (varSize == sizeof(double)) {
      xWrite = &Xdr::doubleWrite; tstr = "DOUBLE";
    } else {
      xWrite = &Xdr::singleWrite; tstr = "FLOAT";
    }
    fprintf(file, "%s\n%d, %d\n3\nTIME\n", title.c_str(), numpart, -numvar);
    fprintf(file, "%e\nXDRD\n0\n%s\n0\n%s^^^\n", time, tstr.c_str(), 
	    varname.c_str());
    return ftell(file);
  }

  long readheader(int &numpart, double &time, std::string &title) {
    char dum[80], tit[80], type[80];
    int  nv;

    fscanf(file, "%s\n%d, %d\n%s\n%s\n", tit, &numpart, &nv, dum, dum);
    fscanf(file, "%lf\n%s\n%s\n%s\n%s\n", &time, dum, dum, type, dum);
    if (strcmp(type, "FLOAT") != 0) {
      if (global::noisy) std::cout << "Reading double precision" << std::endl; 
      xRead = &Xdr::doubleRead; varSize = sizeof(double);
    } else {
      if (global::noisy) std::cout << "Reading single precision" << std::endl; 
      xRead = &Xdr::singleRead; varSize = sizeof(float);
    }
    tstr = type; title = tit;
    if (-nv != numvar) throw ErrorWrongVariableNumber(-nv, name);
    for (int i = 0; i <= numvar; i++) fscanf(file, "%s\n", dum);
    if (noisy) std::cout << "Reading '"<< name << "', " << numpart
			 << " parts., " << numvar << " vars., time = "
			 << time << ", title = " << title << std::endl;
    return ftell(file);
  }

  // Before reading or writing we need to initialize and set the position
  // in the file (the latter is needed in parallel for random access)
  void createXdr() { xdrstdio_create(xdrs, file, xop); }

  void setpos(const long &pos) { fseek(file, pos, SEEK_SET); }

  // There is just one method to read/write data, the xdr_type depends
  // on the 'mode' in the 'open' method.
  void read(ftype *store) {
    for (int i = 0; i < numvar; i++) (*this.*xRead)(xdrs, &store[i]);
  }
  void write(ftype *store) {
    for (int i = 0; i < numvar; i++) (*this.*xWrite)(xdrs, &store[i]);
  }

  void close() {
    xdr_destroy(xdrs);
    fflush(file);
    fclose(file);
  }
};

int Xdr::varSize;

#endif
