#ifndef SPHLATCH_CDAT_WRITER_H
#define SPHLATCH_CDAT_WRITER_H

/*
 *  cdat_writer.h
 *
 *  NOT PARALLEL SAFE!!!
 *
 *
 *  Created by Andreas Reufer on 15.11.07.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include <fstream>
#include <string>

#include <rpc/rpc.h>

#include "typedefs.h"

//#ifdef SPHLATCH_PARALLEL
//#include "communication_manager.h"
//#endif

namespace sphlatch
{

class CDATwriter
{
public:
  void open( std::string _filename,
             std::string _title,
             fType _time,
             size_t _noParts,
             std::vector<std::string> _varNames,
             bool _dblPrec);
  void close();

  void write(size_t _part, size_t _var, fType _val);

  CDATwriter();
  ~CDATwriter();

private:
  std::fstream fout;
  std::fstream::pos_type curPos, startXDRpos;

  stringVectType vars;
  //size_t noParts, noVars, myOffset;
  size_t noVars;

  XDR *xdrs;

  size_t xdrFloatSize, xdrDoubleSize;
  size_t buffsize;
  char *xdr_buff;

  float flt_buff;
  double dbl_buff;

  bool dblPrec;

  std::string readLine( std::fstream& _fin );

//#ifdef SPHLATCH_PARALLEL
//  typedef sphlatch::CommunicationManager commManagerType;
//  commManagerType& CommManager;
//#endif
};

CDATwriter::CDATwriter()
//#ifdef SPHLATCH_PARALLEL
//  : CommManager(commManagerType::instance())
//#endif
{
  xdrFloatSize = 4;
  xdrDoubleSize = 8;

  buffsize = xdrFloatSize;
  dblPrec = false;
}

CDATwriter::~CDATwriter()
{
  delete [] xdr_buff;
  delete xdrs;
}

void CDATwriter::open( std::string _filename,
                       std::string _title,
                       fType _time,
                       size_t _noParts,
                       std::vector<std::string> _varNames,
                       bool _dblPrec)
{
  noVars = _varNames.size();
  dblPrec = _dblPrec;

  fout.open( _filename.c_str(), std::ios::out | std::ios::binary);
  
  fout << _title << "\n";
  fout << _noParts << "," << -static_cast<int>(noVars) << "\n3\n";

  ///
  /// the order of those parameters is strict!
  ///
  fout << "TIME\n" << _time << "\n";
  fout << "XDRD\n0\n";
  
  if ( _dblPrec)
    fout << "DOUBLE\n0\n";
  else
    fout << "FLOAT\n0\n";


  for (size_t i = 0; i < noVars; i++)
  {
    fout << _varNames[i]<< "\n";
  }

  fout << "^^^\n";

  ///
  /// prepare XDR decoding stuff
  ///
  if ( dblPrec )
    buffsize = xdrDoubleSize;
  else
    buffsize = xdrFloatSize;

  xdr_buff = new char[buffsize];
  xdrs = new XDR;

  startXDRpos = fout.tellp();
}

void CDATwriter::close()
{
  fout.close();
}


void CDATwriter::write(size_t _part, size_t _var, fType _value)
{
  const size_t curOffset = ( _part*noVars + _var )*buffsize; 
  curPos = startXDRpos;
  
  //curPos += (myOffset + curOffset);
  curPos += curOffset;
  fout.seekg(curPos);
  
  xdrmem_create(xdrs, xdr_buff, buffsize, XDR_ENCODE);

  if ( buffsize == xdrFloatSize )
  {
    flt_buff = static_cast<float>(_value);
    xdr_float(xdrs, &flt_buff);
  }
  else
  {
    dbl_buff = static_cast<double>(_value);
    xdr_double(xdrs, &dbl_buff);
  }

  fout.write(xdr_buff, buffsize);
}

};


#endif
