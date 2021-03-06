#ifndef SPHLATCH_CDAT_READER_H
#define SPHLATCH_CDAT_READER_H

/*
 *  cdat_reader.h
 *
 *
 *  Created by Andreas Reufer on 15.11.07.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include <fstream>
#include <string>

#include <rpc/rpc.h>
#include <boost/lexical_cast.hpp>

#include "typedefs.h"

#ifdef SPHLATCH_PARALLEL
 #include "communication_manager.h"
#endif

namespace sphlatch
{
class CDATreader
{
public:
   void open(std::string _filename);
   stringVectType getVars();
   attrMapType getAttrs();
   size_t getNoParts();
   void close();

   fType read(size_t _part, size_t _var);

   CDATreader();
   ~CDATreader();

private:
   std::fstream           fin;
   std::fstream::pos_type curPos, startXDRpos;

   stringVectType vars;
   attrMapType    attrs;
   size_t         noParts, noVars, myOffset;

   XDR* xdrs;

   size_t xdrFloatSize, xdrDoubleSize;
   size_t buffsize;
   char   * xdr_buff;

   float  flt_buff;
   double dbl_buff;

   std::string readLine(std::fstream& _fin);

#ifdef SPHLATCH_PARALLEL
   typedef sphlatch::CommunicationManager   commManagerType;
   commManagerType& CommManager;
#endif
};

CDATreader::CDATreader()
#ifdef SPHLATCH_PARALLEL
   : CommManager(commManagerType::instance())
#endif
{
   xdrFloatSize  = 4;
   xdrDoubleSize = 8;

   buffsize = xdrFloatSize;
}

CDATreader::~CDATreader()
{
   delete [] xdr_buff;
   delete xdrs;
}

void CDATreader::open(std::string _filename)
{
   fin.open(_filename.c_str(), std::ios::in | std::ios::binary);
   std::string curLine;
   curLine = readLine(fin); /// reads the name, currently not used

   curLine = readLine(fin);
   size_t noTotParts = boost::lexical_cast<size_t>(
      curLine.substr(0, curLine.find(',', 0)));
   int noVarsRaw = boost::lexical_cast<int>(
      curLine.substr(curLine.find(',', 0) + 2));

   /// noAttrsRaw means, there are also attributes in the header
   if (noVarsRaw < 0)
   {
      size_t noAttrs = boost::lexical_cast<size_t>(readLine(fin));

      for (size_t i = 0; i < noAttrs; i++)
      {
         std::string attrKey = readLine(fin);
         std::string attrVal = readLine(fin);

         if (attrKey == "FLOAT")
         {
            buffsize = xdrFloatSize;
         }
         else if (attrKey == "DOUBLE")
         {
            buffsize = xdrDoubleSize;
         }
         else if (attrKey == "XDRD")
         { }
         else if (attrKey == "DIMENSION")
         { }
         else
         {
            attrs[attrKey] = boost::lexical_cast<fType>(attrVal);
         }
      }
   }
   noVars = abs(noVarsRaw);

   vars.resize(noVars);
   for (size_t i = 0; i < noVars; i++)
   {
      vars[i] = readLine(fin);
   }


   noParts  = noTotParts;
   myOffset = 0;

   if (readLine(fin) != "^^^")
   {
      std::cerr << "file does not contain magic separator ^^^!\n";
      vars.resize(0);
      return;
   }

   ///
   /// prepare XDR decoding stuff
   ///
   xdr_buff = new char[buffsize];
   xdrs     = new XDR;

   startXDRpos = fin.tellp();
}

stringVectType CDATreader::getVars()
{
   return(vars);
}

attrMapType CDATreader::getAttrs()
{
   return(attrs);
}

size_t CDATreader::getNoParts()
{
   return(noParts);
}

void CDATreader::close()
{
   fin.close();
}

fType CDATreader::read(size_t _part, size_t _var)
{
   const size_t curOffset = (_part * noVars + _var) * buffsize;

   curPos  = startXDRpos;
   curPos += (myOffset + curOffset);
   fin.seekg(curPos);

   ///
   /// read desired variable and convert it to native type
   ///
   fin.read(xdr_buff, buffsize);
   xdrmem_create(xdrs, xdr_buff, buffsize, XDR_DECODE);

   if (buffsize == xdrFloatSize)
   {
      xdr_float(xdrs, &flt_buff);
      return(static_cast<fType>(flt_buff));
   }
   else
   {
      xdr_double(xdrs, &dbl_buff);
      return(static_cast<fType>(dbl_buff));
   }
}

std::string CDATreader::readLine(std::fstream& _fin)
{
   const size_t lineBuffsize = 16384;
   std::string  line;
   char         linebuff[lineBuffsize];

   _fin.getline(linebuff, lineBuffsize);
   // remove LF
   line.assign(linebuff, 0, _fin.gcount() - 1);
   return(line);
}
};


#endif
