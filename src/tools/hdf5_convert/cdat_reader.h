#ifndef SPHLATCH_CDAT_READER_H
#define SPHLATCH_CDAT_READER_H

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


std::vector<int> LoadCDAT(std::string FileName);
bool SaveCDAT(std::string FileName, std::vector<int> Attributes);
void SetSinglePrecOut(void);
void SetDoublePrecOut(void);

VarMapType& VM;
MemoryManagerType& MemoryManager;

private:
std::string ReadLine(std::fstream *_fin);
static self_pointer _instance;

size_t sizeOfXDRfloat, sizeOfXDRdouble, datawidthSave;
};

template <typename SimTrait>
typename IOManager<SimTrait>::self_pointer IOManager<SimTrait>::_instance = NULL;

template <typename SimTrait>
typename IOManager<SimTrait>::self_reference IOManager<SimTrait>::Instance(void)
{
  if (_instance != NULL)
    return *_instance;
  else
    {
      _instance = new IOManager;
      return *_instance;
    }
}

template <typename SimTrait>
IOManager<SimTrait>::IOManager(void) : VM(VarMapType::Instance()), MemoryManager(MemoryManagerType::Instance())
{
  sizeOfXDRfloat = 4;
  sizeOfXDRdouble = 8;

  if (sizeof(value_type) == sizeOfXDRfloat)
    {
      datawidthSave = sizeOfXDRfloat;
    }
  else
    {
      datawidthSave = sizeOfXDRdouble;
    }
}

template <typename SimTrait>
IOManager<SimTrait>::~IOManager(void)
{
}

template <typename SimTrait>
std::vector<int> IOManager<SimTrait>::LoadCDAT(std::string FileName)
{
  matrix_reference Data(MemoryManager.Data);
  string_reference Name(MemoryManager.Name);

  using boost::lexical_cast;

#ifndef OOSPH_STANDALONE
  const size_t RANK = MPI::COMM_WORLD.Get_rank();
  const size_t SIZE = MPI::COMM_WORLD.Get_size();
#else
  const size_t RANK = 0;
  const size_t SIZE = 1;
#endif

  std::fstream *fin = new std::fstream;
  std::fstream::pos_type filepos;
  int numtotparticles, numlocparticles, numparams, numattributes, attrmaxindex = 0;
  size_t buffsize, dataoffset;
  size_t datawidthLoad = sizeOfXDRfloat;
  std::vector<int> attributes;


  fin->open(FileName.c_str(), std::ios::in | std::ios::binary);
  if (!fin)
    {
      std::cerr << "Error Loading Data " << FileName << std::endl;
      attributes.resize(0);
      return attributes;
    }

  Name = ReadLine(fin);

  std::string secondline = ReadLine(fin);
  numtotparticles = lexical_cast<int>(secondline.substr(0, secondline.find(',', 0)));
  numattributes = lexical_cast<int>(secondline.substr(secondline.find(' ', 0) + 1));

  if (numattributes < 0)
    {
      numparams = lexical_cast<int>(ReadLine(fin));
      for (int i = 0; i < numparams; i++)
        {
          std::string paramlabel = ReadLine(fin);
          std::string paramvalue = ReadLine(fin);

          if (paramlabel == "FLOAT")
            {
              datawidthLoad = sizeOfXDRfloat;
            }
          else if (paramlabel == "DOUBLE")
            {
              datawidthLoad = sizeOfXDRdouble;
            }
          else if (paramlabel == "XDRD")
            {
            }
          else if (paramlabel == "DIMENSION")
            {
            }
          else
            {
              MemoryManager.SaveParameter(paramlabel, lexical_cast<value_type>(paramvalue), true);
            }
        }
      ;
      numattributes = -numattributes;
    }

  numlocparticles = (lrint(((double)(RANK + 1) / SIZE) * numtotparticles) - lrint(((double)(RANK) / SIZE) * numtotparticles));
  buffsize = numlocparticles * datawidthLoad * numattributes;
  dataoffset = lrint(((double)RANK / SIZE) * numtotparticles) * datawidthLoad * numattributes;
  attributes.resize(numattributes);

  for (int i = 0; i < numattributes; i++)
    {
      std::string attrname = ReadLine(fin);
      attributes[i] = VM.getIndex(attrname);
      if (attributes[i] > attrmaxindex)
        {
          attrmaxindex = attributes[i];
        }
    }

  if (Data.size2() == 0)
    {
      Data.resize(numlocparticles, attrmaxindex + 1);
    }
  else
    {
      Data.resize(numlocparticles, Data.size2());
    }

  if (ReadLine(fin) != "^^^")
    {
      std::cerr << "File does not contain magic separator ^^^, no data loaded\n";
      attributes.resize(0);
      return attributes;
    }

  filepos = fin->tellp();
  filepos += dataoffset;

  for (int i = 0; i < numattributes; i++)
    {
      if (attributes[i] > (int)Data.size2())
        {
          std::cerr << "Data container too small for attribute: " << VM.getName(attributes[i]) << "\n";
          attributes.resize(0);
          return attributes;
        }
    }

  fin->seekg(filepos);
  char *buff = new char[buffsize];
  XDR *xdrs = new XDR;
  xdrmem_create(xdrs, buff, buffsize, XDR_DECODE);

  fin->read(buff, buffsize);

  if (datawidthLoad == sizeOfXDRfloat)
    {
      float *valuep = new float;
      for (int i = 0; i < numlocparticles; i++)
        {
          for (int j = 0; j < numattributes; j++)
            {
              xdr_float(xdrs, valuep);
              if (attributes[j] != -1)
                {
                  Data(i, attributes[j]) = (value_type) * valuep;
                }
            }
        }
    }
  else if (datawidthLoad == sizeOfXDRdouble)
    {
      double *valuep = new double;
      for (int i = 0; i < numlocparticles; i++)
        {
          for (int j = 0; j < numattributes; j++)
            {
              xdr_double(xdrs, valuep);
              if (attributes[j] != -1)
                {
                  Data(i, attributes[j]) = (value_type) * valuep;
                }
            }
        }
    }
  delete buff;
  delete fin;
  delete xdrs;

  // Erase attributes which where not save to Data
  attributes.erase(remove(attributes.begin(), attributes.end(), -1), attributes.end());

  return attributes;
}

template <typename SimTrait>
bool IOManager<SimTrait>::SaveCDAT(std::string FileName, std::vector<int> Attributes)
{
  matrix_reference Data(MemoryManager.Data);
  string_reference Name(MemoryManager.Name);

  using boost::lexical_cast;

#ifndef OOSPH_STANDALONE
  const size_t RANK = MPI::COMM_WORLD.Get_rank();
  const size_t SIZE = MPI::COMM_WORLD.Get_size();
#else
  const size_t RANK = 0;
  const size_t SIZE = 1;
#endif
  int numtotparticles, numlocparticles, numremainparticles[SIZE];
  size_t datasize, buffsize, dataoffset, numattributes;
  std::vector<int> attributes = Attributes;
  std::string header;
  StringSetType ParamSet = MemoryManager.DumpParameter();
  StringSetType::iterator ParamSetIter;

  numattributes = attributes.size();
  numlocparticles = Data.size1();
  buffsize = numlocparticles * datawidthSave * numattributes;

#ifndef OOSPH_STANDALONE
  MPI::COMM_WORLD.Allgather(&numlocparticles, 1, MPI::INT, &numremainparticles, 1, MPI::INT);
#endif

  for (size_t i = 1; i < SIZE; i++)
    {
      numremainparticles[i] += numremainparticles[i - 1];
    }

  if (RANK == 0)
    {
      dataoffset = 0;
    }
  else
    {
      dataoffset = numremainparticles[RANK - 1] * datawidthSave * numattributes;
    }

#ifndef OOSPH_STANDALONE
  numtotparticles = numremainparticles[SIZE - 1];
#else
  numtotparticles = numlocparticles;
#endif

  datasize = numtotparticles * datawidthSave * numattributes;

  header += Name;
  header += "\n";
  header += lexical_cast<std::string>(numtotparticles);
  header += ", -";
  header += lexical_cast<std::string>(numattributes);
  header += "\n";
  header += lexical_cast<std::string>(ParamSet.size() + 3);
  header += "\nDIMENSION\n3\nXDRD\n0\n";
  if (datawidthSave == sizeOfXDRfloat)
    {
      header += "FLOAT\n0\n";
    }
  else
    {
      header += "DOUBLE\n0\n";
    }

  for (ParamSetIter = ParamSet.begin(); ParamSetIter != ParamSet.end(); ++ParamSetIter)
    {
      header += *ParamSetIter;
      header += "\n";
      header += lexical_cast<std::string>(MemoryManager.LoadParameter(*ParamSetIter));
      header += "\n";
    }

  for (size_t i = 0; i < numattributes; i++)
    {
      header += VM.getName(attributes[i]);
      header += "\n";
    }

  header += "^^^\n";

  std::fstream fout;
  std::fstream::pos_type filepos;
  fout.open(FileName.c_str(), std::ios::out | std::ios::binary | std::ios::ate);

  if (!fout) return false;

  // Write an empty file first, so we don't get any trouble with NFS caches
  fout.seekp(header.size() + datasize - 1);
  fout << "\n";
  fout.seekp(0);
  fout.flush();

#ifndef OOSPH_STANDALONE
  MPI::COMM_WORLD.Barrier();
#endif

  if (RANK == 0)
    {
      fout << header;
    }

  filepos = header.size();
  filepos += dataoffset;

  for (size_t i = 0; i < numattributes; i++)
    {
      int attribute = attributes[i];
      if (attribute + 1 > (int)Data.size2())
        {
          std::cerr << "Data array too small! Attribute not found: " << VM.getName(attribute) << "\n";
          return false;
        }
    }
  char *buff = new char[buffsize];
  XDR *xdrs = new XDR;
  xdrmem_create(xdrs, buff, buffsize, XDR_ENCODE);
  if (datawidthSave == sizeOfXDRfloat)
    {
      float *valuep = new float;
      for (int i = 0; i < numlocparticles; i++)
        {
          for (size_t j = 0; j < numattributes; j++)
            {
              *valuep = Data(i, attributes[j]);
              xdr_float(xdrs, valuep);
            }
        }
      delete valuep;
    }
  else           // In every other case, save as doubles
    {
      double *valuep = new double;
      for (int i = 0; i < numlocparticles; i++)
        {
          for (size_t j = 0; j < numattributes; j++)
            {
              *valuep = Data(i, attributes[j]);
              xdr_double(xdrs, valuep);
            }
        }
      delete valuep;
    }
  fout.seekp(filepos);
  fout.write(buff, buffsize);
  fout.close();
  delete buff;

#ifndef OOSPH_STANDALONE
  MPI::COMM_WORLD.Barrier();
#endif

  return true;
}

template <typename SimTrait>
void IOManager<SimTrait>::SetSinglePrecOut(void)
{
  datawidthSave = sizeOfXDRfloat;
}

template <typename SimTrait>
void IOManager<SimTrait>::SetDoublePrecOut(void)
{
  datawidthSave = sizeOfXDRdouble;
}

template <typename SimTrait>
std::string IOManager<SimTrait>::ReadLine(std::fstream *_fin)
{
  const size_t _HEADBUFFSIZE = 512;
  std::string _hdrline;
  char *_buff = new char[_HEADBUFFSIZE];

  _fin->getline(_buff, _HEADBUFFSIZE);
  _hdrline.assign(_buff, 0, _fin->gcount() - 1);       // Removes LF, doesn't work with DOS files (CRLF)
  delete _buff;
  
  return _hdrline;
};
};


#endif
