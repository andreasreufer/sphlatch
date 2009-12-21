#ifndef SPHLATCH_IOMANAGER_H
#define SPHLATCH_IOMANAGER_H

#include <string>
#include <fstream>

#include <rpc/rpc.h>
#include <boost/lexical_cast.hpp>

#include "typedefs.h"
#include "memorymanager.h"
#include "particle.h"

#ifdef SPHLATCH_PARALLEL
#include "communicationmanager.h"
#endif

namespace sphlatch
{
///
/// \brief I/O Manager for reading and writing CDAT files
///
/// \author Andreas Reufer andreas.reufer@space.unibe.ch
///

class IOManager
{
public:

typedef IOManager self_type;
typedef IOManager& self_reference;
typedef IOManager* self_pointer;

typedef sphlatch::ParticleVarMap varMapType;
typedef sphlatch::MemoryManager memoryManagerType;
#ifdef SPHLATCH_PARALLEL
typedef sphlatch::CommunicationManager commManagerType;
#endif

typedef std::set<std::string> StringSetType;

static self_reference instance();

std::vector<int> loadCDAT(std::string _fileName);
bool saveCDAT(std::string _fileName, std::vector<int> _attributes);
void setSinglePrecOut(void);
void setDoublePrecOut(void);

varMapType& VM;
memoryManagerType& MemoryManager;
#ifdef SPHLATCH_PARALLEL
commManagerType& commManager;
#endif

protected:
IOManager(void);
~IOManager(void);

private:
std::string readLine(std::fstream *_fin);
static self_pointer _instance;

size_t sizeOfXDRfloat, sizeOfXDRdouble, datawidthSave;
};

IOManager::self_pointer IOManager::_instance = NULL;

IOManager::self_reference IOManager::instance(void)
{
  if (_instance != NULL)
    return *_instance;
  else
    {
      _instance = new IOManager;
      return *_instance;
    }
}

IOManager::IOManager(void) :
  VM(varMapType::instance()),
  MemoryManager(memoryManagerType::instance())
#ifdef SPHLATCH_PARALLEL
  ,commManager(commManagerType::instance())
#endif
{
  sizeOfXDRfloat = 4;
  sizeOfXDRdouble = 8;

  if (sizeof(fType) == sizeOfXDRfloat)
    {
      datawidthSave = sizeOfXDRfloat;
    }
  else
    {
      datawidthSave = sizeOfXDRdouble;
    }
}

IOManager::~IOManager(void)
{
}

std::vector<int> IOManager::loadCDAT(std::string FileName)
{
  matrixRefType Data(IOManager::MemoryManager.Data);
  std::string&  Name(MemoryManager.simName);

  using boost::lexical_cast;

#ifdef SPHLATCH_PARALLEL
  const size_t myDomain = commManager.getMyDomain();
  const size_t SIZE = commManager.getNoDomains();
#else
  const size_t myDomain = 0;
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

  Name = readLine(fin);

  std::string secondline = readLine(fin);
  numtotparticles = lexical_cast<int>(secondline.substr(0, secondline.find(',', 0)));
  numattributes = lexical_cast<int>(secondline.substr(secondline.find(' ', 0) + 1));

  if (numattributes < 0)
    {
      numparams = lexical_cast<int>(readLine(fin));
      for (int i = 0; i < numparams; i++)
        {
          std::string paramlabel = readLine(fin);
          std::string paramvalue = readLine(fin);

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
              MemoryManager.saveParameter(paramlabel, lexical_cast<fType>(paramvalue), true);
            }
        }
      ;
      numattributes = -numattributes;
    }

  numlocparticles = (lrint(((double)(myDomain + 1) / SIZE) * numtotparticles) - lrint(((double)(myDomain) / SIZE) * numtotparticles));
  buffsize = numlocparticles * datawidthLoad * numattributes;
  dataoffset = lrint(((double)myDomain / SIZE) * numtotparticles) * datawidthLoad * numattributes;
  attributes.resize(numattributes);

  for (int i = 0; i < numattributes; i++)
    {
      std::string attrname = readLine(fin);
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

  if (readLine(fin) != "^^^")
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
                  Data(i, attributes[j]) =
                    static_cast<fType>(*valuep);
                }
            }
        }
      delete valuep;
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
                  Data(i, attributes[j]) =
                    static_cast<fType>(*valuep);
                }
            }
        }
      delete valuep;
    }
  delete buff;
  delete fin;
  delete xdrs;

  // Erase attributes which where not save to Data
  attributes.erase(remove(attributes.begin(), attributes.end(), -1), attributes.end());

  return attributes;
}

bool IOManager::saveCDAT(std::string FileName, std::vector<int> Attributes)
{
  matrixRefType Data(IOManager::MemoryManager.Data);
  std::string&  Name(MemoryManager.simName);

  using boost::lexical_cast;

#ifdef SPHLATCH_PARALLEL
  const size_t myDomain = commManager.getMyDomain();
  const size_t SIZE = commManager.getNoDomains();
#else
  const size_t myDomain = 0;
  const size_t SIZE = 1;
#endif
  int numtotparticles, numlocparticles;
  size_t datasize, buffsize, dataoffset, numattributes;
  std::vector<int> attributes = Attributes;
  std::string header;
  StringSetType ParamSet = MemoryManager.dumpParameter();
  StringSetType::iterator ParamSetIter;

  numattributes = attributes.size();
  numlocparticles = Data.size1();
  buffsize = numlocparticles * datawidthSave * numattributes;

  countsVectType partCounts;
  partCounts.resize(SIZE);

  partCounts[myDomain] = Data.size1();

#ifdef SPHLATCH_PARALLEL
  commManager.sumUpCounts(partCounts);
#endif

  for (size_t i = 1; i < SIZE; i++)
    {
      partCounts[i] += partCounts[i - 1];
    }

  if (myDomain == 0)
    {
      dataoffset = 0;
    }
  else
    {
      dataoffset = partCounts[myDomain - 1] * datawidthSave * numattributes;
    }


#ifdef SPHLATCH_PARALLEL
  numtotparticles = partCounts[ SIZE - 1 ];
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
      header += lexical_cast<std::string>(MemoryManager.loadParameter(*ParamSetIter));
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

#ifdef SPHLATCH_PARALLEL
  commManager.barrier();
#endif

  if (myDomain == 0)
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

#ifdef SPHLATCH_PARALLEL
  commManager.barrier();
#endif

  return true;
}

void IOManager::setSinglePrecOut(void)
{
  datawidthSave = sizeOfXDRfloat;
}

void IOManager::setDoublePrecOut(void)
{
  datawidthSave = sizeOfXDRdouble;
}

std::string IOManager::readLine(std::fstream *_fin)
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
