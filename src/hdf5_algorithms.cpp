#ifndef SPHLATCH_IOMANAGER_H
#define SPHLATCH_IOMANAGER_H

///
/// we use HDF5 >= 1.8 functions
///
#define H5_NO_DEPRECATED_SYMBOLS

#include <hdf5.h>

#include "typedefs.h"

namespace sphlatch
{
///
/// \brief I/O Manager for reading and writing particle dumps with HDF5
///
/// \author Andreas Reufer andreas.reufer@space.unibe.ch
///

class HDF5Simple {
  public:
    HDF5Simple() {};
    ~HDF5Simple() {};

};

void HDF5Simple::savePrimitive(matrixRefType _matr,
                              std::string _name,
                              std::string _outputFile)
{
  hid_t fileHandle = getLocFilehandleRW(_outputFile);
  hid_t accessPropList = H5Pcreate(H5P_DATASET_XFER);
  hid_t rootGroup = H5Gopen(fileHandle, "/", H5P_DEFAULT);

  saveVar(_matr, _name, rootGroup,
          _matr.size1(), _matr.size1(), 0,
          accessPropList);

  H5Gclose(rootGroup);
  H5Pclose(accessPropList);
  H5Fclose(fileHandle);
}

void HDF5Simple::savePrimitive(valvectRefType _valvect,
                              std::string _name,
                              std::string _outputFile)
{
  hid_t fileHandle = getLocFilehandleRW(_outputFile);
  hid_t accessPropList = H5Pcreate(H5P_DATASET_XFER);
  hid_t rootGroup = H5Gopen(fileHandle, "/", H5P_DEFAULT);

  saveVar(_valvect, _name, rootGroup,
          _valvect.size(), _valvect.size(), 0,
          accessPropList);

  H5Gclose(rootGroup);
  H5Pclose(accessPropList);
  H5Fclose(fileHandle);
}

void HDF5Simple::savePrimitive(idvectRefType _idvect,
                              std::string _name,
                              std::string _outputFile)
{
  hid_t fileHandle = getLocFilehandleRW(_outputFile);
  hid_t accessPropList = H5Pcreate(H5P_DATASET_XFER);
  hid_t rootGroup = H5Gopen(fileHandle, "/", H5P_DEFAULT);

  saveVar(_idvect, _name, rootGroup,
          _idvect.size(), _idvect.size(), 0,
          accessPropList);

  H5Gclose(rootGroup);
  H5Pclose(accessPropList);
  H5Fclose(fileHandle);
}

void HDF5Simple::loadPrimitive(valvectRefType _vect,
                              std::string _name,
                              std::string _inputFile)
{
  hsize_t dimsMem[2], dimsFile[2], offset[2];

  offset[0] = 0;
  offset[1] = 0;

  ///
  /// open file, group and dataset
  ///
  hid_t fileHandle = getLocFilehandleRW(_inputFile);
  hid_t accessPropList = H5Pcreate(H5P_DATASET_XFER);
  hid_t rootGroup = H5Gopen(fileHandle, "/", H5P_DEFAULT);
  hid_t dataset = H5Dopen(rootGroup, _name.c_str(), H5P_DEFAULT);

  hid_t fileSpace = H5Dget_space(dataset);
  H5Sget_simple_extent_dims(fileSpace, dimsFile, NULL);

  ///
  /// prepare vector its memspace
  ///
  _vect.resize(dimsFile[0]);
  dimsMem[0] = dimsFile[0];
  dimsMem[1] = 1;
  hid_t memSpace = H5Screate_simple(1, dimsMem, NULL);

  ///
  /// select hyperslab in file and read it
  ///
  H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
                      offset, NULL, dimsMem, NULL);

  H5Dread(dataset, H5floatMemType, memSpace,
          fileSpace, accessPropList, &_vect(0));

  ///
  /// housekeeping
  ///
  H5Sclose(memSpace);
  H5Sclose(fileSpace);
  H5Dclose(dataset);
  H5Gclose(rootGroup);
  H5Pclose(accessPropList);
  H5Fclose(fileHandle);
}

fType HDF5Simple::loadAttribute(std::string _name, std::string _inputFile)
{
  hsize_t dimsMem[1];

  dimsMem[0] = 1;
  fType attrBuff;

  ///
  /// open file, group and dataset
  ///
  hid_t fileHandle = getLocFilehandleRW(_inputFile);
  hid_t accessPropList = H5Pcreate(H5P_DATASET_XFER);
  hid_t rootGroup = H5Gopen(fileHandle, "/", H5P_DEFAULT);

  hid_t curAttr = H5Aopen(rootGroup, _name.c_str(), H5P_DEFAULT);
  H5Aread(curAttr, H5floatMemType, &attrBuff);
  H5Aclose(curAttr);

  ///
  /// housekeeping
  ///
  H5Gclose(rootGroup);
  H5Pclose(accessPropList);
  H5Fclose(fileHandle);

  return attrBuff;
}

void HDF5Simple::saveAttribute(const fType _attr, std::string _name,
                         std::string _inputFile)
{
  hsize_t dimsMem[1];
  dimsMem[0] = 1;

  ///
  /// open file, group
  ///
  hid_t fileHandle = getLocFilehandleRW(_inputFile);
  hid_t accessPropList = H5Pcreate(H5P_DATASET_XFER);
  hid_t rootGroup = H5Gopen(fileHandle, "/", H5P_DEFAULT);
  hid_t memspace = H5Screate_simple(1, dimsMem, NULL);

  ///
  /// if attribute already exists, delete it
  ///
  if ( H5Aexists(rootGroup, _name.c_str()) )
  {
    H5Adelete(rootGroup, _name.c_str());
  }

  ///
  /// write attribute
  ///
  hid_t curAttr = H5Acreate_by_name(rootGroup, ".", _name.c_str(),
                                    H5floatFileType, memspace,
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Awrite(curAttr, H5floatMemType, &_attr);

  ///
  /// housekeeping
  ///
  H5Aclose(curAttr);

  H5Sclose(memspace);
  H5Gclose(rootGroup);
  H5Pclose(accessPropList);
  H5Fclose(fileHandle);
}


void HDF5Simple::saveVar(std::string _name, 
                        hid_t _group,
                        hid_t _writeProplist,
                        hid_t _mtype,
                        hid_t _ftype,
                        const size_t _n, 
                        const size_t _m, 
                        const void* _data)
{
  hsize_t dimsMem[2], dimsFile[2], offset[2];

  size_t noDims = 1;
  if (_m > 1 )
    noDims = 2;

  /// dataset dimensions in memory
  dimsMem[0] = _n;
  dimsMem[1] = _m;
  hid_t memspace = H5Screate_simple(noDims, dimsMem, NULL);
  
  /// file dataset dimensions
  dimsFile[0] = _n;
  dimsFile[1] = _m;
  hid_t filespace = H5Screate_simple(noDims, dimsFile, NULL);

  /// offset
  offset[0] = 0;
  offset[1] = 0;

  /// check if dataset already exists in the right size
  hid_t curDataset;
  if (objectExist(_group, _name))
    {
      curDataset = H5Dopen(_group, _name.c_str(), H5P_DEFAULT);
    }
  else
    {
      curDataset = H5Dcreate(_group, _name.c_str(), _ftype,
                             filespace, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    }

  /// select hyperslab
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                      offset, NULL, dimsMem, NULL);

  /// write the data
  H5Dwrite(curDataset, _mtype, memspace,
           filespace, _writeProplist, _data);

  H5Dclose(curDataset);
  H5Sclose(memspace);
  H5Sclose(filespace);
}


hid_t HDF5Simple::getLocFilehandleRW(std::string _outputFile)
{
  hid_t filePropList = H5Pcreate(H5P_FILE_ACCESS);

  struct stat statBuff;
  hid_t fileHandle;

  if (stat(_outputFile.c_str(), &statBuff) == -1)
    {
      fileHandle = H5Fcreate(_outputFile.c_str(),
                             H5F_ACC_TRUNC, H5P_DEFAULT, filePropList);
    }
  else
    {
      fileHandle = H5Fopen(_outputFile.c_str(), H5F_ACC_RDWR, filePropList);
    }
  H5Pclose(filePropList);
  return fileHandle;
}


/*HDF5Simple::stringListType
HDF5Simple::discoverVars(std::string _inputFile,
                        std::string _stepName)
{
  hid_t fileHandle, filePropList;

  ///
  /// get a new HDF5 file handle
  ///
  filePropList = H5Pcreate(H5P_FILE_ACCESS);
#ifdef SPHLATCH_PARALLEL
  H5Pset_fapl_mpio(filePropList, MPI::COMM_WORLD, MPI::INFO_NULL);
#endif
  fileHandle = H5Fopen(_inputFile.c_str(), H5F_ACC_RDONLY, filePropList);
  H5Pclose(filePropList);

  hid_t curGroup = H5Gopen(fileHandle, _stepName.c_str(), H5P_DEFAULT);

  datasetsInFile.clear();
  H5Lvisit_by_name(fileHandle, _stepName.c_str(), H5_INDEX_CRT_ORDER,
                   H5_ITER_NATIVE, getObsCallback, NULL, H5P_DEFAULT);

  H5Gclose(curGroup);
  H5Fclose(fileHandle);

  return datasetsInFile;
}

HDF5Simple::stringListType
HDF5Simple::discoverVars(std::string _inputFile)
{
  return discoverVars(_inputFile, "/current");
}*/



};
#endif
