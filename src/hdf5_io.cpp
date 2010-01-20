#ifndef SPHLATCH_HDF5_IO_H
#define SPHLATCH_HDF5_IO_H

///
/// we use HDF5 >= 1.8 functions
///
#define H5_NO_DEPRECATED_SYMBOLS

#include <hdf5.h>
#include <sys/stat.h>

#include "typedefs.h"

namespace sphlatch
{
///
/// \brief I/O Manager for reading and writing particle dumps with HDF5
///
/// \author Andreas Reufer andreas.reufer@space.unibe.ch
///

class HDF5File {
public:
   HDF5File(std::string _file);
   ~HDF5File();

   void getDims(std::string _name, size_t& _nx, size_t& _ny);
   void loadPrimitive(std::string _name, valvectType& _v);
   void savePrimitive(std::string _name, valvectType& _v);

   void singlePrecOut();
   void doublePrecOut();

private:
   void saveRaw(std::string  _name,
                hid_t        _group,
                hid_t        _mtype,
                hid_t        _ftype,
                const size_t _n,
                const size_t _m,
                const void   * _data);

   void loadRaw(std::string  _name,
                hid_t        _group,
                hid_t        _mtype,
                void         * _data);

   bool objExist(hid_t _fh, std::string _op);

   hid_t fpl, fh, apl, rgrp;

   hid_t h5mITYPE, h5mFTYPE, h5fITYPE, h5fFTYPE;

   /*  loadAttribute(std::string _name, std::string _inputFile);
    */
};

HDF5File::HDF5File(std::string _file)
{
   fpl = H5Pcreate(H5P_FILE_ACCESS);

   struct stat statBuff;
   if (stat(_file.c_str(), &statBuff) == -1)
      fh = H5Fcreate(_file.c_str(),
                     H5F_ACC_TRUNC, H5P_DEFAULT, fpl);
   else
      fh = H5Fopen(_file.c_str(), H5F_ACC_RDWR, fpl);
   H5Pclose(fpl);

   apl  = H5Pcreate(H5P_DATASET_XFER);
   rgrp = H5Gopen(fh, "/", H5P_DEFAULT);

   // determine datatypes
   if (sizeof(fType) == 8)
      h5mFTYPE = H5T_NATIVE_DOUBLE;
   else
      h5mFTYPE = H5T_NATIVE_FLOAT;

   h5mITYPE = H5T_NATIVE_INT;

   h5fFTYPE = H5T_IEEE_F32LE;
   h5fITYPE = H5T_STD_I32LE;
}

HDF5File::~HDF5File()
{
   H5Gclose(rgrp);
   H5Pclose(apl);
   H5Fclose(fh);
}

void HDF5File::loadPrimitive(std::string _name, valvectType& _v)
{
  loadRaw(_name, rgrp, h5mFTYPE, &_v(0));
}

void HDF5File::savePrimitive(std::string _name, valvectType& _v)
{
  saveRaw(_name, rgrp, h5mFTYPE, h5fFTYPE, _v.size(), 1, &_v(0));
}
   
void HDF5File::saveRaw(std::string  _name,
                       hid_t        _group,
                       hid_t        _mtype,
                       hid_t        _ftype,
                       const size_t _n,
                       const size_t _m,
                       const void   * _data)
{
   hsize_t dimsm[2], dimsf[2], offs[2];

   size_t nodims = 1;
   if (_m > 1)
      nodims = 2;

   /// dataset dimensions in memory
   dimsm[0] = _n;
   dimsm[1] = _m;
   hid_t mspc = H5Screate_simple(nodims, dimsm, NULL);

   /// file dataset dimensions
   dimsf[0] = _n;
   dimsf[1] = _m;
   hid_t fspc = H5Screate_simple(nodims, dimsf, NULL);

   /// offset
   offs[0] = 0;
   offs[1] = 0;

   /// check if dataset already exists in the right size
   hid_t curds;
   if (objExist(_group, _name))
      curds = H5Dopen(_group, _name.c_str(), H5P_DEFAULT);
   else
      curds = H5Dcreate(_group, _name.c_str(), _ftype, fspc, H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);

   /// select hyperslab
   H5Sselect_hyperslab(fspc, H5S_SELECT_SET, offs, NULL, dimsm, NULL);

   /// write the data
   H5Dwrite(curds, _mtype, mspc, fspc, H5P_DEFAULT, _data);

   H5Dclose(curds);
   H5Sclose(mspc);
   H5Sclose(fspc);
}

void HDF5File::loadRaw(std::string  _name,
                       hid_t        _group,
                       hid_t        _mtype,
                       void         * _data)
{
   hsize_t dimsm[2], dimsf[2], offs[2];

   /// offset
   offs[0] = 0;
   offs[1] = 0;
   
   dimsf[0] = 0;
   dimsf[1] = 0;

   // open dataset and create file and mem spaces
   hid_t cuds = H5Dopen(_group, _name.c_str(), H5P_DEFAULT);
   hid_t fspc = H5Dget_space(cuds);
   H5Sget_simple_extent_dims(fspc, dimsf, NULL);
   const size_t nodims = H5Sget_simple_extent_ndims(fspc);

   dimsm[0] = dimsf[0];
   dimsm[1] = dimsf[1];
   hid_t mspc = H5Screate_simple(nodims, dimsm, NULL);

   // select hyperslab in file
   H5Sselect_hyperslab(fspc, H5S_SELECT_SET, offs, NULL, dimsm, NULL);

   // read the dataset
   H5Dread(cuds, _mtype, mspc, fspc, H5P_DEFAULT, _data);

   // housekeeping
   H5Sclose(mspc);
   H5Sclose(fspc);
   H5Dclose(cuds);
}

bool HDF5File::objExist(hid_t _fh, std::string _op)
{
   H5Eset_auto(H5E_DEFAULT, NULL, NULL);

   H5O_info_t oi;
   herr_t     stat;

   stat = H5Oget_info_by_name(_fh, _op.c_str(), &oi, H5P_DEFAULT);
   if (stat != 0)
      return(false);
   else
      return(true);
}

void HDF5File::getDims(std::string _name, size_t& _nx, size_t& _ny)
{
   hid_t curds = H5Dopen(rgrp, _name.c_str(), H5P_DEFAULT);
   hid_t fspc  = H5Dget_space(curds);

   hsize_t dimsf[2];

   dimsf[0] = 0;
   dimsf[1] = 0;

   H5Sget_simple_extent_dims(fspc, dimsf, NULL);

   _nx = dimsf[0];
   _ny = dimsf[1];

   H5Sclose(fspc);
   H5Dclose(curds);
}

void HDF5File::singlePrecOut()
{
  h5fFTYPE = H5T_IEEE_F32LE;
}

void HDF5File::doublePrecOut()
{
  h5fFTYPE = H5T_IEEE_F64LE;
}

/*
   void HDF5File::savePrimitive(idvectRefType _idvect,
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

   void HDF5File::getDims(std::string _name, std::string _inputFile)
   {
   hid_t fileHandle = getLocFilehandleRW(_inputFile);
   hid_t accessPropList = H5Pcreate(H5P_DATASET_XFER);
   hid_t rootGroup = H5Gopen(fileHandle, "/", H5P_DEFAULT);
   hid_t dataset = H5Dopen(rootGroup, _name.c_str(), H5P_DEFAULT);

   hid_t fileSpace = H5Dget_space(dataset);
   H5Sget_simple_extent_dims(fileSpace, dimsFile, NULL);

   }

   void HDF5File::loadPrimitive(valvectRefType _vect,
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
   /// prepare memspace
   ///
   dimsMem = dimsFile;
   hid_t memSpace = H5Screate_simple(2, dimsMem, NULL);

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

   fType HDF5File::loadAttribute(std::string _name, std::string _inputFile)
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

   void HDF5File::saveAttribute(const fType _attr, std::string _name,
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


   void HDF5File::saveVar(std::string _name,
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


   hid_t HDF5File::getLocFilehandleRW(std::string _outputFile)
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
   }*/


/*HDF5File::stringListType
   HDF5File::discoverVars(std::string _inputFile,
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

   HDF5File::stringListType
   HDF5File::discoverVars(std::string _inputFile)
   {
   return discoverVars(_inputFile, "/current");
   }*/
};
#endif
