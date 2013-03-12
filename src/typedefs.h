#ifndef SPHLATCH_TYPEDEFS
#define SPHLATCH_TYPEDEFS

/*
 *  typedefs.h
 *
 *
 *  Created by Andreas Reufer on 15.11.07.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include <vector>
#include <map>
#include <set>
#include <list>

#include <boost/dynamic_bitset.hpp>

#include <cmath>
#include <limits>
#include <valarray>

#define BZ_THREADSAFE
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

#include <boost/numeric/ublas/vector.hpp>

namespace sphlatch {
///
/// the BLAS namespace
///
//namespace blas = boost::numeric::ublas;

///
/// fType: should be a float
///
#ifdef SPHLATCH_SINGLEPREC
typedef float                         fType;
#else
typedef double                        fType;
#endif
typedef fType*                        fPtrType;
typedef fType&                        fRefType;

const fType fTypeInf = std::numeric_limits<fType>::infinity();
const fType fTypeNan = std::numeric_limits<fType>::quiet_NaN();

typedef blitz::TinyVector<fType, 3>   vect3dT;
//typedef blitz::Array<fType, 1>        farrT;
enum dims3
{
   X, Y, Z
};

struct box3dT
{
   vect3dT cen;
   fType   size;

   box3dT& operator=(const box3dT& _rhs)
   {
      cen  = _rhs.cen;
      size = _rhs.size;

      return(* this);
   }

   box3dT operator*(const fType mul)
   {
      box3dT tmp;

      tmp       = *this;
      tmp.size *= mul;
      return(tmp);
   }

   box3dT operator+(const box3dT _rhs)
   {
     const fType rhsize = 0.5*_rhs.size;
     const fType lhsize = 0.5*(this->size);

     const fType xmin = _rhs.cen[0] - rhsize < this->cen[0] - lhsize ?
       _rhs.cen[0] - rhsize : this->cen[0] - lhsize;
     const fType ymin = _rhs.cen[1] - rhsize < this->cen[1] - lhsize ?
       _rhs.cen[1] - rhsize : this->cen[1] - lhsize;
     const fType zmin = _rhs.cen[2] - rhsize < this->cen[2] - lhsize ?
       _rhs.cen[2] - rhsize : this->cen[2] - lhsize;
     
     const fType xmax = _rhs.cen[0] + rhsize > this->cen[0] + lhsize ?
       _rhs.cen[0] + rhsize : this->cen[0] + lhsize;
     const fType ymax = _rhs.cen[1] + rhsize > this->cen[1] + lhsize ?
       _rhs.cen[1] + rhsize : this->cen[1] + lhsize;
     const fType zmax = _rhs.cen[2] + rhsize > this->cen[2] + lhsize ?
       _rhs.cen[2] + rhsize : this->cen[2] + lhsize;

     box3dT nbox;
     nbox.cen = 0.5 * (xmax + xmin), 0.5 * (ymax + ymin), 0.5 * (zmax + zmin);
     nbox.size = std::max(xmax - xmin, std::max(ymax - ymin, zmax - zmin));
     return nbox;
   }
};


///
/// type for particle or tree node IDs
///
typedef int                     idType;
typedef int                     iType;

///
/// a vector of (particle) indices
///
typedef std::vector<size_t>     partsIndexVectType;
typedef std::vector<size_t>&    partsIndexVectRefType;
typedef std::vector<size_t> *   partsIndexVectPtrType;

typedef std::list<size_t>       partsIndexListT;
typedef std::list<size_t>&      partsIndexListRefT;
typedef std::list<size_t> *     partsIndexListPtrT;

///
/// a vector of counts
/// countsType has to be compatible with MPI_INT
///
typedef int                         countsType;
typedef int                         cType;
typedef countsType&                 countsRefType;
typedef countsType*                 countsPtrType;
typedef std::vector<countsType>     countsVectType;
typedef std::vector<countsType>     cVType;
typedef std::vector<countsType>&    countsVectRefType;
typedef std::vector<countsType> *   countsVectPtrType;

typedef std::vector<countsType>     ivectT;
typedef std::vector<countsType>&    ivectRefT;
typedef std::vector<countsType> *   ivectPtrT;

///
/// a vector of particle indices vectors
/// used to map lists of particles to a list of domains
///
typedef std::vector<partsIndexVectType>     domainPartsIndexType;
typedef std::vector<partsIndexVectType>&    domainPartsIndexRefType;
typedef std::vector<partsIndexVectType> *   domainPartsIndexPtrType;

///
/// storage type for bitset
/// the longer, the better
///
typedef unsigned long int                         bitsetBlockType;

///
/// a bit set, useful as a boolean vector
///
typedef boost::dynamic_bitset<bitsetBlockType>    bitsetType;
typedef boost::dynamic_bitset<bitsetBlockType>&   bitsetRefType;

///
/// attribute map
///
typedef std::map<std::string, fType>              attrMapType;
typedef std::map<std::string, fType>              attrMT;
typedef attrMapType&                              attrMapRefType;
typedef attrMapType*                              attrMapPtrType;


///
/// a vector of strings
///
typedef std::vector<std::string>   stringVectType;

///
/// a vector of strings
///
typedef std::set<std::string>      stringSetType;

///
/// a list of strings
///
typedef std::list<std::string>     stringListType;

///
/// vectors of basic data types
///
namespace blas = boost::numeric::ublas;
typedef blas::vector<fType>    fvectT;
typedef blas::vector<idType>   idvectT;
typedef blas::matrix<fType>    fmatrT;
typedef blas::matrix<idType>   idmatrT;
};

#define LEGACY
#ifdef LEGACY
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/matrix_expression.hpp>
 #include <boost/numeric/ublas/matrix_proxy.hpp>
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/vector_expression.hpp>
 #include <boost/numeric/ublas/vector_proxy.hpp>
 #include <boost/numeric/ublas/io.hpp>

namespace sphlatch {
namespace blas = boost::numeric::ublas;

typedef blas::matrix<fType>                     matrixType;
typedef matrixType&                             matrixRefType;
typedef matrixType*                             matrixPtrType;

typedef blas::zero_matrix<fType>                zeromatrixType;

typedef blas::matrix_row<matrixType>            matrixRowType;
typedef blas::matrix_row<matrixType>&           matrixRowRefType;
typedef blas::matrix_row<matrixType> *          matrixRowPtrType;

typedef blas::matrix_column<matrixType>         matrixColumnType;

typedef blas::matrix_range<matrixType>          matrixRangeType;

typedef blas::range                             rangeType;
typedef blas::slice                             sliceType;

typedef blas::matrix_vector_slice<matrixType>   matrixVectorSliceType;
typedef blas::matrix_vector_slice<const matrixType>
constMatrixVectorSliceType;

typedef matrixRowType                           particleRowType;
typedef matrixRowType&                          particleRowRefType;
typedef matrixRowType*                          particleRowPtrType;

typedef matrixColumnType                        quantColumnType;
typedef matrixColumnType&                       quantColumnRefType;
typedef matrixColumnType*                       quantColumnPtrType;

typedef blas::vector<fType>                     valvectType;
typedef blas::vector<fType>&                    valvectRefType;
typedef blas::vector<fType> *                   valvectPtrType;

typedef blas::vector<fType>                     fvectT;
typedef blas::vector<fType>&                    fvectRefT;
typedef blas::vector<fType> *                   fvectPtrT;

typedef blas::zero_vector<fType>                zerovalvectType;

typedef blas::vector_range<valvectType>         valvectRangeType;

typedef idType                                  identType;
typedef blas::vector<identType>                 idvectType;
typedef blas::vector<identType>&                idvectRefType;
typedef blas::vector<identType> *               idvectPtrType;

typedef blas::vector_range<idvectType>          idvectRangeType;

typedef std::set<matrixPtrType>                 matrixPtrSetType;
typedef std::set<valvectPtrType>                valvectPtrSetType;
typedef std::set<idvectPtrType>                 idvectPtrSetType;

typedef std::map<std::string, matrixPtrType>    matrixPtrMapType;
typedef std::map<std::string, valvectPtrType>   valvectPtrMapType;
typedef std::map<std::string, idvectPtrType>    idvectPtrMapType;

struct quantsType
{
   matrixPtrSetType  vects;
   valvectPtrSetType scalars;
   idvectPtrSetType  ints;
};

typedef quantsType&   quantsRefType;
typedef quantsType*   quantPtrType;
};

#endif
#endif
