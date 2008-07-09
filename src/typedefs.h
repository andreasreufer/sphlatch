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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/dynamic_bitset.hpp>

#include <cmath>
#include <limits>
#include <valarray>

namespace sphlatch {
///
/// the BLAS namespace
///
namespace blas = boost::numeric::ublas;

///
/// valueType: should be a float 
///
#ifdef SPHLATCH_SINGLEPREC
typedef float valueType;
#else
typedef double valueType;
#endif
typedef valueType*              valuePtrType;
typedef valueType&              valueRefType;
typedef const valueType&        valueCRefType;

///
/// a valueType matrix
/// please note that the communication manager
/// assumes this type to have continuous storage
///
typedef blas::matrix<valueType>        matrixType;
typedef matrixType&                    matrixRefType;
typedef matrixType*                    matrixPtrType;

typedef blas::zero_matrix<valueType>   zeromatrixType;

///
/// a matrix row
///
typedef blas::matrix_row<matrixType>   matrixRowType;
typedef blas::matrix_row<matrixType>&  matrixRowRefType;
typedef blas::matrix_row<matrixType>*  matrixRowPtrType;

typedef blas::matrix_column<matrixType>   matrixColumnType;

typedef blas::matrix_range<matrixType> matrixRangeType;

///
/// range and slice
///
typedef blas::range rangeType;
typedef blas::slice sliceType;

///
/// matrix vector slices
///
typedef blas::matrix_vector_slice<matrixType> matrixVectorSliceType;
typedef blas::matrix_vector_slice<const matrixType> constMatrixVectorSliceType;

///
/// a matrix row represents a particle in a matrix
///
typedef matrixRowType particleRowType;
typedef matrixRowType& particleRowRefType;
typedef matrixRowType* particleRowPtrType;

///
/// a matrix column represents a quantitiy in a matrix
///
typedef matrixColumnType quantColumnType;
typedef matrixColumnType& quantColumnRefType;
typedef matrixColumnType* quantColumnPtrType;

///
/// a valueType vector
///
typedef blas::vector<valueType> valvectType;
typedef blas::vector<valueType>& valvectRefType;
typedef blas::vector<valueType>* valvectPtrType;

typedef blas::zero_vector<valueType> zerovalvectType;


typedef blas::vector_range<valvectType> valvectRangeType;

///
/// type for particle or tree node IDs
///
typedef int identType;

///
/// a vector of IDs
///
typedef blas::vector<identType> idvectType;
typedef blas::vector<identType>& idvectRefType;
typedef blas::vector<identType>* idvectPtrType;

typedef blas::vector_range<idvectType> idvectRangeType;

///
/// a vector of (particle) indices
///
typedef std::vector<size_t> partsIndexVectType;
typedef std::vector<size_t>& partsIndexVectRefType;
typedef std::vector<size_t>* partsIndexVectPtrType;

///
/// a vector of counts
/// countsType has to be compatible with MPI_INT
///
typedef int countsType;
typedef countsType& countsRefType;
typedef countsType* countsPtrType;
typedef std::vector<countsType> countsVectType;
typedef std::vector<countsType>& countsVectRefType;
typedef std::vector<countsType>* countsVectPtrType;

///
/// a vector of particle indices vectors
/// used to map lists of particles to a list of domains
///
typedef std::vector<partsIndexVectType> domainPartsIndexType;
typedef std::vector<partsIndexVectType>& domainPartsIndexRefType;
typedef std::vector<partsIndexVectType>* domainPartsIndexPtrType;

///
/// storage type for bitset
/// the longer, the better
///
typedef unsigned long int bitsetBlockType;

///
/// a bit set, useful as a boolean vector
///
typedef boost::dynamic_bitset<bitsetBlockType> bitsetType;
typedef boost::dynamic_bitset<bitsetBlockType>& bitsetRefType;

///
/// attribute map
///
typedef std::map< std::string, valueType > attrMapType;
typedef attrMapType& attrMapRefType;
typedef attrMapType* attrMapPtrType;

typedef std::set<matrixPtrType>  matrixPtrSetType;
typedef std::set<valvectPtrType> valvectPtrSetType;
typedef std::set<idvectPtrType>  idvectPtrSetType;

///
/// a list of strings
///
typedef std::list< std::string > stringListType;

///
/// a vector of strings
///
typedef std::vector< std::string > stringVectType;

///
/// a struct with a set of pointers to the
/// three possible containers for physical
/// quantities
///
struct quantsType {
  matrixPtrSetType  vects;
  valvectPtrSetType scalars;
  idvectPtrSetType  ints;
};

typedef quantsType& quantsRefType;
typedef quantsType* quantPtrType;

};

#endif
