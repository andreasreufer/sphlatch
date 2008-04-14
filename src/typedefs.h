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
#include <list>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/dynamic_bitset.hpp>

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

///
/// a valueType matrix
/// please note that the communication manager
/// assumes this type to have continuous storage
///
typedef blas::matrix<valueType>        matrixType;
typedef blas::matrix<valueType>&       matrixRefType;
typedef blas::matrix<valueType>*       matrixPtrType;

///
/// a matrix row
///
typedef blas::matrix_row<matrixType>   matrixRowType;
typedef blas::matrix_row<matrixType>&  matrixRowRefType;
typedef blas::matrix_row<matrixType>*  matrixRowPtrType;

///
/// matrix range and slice
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
/// a valueType vector
///
/// use valarray instead
typedef blas::vector<valueType> valvectType;
typedef blas::vector<valueType>& valvectRefType;
typedef blas::vector<valueType>* valvectPtrType;

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

};

#endif
