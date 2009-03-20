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

/*#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>*/
#include <boost/dynamic_bitset.hpp>

#include <cmath>
#include <limits>
#include <valarray>

#define BZ_THREADSAFE
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

namespace sphlatch {
///
/// the BLAS namespace
///
//namespace blas = boost::numeric::ublas;

///
/// fType: should be a float
///
#ifdef SPHLATCH_SINGLEPREC
typedef float    fType;
#else
typedef double   fType;
#endif
typedef fType*   fPtrType;
typedef fType&   fRefType;

typedef blitz::TinyVector<fType,3> vect3dT;

///
/// type for particle or tree node IDs
///
typedef int                               idType;

///
/// a vector of (particle) indices
///
typedef std::vector<size_t>               partsIndexVectType;
typedef std::vector<size_t>&              partsIndexVectRefType;
typedef std::vector<size_t> *             partsIndexVectPtrType;

typedef std::list<size_t>                 partsIndexListT;
typedef std::list<size_t>&                partsIndexListRefT;
typedef std::list<size_t> *               partsIndexListPtrT;

///
/// a vector of counts
/// countsType has to be compatible with MPI_INT
///
typedef int                         countsType;
typedef countsType&                 countsRefType;
typedef countsType*                 countsPtrType;
typedef std::vector<countsType>     countsVectType;
typedef std::vector<countsType>&    countsVectRefType;
typedef std::vector<countsType> *   countsVectPtrType;

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
/// a struct with a set of pointers to the
/// three possible containers for physical
/// quantities
///
/*struct quantsType
{
   matrixPtrSetType  vects;
   valvectPtrSetType scalars;
   idvectPtrSetType  ints;
};

typedef quantsType&   quantsRefType;
typedef quantsType*   quantPtrType;*/
};

#endif
