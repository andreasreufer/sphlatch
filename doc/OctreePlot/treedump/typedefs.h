/*
 *  typedefs.h
 *  
 *
 *  Created by Andreas Reufer on 15.11.07.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/dynamic_bitset.hpp>

#ifdef SINGLEPREC
typedef float		valueType;
#else
typedef double		valueType;
#endif
typedef valueType*		valuePtrType;
typedef valueType&		valueRefType;

typedef int	identType;

typedef boost::numeric::ublas::matrix<valueType>	matrixType;
typedef boost::numeric::ublas::matrix<valueType>*	matrixPtrType;
typedef boost::numeric::ublas::matrix_row<matrixType>	particleRowType;
typedef boost::numeric::ublas::matrix_row<matrixType>*	particleRowPtr;

typedef unsigned long int bitsetBlockType;
typedef boost::dynamic_bitset<bitsetBlockType> bitsetType;
typedef boost::dynamic_bitset<bitsetBlockType>& bitsetRefType;

