#ifndef BHTREE_H
#define BHTREE_H

/*
 *  bhtree.h
 *
 *  base header file for the SPHLATCH Barnes&Hut tree
 *
 *  Created by Andreas Reufer on 12.11.07.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#ifdef SPHLATCH_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <vector>
#include <stack>
#include <queue>

#include "typedefs.h"

// generic BHtree methods
#include "bhtree_generic.h"

// specialized BHtree methods for 
// desired degree of series expansion
#include "bhtree_monopoles.h"
#include "bhtree_quadrupoles.h"
#include "bhtree_octupoles.h"

#endif
