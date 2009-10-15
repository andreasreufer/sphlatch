#ifndef BHTREE_CONFIG_H
#define BHTREE_CONFIG_H

/*
 *  bhtree_config.h
 *
 *  config header file for the dynamic SPHLATCH Barnes&Hut tree
 *
 *  Created by Andreas Reufer on 15.10.09.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

#include "bhtree_nodes.h"
#include "bhtree_config.h"

namespace sphlatch {
  struct BHTreeConfig 
  {
    const size_t maxDepth;
  };
};
#endif
