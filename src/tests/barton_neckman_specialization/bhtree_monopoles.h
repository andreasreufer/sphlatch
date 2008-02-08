#ifndef BHTREE_MONOPOLES_H
#define BHTREE_MONOPOLES_H

/*
 *  bhtree_monopoles.h
 *
 *
 *  Created by Andreas Reufer on 07.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "bhtree_monopole_node.h"

namespace sphlatch {
class Monopoles : public BHtree<Monopoles>
{

public:

void allocRootNode(void)
{
  rootPtr = new monopoleCellNode;
}

void allocNewChildNode(const size_t _n)
{
}

void calcGravCell()
{
}

private:
protected:

};
};

#endif