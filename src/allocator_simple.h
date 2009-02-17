#ifndef ALLOCATOR_SIMPLE_H
#define ALLOCATOR_SIMPLE_H

/*
 *  allocator_simple.h
 *
 *  generic simple allocator
 *
 *  Created by Andreas Reufer on 03.02.09.
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

namespace sphlatch {
template<class T>
class SimpleAllocator {
public:
   typedef T*   Tptr;

private:

public:
   SimpleAllocator(const size_t _size) { }

   ~SimpleAllocator() { }

   Tptr pop()
   {
      return(new T);
   }

   void push(Tptr _ptr)
   {
      delete _ptr;
   }
};
};

#endif
