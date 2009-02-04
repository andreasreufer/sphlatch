#ifndef ALLOCATOR_FIXEDBUCKET_H
#define ALLOCATOR_FIXEDBUCKET_H

/*
 *  allocator_fixedbucket.h
 *
 *  generic stacked allocator. pushed elements are put back on
 *  the stack, when stack does not contain more than <maxScratch>
 *  elements
 *
 *  Created by Andreas Reufer on 17.12.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

namespace sphlatch {
template<class T, size_t maxScratch>
class FixedBucketAllocator {
public:
   typedef T*   Tptr;

private:
   std::vector<Tptr> ptrsCont;
   Tptr ptrPivot;

public:
   FixedBucketAllocator()
   {
      ptrsCont.reserve(maxScratch);
   }

   ~FixedBucketAllocator()
   {
      const size_t noElems = ptrsCont.size();
      for (size_t i = 0; i < noElems; i++)
         delete ptrsCont[i];
   }

   Tptr pop()
   {
      if (not ptrsCont.empty())
      {
         ptrPivot = ptrsCont.back();
         ptrsCont.pop_back();
         return(ptrPivot);
      }
      else
         return(new T);
   }

   void push(Tptr _ptr)
   {
      if (ptrsCont.size() < maxScratch)
         ptrsCont.push_back(_ptr);
      else
         delete _ptr;
   }
};
};

#endif
