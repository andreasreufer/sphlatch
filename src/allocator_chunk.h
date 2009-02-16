#ifndef ALLOCATOR_CHUNK_H
#define ALLOCATOR_CHUNK_H

/*
 *  allocator_chunk.h
 *
 *  generic chunked allocator
 *
 *  Created by Andreas Reufer on 03.02.09.
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include <vector>
#include <stack>

namespace sphlatch {
template<class T>
class ChunkAllocator {
public:
   typedef T*                                      Tptr;

private:
   class Chunk {
public:
      typedef std::vector<T>                       dataVectT;
      typedef typename dataVectT::iterator         dataVectItrT;
      typedef typename dataVectT::const_iterator   dataVectCItrT;

      Chunk(const size_t _chunkSize) :
         data(_chunkSize),
         curElem(data.begin()),
         fstElem(data.begin()),
         endElem(data.end())
      { }
      ~Chunk() { }

      dataVectT           data;
      dataVectItrT        curElem;
      const dataVectCItrT fstElem;
      const dataVectCItrT endElem;
   };

   typedef Chunk    ChunkT;
   typedef Chunk*   ChunkPtrT;

   std::stack<ChunkPtrT> chunks;
   const size_t          chunkSize;

public:
   ChunkAllocator(size_t _chunkSize) :
      chunkSize(_chunkSize)
   {
      chunks.push(new Chunk(_chunkSize));
   }

   ~ChunkAllocator()
   {
      while (chunks.size() > 0)
      {
         delete chunks.top();
         chunks.pop();
      }
   }

   Tptr pop()
   {
      if (chunks.top()->curElem == chunks.top()->endElem)
         chunks.push(new Chunk(chunkSize));

      return(& * (chunks.top()->curElem++));
   }

   void push(Tptr _ptr)
   {
      if (chunks.top()->curElem == chunks.top()->fstElem)
      {
         delete chunks.top();
         chunks.pop();
      }

      chunks.top()->curElem--;
   }
};
};

#endif
