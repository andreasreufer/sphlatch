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

#include <stack>

namespace sphlatch {
template<class T, size_t chunkSize>
class ChunkAllocator {
public:
   typedef T*   Tptr;

private:
   class Chunk {
public:
      Chunk()
      {
         elemsUsed = 0;
      }

      ~Chunk() { }

      Tptr pop()
      {
         return(&(data[elemsUsed++]));
      }

      void push()
      {
         if (elemsUsed > 0)
         {
            elemsUsed--;
         }
      }

      bool isFull()
      {
         return(elemsUsed >= chunkSize);
      }

      bool isEmpty()
      {
         return(elemsUsed == 0);
      }

private:
      size_t elemsUsed;
      T      data[chunkSize];
   };

   typedef Chunk    ChunkT;
   typedef Chunk*   ChunkPtrT;

   std::stack<ChunkPtrT> chunks;

public:
   ChunkAllocator()
   {
      chunks.push(new Chunk);
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
      if (chunks.top()->isFull())
      {
         chunks.push(new Chunk);
      }

      return(chunks.top()->pop());
   }

   void push(Tpr _ptr)
   {
      chunks.top()->push();

      if (chunks.top()->isEmpty() && (chunks.size() > 1))
      {
         delete chunks.top();
         chunks.pop();
      }
   }
};
};

#endif
