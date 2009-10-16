#ifndef BHTREE_WORKER_TREEDUMP_CPP
#define BHTREE_WORKER_TREEDUMP_CPP

/*
 *  bhtree_worker_treedump.cpp
 *
 *  Created by Andreas Reufer on 14.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include <fstream>
#include "bhtree_worker.cpp"

namespace sphlatch {
class BHTreeDump : public BHTreeWorker {
public:
   BHTreeDump(const treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   BHTreeDump(const BHTreeDump& _dumper) : BHTreeWorker(_dumper) { }
   ~BHTreeDump() { }

public:
   void dotDump(std::string _dotFilename);

   void ptrDump(const std::string _filename);
   void ptrDump(const std::string _filename, const nodePtrT _node);

private:
   std::fstream dumpFile;

   void dotRecursor();
   void ptrRecursor();
};

void BHTreeDump::dotDump(std::string _dotFilename)
{
   dumpFile.open(_dotFilename.c_str(), std::ios::out);

   goRoot();
   dumpFile << "digraph graphname { \n";
   dotRecursor();
   dumpFile << "}\n";
   dumpFile.close();
}

void BHTreeDump::dotRecursor()
{
   if (curPtr->isParticle)
      dumpFile << "P";
   else
      dumpFile << "C";

   dumpFile << abs(curPtr->ident) << " [";
   //if (curPtr->isParticle)
   //   dumpFile << "label=" << curPtr->ident << ",";
   //else
   dumpFile << "label=\"\",";

   if (curPtr->isParticle)
      dumpFile << "shape=circle,color=green";
   else if (curPtr->isCZ)
   {
      dumpFile << "shape=diamond,";
      if (static_cast<czllPtrT>(curPtr)->atBottom)
         dumpFile << "color=darkorange";
      else
         dumpFile << "color=red";
   }
   else
      dumpFile << "shape=box,color=red";

   dumpFile << ",style=filled";

   //if (curPtr->isParticle)
   dumpFile << ",width=0.1,height=0.1];\n";

   /*else
      {
      const fType curSize = 0.1*sqrt( static_cast<fType>( static_cast<gcllPtrT>(curPtr)->noParts ) );
      dumpFile << ",width=" << curSize << ",height=" << curSize << "];\n";
      }*/

   if (curPtr->isParticle == false)
   {
      for (size_t i = 0; i < 8; i++)
      {
         if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
         {
            dumpFile << "C" << abs(curPtr->ident)
                     << " -> ";
            if (static_cast<gcllPtrT>(curPtr)->child[i]->isParticle)
               dumpFile << "P";
            else
               dumpFile << "C";

            dumpFile << abs(static_cast<gcllPtrT>(curPtr)->child[i]->ident)
                     << " \n";
            goChild(i);
            dotRecursor();
            goUp();
         }
      }
   }
}

void BHTreeDump::ptrDump(std::string _filename)
{
   dumpFile.open(_filename.c_str(), std::ios::out);
   goRoot();
   ptrRecursor();
   dumpFile.close();
}

void BHTreeDump::ptrDump(std::string _filename, const nodePtrT _node)
{
   dumpFile.open(_filename.c_str(), std::ios::out);
   curPtr = _node;
   ptrRecursor();
   dumpFile.close();
}

void BHTreeDump::ptrRecursor()
{
   dumpFile << curPtr << " @d" << curPtr->depth << "  "
            << "part?" << curPtr->isParticle << "  "
            << "CZ?" << curPtr->isCZ << "  "
            << "rem?" << curPtr->isRemote << "  "
            << "nSet?" << curPtr->neighSet << "  "
            << "atBo?" << curPtr->atBottom << "  "
            << "setl?" << curPtr->isSettled << "\n";

   dumpFile << "  p -> " << curPtr->parent << "\n";
   dumpFile << "  n -> " << curPtr->next << "\n";

   if (not curPtr->isParticle)
   {
      dumpFile << "  s -> " << static_cast<gcllPtrT>(curPtr)->skip << "\n";

      if (curPtr->isCZ)
      {
         const nodePtrT orphFrst = static_cast<czllPtrT>(curPtr)->orphFrst;

         dumpFile << " of -> " << orphFrst;

         if (orphFrst != NULL)
         {
            curPtr = orphFrst;
            dumpFile << orphFrst << " -> ";
            while (curPtr->next != NULL)
            {
               goNext();
               dumpFile << curPtr << " -> ";
            }
            dumpFile << curPtr->next;
         }
         dumpFile << "\n";
         dumpFile << " ol -> " << static_cast<czllPtrT>(curPtr)->orphLast
                  << "\n";
         dumpFile << " cf -> " << static_cast<czllPtrT>(curPtr)->chldFrst
                  << "\n";
         dumpFile << " cl -> " << static_cast<czllPtrT>(curPtr)->chldLast
                  << "\n";
         dumpFile << " noParts:" << static_cast<czllPtrT>(curPtr)->noParts
                  << "\n";
         dumpFile << " absCost:" << static_cast<czllPtrT>(curPtr)->absCost
                  << "\n";
      }

      for (size_t i = 0; i < 8; i++)
         dumpFile << " c" << i << " -> " <<
         static_cast<gcllPtrT>(curPtr)->child[i] << "\n";

      for (size_t i = 0; i < 8; i++)
      {
         if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
         {
            goChild(i);
            ptrRecursor();
            goUp();
         }
      }
   }
}
};

#endif
