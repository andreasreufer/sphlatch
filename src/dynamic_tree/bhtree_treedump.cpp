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
   void dotDump(std::string _dotFilename, const nodePtrT _node);

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

void BHTreeDump::dotDump(std::string _dotFilename, const nodePtrT _node)
{
   dumpFile.open(_dotFilename.c_str(), std::ios::out);

   curPtr = _node;
   dumpFile << "digraph graphname { \n";
   dotRecursor();
   dumpFile << "}\n";
   dumpFile.close();
}


void BHTreeDump::dotRecursor()
{
   dumpFile << "N" << curPtr << " [";

   //if (curPtr->isParticle)
   //   dumpFile << "label=" << curPtr->ident << ",";
   //else

   if (curPtr->isCZ)
      dumpFile << "label=\"" << curPtr << " (" <<
      static_cast<czllPtrT>(curPtr)->relCost << ")\",";
   else if (curPtr->isParticle)
      //dumpFile << "label=\"" << curPtr->ident << "\",";
      //dumpFile << "label=\"\",";
      dumpFile << "label=\"" << curPtr << "\",";
   else
      //dumpFile << "label=\"\",";
      dumpFile << "label=\"" << curPtr << "\",";

   /*if (curPtr->isParticle)
      dumpFile << "label=\"\",";
      else
      dumpFile << "label=\" " << curPtr << "\",";*/

   //dumpFile << "label=\"" << curPtr << "\",";

   if (curPtr->isParticle)
      //dumpFile << "shape=circle,color=green";
      dumpFile << "shape=box,color=green";
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

   if (not curPtr->isParticle)
   {
      for (size_t i = 0; i < 8; i++)
      {
         if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
         {
            dumpFile << "N" << curPtr << " -> N"
                     << static_cast<gcllPtrT>(curPtr)->child[i] << ";\n";
            goChild(i);
            dotRecursor();
            goUp();
         }
      }
   }

   //if (curPtr->next != NULL)
   //  dumpFile << "N" << curPtr << " -> N" << curPtr->next << " [style=dotted,weight=0];\n";
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
      dumpFile << "  pos: " << static_cast<gcllPtrT>(curPtr)->cen << "\n";
      dumpFile << "  size " << static_cast<gcllPtrT>(curPtr)->clSz << "\n";
      dumpFile << "  s -> " << static_cast<gcllPtrT>(curPtr)->skip << "\n";

      for (size_t i = 0; i < 8; i++)
      {
         dumpFile << " c" << i << " -> "
                  << static_cast<gcllPtrT>(curPtr)->child[i] << "\n";
      }

      if (curPtr->isCZ)
      {
         dumpFile << " noParts:" << static_cast<czllPtrT>(curPtr)->noParts
                  << "\n";
         dumpFile << " relCost:" << static_cast<czllPtrT>(curPtr)->relCost
                  << "\n";
         dumpFile << " ol -> " << static_cast<czllPtrT>(curPtr)->orphLast
                  << "\n";
         dumpFile << " cf -> " << static_cast<czllPtrT>(curPtr)->chldFrst
                  << "\n";
         dumpFile << " cl -> " << static_cast<czllPtrT>(curPtr)->chldLast
                  << "\n";
         dumpFile << " of -> ";

         nodePtrT curOrph = static_cast<czllPtrT>(curPtr)->orphFrst;
         if (curOrph == NULL)
            dumpFile << "0";
         while (curOrph != NULL)
         {
            dumpFile << curOrph << " -> ";
            curOrph = curOrph->next;
         }
         dumpFile << "\n";

         const nodePtrT oldPtr = curPtr;
         curOrph = static_cast<czllPtrT>(curPtr)->orphFrst;
         while (curOrph != NULL)
         {
            curPtr = curOrph;
            ptrRecursor();
            curOrph = curOrph->next;
         }
         curPtr = oldPtr;
      }

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
   else
   {
      dumpFile << "  pos: " << static_cast<pnodPtrT>(curPtr)->pos << "\n";
   }
}
};

#endif
