#ifndef SPHLATCH_COMMUNICATION_MANAGER_H
#define SPHLATCH_COMMUNICATION_MANAGER_H

#include <vector>
#include <list>
#include <mpi.h>

//#include <boost/assign/std/vector.hpp>

namespace sphlatch
{
class CommunicationManager
{
public:
   typedef CommunicationManager    self_type;
   typedef CommunicationManager&   self_reference;
   typedef CommunicationManager*   self_pointer;

   static self_reference instance(void);

   typedef MPI::Datatype           dataType;
   class DataTypeCreator;

private:
   void calcOffsets(domainPartsIndexRefType _indices,
                    size_t                  _totalOffset,
                    countsVectRefType       _sendOffset,
                    countsVectRefType       _partsTo,
                    countsVectRefType       _recvOffset,
                    countsVectRefType       _partsFrom);

   domainPartsIndexType ghostIndices;
   countsVectType       gSendOffsets, gRecvOffsets, gPartsFrom, gPartsTo;

   size_t identSize, floatSize, newNoParts, noSendGhosts;

///
/// various little helper functions
///
///
public:

   void sendBitset(bitsetRefType _bitset, size_t _recvDomain);
   void recvBitset(bitsetRefType _bitset, size_t _sendDomain);

   template<class T> void sendVector(std::vector<T>& _vector,
                                     size_t          _recvDomain);

   template<class T> void recvVector(std::vector<T>& _vector,
                                     size_t          _sendDomain);

   void sumUpCounts(countsVectRefType _indexVect);

   void max(fRefType _val);
   void min(fRefType _val);
   void sum(fRefType _val);
   void sum(vect3dT& _vect);
   void sum(countsRefType _cnt);

   size_t getMyDomain();
   size_t getNoDomains();

   void barrier();
   double wtime();

   size_t domainToMPIrank(size_t _domain);
   size_t MPIrankToDomain(size_t _rank);

protected:
   CommunicationManager(void);
   ~CommunicationManager(void);

private:
   static self_pointer _instance;
   static size_t       commBuffSize;
   size_t myDomain, myRank, noDomains, noRanks;

   std::vector<size_t> domainToRank;
   std::vector<size_t> rankToDomain;

   MPI_Datatype mpiFloatType, mpiCountType;
};

CommunicationManager::self_pointer CommunicationManager::_instance = NULL;

CommunicationManager::self_reference CommunicationManager::instance(void)
{
   if (_instance != NULL)
      return(*_instance);
   else
   {
      _instance = new CommunicationManager;
      return(*_instance);
   }
}

CommunicationManager::CommunicationManager(void) :
   noDomains(MPI::COMM_WORLD.Get_size()),
   noRanks(noDomains),
   domainToRank(noDomains),
   rankToDomain(noDomains),
#ifdef SPHLATCH_SINGLEPREC
   mpiFloatType(MPI::FLOAT),
#else
   mpiFloatType(MPI::DOUBLE),
#endif
   mpiCountType(MPI::INT)
{
   ///
   /// insert mapping here
   ///
   /// right now rank maps directly to the domain number
   ///
   for (size_t i = 0; i < noDomains; i++)
   {
      domainToRank[i] = i;
      rankToDomain[i] = i;
   }

   myRank   = MPI::COMM_WORLD.Get_rank();
   myDomain = rankToDomain[myRank];

   gSendOffsets.resize(noDomains);
   gRecvOffsets.resize(noDomains);
   gPartsFrom.resize(noDomains);
   gPartsTo.resize(noDomains);

   //identSize = sizeof(idType);
   //floatSize = sizeof(fType);
}

CommunicationManager::~CommunicationManager(void)
{ }

///
/// buffer size in bytes, MPI transfers will have on average this size
///
size_t CommunicationManager::commBuffSize = 131072;   //  128 kByte
//size_t CommunicationManager::commBuffSize = 1048576;    // 1024 kByte
//size_t CommunicationManager::commBuffSize = 8388608;  // 8096 kByte

///
/// calculate from a particle index vector and a total offset constant the resulting
/// counts and offsets for sending and receiving
///
/// MPI ranks are used insted of domain numbers for the indices of the offset
/// and count vectors!
///
void CommunicationManager::calcOffsets(domainPartsIndexRefType _indices,
                                       size_t                  _totalOffset,
                                       countsVectRefType       _sendOffset,
                                       countsVectRefType       _partsTo,
                                       countsVectRefType       _recvOffset,
                                       countsVectRefType       _partsFrom)
{
   ///
   /// determine number of particles from local rank to remote ranks
   ///
   _partsTo.resize(noRanks);
   _partsFrom.resize(noRanks);
   for (size_t i = 0; i < noDomains; i++)
   {
      _partsTo[domainToRank[i]] = _indices[i].size();
   }

   ///
   /// determine number of particles from remote ranks to local ranks
   ///
   barrier();
   MPI::COMM_WORLD.Alltoall(&_partsTo[0], 1, MPI::INT,
                            &_partsFrom[0], 1, MPI::INT);
   barrier();

   ///
   /// determine send and receive offsets. total offset is not added
   /// for send offset, as there is no total offset in the send buffer.
   ///
   _sendOffset[0] = 0;
   _recvOffset[0] = _totalOffset;
   for (size_t i = 1; i < noRanks; i++)
   {
      _sendOffset[i] = _sendOffset[i - 1] + _partsTo[i - 1];
      _recvOffset[i] = _recvOffset[i - 1] + _partsFrom[i - 1];
   }
}

void CommunicationManager::sendBitset(bitsetRefType _bitset, size_t _recvRank)
{
   /// noBlock = no. of bits / no. of bits per bitsetBlockType
   size_t noBlocks = lrint(
      ceil(static_cast<double>(_bitset.size())
           / static_cast<double>(sizeof(bitsetBlockType) * 8.)));

   ///
   /// allocate buffer and fill it
   ///
   std::vector<bitsetBlockType> sendBuff(noBlocks);
   boost::to_block_range(_bitset, sendBuff.begin());

   ///
   /// send it
   ///
   sendVector<bitsetBlockType>(sendBuff, _recvRank);
}

void CommunicationManager::recvBitset(bitsetRefType _bitset, size_t _sendDomain)
{
   /// noBlock = no. of bits / no. of bits per bitsetBlockType
   size_t noBlocks = lrint(
      ceil(static_cast<double>(_bitset.size())
           / static_cast<double>(sizeof(bitsetBlockType) * 8.)));

   ///
   /// allocate buffer
   ///
   std::vector<bitsetBlockType> recvBuff(noBlocks);

   ///
   /// send it
   ///
   recvVector<bitsetBlockType>(recvBuff, _sendDomain);

   ///
   /// copy buffer to bitset
   ///
   boost::from_block_range(recvBuff.begin(), recvBuff.end(), _bitset);
}

void CommunicationManager::sumUpCounts(countsVectRefType _indexVect)
{
   static size_t noRemElems, noCurElems, round;
   const size_t  noBuffElems = commBuffSize / sizeof(countsType);

   static countsVectType rootBuff;

   noRemElems = _indexVect.size();
   round      = 0;

   while (noRemElems > 0)
   {
      noCurElems = std::min(noRemElems, noBuffElems);
      ///
      /// check whether this MPI implementation supports
      /// the faster MPI_IN_PLACE method
      ///
#ifdef MPI_IN_PLACE
      ///
      /// the vector cast is ugly, but actually legal:
      /// http://www.parashift.com/c++-faq-lite/containers.html#faq-34.3
      ///
      MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE,
                                &_indexVect[round * noBuffElems],
                                noCurElems, MPI::INT, MPI::SUM);
#else
      ///
      /// slower without MPI_IN_PLACE
      ///
      static countsVectType recvBuffer;
      recvBuffer.resize(noCurElems);
      MPI::COMM_WORLD.Allreduce(&_indexVect[round * noBuffElems],
                                &recvBuffer[0],
                                noCurElems, MPI::INT, MPI::SUM);
      const size_t offset = round * noBuffElems;
      for (size_t i = 0; i < noCurElems; i++)
      {
         _indexVect[i + offset] = recvBuffer[i];
      }
#endif
      noRemElems -= noCurElems;
      round++;
   }
}

void CommunicationManager::max(fRefType _val)
{
   static double doubleBuff;

   doubleBuff = static_cast<double>(_val);

#ifdef MPI_IN_PLACE
   MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &doubleBuff, 1,
                             MPI::DOUBLE, MPI::MAX);
#else
   static double recvDoubleBuff;
   MPI::COMM_WORLD.Allreduce(&doubleBuff, &recvDoubleBuff, 1,
                             MPI::DOUBLE, MPI::MAX);
   doubleBuff = recvDoubleBuff;
#endif
   _val = static_cast<fType>(doubleBuff);
}

void CommunicationManager::min(fRefType _val)
{
   static double doubleBuff;

   doubleBuff = static_cast<double>(_val);

#ifdef MPI_IN_PLACE
   MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &doubleBuff, 1,
                             MPI::DOUBLE, MPI::MIN);
#else
   static double recvDoubleBuff;
   MPI::COMM_WORLD.Allreduce(&doubleBuff, &recvDoubleBuff, 1,
                             MPI::DOUBLE, MPI::MIN);
   doubleBuff = recvDoubleBuff;
#endif
   _val = static_cast<fType>(doubleBuff);
}

void CommunicationManager::sum(fRefType _val)
{
   static double doubleBuff;

   doubleBuff = static_cast<double>(_val);

#ifdef MPI_IN_PLACE
   MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &doubleBuff, 1,
                             MPI::DOUBLE, MPI::SUM);
#else
   static double recvDoubleBuff;
   MPI::COMM_WORLD.Allreduce(&doubleBuff, &recvDoubleBuff, 1,
                             MPI::DOUBLE, MPI::SUM);
   doubleBuff = recvDoubleBuff;
#endif
   _val = static_cast<fType>(doubleBuff);
}

void CommunicationManager::sum(vect3dT& _vect)
{
   const size_t vectLength = 3;

#ifdef MPI_IN_PLACE
   MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &_vect(0), vectLength,
                             mpiFloatType, MPI::SUM);
#else
   vect3dT recvBuff(_vect);
   MPI::COMM_WORLD.Allreduce(&_vect(0), &recvBuff(0), vectLength,
                             mpiFloatType, MPI::SUM);
   _vect = recvBuff;
#endif
}

void CommunicationManager::sum(countsRefType _cnt)
{
#ifdef MPI_IN_PLACE
   MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &_cnt, 1,
                             mpiCountType, MPI::SUM);
#else
   idType countBuff;
   MPI::COMM_WORLD.Allreduce(&_cnt, &countBuff, 1,
                             mpiCountType, MPI::SUM);
   _cnt = countBuff;
#endif
}

template<class T>
void CommunicationManager::sendVector(std::vector<T>& _vector,
                                      size_t          _recvDomain)
{
   static size_t noRemElems, noCurElems, noCurBytes, round;
   const size_t  noBuffElems = commBuffSize / sizeof(T);

   const size_t recvRank = domainToRank[_recvDomain];

   ///
   /// no checks are performed whether the vector
   /// on the receiving side actually matches the amount
   /// of incoming data!
   ///
   noRemElems = _vector.size();
   round      = 0;

   while (noRemElems > 0)
   {
      noCurElems = std::min(noRemElems, noBuffElems);
      noCurBytes = noCurElems * sizeof(T);

      ///
      /// this pointer arithmetics stuff assumes continous
      /// storage for std::vector<T>
      ///
      /// does not work for std::vector<bool> !
      ///
      MPI::COMM_WORLD.Send((&_vector[0] + round * noBuffElems),
                           noCurBytes, MPI_BYTE, recvRank,
                           myRank + 255 + round);
      noRemElems -= noCurElems;
      round++;
   }
}

template<class T>
void CommunicationManager::recvVector(std::vector<T>& _vector,
                                      size_t          _sendDomain)
{
   static size_t noRemElems, noCurElems, noCurBytes, round;
   const size_t  noBuffElems = commBuffSize / sizeof(T);

   const size_t sendRank = domainToRank[_sendDomain];

   ///
   /// no checks are performed whether the vector
   /// on the receiving side actually matches the amount
   /// of incoming data!
   ///
   noRemElems = _vector.size();
   round      = 0;

   while (noRemElems > 0)
   {
      noCurElems = std::min(noRemElems, noBuffElems);
      noCurBytes = noCurElems * sizeof(T);

      ///
      /// this pointer arithmetics stuff assumes continous
      /// storage for std::vector<T>
      ///
      /// does not work for std::vector<bool> !
      ///
      MPI::COMM_WORLD.Recv((&_vector[0] + round * noBuffElems),
                           noCurBytes, MPI_BYTE, sendRank,
                           sendRank + 255 + round);
      noRemElems -= noCurElems;
      round++;
   }
}

size_t CommunicationManager::getMyDomain()
{
   return(myDomain);
}

size_t CommunicationManager::getNoDomains()
{
   return(noDomains);
}

void CommunicationManager::barrier()
{
   MPI::COMM_WORLD.Barrier();
}

double CommunicationManager::wtime()
{
   return(MPI::Wtime());
}

size_t CommunicationManager::domainToMPIrank(size_t _domain)
{
   return(domainToRank[_domain]);
}

size_t CommunicationManager::MPIrankToDomain(size_t _rank)
{
   return(rankToDomain[_rank]);
}

class CommunicationManager::DataTypeCreator {
public:

   DataTypeCreator(void* _basePtr) :
      baseAddr(MPI::Get_address(_basePtr)),
#ifdef SPHLATCH_SINGLEPREC
      MPIfType(MPI::FLOAT)
#else
      MPIfType(MPI::DOUBLE)
#endif
   { };
   DataTypeCreator() {};
   ~DataTypeCreator() {};
   void operator+=(vect3dT* _vectPtr)
   {
      elems.push_back(quart(MPI::Get_address(_vectPtr) - baseAddr,
                            MPIfType, 3, sizeof(vect3dT)));
   }

   void operator+=(fType* _fPtr)
   {
      elems.push_back(quart(MPI::Get_address(_fPtr) - baseAddr,
                            MPIfType, 1, sizeof(fType)));
   }

   void operator+=(countsType* _cPtr)
   {
      elems.push_back(quart(MPI::Get_address(_cPtr) - baseAddr,
                            MPI::INT, 1, sizeof(countsType)));
   }

   void operator()(void* _basePtr)
   {
     elems.clear();
     baseAddr = MPI::Get_address(_basePtr);
   }

   MPI::Datatype finalize()
   {
      // sort elements after displacement
      elems.sort();
      const size_t noElems = elems.size() + 2;

      int           blength[noElems];
      MPI::Aint     displac[noElems];
      MPI::Datatype dtypes[noElems];

      elemsListT::iterator       elemItr = elems.begin();
      elemsListT::const_iterator elemEnd = elems.end();

      displac[0] = 0;
      dtypes[0]  = MPI::LB;
      blength[0] = 1;

      size_t count = 1;
      while (elemItr != elemEnd)
      {
         blength[count] = elemItr->blength;
         displac[count] = elemItr->displ;
         dtypes[count]  = elemItr->dtype;
         count++;
         elemItr++;
      }

      if ( count > 1 )
        displac[count] = displac[count - 1] + (elemItr--)->size;
      else
        displac[count] = 0;
      dtypes[count]  = MPI::UB;
      blength[count] = 1;

      MPI::Datatype retDatatype = 
        MPI::Datatype::Create_struct(noElems, blength, displac, dtypes);
      retDatatype.Commit();
      return(retDatatype);
   }

private:
   class quart {
public:
      quart(
         MPI::Aint     _displ,
         MPI::Datatype _dtype,
         int           _blength,
         size_t        _size
         ) :
         displ(_displ),
         dtype(_dtype),
         blength(_blength),
         size(_size)
      { }
      ~quart() { }

      bool operator<(const quart& _rhs)
      {
         return(displ < _rhs.displ);
      }

      MPI::Aint     displ;
      MPI::Datatype dtype;
      int           blength;
      size_t        size;
   };

private:
   MPI::Aint baseAddr;
   typedef std::list<quart>   elemsListT;
   elemsListT    elems;
   MPI::Datatype MPIfType;
};
};

#endif
