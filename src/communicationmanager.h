#ifndef SPHLATCH_COMMUNICATION_MANAGER_H
#define SPHLATCH_COMMUNICATION_MANAGER_H

#include <vector>
#include <queue>
#include <mpi.h>

#include <boost/assign/std/vector.hpp>

namespace sphlatch
{
class CommunicationManager
{
public:
typedef CommunicationManager self_type;
typedef CommunicationManager& self_reference;
typedef CommunicationManager* self_pointer;

//typedef MPI_SingleCom < typename SimulationTrait::value_type > single_com_type;

static self_reference instance(void);

void exchange(matrixRefType _dataSource,
              domainPartsIndexRefType domainIndices,
              matrixRefType _dataTarget);

void sendBitSet(bitsetRefType _bitset, size_t _recvRank);
void recvBitSet(bitsetRefType _bitset, size_t _sendRank);

void sendMatrix(matrixRefType _bitset, size_t _recvRank);
void recvMatrix(matrixRefType _bitset, size_t _sendRank);

void sumUpCounts(countsVectRefType _indexVect);

protected:

CommunicationManager(void);
~CommunicationManager(void);

private:

static self_pointer _instance;
//static single_com_type SingleCom;

static size_t commBuffSize;

std::vector<matrixType> composeBuffers;
std::vector<matrixType> recvBuffers;
matrixType sendBuffer;

std::vector<size_t> sendNoPart;
std::vector<size_t> recvNoPart;

std::queue< std::pair<size_t, std::vector<size_t> > > sendQueue;
};

CommunicationManager::self_pointer CommunicationManager::_instance = NULL;

CommunicationManager::self_reference CommunicationManager::instance(void)
{
  if (_instance != NULL)
    return *_instance;
  else
    {
      _instance = new CommunicationManager;
      return *_instance;
    }
}

//typename CommunicationManager<SimTraitType>::single_com_type CommunicationManager<SimTraitType>::SingleCom;

CommunicationManager::CommunicationManager(void)
{
}

CommunicationManager::~CommunicationManager(void)
{
}

///
/// buffer size in bytes, MPI transfers will have on average this size
///
size_t CommunicationManager::commBuffSize = 131072; //  128 kByte
//size_t CommunicationManager::commBuffSize = 1048576; // 1024 kByte

/*
   template <int N>
   typename CommunicationManager<SimTraitType>::value_type CommunicationManager<SimTraitType>::Min(void)
   {
   static const matrix_column Column(Data, N);
   value_type local = *std::min_element(Column.begin(), Column.end());
   value_type global;

   SingleCom.min(local, global);
   return global;
   };

   template <int N>
   typename CommunicationManager<SimTraitType>::value_type
   CommunicationManager<SimTraitType>::Max(void)
   {
   static const matrix_column Column(Data, N);
   value_type local = *std::max_element(Column.begin(), Column.end());
   value_type global;

   SingleCom.max(local, global);
   return global;
   }
 */

void CommunicationManager::exchange(matrixRefType _dataSource,
                                    domainPartsIndexRefType domainIndices,
                                    matrixRefType _dataTarget)
{
  const size_t RANK = MPI::COMM_WORLD.Get_rank();
  const size_t SIZE = MPI::COMM_WORLD.Get_size();

  /// Only common range will be exchanged,
  /// _dataTarget will be resized to that size
  const size_t partSize = std::min(_dataSource.size2(), _dataTarget.size2());

  /// number of particles in a comm buffer
  const size_t noBuffParts = lrint(floor(commBuffSize / (partSize * sizeof(valueType))));
  const size_t commBuffByteEffSize = noBuffParts * partSize * sizeof(valueType);

  /// Slice of matrix which will be copied:
  /// first index is 0, iterate through
  /// every index until partSize */
  sliceType partSlice(0, 1, partSize);

  /// Exchange knowledge about who is sending whom how much
  sendNoPart.resize(SIZE);
  recvNoPart.resize(SIZE);

  assert(domainIndices.size() == SIZE);
  for (size_t i = 0; i < SIZE; i++)
    {
      sendNoPart[i] = domainIndices[i].size();
    }

  MPI::COMM_WORLD.Alltoall(&sendNoPart[0], 1, MPI::UNSIGNED_LONG,
                           &recvNoPart[0], 1, MPI::UNSIGNED_LONG);

  /// Fill up send queue
  for (size_t t = 1; t < SIZE; t++)
    {
      const size_t curSendRank = (RANK + t) % SIZE;
      const size_t noSendParts = domainIndices[curSendRank].size();

      using namespace boost::assign;
      std::vector<size_t> curIndices;

      for (size_t i = 0; i < noSendParts; i++)
        {
          // terribly slow to assign
          // make static
          curIndices += domainIndices[curSendRank][i];
          if (curIndices.size() == noBuffParts || i == noSendParts - 1)
            {
              std::pair<size_t, std::vector<size_t> > rankIndicesPair;
              rankIndicesPair.first = curSendRank;
              rankIndicesPair.second = curIndices;
              sendQueue.push(rankIndicesPair);
              curIndices.resize(0);
            }
        }
    }

  /// Prepare composeBuffers
  composeBuffers.resize(SIZE);
  for (size_t i = 0; i < SIZE; i++)
    {
      composeBuffers[i].resize(recvNoPart[i], partSize);
    }

  /// Copy particles staying on node to ComposeBuffer
  for (size_t i = 0; i < domainIndices[RANK].size(); i++)
    {
      matrixRowType(composeBuffers[RANK], i) =
        constMatrixVectorSliceType(_dataSource,
                                   sliceType(domainIndices[RANK][i], 0, partSize),
                                   partSlice);
    }
  recvNoPart[RANK] = 0;

  /// Prepare Communication
  std::vector<MPI::Request> recvReq;
  recvReq.resize(SIZE);
  recvBuffers.resize(SIZE);

  for (size_t i = 0; i < SIZE; i++)
    {
      if (recvNoPart[i] > 0)
        {
          (recvBuffers[i]).resize(noBuffParts, partSize);
          recvReq[i] = MPI::COMM_WORLD.Irecv(&recvBuffers[i](0, 0),
                                             commBuffByteEffSize, MPI::BYTE, i, i);
        }
    }

  sendBuffer.resize(noBuffParts, partSize);

  /// Start Communication!
  bool commFinished = false;

  std::vector<size_t> noPartsRecvd(SIZE);
  for (size_t i = 0; i < SIZE; i++)
    {
      noPartsRecvd[i] = 0;
    }

  while (!commFinished)
    {
      /// Assume we have finished, prove otherwise
      commFinished = true;

      /// Check every receiver for arrived particles
      for (size_t i = 0; i < SIZE; i++)
        {
          /// Still expecting particles and buffer is filled
          if (recvNoPart[i] > 0 && recvReq[i].Test())
            {
              /// Flush the receiving buffer...
              size_t noFlushs = std::min(noBuffParts, recvNoPart[i]);
              for (size_t j = 0; j < noFlushs; j++)
                {
                  matrixRowType(composeBuffers[i], noPartsRecvd[i]) =
                    matrixRowType(recvBuffers[i], j);
                  noPartsRecvd[i]++;
                }
              recvNoPart[i] -= noFlushs;

              /// ... and setup new Receiver if there's more to come
              if (recvNoPart[i] > 0)
                {
                  recvReq[i] = MPI::COMM_WORLD.Irecv(&recvBuffers[i](0, 0),
                                                     commBuffByteEffSize, MPI::BYTE, i, i);
                }
            }
        }

      /// If sendQueue not empty send something
      if (sendQueue.size())
        {
          /// Pick particles to be sent and put them into the sending buffer
          // make this stuff static!
          std::pair<size_t, std::vector<size_t> > sendPair = sendQueue.front();
          sendQueue.pop();
          const size_t sendRank = sendPair.first;
          std::vector<size_t> sendIndices = sendPair.second;
          const size_t noSends = sendIndices.size();

          for (size_t i = 0; i < noSends; i++)
            {
              matrixRowType(sendBuffer, i) =
                constMatrixVectorSliceType(_dataSource,
                                           sliceType(sendIndices[i], 0, partSize), partSlice);
            }

          /// Fire the buffer! (Synchronous send)
          MPI::COMM_WORLD.Ssend(&sendBuffer(0, 0), commBuffByteEffSize, MPI::BYTE, sendRank, RANK);
          commFinished = false;
        }

      /// Check whether communication has finished
      for (size_t i = 0; i < SIZE; i++)
        {
          if (recvNoPart[i] > 0)
            {
              commFinished = false;
            }
        }
    }

  /// Determine number of particles in ComposeBuffer
  size_t noCompParts = 0;
  for (size_t i = 0; i < SIZE; i++)
    {
      noCompParts += composeBuffers[i].size1();
    }

  /// Prepare target matrix
  _dataTarget.resize(noCompParts, partSize);

  /// Flush composeBuffers
  size_t partCounter = 0;
  for (size_t i = 0; i < SIZE; i++)
    {
      const size_t noCurBuffParts = composeBuffers[i].size1();
      for (size_t j = 0; j < noCurBuffParts; j++)
        {
          matrixRowType(_dataTarget, partCounter) = matrixRowType(composeBuffers[i], j);
          partCounter++;
        }
      composeBuffers[i].resize(0, 0);
    }

  /// Be polite and wait for the others (it's useful too by the way)
  MPI::COMM_WORLD.Barrier();

  return;
}

void CommunicationManager::sumUpCounts(countsVectRefType _indexVect)
{
  static size_t noRemElems, noCurElems, round;
  const size_t noBuffElems = commBuffSize / sizeof(countsType);

  static countsVectType rootBuff;

  noRemElems = _indexVect.size();
  round = 0;

  while (noRemElems > 0)
    {
      noCurElems = std::min(noRemElems, noBuffElems);
      ///
      /// ugly, but actually legal:
      /// http://www.parashift.com/c++-faq-lite/containers.html#faq-34.3
      ///
      MPI::COMM_WORLD.Allreduce(&_indexVect[round*noBuffElems],
                                &_indexVect[round*noBuffElems],
                                noCurElems, MPI::INT, MPI::SUM);
      noRemElems -= noCurElems;
      round++;
    }
}

void CommunicationManager::sendBitSet(bitsetRefType _bitset, size_t _recvRank)
{
}
void CommunicationManager::recvBitSet(bitsetRefType _bitset, size_t _sendRank)
{
}

void CommunicationManager::sendMatrix(matrixRefType _bitset, size_t _recvRank)
{
}
void CommunicationManager::recvMatrix(matrixRefType _bitset, size_t _sendRank)
{
}


/*template <typename T>
   struct MPI_SingleCom {
   void min(T& t1, T& t2)
   {
    MPI::COMM_WORLD.Reduce(&t1, &t2, sizeof(T), MPI::BYTE, MPI::MIN, 0);
    MPI::COMM_WORLD.Bcast(&t1, sizeof(T), MPI::BYTE, 0);
   }
   void max(T& t1, T& t2)
   {
    MPI::COMM_WORLD.Reduce(&t1, &t2, sizeof(T), MPI::BYTE, MPI::MAX, 0);
    MPI::COMM_WORLD.Bcast(&t2, sizeof(T), MPI::BYTE, 0);
   }
   }
   ;

   template <>
   struct MPI_SingleCom<float>{
   void min(float& t1, float& t2)
   {
    MPI::COMM_WORLD.Reduce(&t1, &t2, 1, MPI::FLOAT, MPI::MIN, 0);
    MPI::COMM_WORLD.Bcast(&t2, 1, MPI::FLOAT, 0);
   }

   void max(float& t1, float& t2)
   {
    MPI::COMM_WORLD.Reduce(&t1, &t2, 1, MPI::FLOAT, MPI::MAX, 0);
    MPI::COMM_WORLD.Bcast(&t2, 1, MPI::FLOAT, 0);
   }
   };

   template <>
   struct MPI_SingleCom<double>{
   void min(double& t1, double& t2)
   {
    MPI::COMM_WORLD.Reduce(&t1, &t2, 1, MPI::DOUBLE, MPI::MIN, 0);
    MPI::COMM_WORLD.Bcast(&t2, 1, MPI::DOUBLE, 0);
   }

   void max(double& t1, double& t2)
   {
    MPI::COMM_WORLD.Reduce(&t1, &t2, 1, MPI::DOUBLE, MPI::MAX, 0);
    MPI::COMM_WORLD.Bcast(&t2, 1, MPI::DOUBLE, 0);
   }
   };*/
}
;

#endif
