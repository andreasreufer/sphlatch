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

static self_reference instance(void);

void exchange(matrixRefType _dataSource,
              domainPartsIndexRefType domainIndices,
              matrixRefType _dataTarget);

void sendMatrix(matrixRefType _matrix, size_t _recvDomain);
void recvMatrix(matrixRefType _matrix, size_t _sendDomain);

void sendBitset(bitsetRefType _bitset, size_t _recvDomain);
void recvBitset(bitsetRefType _bitset, size_t _sendDomain);

void sumUpCounts(countsVectRefType _indexVect);

void max(valueRefType _val);
void min(valueRefType _val);
void sum(valueRefType _val);

template<class T> void sendVector(std::vector<T>& _vector, size_t _recvDomain);
template<class T> void recvVector(std::vector<T>& _vector, size_t _sendDomain);

size_t getMyDomain();
size_t getNoDomains();

size_t domainToMPIrank(size_t _domain);
size_t MPIrankToDomain(size_t _rank);

protected:
CommunicationManager(void);
~CommunicationManager(void);

private:
static self_pointer _instance;
static size_t commBuffSize;
size_t myDomain, noDomains;

std::vector<matrixType> composeBuffers;
std::vector<matrixType> recvBuffers;
matrixType sendBuffer;

std::vector<size_t> sendNoPart;
std::vector<size_t> recvNoPart;

std::vector<size_t> domainToRank;
std::vector<size_t> rankToDomain;

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

CommunicationManager::CommunicationManager(void)
{
  noDomains = MPI::COMM_WORLD.Get_size();
  domainToRank.resize( noDomains );
  rankToDomain.resize( noDomains );
  
  ///
  /// insert mapping here
  ///
  for (size_t i = 0; i < noDomains; i++ )
  {
    domainToRank[i] = i;
    rankToDomain[i] = i;
  }

  myDomain = rankToDomain[ MPI::COMM_WORLD.Get_rank() ];
}

CommunicationManager::~CommunicationManager(void)
{
}

///
/// buffer size in bytes, MPI transfers will have on average this size
///
//size_t CommunicationManager::commBuffSize = 131072; //  128 kByte
size_t CommunicationManager::commBuffSize = 1048576; // 1024 kByte

void CommunicationManager::exchange(matrixRefType _dataSource,
                                    domainPartsIndexRefType domainIndices,
                                    matrixRefType _dataTarget)
{
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
  sendNoPart.resize(noDomains);
  recvNoPart.resize(noDomains);

  assert(domainIndices.size() == noDomains);
  for (size_t i = 0; i < noDomains; i++)
    {
      sendNoPart[i] = domainIndices[i].size();
    }

  MPI::COMM_WORLD.Alltoall(&sendNoPart[0], 1, MPI::UNSIGNED_LONG,
                           &recvNoPart[0], 1, MPI::UNSIGNED_LONG);

  /// Fill up send queue
  for (size_t t = 1; t < noDomains; t++)
    {
      const size_t curSendDomain = (myDomain + t) % noDomains;
      const size_t noSendParts = domainIndices[curSendDomain].size();

      using namespace boost::assign;
      std::vector<size_t> curIndices;

      for (size_t i = 0; i < noSendParts; i++)
        {
          // terribly slow to assign
          // make static
          curIndices += domainIndices[curSendDomain][i];
          if (curIndices.size() == noBuffParts || i == noSendParts - 1)
            {
              std::pair<size_t, std::vector<size_t> > domainIndicesPair;
              domainIndicesPair.first = curSendDomain;
              domainIndicesPair.second = curIndices;
              sendQueue.push(domainIndicesPair);
              curIndices.resize(0);
            }
        }
    }

  /// Prepare composeBuffers
  composeBuffers.resize(noDomains);
  for (size_t i = 0; i < noDomains; i++)
    {
      composeBuffers[i].resize(recvNoPart[i], partSize);
    }

  /// Copy particles staying on node to ComposeBuffer
  for (size_t i = 0; i < domainIndices[myDomain].size(); i++)
    {
      matrixRowType(composeBuffers[myDomain], i) =
        constMatrixVectorSliceType(_dataSource,
                                   sliceType(domainIndices[myDomain][i], 0, partSize),
                                   partSlice);
    }
  recvNoPart[myDomain] = 0;

  /// Prepare Communication
  std::vector<MPI::Request> recvReq;
  recvReq.resize(noDomains);
  recvBuffers.resize(noDomains);

  for (size_t i = 0; i < noDomains; i++)
    {
      if (recvNoPart[i] > 0)
        {
          const size_t recvRank = domainToRank[i];
          (recvBuffers[i]).resize(noBuffParts, partSize);
          recvReq[i] = MPI::COMM_WORLD.Irecv(&recvBuffers[i](0, 0),
                                             commBuffByteEffSize, MPI::BYTE, recvRank, recvRank);
        }
    }

  sendBuffer.resize(noBuffParts, partSize);

  /// Start Communication!
  bool commFinished = false;

  std::vector<size_t> noPartsRecvd(noDomains);
  for (size_t i = 0; i < noDomains; i++)
    {
      noPartsRecvd[i] = 0;
    }

  while (!commFinished)
    {
      /// Assume we have finished, prove otherwise
      commFinished = true;

      /// Check every receiver for arrived particles
      for (size_t i = 0; i < noDomains; i++)
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
                  const size_t recvRank = domainToRank[i];
                  recvReq[i] = MPI::COMM_WORLD.Irecv(&recvBuffers[i](0, 0),
                                                     commBuffByteEffSize, MPI::BYTE, recvRank, recvRank);
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
          const size_t sendRank = domainToRank[ sendPair.first ];
          const size_t myRank = domainToRank[myDomain];
          std::vector<size_t> sendIndices = sendPair.second;
          const size_t noSends = sendIndices.size();

          for (size_t i = 0; i < noSends; i++)
            {
              matrixRowType(sendBuffer, i) =
                constMatrixVectorSliceType(_dataSource,
                                           sliceType(sendIndices[i], 0, partSize), partSlice);
            }

          /// Fire the buffer! (Synchronous send)
          MPI::COMM_WORLD.Ssend(&sendBuffer(0, 0), commBuffByteEffSize, MPI::BYTE, sendRank, myRank);
          commFinished = false;
        }

      /// Check whether communication has finished
      for (size_t i = 0; i < noDomains; i++)
        {
          if (recvNoPart[i] > 0)
            {
              commFinished = false;
            }
        }
    }

  /// Determine number of particles in ComposeBuffer
  size_t noCompParts = 0;
  for (size_t i = 0; i < noDomains; i++)
    {
      noCompParts += composeBuffers[i].size1();
    }

  /// Prepare target matrix
  _dataTarget.resize(noCompParts, partSize);

  /// Flush composeBuffers
  size_t partCounter = 0;
  for (size_t i = 0; i < noDomains; i++)
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

void CommunicationManager::sendMatrix(matrixRefType _matrix, size_t _recvDomain)
{
  static size_t noRemElems, noCurElems, noCurBytes, round;
  const size_t noBuffElems = commBuffSize / sizeof(valueType);

  const size_t recvRank = domainToRank[_recvDomain];
  const size_t myRank = domainToRank[myDomain];

  ///
  /// no checks are performed whether the matrix
  /// on the receiving side actually matches the amount
  /// of incoming data!
  ///
  noRemElems = _matrix.size1()*_matrix.size2();
  round = 0;

  while (noRemElems > 0)
    {
      noCurElems = std::min(noRemElems, noBuffElems);
      noCurBytes = noCurElems * sizeof(valueType);

      ///
      /// this pointer arithmetics stuff assumes continous
      /// storage for matrixType
      ///
      MPI::COMM_WORLD.Send((&_matrix(0,0) + round*noBuffElems),
                           noCurBytes, MPI_BYTE, recvRank, myRank + round);
      noRemElems -= noCurElems;
      round++;
    }
}

void CommunicationManager::recvMatrix(matrixRefType _matrix, size_t _sendDomain)
{
  static size_t noRemElems, noCurElems, noCurBytes, round;
  const size_t noBuffElems = commBuffSize / sizeof(valueType);

  const size_t sendRank = domainToRank[_sendDomain];
  ///
  /// no checks are performed whether the matrix
  /// on the receiving side actually matches the amount
  /// of incoming data!
  ///
  noRemElems = _matrix.size1()*_matrix.size2();
  round = 0;

  while (noRemElems > 0)
    {
      noCurElems = std::min(noRemElems, noBuffElems);
      noCurBytes = noCurElems * sizeof(valueType);

      ///
      /// this pointer arithmetics stuff assumes continous
      /// storage for matrixType
      ///
      MPI::COMM_WORLD.Recv((&_matrix(0,0) + round*noBuffElems),
                           noCurBytes, MPI_BYTE, sendRank, sendRank + round);
      noRemElems -= noCurElems;
      round++;
    }
}

void CommunicationManager::sendBitset(bitsetRefType _bitset, size_t _recvRank)
{
  /// noBlock = no. of bits / no. of bits per bitsetBlockType
  size_t noBlocks = lrint(
                      ceil(static_cast<double>(_bitset.size())
                      / static_cast<double>(sizeof(bitsetBlockType) * 8.)) );

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
                      / static_cast<double>(sizeof(bitsetBlockType) * 8.)) );

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
  const size_t noBuffElems = commBuffSize / sizeof(countsType);

  static countsVectType rootBuff;

  noRemElems = _indexVect.size();
  round = 0;

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
      MPI::COMM_WORLD.Allreduce( MPI_IN_PLACE ,
                                &_indexVect[round*noBuffElems],
                                noCurElems, MPI::INT, MPI::SUM);
#else
      ///
      /// slower without MPI_IN_PLACE
      ///
      static countsVectType recvBuffer;
      recvBuffer.resize(noCurElems);
      MPI::COMM_WORLD.Allreduce(&_indexVect[round*noBuffElems],
                                &recvBuffer[0],
                                noCurElems, MPI::INT, MPI::SUM);
      const size_t offset = round*noBuffElems;
      for (size_t i = 0; i < noCurElems; i++)
      {
        _indexVect[i + offset] = recvBuffer[i];
      }
#endif
      noRemElems -= noCurElems;
      round++;
    }
}

void CommunicationManager::max(valueRefType _val)
{
  static double doubleBuff;
  doubleBuff = static_cast<double>(_val);

#ifdef MPI_IN_PLACE
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE , &doubleBuff, 1, MPI::DOUBLE, MPI::MAX);
#else
  static double recvDoubleBuff;
  MPI::COMM_WORLD.Allreduce(&doubleBuff, &recvDoubleBuff, 1, MPI::DOUBLE, MPI::MAX);
  doubleBuff = recvDoubleBuff;
#endif
  _val = static_cast<valueType>( doubleBuff );
}

void CommunicationManager::min(valueRefType _val)
{
  static double doubleBuff;
  doubleBuff = static_cast<double>(_val);

#ifdef MPI_IN_PLACE
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE , &doubleBuff, 1, MPI::DOUBLE, MPI::MIN);
#else
  static double recvDoubleBuff;
  MPI::COMM_WORLD.Allreduce(&doubleBuff, &recvDoubleBuff, 1, MPI::DOUBLE, MPI::MIN);
  doubleBuff = recvDoubleBuff;
#endif
  _val = static_cast<valueType>( doubleBuff );
}

void CommunicationManager::sum(valueRefType _val)
{
  static double doubleBuff;
  doubleBuff = static_cast<double>(_val);

#ifdef MPI_IN_PLACE
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE , &doubleBuff, 1, MPI::DOUBLE, MPI::SUM);
#else
  static double recvDoubleBuff;
  MPI::COMM_WORLD.Allreduce(&doubleBuff, &recvDoubleBuff, 1, MPI::DOUBLE, MPI::SUM);
  doubleBuff = recvDoubleBuff;
#endif
  _val = static_cast<valueType>( doubleBuff );
}

template<class T>
void CommunicationManager::sendVector(std::vector<T>& _vector, size_t _recvDomain)
{
  static size_t noRemElems, noCurElems, noCurBytes, round;
  const size_t noBuffElems = commBuffSize / sizeof(T);

  const size_t recvRank = domainToRank[_recvDomain];
  const size_t myRank = domainToRank[myDomain];

  ///
  /// no checks are performed whether the vector
  /// on the receiving side actually matches the amount
  /// of incoming data!
  ///
  noRemElems = _vector.size();
  round = 0;

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
      MPI::COMM_WORLD.Send((&_vector[0] + round*noBuffElems),
                           noCurBytes, MPI_BYTE, recvRank, myRank + 255 + round);
      noRemElems -= noCurElems;
      round++;
    }
}

template<class T>
void CommunicationManager::recvVector(std::vector<T>& _vector, size_t _sendDomain)
{
  static size_t noRemElems, noCurElems, noCurBytes, round;
  const size_t noBuffElems = commBuffSize / sizeof(T);

  const size_t sendRank = domainToRank[_sendDomain];
  ///
  /// no checks are performed whether the vector
  /// on the receiving side actually matches the amount
  /// of incoming data!
  ///
  noRemElems = _vector.size();
  round = 0;

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
      MPI::COMM_WORLD.Recv((&_vector[0] + round*noBuffElems),
                           noCurBytes, MPI_BYTE, sendRank, sendRank + 255 + round);
      noRemElems -= noCurElems;
      round++;
    }
}

size_t CommunicationManager::getMyDomain()
{
  return myDomain;
}
size_t CommunicationManager::getNoDomains()
{
  return noDomains;
}

size_t CommunicationManager::domainToMPIrank(size_t _domain)
{
  return domainToRank[_domain];
}
size_t CommunicationManager::MPIrankToDomain(size_t _rank)
{
  return rankToDomain[_rank];
}

};

#endif
