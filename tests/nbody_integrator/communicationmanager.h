#ifndef OOSPH_COMMUNICATION_MANAGER_H
#define OOSPH_COMMUNICATION_MANAGER_H

#include <cstdlib>

#include <vector>
#include <queue>
#include <mpi.h>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/assign/std/vector.hpp>

#include "mpi_help.h"
#include "memorymanager.h"

namespace oosph
{

namespace mpl = boost::mpl;

template <typename SimTraitType>
class CommunicationManager
{
public:
    typedef CommunicationManager  self_type;
    typedef CommunicationManager& self_reference;
    typedef CommunicationManager* self_pointer;

    typedef SimTraitType SimulationTrait;

    typedef typename SimulationTrait::matrix_type matrix_type;
    typedef typename SimulationTrait::matrix_reference matrix_reference;
    typedef typename SimulationTrait::matrix_pointer matrix_pointer;

    typedef typename SimulationTrait::matrix_column matrix_column;
    typedef typename SimulationTrait::matrix_row matrix_row;
    
    typedef typename SimulationTrait::matrix_vector_slice matrix_vector_slice;
    typedef typename SimulationTrait::const_matrix_vector_slice const_matrix_vector_slice;
    typedef typename SimulationTrait::slice slice;

    typedef typename SimulationTrait::value_type value_type;
    typedef typename SimulationTrait::value_reference value_reference;
    typedef typename SimulationTrait::value_pointer value_pointer;

    typedef typename SimulationTrait::matrix_type buffer_type;
    typedef typename SimulationTrait::matrix_reference buffer_reference;
    typedef typename SimulationTrait::matrix_pointer buffer_pointer;
    
    typedef MPI_SingleCom<typename SimulationTrait::value_type> single_com_type;

    typedef typename SimTraitType::index_vector_type index_vector_type;
    typedef typename SimTraitType::index_vector_reference index_vector_reference;
    typedef typename SimTraitType::index_vector_pointer index_vector_pointer;

    typedef typename SimTraitType::index_list_type index_list_type;
    typedef typename SimTraitType::index_list_reference index_list_reference;
    typedef typename SimTraitType::index_list_pointer index_list_pointer;

    typedef std::vector<index_vector_type> domain_index_vector_type;
    typedef std::vector<index_vector_type>& domain_index_vector_reference;
    typedef std::vector<index_vector_type>* domain_index_vector_pointer;

    typedef std::vector<index_list_type> domain_index_list_type;
    typedef std::vector<index_list_type>& domain_index_list_reference;
    typedef std::vector<index_list_type>* domain_index_list_pointer;

    typedef MemoryManager<SimulationTrait> memory_manager_type;
    typedef MemoryManager<SimulationTrait>& memory_manager_reference;
    typedef MemoryManager<SimulationTrait>* memory_manager_pointer;

    template <int N> static value_type Min(void);
    template <int N> static value_type Max(void);

    static self_reference Instance(void);

    void Exchange(matrix_reference DataSource, domain_index_vector_reference Index, matrix_reference DataTarget);

protected:

    CommunicationManager(void);
    ~CommunicationManager(void);

private:

    static self_pointer _instance;
    static single_com_type SingleCom;
    static memory_manager_reference MemManager;
    static matrix_reference Data;
    static matrix_reference GData;

    static size_t CommBuffSize;

    std::vector<buffer_type> ComposeBuffers;
    std::vector<buffer_type> RecvBuffs;
    buffer_type SendBuff;

    std::vector<size_t> SendNoPart;
    std::vector<size_t> RecvNoPart;

    std::queue< std::pair<size_t, std::vector<size_t> > > SendQueue;
};

template <typename SimTraitType>
typename CommunicationManager<SimTraitType>::self_pointer CommunicationManager<SimTraitType>::_instance = NULL;

template <typename SimTraitType>
typename CommunicationManager<SimTraitType>::self_reference CommunicationManager<SimTraitType>::Instance(void)
{
    if (_instance != NULL)
        return *_instance;
    else
    {
        _instance = new CommunicationManager;
        return *_instance;
    }
}

template <typename SimTraitType>
typename CommunicationManager<SimTraitType>::single_com_type CommunicationManager<SimTraitType>::SingleCom;

template <typename SimTraitType>
typename CommunicationManager<SimTraitType>::memory_manager_reference
CommunicationManager<SimTraitType>::MemManager(memory_manager_type::Instance() );

template <typename SimTraitType>
typename CommunicationManager<SimTraitType>::matrix_reference
CommunicationManager<SimTraitType>::Data(MemManager.Data);

template <typename SimTraitType>
typename CommunicationManager<SimTraitType>::matrix_reference
CommunicationManager<SimTraitType>::GData(MemManager.Data);

template <typename SimTraitType>
size_t CommunicationManager<SimTraitType>::CommBuffSize = 128;	

template <typename SimTraitType>
CommunicationManager<SimTraitType>::CommunicationManager(void)
{}

template <typename SimTraitType>
CommunicationManager<SimTraitType>::~CommunicationManager(void)
{}

template <typename SimTraitType>
template <int N>
typename CommunicationManager<SimTraitType>::value_type CommunicationManager<SimTraitType>::Min(void)
{
    static const matrix_column Column(Data, N);
    value_type local = *std::min_element(Column.begin(), Column.end() );
    value_type global;
    SingleCom.min(local, global);
    return global;
};

template <typename SimTraitType>
template <int N>
typename CommunicationManager<SimTraitType>::value_type
CommunicationManager<SimTraitType>::Max(void)
{
    static const matrix_column Column(Data, N);
    value_type local = *std::max_element(Column.begin(), Column.end() );
    value_type global;
    SingleCom.max(local, global);
    return global;
}

template <typename SimTraitType>
void CommunicationManager<SimTraitType>::Exchange(matrix_reference DataSource, domain_index_vector_reference Index, matrix_reference DataTarget)
{
    const size_t RANK = MPI::COMM_WORLD.Get_rank();
    const size_t SIZE = MPI::COMM_WORLD.Get_size();

    /** Only common range will be exchanged, 
     * DataTarget will be resized to that size */
    const size_t PARTSIZE = std::min(DataSource.size2(), DataTarget.size2());
    const size_t CommBuffByteSize = CommBuffSize * PARTSIZE * sizeof(value_type);

    /** Slice of matrix which will be copied:
     * first index is 0, iterate through
     * every index until PARTSIZE */
    slice PartSlice(0,1,PARTSIZE);

    /** Exchange knowledge about who is sending whom how much */
    SendNoPart.resize(SIZE);
    RecvNoPart.resize(SIZE);

    assert(Index.size() == SIZE);
    for (size_t i = 0; i < SIZE; i++) {
	    SendNoPart[i] = Index[i].size();
    }

    MPI::COMM_WORLD.Alltoall(&SendNoPart[0], 1, MPI::UNSIGNED_LONG, &RecvNoPart[0], 1, MPI::UNSIGNED_LONG);

    /** Fill up send queue */
    for (size_t t = 1; t < SIZE; t++) {
	     const size_t CurSendRank = ( RANK + t ) % SIZE;
	     const size_t NoSendParts = Index[CurSendRank].size();
	     
	     using namespace boost::assign;
	     std::vector<size_t> Indices;

	     for (size_t i = 0; i < NoSendParts; i++) {
		     Indices += Index[CurSendRank][i];
		     if ( Indices.size() == CommBuffSize || i == NoSendParts - 1 ) {
			     std::pair<size_t, std::vector<size_t> > RankIndicesPair;
			     RankIndicesPair.first = CurSendRank;
			     RankIndicesPair.second = Indices;
			     SendQueue.push(RankIndicesPair);
			     Indices.resize(0);
		     }
	     }
    }

    /** Prepare ComposeBuffers */
    ComposeBuffers.resize( SIZE );
    for (size_t i = 0; i < SIZE; i++) {
	    ComposeBuffers[i].resize( RecvNoPart[i], PARTSIZE );
    }

    /** Copy particles staying on node to ComposeBuffer */
    for (size_t i = 0; i < Index[RANK].size(); i++) {
	    matrix_row(ComposeBuffers[RANK], i) = const_matrix_vector_slice(DataSource, slice(Index[RANK][i], 0, PARTSIZE), PartSlice);

    }
    RecvNoPart[RANK] = 0;

    /** Prepare Communication */
    std::vector<MPI::Request> RecvReq;
    RecvReq.resize(SIZE);
    RecvBuffs.resize(SIZE);

    for (size_t i = 0; i < SIZE; i++) {
	    if ( RecvNoPart[i] > 0 ) {
		    (RecvBuffs[i]).resize( CommBuffSize, PARTSIZE );
		    RecvReq[i] = MPI::COMM_WORLD.Irecv(&RecvBuffs[i](0, 0), CommBuffByteSize, MPI::BYTE, i, i);
	    }
    }

    SendBuff.resize(CommBuffSize, PARTSIZE);
   
    /** Start Communication! */
    bool CommFinished = false;
    
    std::vector<size_t> PartsRecvd(SIZE);
    for (size_t i = 0; i < SIZE; i++) {
	    PartsRecvd[i] = 0;
    }

    while (!CommFinished) {
	    /* Assume we have finished, 
	     * prove otherwise */
	    CommFinished = true;

	    /** Check every receiver for arrived particles */
	    for (size_t i = 0; i < SIZE; i++) {
		    /** Still expecting particles and buffer is filled */
		    if ( RecvNoPart[i] > 0 && RecvReq[i].Test() ) {
			    /** Flush the receiving buffer... */
			    size_t NoFlushs = std::min(CommBuffSize, RecvNoPart[i]);
			    for (size_t j = 0; j < NoFlushs; j++) {
				    matrix_row(ComposeBuffers[i], PartsRecvd[i]) = matrix_row(RecvBuffs[i], j);
				    PartsRecvd[i]++;
			    }
			    RecvNoPart[i] -= NoFlushs;

			    /** .. and setup new Receiver if there's more expected */
			    if ( RecvNoPart[i] > 0 ) {
				    RecvReq[i] = MPI::COMM_WORLD.Irecv(&RecvBuffs[i](0, 0), CommBuffByteSize, MPI::BYTE, i, i);
			    }
		    }
	    }

	    /** If SendQueue not empty send something */
	    if ( SendQueue.size() ) {
		    /** Pick particles to be sent and put them into the sending buffer */
		    std::pair<size_t, std::vector<size_t> > SendPair = SendQueue.front();
		    SendQueue.pop();
		    const size_t SendRank = SendPair.first;
		    std::vector<size_t> SendIndices = SendPair.second;
		    const size_t NoSends = SendIndices.size();

		    for (size_t i = 0; i < NoSends; i++) {
			    matrix_row(SendBuff, i) = const_matrix_vector_slice(DataSource, slice(SendIndices[i], 0, PARTSIZE), PartSlice);
		    }

		    /** Fire the buffer! (Synchronous send) */
		    MPI::COMM_WORLD.Ssend(&SendBuff(0, 0), CommBuffByteSize, MPI::BYTE, SendRank, RANK);
		    CommFinished = false;
	    }
	    
	    /** Check whether communication has finished */
	    for (size_t i = 0; i < SIZE; i++) {
		    if ( RecvNoPart[i] > 0 ) {
			    CommFinished = false;
		    }
	    }
    }

    /** Determine number of particles in ComposeBuffer */
    size_t CompParts = 0;
    for (size_t i = 0; i < SIZE; i++) {
	    CompParts += ComposeBuffers[i].size1();
    }

    /** Prepare target matrix */
    DataTarget.resize( CompParts, PARTSIZE );

    /** Flush ComposeBuffers */
    size_t PartCounter = 0;
    for (size_t i = 0; i < SIZE; i++) {
	   const size_t BuffParts = ComposeBuffers[i].size1();
	   for (size_t j = 0; j < BuffParts; j++) {
		   matrix_row(DataTarget, PartCounter) = matrix_row(ComposeBuffers[i], j);
		   PartCounter++;
	   }
	   ComposeBuffers[i].resize(0,0);
    }
    
    /** Be polite and wait for the others
     *  (it's useful too by the way) */
    MPI::COMM_WORLD.Barrier();
    
    return;
}

}
;

#endif
