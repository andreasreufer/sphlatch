#ifndef OOSPH_MEMORYMANAGER_H
#define OOSPH_MEMORYMANAGER_H

#include <map>
#include <set>

#include "particle.h"
#include "simulation_trait.h"

#include <boost/assign/std/vector.hpp>

namespace oosph
{

namespace num = boost::numeric::ublas;

template <typename SimTraitType>
class MemoryManager
{
public:

    typedef MemoryManager self_type;
    typedef MemoryManager& self_reference;
    typedef MemoryManager* self_pointer;

    typedef SimTraitType SimulationTrait;

    typedef typename SimTraitType::matrix_type matrix_type;
    typedef typename SimTraitType::matrix_reference matrix_reference;
    typedef typename SimTraitType::matrix_pointer matrix_pointer;
    typedef typename SimTraitType::matrix_row matrix_row;

    typedef typename SimTraitType::value_type value_type;
    typedef typename SimTraitType::vector_type vector_type;
    typedef typename SimTraitType::string_type string_type;
    
    typedef typename SimTraitType::index_vector_type index_vector_type;

    typedef std::map<std::string, value_type> StrValMap;

    matrix_type Data;
    matrix_type GData;
    matrix_type Calc;

    string_type Name;

    const size_t size(void);
    const size_t gsize(void);
    const size_t particle_size(void);

    /** PopParticle pops the particle with
     *  a given row index or particle ID */ 
    vector_type PopParticle(value_type ParticleID);
    vector_type PopParticle(size_t ParticleIndex);

    /** PopParticles pops particles with
     *  given row indeces or particle ID */ 
    matrix_type PopParticles(value_type ParticleID);
    matrix_type PopParticles(index_vector_type ParticleIndices);
    
    /** DelParticle(s) deletes one or more 
     *  particle with given row indices */
    void DelParticle(size_t ParticleIndex);

    void DelParticles(index_vector_type ParticleIndices);

    /** PushParticle() pushes a particle
     *  to the end of Data and returns the row index */
    size_t PushParticle(vector_type Particle);

    bool SaveParameter(std::string ParamName, value_type ParamValue, bool QueueParam);
    value_type LoadParameter(std::string ParamName);
    std::set<std::string> DumpParameter();

    size_t Usage(void);

    static self_reference Instance(void);

protected:

    MemoryManager(void);
    ~MemoryManager(void);

private:

    static self_pointer _instance;
    StrValMap ParamMap;
    std::set<std::string> ParamSave;

};

template <typename SimTraitType>
typename MemoryManager<SimTraitType>::self_pointer MemoryManager<SimTraitType>::_instance = NULL;

template <typename SimTraitType>
typename MemoryManager<SimTraitType>::self_reference MemoryManager<SimTraitType>::Instance(void)
{
    if (_instance != NULL)
        return *_instance;
    else
    {
        _instance = new MemoryManager;
        return *_instance;
    }
}

template <typename SimTraitType>
MemoryManager<SimTraitType>::MemoryManager(void)
{}

template <typename SimTraitType>
MemoryManager<SimTraitType>::~MemoryManager(void)
{}

template <typename SimTraitType>
inline size_t MemoryManager<SimTraitType>::Usage(void)
{
  return sizeof(Data) + sizeof(GData);
}

template <typename SimTraitType>
inline const size_t MemoryManager<SimTraitType>::size(void)
{
  return Data.size1();
}

template <typename SimTraitType>
inline const size_t MemoryManager<SimTraitType>::gsize(void)
{
  return GData.size1();
}

template <typename SimTraitType>
inline const size_t MemoryManager<SimTraitType>::particle_size(void)
{
  return Data.size2();
}

template <typename SimTraitType>
typename SimTraitType::value_type
MemoryManager<SimTraitType>::LoadParameter(std::string ParamName)
{
  if (ParamMap.count(ParamName) > 0) {
	  return ParamMap[ParamName];
  } else {
	  return std::numeric_limits<value_type>::quiet_NaN();
  }
}

template <typename SimTraitType>
bool MemoryManager<SimTraitType>::SaveParameter(std::string ParamName, value_type ParamValue, bool QueueParam)
{
  ParamMap[ParamName] = ParamValue;
  if (QueueParam == true) {
	  ParamSave.insert(ParamName);
  }
  return true;
}

template <typename SimTraitType>
std::set<std::string> MemoryManager<SimTraitType>::DumpParameter()
{
  return ParamSave;
}

template <typename SimTraitType>
typename SimTraitType::vector_type
MemoryManager<SimTraitType>::PopParticle(value_type ParticleID) {
	
	const size_t NoParticles = Data.size1();
	bool PartFound = false;
	size_t i = 0;
	vector_type PopParticle;

	if ( NoParticles > 0 ) {
		while( !PartFound && i < NoParticles ) {
			if ( Data(i, ID) == ParticleID ) {
				PartFound = true;

				matrix_row MatchedParticle(Data, i);
				matrix_row LastParticle(Data, NoParticles - 1);
				PopParticle = MatchedParticle;

				LastParticle.swap(MatchedParticle);
				Data.resize( NoParticles - 1, Data.size2(), true);
			}
			i++;
		}
	}

	return PopParticle;

}

template <typename SimTraitType>
typename SimTraitType::vector_type
MemoryManager<SimTraitType>::PopParticle(size_t ParticleIndex) {
	
	const size_t NoParticles = Data.size1();
	vector_type PopParticle;

	if ( ParticleIndex < NoParticles ) {
		matrix_row MatchedParticle(Data, ParticleIndex);
		matrix_row LastParticle(Data, NoParticles - 1);
		PopParticle = MatchedParticle;

		LastParticle.swap(MatchedParticle);
		Data.resize( NoParticles - 1, Data.size2(), true);
	}

	return PopParticle;
}

template <typename SimTraitType>
typename SimTraitType::matrix_type
MemoryManager<SimTraitType>::PopParticles(index_vector_type ParticleIndices) {
	
	const size_t NoParticles = Data.size1();
	const size_t NoVars = Data.size2();
	const size_t NoPopParticles = ParticleIndices.size();

	matrix_type PoppedData;
	PoppedData.resize(NoPopParticles, NoVars);

	size_t ParticlesPopped = 0;

	std::sort( ParticleIndices.begin(), ParticleIndices.end(), std::greater<size_t>() );
	
	for (size_t i = 0; i < NoPopParticles; i++) {
		const size_t ParticleIndex = ParticleIndices[i];
		if ( ParticleIndex < NoParticles ) {
			matrix_row MatchedParticle(Data, ParticleIndex);
			matrix_row(PoppedData, ParticlesPopped) = MatchedParticle;
			matrix_row LastParticle(Data, NoParticles - 1 - ParticlesPopped);
			LastParticle.swap(MatchedParticle);
			ParticlesPopped++;
		}
	}

	Data.resize( NoParticles - ParticlesPopped, NoVars, true);
	PoppedData.resize( ParticlesPopped, NoVars, true);

	return PoppedData;
}


template <typename SimTraitType>
typename SimTraitType::matrix_type
MemoryManager<SimTraitType>::PopParticles(value_type ParticleID) {
	
	using namespace boost::assign;

	const size_t NoParticles = Data.size1();
	const size_t NoVars = Data.size2();
	
	index_vector_type ParticleIndices;

	for (int i = NoParticles - 1; i >= 0; --i) {
		if ( lrint( Data(i, ID) ) == lrint(ParticleID) ) {
			ParticleIndices += i;
		}
	}

	const size_t NoPopParticles = ParticleIndices.size();

	matrix_type PoppedData;
	PoppedData.resize(NoPopParticles, NoVars);

	for (size_t i = 0; i < NoPopParticles; i++) {
		const size_t ParticleIndex = ParticleIndices[i];
		
		matrix_row MatchedParticle(Data, ParticleIndex);
		matrix_row(PoppedData, i) = MatchedParticle;
		matrix_row LastParticle(Data, NoParticles - 1 - i);
		LastParticle.swap(MatchedParticle);
	}

	Data.resize( NoParticles - NoPopParticles, NoVars, true);
	PoppedData.resize( NoPopParticles, NoVars, true);

	return PoppedData;
}


template <typename SimTraitType>
void MemoryManager<SimTraitType>::DelParticle(size_t ParticleIndex) {
	
	const size_t NoParticles = Data.size1();
	
	if ( ParticleIndex < NoParticles ) {
		matrix_row MatchedParticle(Data, ParticleIndex);
		matrix_row LastParticle(Data, NoParticles - 1);
		LastParticle.swap(MatchedParticle);
		Data.resize( NoParticles - 1, Data.size2(), true);
	}

}

template <typename SimTraitType>
void MemoryManager<SimTraitType>::DelParticles(index_vector_type ParticleIndices) {

	const size_t NoParticles = Data.size1();
	const size_t NoDelParticles = ParticleIndices.size();
	size_t ParticlesDeleted = 0;

	std::sort( ParticleIndices.begin(), ParticleIndices.end(), std::greater<size_t>() );
	
	for (size_t i = 0; i < NoDelParticles; i++) {
		const size_t ParticleIndex = ParticleIndices[i];
		if ( ParticleIndex < NoParticles ) {
			matrix_row MatchedParticle(Data, ParticleIndex);
			matrix_row LastParticle(Data, NoParticles - 1 - ParticlesDeleted);
			LastParticle.swap(MatchedParticle);
			ParticlesDeleted++;
		}
	}
	Data.resize( NoParticles - ParticlesDeleted, Data.size2(), true);

}

template <typename SimTraitType>
size_t MemoryManager<SimTraitType>::PushParticle(vector_type Particle) {
	
	const size_t InsertionIndex = Data.size1();
	const size_t PARTSIZE = Data.size2();
	
	Data.resize( InsertionIndex + 1, PARTSIZE, true);
	Particle.resize( PARTSIZE );
	
	matrix_row (Data, InsertionIndex) = Particle;

	return InsertionIndex;
}

};

#endif
