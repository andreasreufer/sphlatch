#ifndef OOSPH_LOADBALANCE_H
#define OOSPH_LOADBALANCE_H

#include <cstdlib>
#include <vector>
#include <list>

namespace oosph
{

template <typename T_LeafType, typename SimTraitType>
class LoadBalance
{
public:

    typedef LoadBalance self_type;
    typedef LoadBalance& self_reference;
    typedef LoadBalance* self_pointer;

    typedef T_LeafType leaf_type;
    typedef T_LeafType& leaf_reference;
    typedef T_LeafType* leaf_pointer;

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

    LoadBalance(void);
    ~LoadBalance(void);

    leaf_reference AsLeaf(void);

    static domain_index_vector_reference CreateDomainIndexVector(void);
    static domain_index_list_reference   CreateDomainIndexList(void);

    static domain_index_vector_reference CreateDomainGhostIndexVector(void);
    static domain_index_list_reference   CreateDomainGhostIndexList(void);

    static domain_index_vector_type div; // domain index vector
    static domain_index_list_type   dil; // domain index list

    static domain_index_vector_type gdiv; // ghost domain index vector;
    static domain_index_list_type   gdil; // ghost domain index list;

protected:


private:

};

template <typename T_LeafType, typename SimTraitType>
LoadBalance<T_LeafType, SimTraitType>::LoadBalance(void)
{
}

template <typename T_LeafType, typename SimTraitType>
LoadBalance<T_LeafType, SimTraitType>::~LoadBalance(void)
{}

template <typename T_LeafType, typename SimTraitType>
inline typename LoadBalance<T_LeafType, SimTraitType>::leaf_reference LoadBalance<T_LeafType, SimTraitType>::AsLeaf(void)
{
    return static_cast<leaf_reference>(*this);
}

template <typename T_LeafType, typename SimTraitType>
inline typename LoadBalance<T_LeafType, SimTraitType>::domain_index_vector_reference
LoadBalance<T_LeafType, SimTraitType>::CreateDomainIndexVector(void)
{
    std::cout << "parent: " << std::endl;
    AsLeaf().CreateDomainIndexVector();

    return div;
}

template <typename T_LeafType, typename SimTraitType>
inline typename LoadBalance<T_LeafType, SimTraitType>::domain_index_list_reference
LoadBalance<T_LeafType, SimTraitType>::CreateDomainIndexList(void)
{
    return T_LeafType::CreateDomainIndexList();
}

template <typename T_LeafType, typename SimTraitType>
inline typename LoadBalance<T_LeafType, SimTraitType>::domain_index_vector_reference
LoadBalance<T_LeafType, SimTraitType>::CreateDomainGhostIndexVector(void)
{
    return T_LeafType::CreateDomainGhostIndexVector();
}

template <typename T_LeafType, typename SimTraitType>
inline typename LoadBalance<T_LeafType, SimTraitType>::domain_index_list_reference
LoadBalance<T_LeafType, SimTraitType>::CreateDomainGhostIndexList(void)
{
    return T_LeafType::CreatDomainGhostIndexList();
}

template <typename T_LeafType, typename SimTraitType>
typename LoadBalance<T_LeafType, SimTraitType>::domain_index_vector_type
LoadBalance<T_LeafType, SimTraitType>::div;

template <typename T_LeafType, typename SimTraitType>
typename LoadBalance<T_LeafType, SimTraitType>::domain_index_list_type
LoadBalance<T_LeafType, SimTraitType>::dil;

template <typename T_LeafType, typename SimTraitType>
typename LoadBalance<T_LeafType, SimTraitType>::domain_index_vector_type
LoadBalance<T_LeafType, SimTraitType>::gdiv;

template <typename T_LeafType, typename SimTraitType>
typename LoadBalance<T_LeafType, SimTraitType>::domain_index_list_type
LoadBalance<T_LeafType, SimTraitType>::gdil;

};

#endif
