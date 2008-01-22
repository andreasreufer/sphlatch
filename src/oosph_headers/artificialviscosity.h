#ifndef OOSPHARTIFICIALVISCOSITY_H
#define OOSPHARTIFICIALVISCOSITY_H

namespace oosph
{

/**
 * \brief Artificial Viscosity Parent Class
 * 
 * This class uses the Barton Nackman Trick for static Polymorphism
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 */



template <typename T_LeafType, typename Sim_Trait>
class ArtificialViscosity
{
public:

    typedef ArtificialViscosity  self_type;
    typedef ArtificialViscosity& self_reference;
    typedef ArtificialViscosity* self_pointer;

    typedef T_LeafType  leaf_type;
    typedef T_LeafType& leaf_reference;
    typedef T_LeafType* leaf_pointer;

    typedef typename Sim_Trait::value_type value_type;
    typedef typename Sim_Trait::vector_type vector_type;
    typedef typename Sim_Trait::matrix_row matrix_row;

    ArtificialViscosity(void);
    ~ArtificialViscosity(void);

    leaf_reference AsLeaf(void);

    static void PreLocal(void);
    static value_type Local(const matrix_row& i, const matrix_row& j);

protected:

private:

};

/** \brief Standart Constructor */
template <typename T_LeafType, typename Sim_Trait>
ArtificialViscosity<T_LeafType, Sim_Trait>::ArtificialViscosity(void)
{}

/** \brief Standart Destructor */
template <typename T_LeafType, typename Sim_Trait>
ArtificialViscosity<T_LeafType, Sim_Trait>::~ArtificialViscosity(void)
{}

/** \brief Method to return the child_type reference.
 * 
 * This is used for the Barton Nackman Trick
*/
template <typename T_LeafType, typename Sim_Trait>
typename ArtificialViscosity<T_LeafType, Sim_Trait>::leaf_reference ArtificialViscosity<T_LeafType, Sim_Trait>::AsLeaf(void)
{
    return static_cast<leaf_reference>(*this);
}

/** \brief Method to Call the PreLocal Routin of the Child Class */
template <typename T_LeafType, typename Sim_Trait>
inline void ArtificialViscosity<T_LeafType, Sim_Trait>::PreLocal(void)
{
    return leaf_type::PreLocal();
}

/** \brief Method to Class the Local Routin of the Child Class*/
template <typename T_LeafType, typename Sim_Trait>
inline typename ArtificialViscosity<T_LeafType, Sim_Trait>::value_type ArtificialViscosity<T_LeafType, Sim_Trait>::Local(const matrix_row& i, const matrix_row& j)
{

    assert(i.index() < i.data().size1() );
    assert(j.index() < j.data().size1() );

    return leaf_type::Local(i, j);
}

}

#endif
