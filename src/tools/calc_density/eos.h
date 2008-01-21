#ifndef OOSPHEOS_H
#define OOSPHEOS_H

namespace oosph
{

/**
 * \brief Kernel Parent Class
 * 
 * \author Pascal Bauer pbauer@phim.unibe.ch
 */

template <typename T_LeafType, typename Sim_Trait>
class EOS
{
public:

    typedef EOS  self_type;
    typedef EOS& self_reference;
    typedef EOS* self_pointer;

    typedef T_LeafType  leaf_type;
    typedef T_LeafType& leaf_reference;
    typedef T_LeafType* leaf_pointer;

    typedef Sim_Trait simulation_trait;

    typedef typename simulation_trait::value_type value_type;
    typedef typename simulation_trait::matrix_column matrix_column;


    EOS(void);
    ~EOS(void);

    static void Pressure(void);

protected:

private:
};

template <typename T_LeafType, typename Sim_Trait>
EOS<T_LeafType, Sim_Trait>::EOS(void)
{}

template <typename T_LeafType, typename Sim_Trait>
EOS<T_LeafType, Sim_Trait>::~EOS(void)
{}

template <typename T_LeafType, typename Sim_Trait>
inline void EOS<T_LeafType, Sim_Trait>::Pressure(void)
{
    return leaf_type::Pressure();
}

}

#endif
