#ifndef OOSPHKERNEL_H
#define OOSPHKERNEL_H

namespace oosph
{

  /**
   * \brief Kernel Parent Class
   * 
   * \author Pascal Bauer pbauer@phim.unibe.ch
   * 
   */
template <typename T_LeafType, typename Sim_Trait>
class Kernel
{
public:

    typedef Sim_Trait sph_trait;

    typedef Kernel  self_type;
    typedef Kernel& self_reference;
    typedef Kernel* self_pointer;

    typedef T_LeafType  leaf_type;
    typedef T_LeafType& leaf_reference;
    typedef T_LeafType* leaf_pointer;
    
    typedef typename Sim_Trait::value_type value_type;
    typedef typename Sim_Trait::vector_type vector_type;
    typedef typename Sim_Trait::matrix_row_range matrix_row_range;
    typedef typename Sim_Trait::const_matrix_row_range const_matrix_row_range;

    typedef typename Sim_Trait::matrix_row matrix_row;
    
    typedef typename Sim_Trait::zero_vector zero_vector;

    Kernel(void);
    ~Kernel(void);
    
    leaf_reference AsLeaf(void);

    static value_type K( const matrix_row_range& i, const matrix_row_range& j);
    static vector_type Grad( const matrix_row_range& i, const matrix_row_range& j);

protected:

private:

};

template <typename T_LeafType, typename Sim_Trait>
Kernel<T_LeafType, Sim_Trait>::Kernel(void)
{}

template <typename T_LeafType, typename Sim_Trait>
Kernel<T_LeafType, Sim_Trait>::~Kernel(void)
{}

template <typename T_LeafType, typename Sim_Trait>
    inline typename Kernel<T_LeafType, Sim_Trait>::value_type Kernel<T_LeafType, Sim_Trait>::K(const matrix_row_range& i, const matrix_row_range& j)
{
  return leaf_type::K(i, j);
};

template <typename T_LeafType, typename Sim_Trait>
    inline typename Kernel<T_LeafType, Sim_Trait>::vector_type Kernel<T_LeafType, Sim_Trait>::Grad( const matrix_row_range& i, const matrix_row_range& j)
{
  return leaf_type::Grad(i, j);
}

template <typename T_LeafType, typename Sim_Trait>
    inline typename Kernel<T_LeafType, Sim_Trait>::leaf_reference Kernel<T_LeafType, Sim_Trait>::AsLeaf(void)
{
  return static_cast<leaf_reference>(*this);
}

}

#endif
