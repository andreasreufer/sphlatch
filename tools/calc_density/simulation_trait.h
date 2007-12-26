#ifndef OOSPH_SIMULATION_TRAIT_H
#define OOSPH_SIMULATION_TRAIT_H

#include <cstdlib>

#include <vector>
#include <list>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/conversion/converter.hpp>

namespace oosph
{

namespace num = boost::numeric::ublas;

/**

\brief This class provides a data type forwarding for the most important datatypes.
\author Pascal Bauer

*/

#ifdef OOSPH_SINGLE_PRECISION
template <typename ValueType = float, template <typename T> class MatrixType = num::matrix>
#else
template <typename ValueType = double, template <typename T> class MatrixType = num::matrix>
#endif
class SimulationTrait
{
public:

    typedef ValueType  value_type;
    typedef ValueType& value_reference;
    typedef ValueType* value_pointer;

    typedef num::vector<value_type> vector_type;
    typedef num::vector<value_type>& vector_reference;
    typedef num::vector<value_type>* vector_pointer;
    
    typedef num::vector_range<vector_type> vector_range;

    typedef num::zero_vector<value_type> zero_vector;

    typedef MatrixType<value_type> matrix_type;
    typedef MatrixType<value_type>& matrix_reference;
    typedef MatrixType<value_type>* matrix_pointer;

    typedef num::zero_matrix<value_type> zero_matrix;

    typedef num::matrix_row<matrix_type> matrix_row;
    typedef const num::matrix_row<matrix_type> const_matrix_row;
    typedef num::matrix_column<matrix_type> matrix_column;
    typedef const num::matrix_column<matrix_type> const_matrix_column;

    typedef num::vector_range<const matrix_row> const_matrix_row_range;

    typedef num::vector_range<matrix_row> matrix_row_range;
    typedef num::vector_range<matrix_column> matrix_column_range;

    typedef num::matrix_vector_slice<matrix_type> matrix_vector_slice;
    typedef num::matrix_vector_slice<const matrix_type> const_matrix_vector_slice;

    typedef num::range range;
    typedef num::slice slice;

    typedef std::vector<size_t> index_vector_type;
    typedef std::vector<size_t>& index_vector_reference;
    typedef std::vector<size_t>* index_vector_pointer;

    typedef std::list<size_t> index_list_type;
    typedef std::list<size_t>& index_list_reference;
    typedef std::list<size_t>* index_list_pointer;

    typedef std::pair<bool, size_t> particle_index_type;
    typedef std::pair<bool, size_t>& particle_index_reference;
    typedef std::pair<bool, size_t>* particle_index_pointer;

    typedef std::vector<particle_index_type> particle_index_vector_type;
    typedef std::vector<particle_index_type>& particle_index_vector_reference;
    typedef std::vector<particle_index_type>* particle_index_vector_pointer;

    typedef std::list<size_t>  particle_index_list_type;
    typedef std::list<size_t>& particle_index_list_reference;
    typedef std::list<size_t>* particle_index_list_pointer;

    typedef std::vector<matrix_row>  particle_container_type;
    typedef std::vector<matrix_row>& particle_container_reference;
    typedef std::vector<matrix_row>* particle_container_pointer;

    typedef std::string string_type;
    typedef std::string& string_reference;
    typedef std::string* string_pointer;

    typedef boost::numeric::converter<int, value_type> converter_type;

protected:
private:
};

};

#endif
