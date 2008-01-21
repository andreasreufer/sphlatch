#ifndef OOSPH_DOMAINMAP_H
#define OOSPH_DOMAINMAP_H

#include <cstdlib>
#include <vector>

#include <mpi.h>

namespace oosph
{

class DomainMap
{
public:

    typedef DomainMap self_type;
    typedef DomainMap& self_reference;
    typedef DomainMap* self_pointer;

    typedef std::vector<size_t> domain_map_type;
    typedef std::vector<size_t>& domain_map_reference;
    typedef std::vector<size_t>* domain_map_pointer;

    typedef domain_map_type::iterator iterator;
    typedef domain_map_type::const_iterator const_iterator;

    DomainMap(void)
    {}

    DomainMap(const self_reference src)
    {
        domain_map = src.domain_map;
    }

    self_reference operator=(const self_reference src)
    {
        domain_map = src.domain_map;
        return *this;
    }

    iterator begin(void)
    {
        return domain_map.begin();
    }

    iterator end(void)
    {
        return domain_map.end();
    }

    void clear(void)
    {
        domain_map.clear();
    }

    void resize(const size_t size)
    {
        domain_map.resize(size);
    }

    size_t size(void)
    {
        return domain_map.size();
    }

    size_t& operator[](const size_t index)
    {
        size_t& ret = domain_map[index];
        assert(ret < (size_t)MPI::COMM_WORLD.Get_size() );
        return ret;
    }



protected:
    domain_map_type domain_map;

private:

};

};

#endif
