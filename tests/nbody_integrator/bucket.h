#ifndef OOSPH_BUCKET_H
#define OOSPH_BUCKET_H

#include <vector>

#include <list>
namespace oosph
{

template <typename T, template <typename T> class C = std::list >
class Bucket
{
public:

    typedef Bucket self_type;
    typedef Bucket& self_reference;
    typedef Bucket* self_pointer;

    typedef T value_type;
    typedef T& value_reference;
    typedef T* value_pointer;

    typedef C<T> cell_type;
    typedef C<T>& cell_reference;
    typedef C<T>* cell_pointer;

    typedef std::vector<cell_type> container_type;
    typedef std::vector<cell_type>& container_reference;
    typedef std::vector<cell_type>* container_pointer;

    typedef typename container_type::iterator iterator;

    iterator begin(void);
    iterator end(void);

    Bucket(void);
    Bucket(const self_reference src);
    ~Bucket(void);

    void clear();
    void resize(const size_t i);
    size_t size(void);

    self_reference operator=(const self_reference src);
    cell_reference operator()(const size_t index);
    cell_reference operator [](const size_t index);

protected:

private:

    container_type bucket_data;

}
;

template <typename T, template <typename T> class C>
Bucket<T, C>::Bucket(void)
{}

template <typename T, template <typename T> class C>
Bucket<T, C>::~Bucket(void)
{}

template <typename T, template <typename T> class C>
Bucket<T, C>::Bucket(const self_reference src)
{
    bucket_data = src.bucket_data;
}

template <typename T, template <typename T> class C>
typename Bucket<T, C>::self_reference Bucket<T, C>::operator=(const self_reference src)
{
    bucket_data = src.bucket_data;
    return *this;
}

template <typename T, template <typename T> class C>
void Bucket<T, C>::clear()
{
    bucket_data.clear();
}

template <typename T, template <typename T> class C>
void Bucket<T, C>::resize(const size_t i)
{
    bucket_data.resize(i);
}

template <typename T, template <typename T> class C>
inline typename Bucket<T, C>::cell_reference
Bucket<T, C>::operator()(const size_t index)
{
    return bucket_data[index];
}

template <typename T, template <typename T> class C>
inline typename Bucket<T, C>::cell_reference
Bucket<T, C>::operator[](const size_t index)
{
    return bucket_data[index];
}

template <typename T, template <typename T> class C>
inline typename Bucket<T, C>::iterator Bucket<T, C>::begin(void)
{
  return bucket_data.begin();
}

template <typename T, template <typename T> class C>
inline typename Bucket<T, C>::iterator Bucket<T, C>::end(void)
{
  return bucket_data.end();
}

template <typename T, template <typename T> class C>
inline size_t Bucket<T, C>::size(void)
{
  return bucket_data.size();
}

}
;

#endif
