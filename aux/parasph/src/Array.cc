#ifndef ARRAY_CC
#define ARRAY_CC

#include <assert.h>

template<class T>
class Array {
private:
  T   *element;
  int   max, size;
  
public:
  Array() { max = size = 0; element = new T[max]; }
  Array(const int &_max) { max = _max; size = 0; element = new T[max]; }

  ~Array() { delete [] element; }

  const int &getSize() const { return size; }

  // This function is the reason for bothering with my own array class:
  // We need to have writing access to elements of the array, even if they do
  // not exist yet, and we make sure that they all are sequential
  T &operator[](const int &n) { 
    if (n > 0 && n >= size) 
      std::cout << "array[" << n << "] does not exist " << size << std::endl;
    assert(n == 0 || (0 < n && n < size)); return element[n];
  }

  void setSize(const int &newSize) { 
    if (newSize > max) { 
      delete [] element;
      max = size = newSize;
      element = new T[max];
    } else size = newSize;
  }

  /*void allocate(const int &s) {
    if (s > max) {
      T *newelem = new T[s];
      for (int i = 0; i < size; i++) newelem[i] = element[i];
      delete [] element;
      element = newelem;
      max = s;
    }
  }*/

  void append(const T &item) { 
    if (size >= max) 
      std::cout << global::rank << " array max reached " << max << std::endl;
    assert(size < max); element[size++] = item; 
  }

  void clear(const T &z) { for (int i = 0; i < size; i++) element[i] = z; }

  void grow(const int &s) {
    if (size + s > max) {
      T *newelem = new T[size+s];
      for (int i = 0; i < size; i++) newelem[i] = element[i];
      delete [] element;
      element = newelem;
      max = size + s;
    }
  }

  void growTo(const int s) { grow(s - size); }
};

#endif
