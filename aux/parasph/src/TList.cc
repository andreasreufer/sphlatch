/******************************************************************************
 * ParaSPH -- Version 24.11.2003                                              *
 *----------------------------------------------------------------------------*
 * File:      TList.cc                                                        *
 * Purpose:   This is a template list class, ie. a list that contains objects *
 *            of choosable datatype. The list can grow and one can access its *
 *            elements directly. This is kind of a compromise between speed   *
 *            and object orientation. Sorry for this, but speed is more       *
 *            important here. 'size' is the number of elements that are in    *
 *            the list, 'maxSize' is the size that is allocated.              *
 *****************************************************************************/
#ifndef LIST_CC
#define LIST_CC

template <class T>
class TList {
private:
  T   *element;       // array of listelements
  int  size, maxSize;

public:
  TList() { size = 0; maxSize = 0; element = new T[maxSize]; }

  ~TList() { delete [] element; }
  
  // Get element number n
  T &operator[](const int &n) const { 
    assert(n == 0 || 0 < n && n < size); return element[n]; 
  }

  // Get size of list
  const int &getSize() const { return size; }

  // Change size of list (the calling class must read in s particles
  // into element[0..size-1])
  void setSize(const int &s) {
    if (s > maxSize) {
      size = s; maxSize = s;
      delete [] element;
      element = new T[size];
    } else size = s;
  }

  // Inserts an element into the list. You have to care about memory
  // management
  //int append(T *item) { return append(*item); }

  int append(const T &item) {
    assert(size < maxSize);
    element[size++] = item;
    return size-1;
  }

  // Same as append, but if 'maxSize' is too small, the list will grow.
  // Can be time consuming if size changes often.
  int safeAppend(const T &item) {
    if (size >= maxSize) grow(Max(10, maxSize >> 5));
    return append(item);
  }

  // Delete an element and fill the cavity with the last element of the list
  //void del(const int &n) {
    //assert(size > 0 && n >= 0 && n < size);
    //element[n] = element[--size];
  //}

  // Allocate more memory (will be filled using append())
  void grow(const int &s) {
    if (size + s > maxSize) {
      T *newelem = new T[size+s];
      for (int i = 0; i < size; i++) newelem[i] = element[i];
      delete [] element;
      element = newelem;
      maxSize = size + s;
    }
  }
};

#endif
