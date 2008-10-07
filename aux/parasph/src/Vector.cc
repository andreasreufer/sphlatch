/******************************************************************************
 * ParaSPH -- Version 24.11.2003                                              *
 *----------------------------------------------------------------------------*
 * File:      Vector.cc                                                       *
 * Purpose:   This class has been written due to performance reason and since *
 *            there are some methods missing in the STL vector class.         *
 *****************************************************************************/
#ifndef VECTOR_CC
#define VECTOR_CC

#include <iostream>

#include "Def.cc"

template <class T>
class Vector {
private:
  T x[3];

public:
  Vector() { set(0, 0, 0); }
  
  Vector(const T &x0, const T &x1, const T &x2) { set(x0, x1, x2); }

  void set(const T &x0, const T &x1, const T &x2) {
    x[0] = x0; x[1] = x1; x[2] = x2;
  }

  const T &operator[](const int &i) const { return x[i]; }
  T &operator[](const int &i) { return x[i]; }

  Vector<T> operator+(const Vector<T> &other) const {
    return Vector<T>(x[0]+other.x[0], x[1]+other.x[1], x[2]+other.x[2]);
  }

  Vector<T> operator+(const T &num) const {
    return Vector<T>(x[0]+num, x[1]+num, x[2]+num);
  }

  Vector<T> operator+=(const Vector<T> &other) {
    x[0] += other.x[0]; x[1] += other.x[1]; x[2] += other.x[2];
    return *this;
  }

  Vector<T> operator-(const Vector<T> &other) const {
    return Vector<T>(x[0]-other.x[0], x[1]-other.x[1], x[2]-other.x[2]);
  }
  
  Vector<T> operator-(const T &num) const {
    return Vector<T>(x[0]-num, x[1]-num, x[2]-num);
  }

  Vector<T> operator-=(const Vector<T> &other) {
    x[0] -= other.x[0]; x[1] -= other.x[1]; x[2] -= other.x[2];
    return *this;
  }

  Vector<T> operator*(const T &num) const {
    return Vector<T>(x[0]*num, x[1]*num, x[2]*num);
  }

  T operator*(const Vector<T> &other) const {
    return x[0]*other.x[0] + x[1]*other.x[1] + x[2]*other.x[2];
  }

  Vector<T> operator*=(const T &num) {
    x[0] *= num; x[1] *= num; x[2] *= num;
    return *this;
  }

  Vector<T> operator>>(const int &n) const {
    return Vector<T>(x[0]>>n, x[1]>>n, x[2]>>n);
  }
  Vector<T> operator<<(const int &n) const {
    return Vector<T>(x[0]<<n, x[1]<<n, x[2]<<n);
  }

  T len2() const { return x[0]*x[0]+x[1]*x[1]+x[2]*x[2]; }
  T len() const  { return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); }

  friend std::ostream &operator<<(std::ostream &out, const Vector<T> &v) {
    out << '(' << v.x[0] << ", " << v.x[1] << ", " << v.x[2] << ')';
    return out;
  }
};

template <class T>
class Range{
public:
  Vector<T> min, max;

  void set(Vector<T> _min, Vector<T> _max) { min = _min; max = _max; }

  friend std::ostream &operator<<(std::ostream &out, const Range<T> &r) {
    out << r.min << '-' << r.max;
    return out;
  }
};

Vector<ftype> maxVec(1.e30, 1.e30, 1.e30), minVec(-1.e30, -1.e30, -1.e30);

#endif
