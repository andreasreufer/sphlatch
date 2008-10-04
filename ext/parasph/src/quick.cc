#ifndef QUICK_CC
#define QUICK_CC

template <class T>
class Stack{
private:
  T   ss[100];
  int sp;
public:
  Stack()        { sp = 0; }
  void push(T a) { assert(sp < 100); ss[sp++] = a; }
  T    pop()     { assert(sp > 0);   return ss[--sp]; }
  bool empty()   { return (sp == 0); }
};

template <class T>
void quick(T *ind, const int &size) {
  int        i, j, l, m, r;
  T          piv;
  Stack<int> s;
  
  l = 0; r = size-1;
  s.push(l); s.push(r);
  do {
    if (r > l) {
      m = (r+l) >> 1;
      if (ind[m] < ind[l]) std::swap(ind[m], ind[l]);
      if (ind[r] < ind[m]) {
	std::swap(ind[r], ind[m]);
	if (ind[m] < ind[l]) std::swap(ind[m], ind[l]);
      }
      piv = ind[m];
      
      i = l+1; j = r-1;
      do {
	while (ind[i] < piv) i++;
	while (piv < ind[j]) j--;
	if (i < j) {
	  std::swap(ind[i], ind[j]);
	  i++; j--;
	} else if (i == j) {
	  i++; j--;
	  break;
	}
      } while (i <= j);
      
      if ((j-l) > (r-i)) { s.push(l); s.push(j); l = i; }
      else               { s.push(i); s.push(r); r = j; }
    } else {
      r = s.pop(); l = s.pop();
    }
  } while (!s.empty());
}

template <class T>
void quick(T *list, int *index, const int &size) {
  int        i, j, l, m, r;
  T          piv;
  Stack<int> s;
  
  l = 0; r = size-1;
  s.push(l); s.push(r);
  do {
    if (r > l) {
      m = (r+l) >> 1;
      if (list[index[m]] < list[index[l]]) std::swap(index[m], index[l]);
      if (list[index[r]] < list[index[m]]) {
	std::swap(index[r], index[m]);
	if (list[index[m]] < list[index[l]]) std::swap(index[m], index[l]);
      }
      piv = list[index[m]];

      i = l+1; j = r-1;
      do {
	while (list[index[i]] < piv) i++;
	while (piv < list[index[j]]) j--;
	if (i < j) {
	  std::swap(index[i], index[j]);
	  i++; j--;
	} else if (i == j) {
	  i++; j--;
	  break;
	}
      } while (i <= j);
      
      if ((j-l) > (r-i)) { s.push(l); s.push(j); l = i; }
      else               { s.push(i); s.push(r); r = j; }
    } else {
      r = s.pop(); l = s.pop();
    }
  } while (!s.empty());
}

template <class T>
void quickStart(T *list, int *index, const int &size) {
  for (int i = 0; i < size; i++) index[i] = i;
  quick(list, index, size);
}

#endif
