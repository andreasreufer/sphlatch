// Input:  Particle.pos, Particle.work
// Output: Particle.proc

#ifndef DOMAINDECOMP_CC
#define DOMAINDECOMP_CC

#include "Array.cc"
#include "Debug.cc"
#include "Particle.cc"
#include "TList.cc"

class DomainDecomp {
private:
  struct posVal {
    Vector<ftype> pos; 
    float         val;
    int           nr;
  };

  Array<Particle>  buffer, oversize;
  Array<posVal>    index;
  TList<Particle> &pList;  

  void optimizeArray() {
    int k = 0, p = 0, s = 0, with;
    for (with = global::rank; p < global::npro-1; p++, s++) {
      for (int i = k; i < pList.getSize(); i++)
	if (pList[i].getProc() == with) 
	  std::swap(pList[i], pList[k++]);
      with = global::comm.get(s); 
      if (with == global::rank) with = global::comm.get(++s);
    }
  }

  void pSwapParticles(Particle *from, Particle *to, const int &nr, 
	 	      const int &with) {    
    buffer.setSize(nr);
    MCW.Sendrecv(from,       nr*sizeof(Particle), MPI::BYTE, with, 2*with+1,
	 	 &buffer[0], nr*sizeof(Particle), MPI::BYTE, with, 
		 2*global::rank+1);
    for (int i = 0; i < nr; i++) to[i] = buffer[i];
  }

  void pComputeTraffic(int *toSend, int *toRecv) {
    for (int i = 0; i < global::npro; i++) toSend[i] = 0;
    for (int p = 0; p < pList.getSize(); p++) 
      toSend[pList[p].getProc()]++;

    for (int i = 0; i < global::npro; i++) 
      MCW.Sendrecv(toSend + i, 1, MPI::INT, i, 2*i+1,
		   toRecv + i, 1, MPI::INT, i, 2*global::rank+1);
  }

  void pMinimizeTraffic() {
    DEBUG("pMinimizeTraffic()", "");
    bool stillIn = true;
    int *allMax = new int[global::npro], j, *map = new int[global::npro];
    int  max, p, pos, tot;
    int *toSend = new int[global::npro], *toRecv = new int[global::npro];

    pComputeTraffic(toSend, toRecv);

    for (j = 0; j < global::npro; j++) {
      pos = maxLoc(toSend, global::npro);
      if (stillIn) max = toSend[pos]; else max = -1;
      MCW.Allgather(&max, 1, MPI::INT, allMax, 1, MPI::INT);
      tot = maxLoc(allMax, global::npro);
      if (global::rank == tot) stillIn = false;
      MCW.Bcast(&pos, 1, MPI::INT, tot);
      map[pos] = tot; toSend[pos] = -1;
    }

    for (p = 0; p < pList.getSize(); p++) 
      pList[p].setProc(map[pList[p].getProc()]);
  }

  void pTidyUp() {
    int  excess, maxOver, min, n, newSize, sum, with;
    int *toSend = new int[global::npro], *toRecv = new int[global::npro];

    optimizeArray();

    pComputeTraffic(toSend, toRecv);

    maxOver = 0; newSize = toRecv[global::rank];
    for (n = 0, sum = 0; n < global::comm.getRounds(); n++) {
      with    = global::comm.get(n);
      sum     = Max(sum + toRecv[with] - toSend[with], 0); 
      maxOver = Max(sum, maxOver);
      if (global::rank != with) newSize += toRecv[with];
    }

    Particle *from = &pList[0] + toRecv[global::rank], *to = from;
    oversize.setSize(maxOver);
    Particle *over = &oversize[0];

    for (n = 0; n < global::comm.getRounds(); n++) {
      with = global::comm.get(n);
      if (with != global::rank) {
	excess = toSend[with] - toRecv[with];
	min    = Min(toSend[with], toRecv[with]);

	pSwapParticles(from, to, min, with);
	from += min; to += min;

	if (excess > 0) {
	  MCW.Send(from, excess*sizeof(Particle), MPI::BYTE, with, 2*with+1);
	  from += excess;
	}
	if (excess < 0) {
	  MCW.Recv(over, -excess*sizeof(Particle), MPI::BYTE, with, 
		   2*global::rank+1);
	  over -= excess; 
	}
	
	for (; over > &oversize[0] && to < from ; ) *(to++) = *(--over);
      }
    }
    if (over > &oversize[0]) pList.grow(over - &oversize[0]);
    for (;over > &oversize[0];) pList.append(*(--over)); 
    pList.setSize(newSize);
  }

  ftype pMedian(posVal *ind, const int &size, float split, const int &coord, 
		const int &depth) {
    DEBUG("pMedian()", "sifii", "size, split, coord, depth = ", size, split, 
	  coord, depth);
    int   i, j;
    float work, tot;
    ftype myPivot[2], totPivot[2];
  
    if (size < 1) { myPivot[0] = 0.; myPivot[1] = 0.; }
    else { 
      DMORE("sff", "testing ", ind[0].pos[coord], ind[size-1].pos[coord]);
      myPivot[0] = .5*(ind[0].pos[coord] + ind[size-1].pos[coord]); 
      myPivot[1] = 1.;
    }
    MCW.Allreduce(myPivot, totPivot, 2, MPI_ftype, MPI_SUM);
    DMORE("sff", "totPivot[0, 1] = ", totPivot[0], totPivot[1]);

    if (totPivot[1] == 0.) std::cout << "Division by Zero!!! Uaaahhh ...\n";
    totPivot[0] /= totPivot[1];

    for (i = 0, j = size-1, work = 0.; i <= j/* && i < size*/; ) {
      while (ind[i].pos[coord] <  *totPivot && i <  size) work += ind[i++].val;
      while (j >= 0 && ind[j].pos[coord] >= *totPivot)    j--;
      if (i < j) { work += ind[j].val; std::swap(ind[i++], ind[j--]); }
    }

    MCW.Allreduce(&work, &tot, 1, MPI::FLOAT, MPI::SUM);

    if (tot == 0 || depth > 10) return *totPivot;
    if (tot < split) return pMedian(ind+i, size-i, split-tot, coord, depth+1);
    else             return pMedian(ind,   i,      split,     coord, depth+1);
    DMORE("s", "exit pMedian");
  }

  void recDD(posVal *ind, const int &size, const int &coord, const int &proc) {
    if (proc <= 1) return;
    DEBUG("recDD()", "siii", "size coord proc = ", size, coord, proc); 
    
    int   i, left = proc >> 1;
    float split = (float)left/(float)proc, totalWork = 0., work = 0;

    for (i = 0; i < size; i++) work += ind[i].val;
    MCW.Allreduce(&work, &totalWork, 1, MPI::FLOAT, MPI::SUM);
    DMORE("sff", "work, totalWork = ", work, totalWork);

    ftype med = pMedian(ind, size, split*totalWork, coord, 0);
    for (i = size - 1; i >= 0 && ind[i].pos[coord] > med; i--)
      pList[ind[i].nr].addToProc(left);
    
    recDD(ind,     i+1,      (coord+1)%2, left);
    recDD(ind+i+1, size-i-1, (coord+1)%2, proc-left); 
  }

public:
  DomainDecomp(TList<Particle> &_pList) : pList(_pList) {}

  void pDomainDecomp() {
    DEBUG("pDomainDecomp()", "");
    if (global::npro == 1) return;

    int size = pList.getSize(); 
    index.setSize(size);
    DMORE("si", "index size = ", size);

    for (int p = 0; p < size; p++) {
      index[p].pos = pList[p].getPos();
      index[p].val = pList[p].getWork();
      index[p].nr  = p;
      pList[p].setProc(0);
    }

    recDD(&index[0], size, 0, global::npro);

    pMinimizeTraffic();
    pTidyUp();
    //std::cout << "Rank " << global::rank << ": " << pList.getSize() 
	//      << " particles" << std::endl;
  }
};

#endif
