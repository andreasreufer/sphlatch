#ifndef BHTREE_INTEGRATOR_VERLET_H
#define BHTREE_INTEGRATOR_VERLET_H

/*
 *  integrator_verlet.h
 *
 *
 *  Created by Andreas Reufer on 25.05.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "integrator_generic.h"

namespace sphlatch {

/*template<class T>
class VerletIntegrator
{
  public:
    VerletIntegrator(T& _var, T& _dVar, T& _ddVar)
    {
      var   = *_var;
      dVar  = *_dVar;
      ddVar = *_ddVar;
    }
    void resize() {};
    void bootstrap();
    void integrate();
  public:
    T oVar;
    T* var, dVar, ddVar;
};

template<>
void VerletIntegrator<matrixType>::resize()
{
  oVar.resize( oVar.size1(), var->size2() );
}
*/

///
/// template polymorphism is not allow, so we will
/// use a non-templated hierarchy here
///
/// virtual functions are a sin, but in this case
/// they are not called very often
///
class VerletGeneric {
  public:
    VerletGeneric(void) {};
    virtual ~VerletGeneric() {};

  public:
    virtual void bootstrap(void) {};
    virtual void integrate(valueType _dt) {};
};

class VerletVectO2 : public VerletGeneric {
  public:
    VerletVectO2(size_t _i)
    {
      std::cout << "new VerletVectO2 with arg " << _i << "\n";
    };
    ~VerletVectO2() {};

  public:
    void bootstrap(void)
    {
      std::cout << "VerletVectO2 bootstrapped!\n";
    }

    void integrate(valueType _dt)
    {
    }
  private:
    matrixType oVar;
    matrixPtrType var, dVar, ddVar;

};

class VerletScalO2 : public VerletGeneric {
  public:
    VerletScalO2() {};
    ~VerletScalO2() {};

  public:
    void bootstrap(void)
    {
      std::cout << "VerletScalO2 bootstrapped!\n";
    }

    void integrate(valueType _dt)
    {
      std::cout << " integrate ... " << _dt << "\n";
    }
};



class Verlet : public MetaIntegrator<Verlet>
{
  public:
    ///
    /// only second order integration is implemented
    ///
    /// make template out of it
    //template<class T>
    void regIntegration( valvectRefType _var,
                         valvectRefType _devVar,
                         valvectRefType _ddevVar )
    {
      std::vector<VerletGeneric*> integrators;
      integrators.push_back( new VerletVectO2(1) );

      integrators[0]->bootstrap();
      
      //integrators.push_back( new VerletScalO2(2) );

      delete integrators[0];
      //delete integrators[1];
    }
    /*void regIntegration( T& _var,
                         T& _dVar,
                         T& _ddVar )
    {
      VerletIntegrator<T>* Int = new VerletIntegrator<T>(_var, _dVar, _ddVar);
      std::vector< VerletIntegrator<T>* > vect;
    }*/
};

};

#endif
