#ifndef SPHLATCH_PARTICLE_H
#define SPHLATCH_PARTICLE_H
#include <map>

namespace sphlatch {
enum ParticleIndex { ID, X, Y, Z, VX, VY, VZ, AX, AY, AZ, M, H,
                     DHDT, RHO, E, P, POW, DIV_V, ROTX_V, ROTY_V, ROTZ_V,
                     Q, GRAVEPS, MAXAVMU, ALPHA, DALPHADT, NONEIGH, SIZE };

class ParticleVarMap {
public:
typedef ParticleVarMap self_type;
typedef ParticleVarMap& self_reference;
typedef ParticleVarMap* self_pointer;

typedef std::multimap<size_t, std::string> mapType;
static self_reference instance();

std::string getName(int _index);
int getIndex(std::string _name);

protected:
ParticleVarMap(void);
~ParticleVarMap(void);

private:
mapType mapping;
static self_pointer _instance;
};

ParticleVarMap::self_pointer ParticleVarMap::_instance = NULL;

ParticleVarMap::self_reference ParticleVarMap::instance(void)
{
  if (_instance != NULL)
    return *_instance;
  else
    {
      _instance = new ParticleVarMap;
      return *_instance;
    }
}

ParticleVarMap::ParticleVarMap(void)
{
  // Define string <-> variable mapping HERE!
  mapping.insert(mapType::value_type(ID, "id"));
  mapping.insert(mapType::value_type(X, "x"));
  mapping.insert(mapType::value_type(Y, "y"));
  mapping.insert(mapType::value_type(Z, "z"));
  mapping.insert(mapType::value_type(VX, "vx"));
  mapping.insert(mapType::value_type(VY, "vy"));
  mapping.insert(mapType::value_type(VZ, "vz"));
  mapping.insert(mapType::value_type(AX, "ax"));
  mapping.insert(mapType::value_type(AY, "ay"));
  mapping.insert(mapType::value_type(AZ, "az"));
  mapping.insert(mapType::value_type(M, "m"));
  mapping.insert(mapType::value_type(M, "mass"));
  mapping.insert(mapType::value_type(M, "pmass"));
  mapping.insert(mapType::value_type(H, "h"));
  mapping.insert(mapType::value_type(DHDT, "dhdt"));
  mapping.insert(mapType::value_type(RHO, "rho"));
  mapping.insert(mapType::value_type(E, "u"));
  mapping.insert(mapType::value_type(P, "p"));
  mapping.insert(mapType::value_type(POW, "pow"));
  mapping.insert(mapType::value_type(DIV_V, "divv"));
  mapping.insert(mapType::value_type(ROTX_V, "rotxv"));
  mapping.insert(mapType::value_type(ROTY_V, "rotyv"));
  mapping.insert(mapType::value_type(ROTZ_V, "rotzv"));
  mapping.insert(mapType::value_type(Q, "q"));
  mapping.insert(mapType::value_type(GRAVEPS, "graveps"));
  mapping.insert(mapType::value_type(MAXAVMU, "maxavmu"));
  mapping.insert(mapType::value_type(ALPHA, "alpha"));
  mapping.insert(mapType::value_type(DALPHADT, "dalphadt"));
  mapping.insert(mapType::value_type(NONEIGH, "noneigh"));
}

ParticleVarMap::~ParticleVarMap(void)
{
}

std::string ParticleVarMap::getName(int _index)
{
  return (mapping.lower_bound(_index))->second;
};

int ParticleVarMap::getIndex(std::string _name)
{
  int _returnval = -1;

  for (mapType::iterator pos = mapping.begin(); pos != mapping.end(); ++pos)
    {
      if (pos->second == _name)
        {
          _returnval = pos->first;
        }
    }
  return _returnval;
};
};

#endif
