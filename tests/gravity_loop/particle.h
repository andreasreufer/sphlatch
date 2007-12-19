#ifndef OOSPHPARTICLE_H
#define OOSPHPARTICLE_H
#include <map>

namespace oosph {

enum Particle {ID, X, Y, Z, VX, VY, VZ, AX, AY, AZ, M, H, DHDT, RHO, E, P, POW, DIV_V, ROTX_V, ROTY_V, ROTZ_V, Q, GRAVEPS, MAXAVMU, ALPHA, DALPHADT, NONEIGH, OX, OY, OZ, OVX, OVY, OVZ, OAX, OAY, OAZ, PX, PY, PZ, PVX, PVY, PVZ, PAX, PAY, PAZ, OE, PE, OPOW, PPOW, OH, PH, ODHDT, PDHDT, OALPHA, PALPHA, ODALPHADT, PDALPHADT, SIZE};

template <typename SimTrait>
class ParticleVarMap {

	public:
		typedef ParticleVarMap self_type;
		typedef ParticleVarMap& self_reference;
		typedef ParticleVarMap* self_pointer;
		
		typedef std::multimap<size_t,std::string> MapType;
		static self_reference Instance();

		std::string getName(int);
		int getIndex(std::string);

	protected:
		ParticleVarMap(void);
		~ParticleVarMap(void);
	
	private:
		MapType mapping;
		static self_pointer _instance;
};

template <typename SimTrait>
typename ParticleVarMap<SimTrait>::self_pointer ParticleVarMap<SimTrait>::_instance = NULL;

template <typename SimTrait>
typename ParticleVarMap<SimTrait>::self_reference ParticleVarMap<SimTrait>::Instance(void)
{
	if (_instance != NULL)
		return *_instance;
	else
	{
		_instance = new ParticleVarMap;
		return *_instance;
	}
}

template <typename SimTrait>
ParticleVarMap<SimTrait>::ParticleVarMap(void)
{
	// Define string <-> variable mapping HERE!
	mapping.insert(MapType::value_type(ID,"id"));
	mapping.insert(MapType::value_type(X,"x"));
	mapping.insert(MapType::value_type(Y,"y"));
	mapping.insert(MapType::value_type(Z,"z"));
	mapping.insert(MapType::value_type(VX,"vx"));
	mapping.insert(MapType::value_type(VY,"vy"));
	mapping.insert(MapType::value_type(VZ,"vz"));
	mapping.insert(MapType::value_type(AX,"ax"));
	mapping.insert(MapType::value_type(AY,"ay"));
	mapping.insert(MapType::value_type(AZ,"az"));
	mapping.insert(MapType::value_type(M,"m"));
	mapping.insert(MapType::value_type(M,"mass"));
	mapping.insert(MapType::value_type(M,"pmass"));
	mapping.insert(MapType::value_type(H,"h"));
	mapping.insert(MapType::value_type(OH,"oh"));
	mapping.insert(MapType::value_type(PH,"ph"));
	mapping.insert(MapType::value_type(DHDT,"dhdt"));
	mapping.insert(MapType::value_type(RHO,"rho"));
	mapping.insert(MapType::value_type(E,"u"));
	mapping.insert(MapType::value_type(OE,"ou"));
	mapping.insert(MapType::value_type(PE,"pu"));
	mapping.insert(MapType::value_type(P,"p"));
	mapping.insert(MapType::value_type(POW,"pow"));
	mapping.insert(MapType::value_type(DIV_V,"divv"));
	mapping.insert(MapType::value_type(ROTX_V,"rotxv"));
	mapping.insert(MapType::value_type(ROTY_V,"rotyv"));
	mapping.insert(MapType::value_type(ROTZ_V,"rotzv"));
	mapping.insert(MapType::value_type(Q,"q"));
	mapping.insert(MapType::value_type(GRAVEPS,"graveps"));
	mapping.insert(MapType::value_type(MAXAVMU,"maxavmu"));
	mapping.insert(MapType::value_type(ALPHA,"alpha"));
	mapping.insert(MapType::value_type(DALPHADT,"dalphadt"));
	mapping.insert(MapType::value_type(NONEIGH,"noneigh"));
	mapping.insert(MapType::value_type(OX,"ox"));
	mapping.insert(MapType::value_type(OY,"oy"));
	mapping.insert(MapType::value_type(OZ,"oz"));
	mapping.insert(MapType::value_type(PX,"px"));
	mapping.insert(MapType::value_type(PY,"py"));
	mapping.insert(MapType::value_type(PZ,"pz"));
	mapping.insert(MapType::value_type(OVX,"ovx"));
	mapping.insert(MapType::value_type(OVY,"ovy"));
	mapping.insert(MapType::value_type(OVZ,"ovz"));
	mapping.insert(MapType::value_type(PVX,"pvx"));
	mapping.insert(MapType::value_type(PVY,"pvy"));
	mapping.insert(MapType::value_type(PVZ,"pvz"));
	mapping.insert(MapType::value_type(OAX,"oax"));
	mapping.insert(MapType::value_type(OAY,"oay"));
	mapping.insert(MapType::value_type(OAZ,"oaz"));
	mapping.insert(MapType::value_type(PAX,"pax"));
	mapping.insert(MapType::value_type(PAY,"pay"));
	mapping.insert(MapType::value_type(PAZ,"paz"));
	mapping.insert(MapType::value_type(OE,"ou"));
	mapping.insert(MapType::value_type(PE,"pu"));
	mapping.insert(MapType::value_type(OPOW,"opow"));
	mapping.insert(MapType::value_type(PPOW,"ppow"));
	mapping.insert(MapType::value_type(OH,"oh"));
	mapping.insert(MapType::value_type(PH,"ph"));
	mapping.insert(MapType::value_type(ODHDT,"odhdt"));
	mapping.insert(MapType::value_type(PDHDT,"pdhdt"));
	mapping.insert(MapType::value_type(OALPHA,"oalpha"));
	mapping.insert(MapType::value_type(PALPHA,"palpha"));
	mapping.insert(MapType::value_type(ODALPHADT,"odalphadt"));
	mapping.insert(MapType::value_type(PDALPHADT,"pdalphadt"));
}

template <typename SimTrait>
ParticleVarMap<SimTrait>::~ParticleVarMap(void)
{}

template <typename SimTrait>
std::string ParticleVarMap<SimTrait>::getName(int _index) {
	return (mapping.lower_bound(_index))->second;
};

template <typename SimTrait>
int ParticleVarMap<SimTrait>::getIndex(std::string _name) {
	int _returnval = -1;
	for (MapType::iterator pos = mapping.begin(); pos != mapping.end(); ++pos) {
		if (pos->second == _name) { _returnval = pos->first;}
	}
	return _returnval;
};

};

#endif
