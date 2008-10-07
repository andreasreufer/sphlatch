#include <iostream>

#ifdef SOLID
  #ifndef SPH
    #define SPH
  #endif
#endif

#include "Def.cc"
#include "Cracks.cc"
#include "Particle.cc"

class Init {
public:
  int planet(Particle *part, const ftype &radius, const int &num) {
    ftype         d = .5527*pow(num, 0.33333), h = radius / d;
    ftype         layer, row, sx = d, sy = d/0.82, sz = d/0.82, x, y, z;
    int           n = 0;
    Vector<ftype> X, V;
    
    for (z = -sz, layer = 0.0; z <= sz; z++, layer++) {
      if (layer > 1.0) layer = -1.0;
      for (y = -sy, row = 0.0; y <= sy; y++, row = 1.0 - row) {
	for (x = -sx; x <= sx; x++) {
	  X.set( h * (x + row/2.0),
		 -h * 1.73205080757 * (y/2.0 + layer/3.0),
		 h * 0.81649658093 * z);
	  part[n].pos  = X;
	  part[n].v    = Vector<ftype>(0., 0., 0.);
	  part[n].h    = h;
	  part[n].mass = 2.7 * h*h*h / 1.4142136;
	  part[n].id   = (ftype)n;
#ifdef SPH
	  part[n].mat  = 1.;
	  part[n].rho  = 2.7;
	  part[n].u    = 1.e5;
#endif	  
	  if (X.len() <= radius) n++; 
	}
      }
    }
    return n;
  }
};  

int main(int argc, char *argv[]) {
  double           time = 0.0;
  Init             init;
  int              n1, n2, sum;
  Material<eosMat> eosTab;
  Particle        *part;
  std::string      crackFile, eosDataFile, path, xdrFile;

  if (argc < 4) { 
    std::cerr << "Usage: " << argv[0] << " filename radius #p\n" << std::endl; 
    quit(1);
  }

  try {
    Manager man(argv[1], true);
    Particle::init(man);

    Xdr output(Particle::numVar, true);

    man.getValue("simulation.path", path);
#ifdef SPH                                                                    
    man.getValue("input.eosDataFile", eosDataFile);
    eosDataFile.insert(0, path + "/");
    eosTab.read(true, eosDataFile);
#endif

    sum  = atol(argv[3]);
    // This is an estimate. It's not possible to build a body with exactly
    // 'sum' particles.
    part = new Particle[sum*11/10];
    sum = init.planet(part, atof(argv[2]), sum);

    std::cout << "Setting up " << sum << " particles" << std::endl;
    
    // Here is the way you write a xdr file. It's a little complicated, but
    // the method can be used in parallel as well.
    man.getValue("input.xdrFile", xdrFile);
    xdrFile.insert(0, path + "/");
    output.open(xdrFile, 'w');
    long pos = output.writeheader(Particle::varName, sum, time, "Test");
    output.createXdr();
    output.setpos(pos);
    for (int i = 0; i < sum; i++) part[i].write(output);
    output.close();

    quit(0);
  }

  catch (std::bad_alloc) {
    std::cerr << "Error: no memory left for new!" << std::endl;
  }
  catch (Eos::ErrorNoNumber) {
    std::cerr << "Eos::Error: No EOS number indicated!" << std::endl;
  }
  catch (Manager::ErrorFileNotFound &e) {
    std::cerr << "Manager::Error: unable to open config file '"
	      << e.file << "'!" << std::endl;
  }
  catch (Manager::ErrorVarNotFound &e) {
    std::cerr << "Manager::Error: variable '" << e.name << "' not "
	      << "found in config file '" << e.file << "'!" << std::endl;
  }
  catch (Material<eosMat>::ErrorFileNotFound &e) {
    std::cerr << "Material::Error: could not find file '" << e.file << "'!" 
	      << std::endl;
  }
  catch (Material<eosMat>::ErrorFileCorrupt &e) {
    std::cerr << "Material::Error: file '" << e.file << "' may be corrupt!" 
	      << std::endl;
  }
  catch (Xdr::ErrorFileNotOpen &e) {
    std::cerr << "Error: unable to open xdr file '" << e.file << "'!"
	      << std::endl;
  }
 quit(1);
}