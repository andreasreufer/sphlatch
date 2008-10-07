#include <iostream>

#include "Def.cc"
#include "Manager.cc"
#include "Particle.cc"

class Init {
 public:
  static void init(Particle *part) {
    ftype D = 100., d = 5.;
    
    ftype v = sqrt(2./d), V = sqrt(4./D);

    part[0].pos.set(-D-d, 0., 0.);
    part[0].v.set(0., +V+v, 0.);
    
    part[1].pos.set(-D+d, 0., 0.);
    part[1].v.set(0., +V-v, 0.);
    
    part[2].pos.set(+D-d, 0., 0.);
    part[2].v.set(0., -V+v, 0.);
    
    part[3].pos.set(+D+d, 0., 0.);
    part[3].v.set(0., -V-v, 0.);
    
    for (int i = 0; i < 4; i++) {
      part[i].h = 1.;
      part[i].mass = 1.;
      part[i].id = i;
    }
  }

  static int planet(Particle *part, const ftype &radius, const int &num) {
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

          if (X.len() <= radius) n++;
        }
      }
    }
    return n;
  }

};

int main(int argc, char *argv[]) {
  double           time = 0.0;
  int              sum;
  Particle        *part;
  std::string      path, xdrFile;

  if (argc < 2) { 
    std::cerr << "Usage: " << argv[0] << " filename\n" << std::endl; 
    quit(1);
  }

  try {
    Manager man(argv[1], true);
    Particle::init(man);

    Xdr output(Particle::numVar, true);

    man.getValue("simulation.path", path);

    sum  = 4;
    part = new Particle[sum];
    Init::init(part);

    std::cout << "Setting up " << sum << " particles" << std::endl;
    
    man.getValue("input.xdrFile", xdrFile);
    xdrFile.insert(0, path + "/");
    output.open(xdrFile, 'w');
    long pos = output.writeheader(Particle::varName, sum, time, "grav-test");
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
