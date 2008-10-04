#include <iostream>
#include <mpi.h>

#define MPI_OK

#ifdef SOLID
  #ifndef SPH
    #define SPH
  #endif
#endif

#include "Def.cc"
#include "PProcessor.cc"

int main(int argc, char *argv[]) {
  MPI::Init(argc, argv);
  global::npro  = MCW.Get_size();
  global::rank  = MCW.Get_rank();
  global::noisy = (global::rank == 0);
  
  listSwitches(global::noisy);

  if (argc < 2) { 
    if (global::noisy) std::cerr << "Usage: " << argv[0] << " dirname\n" 
				 << std::endl;
    quit(global::noisy);
  }

  try {
    checkbuild(global::rank, global::noisy);    
    Processor pro(argv[1]);
    pro.integrator();
    quit(0);
  }

  // If adding (definitive) error messages to the code, they should be added
  // here ... makes the code a little more consistent :-)
  catch (std::bad_alloc) {
    if (global::noisy) 
      std::cerr << "std::Error: no memory left for new!" << std::endl;
  }
  catch (ErrorVaryingBuildNumbers) {
    if (global::noisy) 
      std::cerr << "def::Error: build numbers do not match!" << std::endl;
  }
/*
  catch (Cracks::ErrorFileNotOpen &e) {
    if (global::noisy) std::cerr << "Cracks::Error: unable to open crack file '" 
			 << e.file << "'!" << std::endl;
  }
*/
  catch (Eos::ErrorNoNumber) {
    if (global::noisy) 
      std::cerr << "Eos::Error: No EOS number indicated!" << std::endl;
  }
  catch (Manager::ErrorFileNotFound &e) {
    if (global::noisy) 
      std::cerr << "Manager::Error: unable to open config file '"
		<< e.file << "'!" << std::endl;
  }
  catch (Manager::ErrorVarNotFound &e) {
    if (global::noisy) 
      std::cerr << "Manager::Error: variable '" << e.name << "' not "
		<< "found in config file '" << e.file << "'!" << std::endl;
  }
  catch (Material<eosMat>::ErrorFileNotFound &e) {
    if (global::noisy) 
      std::cerr << "Material::Error: could not find eos file '"
		<< e.file << "'!" << std::endl;
  }
  catch (Material<eosMat>::ErrorFileCorrupt &e) {
    if (global::noisy) 
      std::cerr << "Material::Error: eos file '" << e.file
		<< "' may be corrupt!" << std::endl;
  }
  catch (Material<Mat>::ErrorFileNotFound &e) {
    if (global::noisy) 
      std::cerr << "Material::Error: could not find mat file '"
		<< e.file << "'!" << std::endl;
  }
  catch (Material<Mat>::ErrorFileCorrupt &e) {
    if (global::noisy) 
      std::cerr << "Material::Error: mat file '" << e.file
		<< "' may be corrupt!" << std::endl;
  }
#ifdef SPH
  catch (Particle::ErrorHTooBig &e) {
    if (global::noisy) 
      std::cerr << "Particle::Error: h > hMax (" << e.h << " > "
		<< /*Particle::hMax <<*/ ") -> neighbourlist corrupted!"
		<< std::endl;
  }
#endif
  catch (Xdr::ErrorFileNotOpen &e) {
    if (global::noisy) 
      std::cerr << "Xdr::Error: unable to open xdr file '" << e.file
		<< "'!" << std::endl;
  }
  catch (Xdr::ErrorWrongVariableNumber &e) {
    if (global::noisy) 
      std::cerr << "Xdr::Error: wrong number of variables (" << e.num
		<< " in xdr file '" << e.file << "' instead of " 
		<< Particle::numVar << "! Recompile the "
		<< "code with appropriate switches." << std::endl;
  }
  quit(global::noisy);
}
