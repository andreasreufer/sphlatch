/******************************************************************************
 * ParaSPH -- Version 19.01.2004                                              *
 *----------------------------------------------------------------------------*
 * File:      Manager.cc                                                      *
 * Purpose:   This class allows you to read and write configuration files.    *
 *            The advantage in using this class is that you do not have to    *
 *            'hardwire' variables into the code. For example: you want to    *
 *            access a new variable called 'output.title', just write the     *
 *            variable and its value into the '.config' file                  *
 *            (eg. output.title = "this_is_a_title") and call                 *
 *            man.getValue("output.title", title) within the code and the     *
 *            variable 'title' will contain the value from the '.config' file.*
 *            Confused? Never mind, just work with the program and you'll see.*
 *                                                                            *
 * Input:     The actual version of the code expects at least the following   *
 *            variables to be set: simulation.[eos, gravH, heatCon, theta,    *
 *            XSPH], input.[crackFile, eosDataFile, matDataFile, xdrFile],    *
 *            output.saveTime (several possible), physics.G, SPH.[avAlpha,    *
 *            avBeta, courant, damping, dtFactor, dtMax], thres.[dmMin, hMin, *
 *            rhoMin, sMin, uMin, tiny].                                      *
 *            Not everything is needed in every possible combination of       *
 *            compiler switches. If some variable is missing the program will *
 *            stop and inform you about it.                                   *
 *                                                                            *
 * Modifications:                                                             *
 *            - 12.12.2003 WB                                                 *
 *              changed function getMinValueGreater to solve problem with     *
 *              time input format                                             *
 *****************************************************************************/
#ifndef MANAGER_CC
#define MANAGER_CC

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Def.cc"

// This class is needed by class 'Manager' to store all the variables read
// from the .config file
class Variable {
private:
  std::string name, value;
  static u_int max;

public:
  Variable() {}

  Variable(const std::string &_name, const std::string &_value) {
    name = _name; value = _value;
    if (name.size() > max) max = name.size();
  }

  bool read(std::string &line) {
    if (trim(line).size() < 1) return false;
    if (line.at(0) == '#')     return false;

    u_long pos = line.find('=');
    name = trim(line.substr(0, pos));
    if (name.size() > max) max = name.size();
    if (pos != std::string::npos) value = trim(line.substr(pos+1));
    return true;
  }

  bool operator>(const Variable &o) const { 
   return (name > o.name || (name == o.name && value > o.value)); 
  }
  bool operator<(const Variable &o) const {
    return (name < o.name || (name == o.name && value < o.value)); 
  } 

  std::string getName()  const { return name; }

  void getValue(bool   &res) {
    res = (toUpper(value).find("TRUE") != std::string::npos);
  }
  void getValue(int    &res) { res = atoi(value.c_str()); }
  void getValue(float  &res) { res = atof(value.c_str()); }
  void getValue(double &res) { res = atof(value.c_str()); }
  void getValue(std::string &res) { res = value; }

  void setName (const std::string &_name)  { name  = _name; }
  void setValue(const std::string &_value) { value = _value;}
  
  void print(std::ostream &out) { 
    out << std::setw(max) << setiosflags(std::ios::left) << name << " = "
	<< value << std::endl; 
  }
};

u_int Variable::max;



class Manager {
private:
  std::string           configPath, path, simName;
  std::vector<Variable> var;

public:
  class ErrorFileNotFound {
  public:
    std::string file;
    ErrorFileNotFound(const std::string &_file) { file = _file; }
  };
  class ErrorVarNotFound {
  public:
    std::string file, name;
    ErrorVarNotFound(const std::string &_name, const std::string &_file) { 
      name = _name; file = _file; 
    }
  };

  // Important! The Manager needs to be initialized before use. 'pathName'
  // should be a valid ParaSPH path. See ParaSPH.cc for details!
  // If 'first' is true the program tries to load configuration data from
  // file '.config' (otherwise 'config'). 
  //void init(const std::string &pathName, const bool &first) {
  Manager(const std::string &pathName, const bool &first) {
    int len;

    for (len = pathName.size(); pathName[len-1] == '/'; len--);
    path       = pathName.substr(0, len);
    simName    = path.substr(path.rfind("/")+1);
    configPath = path + (first ? "/.config" : "/config");

    readConfig();

    configPath = path + "/config";

    var.push_back(Variable("simulation.path", path));
    var.push_back(Variable("simulation.name", simName));
  }

  void readConfig() {
    std::ifstream file(configPath.c_str());
    std::string   line;
    Variable      v;

    if (!file.good()) throw ErrorFileNotFound(configPath);
    do {
      getline(file, line);
      if (!file.eof() && v.read(line)) var.push_back(v);
    } while (!file.eof());
    file.close();

    if (global::noisy) 
      std::cout << "Reading config file '" << configPath << "', " 
		<< var.size() << " variables read" << std::endl;
  }

  void saveConfig() {
    std::ofstream file(configPath.c_str());
    
    if (!file.good()) throw ErrorFileNotFound(configPath);
    if (global::noisy) std::cout << "Writing config file '" << configPath 
				 << "'" << std::endl;
    std::sort(var.begin(), var.end());
    file << "# Autom. created config file, manual changes may apply\n\n";
    for (u_int i = 0; i < var.size(); i++) var[i].print(file);
    file.close();
  }

  Variable *find(const std::string &name) {
    for (u_int i = 0; i < var.size(); i++)
      if (var[i].getName().find(name) != std::string::npos) return &var[i]; 
    throw ErrorVarNotFound(name, configPath);
  }

  // Get a variables value. If the variable does not exist, an error is 
  // thrown and the program stops.
  template <class T>
  void getValue(const std::string &name, T &value) { 
    find(name)->getValue(value);
  }

  // Every variable can hold several values. Here you can find the value that
  // is greater than a certain minimum 'min'. Is used for outputs at varying 
  // times. The function returns true, if such a value has been found, or
  // false otherwise.
  template <class T>
  bool getMinValueGreater(const std::string &name, const T &min, T &value) {
    T help;
    value = 1.e30;
    for (uint i = 0; i < var.size(); i++)
      if (var[i].getName().find(name) != std::string::npos) {
      var[i].getValue(help);
      if (help > min && help < value) value = help;
    }
    if (value == 1.e30) return false; else return true;
  }

  // Add a value to an existing variable
  void addValue(const std::string &name, const std::string &value) {
    var.push_back(Variable(name, value));
  }

  // Set the value of a variable
  void setValue(const std::string &name, const std::string &value) {
    try { Variable *v = find(name); v->setValue(value); }
    catch (Manager::ErrorVarNotFound) { addValue(name, value); }
  }
};

#endif
