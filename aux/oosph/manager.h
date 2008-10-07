#ifndef OOSPH_MANAGER_H
#define OOSPH_MANAGER_H

#include "memorymanager.h"
#include "iomanager.h"
//#include "protocolmanager.h"
#include "communicationmanager.h"

namespace oosph
{

template <typename SimTrait>
class Manager
{
public:

    typedef Manager self_type;
    typedef Manager& self_reference;
    typedef Manager* self_pointer;

    typedef IOManager<SimTrait> io_manager_type;
    typedef IOManager<SimTrait>& io_manager_reference;
    typedef IOManager<SimTrait>* io_manager_pointer;

    typedef MemoryManager<SimTrait> mem_manager_type;
    typedef MemoryManager<SimTrait>& mem_manager_reference;
    typedef MemoryManager<SimTrait>* mem_manager_pointer;

    /*typedef ProtocolManager<SimTrait> protocol_manager_type;
    typedef ProtocolManager<SimTrait>& protocol_manager_reference;
    typedef ProtocolManager<SimTrait>* protocol_manager_pointer;*/

    typedef CommunicationManager<SimTrait> communication_manager_type;
    typedef CommunicationManager<SimTrait>& communication_manager_reference;
    typedef CommunicationManager<SimTrait>* communication_manager_pointer;

    static self_reference Instance(void);

    static io_manager_reference IO;
    static mem_manager_reference MEM;
    //static protocol_manager_reference PROTOCOL;
    static communication_manager_reference COM;

protected:

    Manager(void);
    ~Manager(void);

private:

    static self_pointer _instance;

};

template <typename SimTrait>
typename Manager<SimTrait>::self_pointer Manager<SimTrait>::_instance = NULL;

template <typename SimTrait>
typename Manager<SimTrait>::self_reference Manager<SimTrait>::Instance(void)
{
    if (_instance != NULL)
        return *_instance;
    else
    {
        _instance = new Manager;
        return *_instance;
    }
}

template <typename SimTrait>
typename Manager<SimTrait>::io_manager_reference Manager<SimTrait>::IO(io_manager_type::Instance() );

template <typename SimTrait>
typename Manager<SimTrait>::mem_manager_reference
Manager<SimTrait>::MEM(mem_manager_type::Instance() );

/*template <typename SimTrait>
typename Manager<SimTrait>::protocol_manager_reference
Manager<SimTrait>::PROTOCOL(protocol_manager_type::Instance() );*/

template <typename SimTrait>
typename Manager<SimTrait>::communication_manager_reference
Manager<SimTrait>::COM(communication_manager_type::Instance() );

template <typename SimTrait>
Manager<SimTrait>::Manager(void)
{
    MEM;
    COM;
    //PROTOCOL;
    IO;
}

template <typename SimTrait>
Manager<SimTrait>::~Manager(void)
{

}



};

#endif
