#include "berylliumhf.h"

BerylliumHF::BerylliumHF() :
    AtomicOrbitals::AtomicOrbitals(4, 4, 4)
{
}

std::string BerylliumHF::str()
{
    return std::string("Beryllium");
}


