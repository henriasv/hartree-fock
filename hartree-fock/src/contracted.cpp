#include "contracted.h"

Contracted::Contracted()
{
}

std::vector<Primitive> Contracted::primitives()
{
    return m_primitives;
}

void Contracted::addPrimitive(Primitive primitive)
{
    m_primitives.push_back(primitive);
}
