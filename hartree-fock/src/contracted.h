#ifndef CONTRACTED_H
#define CONTRACTED_H
#include <vector>
#include <primitive.h>

class Contracted
{
public:
    Contracted();
    std::vector<Primitive> primitives();
    void addPrimitive(Primitive primitive);

private:
    std::vector<Primitive> m_primitives;
};

#endif // CONTRACTED_H
