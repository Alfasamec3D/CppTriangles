#include "geom.h"

bool Test_vectprod()
{
    Vector vector1{1,2,3};
    Vector vector2{-5,3,-2};
    Vector res{-13,-13,13};
    return ((vectprod(vector1,vector2)==res));
}

int main()
{
    if (Test_vectprod())
    {
        return 0;
    }

    return 1;
}