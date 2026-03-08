#include "geom.h"

bool Test_scalprod()
{
    Point point1{1,2,3};
    Point point2{-5,3,-2};
    double res {-5};
    return (eps_eq(scalprod(point1, point2),res));
}

int main()
{
    if (Test_scalprod())
    {
        return 0;
    }

    return 1;
}