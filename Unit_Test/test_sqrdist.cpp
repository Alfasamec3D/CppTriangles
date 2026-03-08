#include "geom.h"

bool Test_sqrdist()
{
    Point point1{1,2,3};
    Point point2{-5,3,-2};
    double res{62};
    return (eps_eq(sqrdist(point1, point2),res));
}

int main()
{
    if (Test_sqrdist())
    {
        return 0;
    }

    return 1;
}