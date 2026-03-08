#include "geom.h"

bool Test_comp_ps()
{
    Point point {2,1,3};
    Point edge1{5,5,8};
    Point edge2{-4,-7,-7};
    Segment segment{edge1, edge2};
    bool res{true};
    return (compare(point,segment)==res);
}

int main()
{
    if (Test_comp_ps())
    {
        return 0;
    }

    return 1;
}