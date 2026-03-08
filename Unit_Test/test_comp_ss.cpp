#include "geom.h"

bool Test_comp_ss()
{
    Point edge11{5,5,8};
    Point edge12{-4,-7,-7};
    Point edge21{14,5,17};
    Point edge22{-16,-5,-18};
    Segment segment1{edge11, edge12};
    Segment segment2{edge21, edge22};
    bool res{true};
    return (compare(segment1,segment2)==res);
}

int main()
{
    if (Test_comp_ss())
    {
        return 0;
    }

    return 1;
}