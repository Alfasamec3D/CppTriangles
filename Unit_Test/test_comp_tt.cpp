#include "geom.h"

bool Test_comp_tt()
{
    Point apex11{5,2,6};
    Point apex12{-2, 2, -4};
    Point apex13{0, -2, 4};
    Point apex21{5,5,8};
    Point apex22{-4,-7,-7};
    Point apex23{3,9,-6};
    Triangle triangle1{apex11, apex12, apex13};
    Triangle triangle2{apex21,apex22,apex23};
    bool res{true};
    return (compare(triangle1,triangle2)==res);
}

int main()
{
    if (Test_comp_tt())
    {
        return 0;
    }

    return 1;
}