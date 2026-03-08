#include "geom.h"

bool Test_comp_pt()
{
    Point point {2,1,3};
    Point apex1 {5,2,6};
    Point apex2{-2, 2, -4};
    Point apex3{0, -2, 4};
    Triangle triangle{apex1, apex2, apex3};
    bool res{true};
    return (compare(point,triangle)==res);
}

int main()
{
    if (Test_comp_pt())
    {
        return 0;
    }

    return 1;
}