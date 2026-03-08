#include "geom.h"

bool Test_tr_segm()
{
    Point point1{1,-9,17};
    Point point2{-3,-1,7};
    Point point3{-5,3,2};
    Triangle triangle{point1, point2, point3};
    Segment res{point1,point3};
    return (triangle.segm()==res);
}

int main()
{
    if (Test_tr_segm())
    {
        return 0;
    }

    return 1;
}