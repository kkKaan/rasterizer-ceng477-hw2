#ifndef __VEC4_H__
#define __VEC4_H__

#define NO_COLOR -1

#include "Scaling.h"
#include "Translation.h"
#include "Matrix4.h"

class Vec4
{
public:
    double x, y, z, t;
    int colorId;

    Vec4();
    Vec4(double x, double y, double z, double t);
    Vec4(double x, double y, double z, double t, int colorId);
    Vec4(const Vec4 &other);

    double getNthComponent(int n);

    friend std::ostream &operator<<(std::ostream &os, const Vec4 &v);

    Vec4 translateVec4(Vec4 v, Translation t);
    Vec4 scaleVec4(Vec4 v, Scaling s);
    // Vec4 multiplyMatrixVec4(Matrix4 m, Vec4 v);
};

#endif