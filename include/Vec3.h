#ifndef __VEC3_H__
#define __VEC3_H__

#define NO_COLOR -1
#define M_PI 3.14159265358979323846

#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"

#include <cmath>

// #include "Helpers.h"

class Vec3
{
public:
    double x, y, z;
    int colorId;

    Vec3();
    Vec3(double x, double y, double z);
    Vec3(double x, double y, double z, int colorId);
    Vec3(const Vec3& other);

    double getNthComponent(int n);

    friend std::ostream& operator<<(std::ostream& os, const Vec3& v);
    Vec3 operator-(const Vec3& v);
};

#endif