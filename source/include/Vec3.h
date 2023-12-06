#ifndef __VEC3_H__
#define __VEC3_H__
#define NO_COLOR -1
#include <cmath>
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Scene.h"
#include "Triangle.h"
#include "Matrix4.h"

class Vec3
{
public:
    double x, y, z;
    int colorId;

    Vec3();
    Vec3(double x, double y, double z);
    Vec3(double x, double y, double z, int colorId);
    Vec3(const Vec3 &other);

    double getNthComponent(int n);

    friend std::ostream &operator<<(std::ostream &os, const Vec3 &v);

    Vec3 translateVec3(Vec3 v, Translation t);
    Vec3 scaleVec3(Vec3 v, Scaling s);
    Vec3 rotateVec3(Vec3 v, Rotation r);
};

#endif