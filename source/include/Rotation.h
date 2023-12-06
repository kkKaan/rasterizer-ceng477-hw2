#ifndef __ROTATION_H__
#define __ROTATION_H__
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Triangle.h"
#include "Matrix4.h"

class Rotation
{
public:
    int rotationId;
    double angle, ux, uy, uz;

    Rotation();
    Rotation(int rotationId, double angle, double x, double y, double z);
    friend std::ostream &operator<<(std::ostream &os, const Rotation &r);
};

#endif