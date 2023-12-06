#ifndef __SCALING_H__
#define __SCALING_H__
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Triangle.h"
#include "Matrix4.h"

class Scaling
{
public:
    int scalingId;
    double sx, sy, sz;

    Scaling();
    Scaling(int scalingId, double sx, double sy, double sz);
    friend std::ostream &operator<<(std::ostream &os, const Scaling &s);
};

#endif