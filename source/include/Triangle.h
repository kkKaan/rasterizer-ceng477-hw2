#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"

class Triangle
{
public:
    int vertexIds[3];

    Triangle();
    Triangle(int vid1, int vid2, int vid3);
    Triangle(const Triangle &other);
    friend std::ostream &operator<<(std::ostream &os, const Triangle &t);

    Triangle transformTriangleToCameraSpace(Mesh mesh, Camera camera);
};

#endif