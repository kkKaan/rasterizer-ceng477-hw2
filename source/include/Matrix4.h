#ifndef __MATRIX4_H__
#define __MATRIX4_H__
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

class Matrix4
{
public:
    double values[4][4];

    Matrix4();
    Matrix4(double values[4][4]);
    Matrix4(const Matrix4 &other);
    friend std::ostream &operator<<(std::ostream &os, const Matrix4 &m);

    Matrix4 operator*(const Matrix4 &other);
    // static Matrix4 createOrthographicProjectionMatrix(Camera *camera);
};

#endif