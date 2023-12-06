#ifndef __TRANSLATION_H__
#define __TRANSLATION_H__
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Camera.h"
#include "Mesh.h"
#include "Triangle.h"
#include "Matrix4.h"

class Translation
{
public:
    int translationId;
    double tx, ty, tz;

    Translation();
    Translation(int translationId, double tx, double ty, double tz);
    friend std::ostream &operator<<(std::ostream &os, const Translation &t);
};

#endif