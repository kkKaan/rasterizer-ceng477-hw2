#include <iomanip>

#include "../include/Vec4.h"

Vec4::Vec4()
{
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
    this->t = 0.0;
    this->colorId = NO_COLOR;
}

Vec4::Vec4(double x, double y, double z, double t)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->t = t;
    this->colorId = NO_COLOR;
}

Vec4::Vec4(double x, double y, double z, double t, int colorId)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->t = t;
    this->colorId = colorId;
}

Vec4::Vec4(const Vec4 &other)
{
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->t = other.t;
    this->colorId = other.colorId;
}

double Vec4::getNthComponent(int n)
{
    switch (n)
    {
    case 0:
        return this->x;

    case 1:
        return this->y;

    case 2:
        return this->z;

    case 3:
    default:
        return this->t;
    }
}

std::ostream &operator<<(std::ostream &os, const Vec4 &v)
{
    os << std::fixed << std::setprecision(6) << "[" << v.x << ", " << v.y << ", " << v.z << ", " << v.t << "]";
    return os;
}

Vec4 Vec4::translateVec4(Vec4 v, Translation t)
{
    return Vec4(v.x + t.tx, v.y + t.ty, v.z + t.tz, v.t);
}

Vec4 Vec4::scaleVec4(Vec4 v, Scaling s)
{
    return Vec4(v.x * s.sx, v.y * s.sy, v.z * s.sz, v.t);
}

// Vec4 Vec4::multiplyMatrixVec4(Matrix4 m, Vec4 v)
// {
//     Vec4 result;

//     result.x = v.x * m.values[0][0] + v.y * m.values[1][0] + v.z * m.values[2][0] + v.t * m.values[3][0];
//     result.y = v.x * m.values[0][1] + v.y * m.values[1][1] + v.z * m.values[2][1] + v.t * m.values[3][1];
//     result.z = v.x * m.values[0][2] + v.y * m.values[1][2] + v.z * m.values[2][2] + v.t * m.values[3][2];
//     result.t = v.x * m.values[0][3] + v.y * m.values[1][3] + v.z * m.values[2][3] + v.t * m.values[3][3];

//     return result;
// }
