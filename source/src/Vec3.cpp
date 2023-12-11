#include <iomanip>
#include "../include/Vec3.h"

Vec3::Vec3()
{
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
    this->colorId = NO_COLOR;
}

Vec3::Vec3(double x, double y, double z)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->colorId = NO_COLOR;
}

Vec3::Vec3(double x, double y, double z, int colorId)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->colorId = colorId;
}

Vec3::Vec3(const Vec3 &other)
{
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->colorId = other.colorId;
}

double Vec3::getNthComponent(int n)
{
    switch (n)
    {
    case 0:
        return this->x;

    case 1:
        return this->y;

    case 2:
    default:
        return this->z;
    }
}

std::ostream &operator<<(std::ostream &os, const Vec3 &v)
{
    os << std::fixed << std::setprecision(6) << "[" << v.x << ", " << v.y << ", " << v.z << "]";
    return os;
}

Vec3 Vec3::operator-(const Vec3 &v)
{
    return Vec3(-v.x, -v.y, -v.z);
}

Vec3 Vec3::translateVec3(Vec3 v, Translation t)
{
    return Vec3(v.x + t.tx, v.y + t.ty, v.z + t.tz);
}

Vec3 Vec3::scaleVec3(Vec3 v, Scaling s)
{
    return Vec3(v.x * s.sx, v.y * s.sy, v.z * s.sz);
}

Vec3 Vec3::rotateVec3(Vec3 v, Rotation r)
{
    // Rodrigues' rotation formula
    double x = v.x;
    double y = v.y;
    double z = v.z;

    double cosTheta = cos(r.angle);
    double sinTheta = sin(r.angle);

    double xPrime = (cosTheta + (1 - cosTheta) * pow(r.ux, 2)) * x + ((1 - cosTheta) * r.ux * r.uy - r.uz * sinTheta) * y + ((1 - cosTheta) * r.ux * r.uz + r.uy * sinTheta) * z;
    double yPrime = ((1 - cosTheta) * r.ux * r.uy + r.uz * sinTheta) * x + (cosTheta + (1 - cosTheta) * pow(r.uy, 2)) * y + ((1 - cosTheta) * r.uy * r.uz - r.ux * sinTheta) * z;
    double zPrime = ((1 - cosTheta) * r.ux * r.uz - r.uy * sinTheta) * x + ((1 - cosTheta) * r.uy * r.uz + r.ux * sinTheta) * y + (cosTheta + (1 - cosTheta) * pow(r.uz, 2)) * z;

    return Vec3(xPrime, yPrime, zPrime);
}