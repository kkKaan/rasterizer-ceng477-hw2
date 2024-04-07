#include "../include/Color.h"

#include <iomanip>

Color::Color()
{
    this->r = 0;
    this->g = 0;
    this->b = 0;
}

Color::Color(double r, double g, double b)
{
    this->r = r;
    this->g = g;
    this->b = b;
}

Color::Color(const Color& other)
{
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
}

std::ostream& operator<<(std::ostream& os, const Color& c)
{
    os << std::fixed << std::setprecision(0) << "rgb(" << c.r << ", " << c.g << ", " << c.b << ")";
    return os;
}

Color Color::operator+(const Color& c)
{
    return Color(std::min(this->r + c.r, 255.0), std::min(this->g + c.g, 255.0), std::min(this->b + c.b, 255.0));
}

Color Color::operator*(const double& c)
{
    return Color(std::min(this->r * c, 255.0), std::min(this->g * c, 255.0), std::min(this->b * c, 255.0));
}

Color Color::operator-(const Color& c)
{
    return Color(std::max(this->r - c.r, 0.0), std::max(this->g - c.g, 0.0), std::max(this->b - c.b, 0.0));
}

Color Color::operator/(const double& d)
{
    return Color(this->r / d, this->g / d, this->b / d);
}