#include <iomanip>
#include "../include/Matrix4.h"

Matrix4::Matrix4()
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            this->values[i][j] = 0;
        }
    }
}

Matrix4::Matrix4(double values[4][4])
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            this->values[i][j] = values[i][j];
        }
    }
}

Matrix4::Matrix4(const Matrix4 &other)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            this->values[i][j] = other.values[i][j];
        }
    }
}

std::ostream &operator<<(std::ostream &os, const Matrix4 &m)
{

    os << std::fixed << std::setprecision(6) << "|" << m.values[0][0] << "|" << m.values[0][1] << "|" << m.values[0][2] << "|" << m.values[0][3] << "|"
       << std::endl
       << "|" << m.values[1][0] << "|" << m.values[1][1] << "|" << m.values[1][2] << "|" << m.values[1][3] << "|"
       << std::endl
       << "|" << m.values[2][0] << "|" << m.values[2][1] << "|" << m.values[2][2] << "|" << m.values[2][3] << "|"
       << std::endl
       << "|" << m.values[3][0] << "|" << m.values[3][1] << "|" << m.values[3][2] << "|" << m.values[3][3] << "|";

    return os;
}

Matrix4 Matrix4::operator*(const Matrix4 &other)
{
    Matrix4 result;

    for (int i = 0; i < 4; i++) // row
    {
        for (int j = 0; j < 4; j++) // column
        {
            for (int k = 0; k < 4; k++) // index
            {
                result.values[i][j] += this->values[i][k] * other.values[k][j];
            }
        }
    }
    return result;
}

Matrix4 Matrix4::createOrthographicProjectionMatrix(Camera *camera)
{
    return Matrix4({
        {2 / (camera->right - camera->left), 0, 0, -(camera->right + camera->left) / (camera->right - camera->left)},
        {0, 2 / (camera->top - camera->bottom), 0, -(camera->top + camera->bottom) / (camera->top - camera->bottom)},
        {0, 0, -2 / (camera->far - camera->near), -(camera->far + camera->near) / (camera->far - camera->near)},
        {0, 0, 0, 1}
    });
}

Matrix4 createPerspectiveProjectionMatrix(Camera *camera)
{
    double aspectRatio = static_cast<double>(camera->horRes) / camera->verRes;
    double fovyRadians = 2 * atan((camera->top - camera->bottom) / (2 * camera->near)); // Field of view in radians
    double f = 1.0 / tan(fovyRadians / 2); // Focal length

    return Matrix4({
        {f / aspectRatio, 0, 0, 0},
        {0, f, 0, 0},
        {0, 0, (camera->far + camera->near) / (camera->near - camera->far), 2 * camera->far * camera->near / (camera->near - camera->far)},
        {0, 0, -1, 0}
    });
}
