#ifndef _SCENE_H_
#define _SCENE_H_

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "Camera.h"
#include "Color.h"
#include "Helpers.h"
#include "Matrix4.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "tinyxml2.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color> > image;
	std::vector<std::vector<double> > depth;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	Scene(const char *xmlPath);

	void assignColorToPixel(int i, int j, Color c);
	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName);
	void forwardRenderingPipeline(Camera *camera);

	void applyModelTransformation(Mesh *mesh, Triangle& triangle);
	void applyCameraTransformation(Camera *camera, Triangle& triangle);
	void applyViewportTransformation(Camera *camera, Triangle& triangle);

	// uint8_t computeOutcode(Camera *camera, double x, double y);
	// void clipLine(Camera *camera, Vec3 &p0, Vec3 &p1);
	// void clipTriangle(Camera *camera, Triangle& triangle);

	bool liangBarskyClip(Camera *camera, Vec3 &p0, Vec3 &p1);
	void drawLine(Camera *camera, Vec3 *v1, Vec3 *v2, Color *c1, Color *c2);

	bool clipping(Camera& camera, Vec4 &vec0, Vec4 &vec1);
	void midpoint1(Vec4 &vec1, Vec4 &vec2);
	void midpoint2(Vec4 &vec1, Vec4 &vec2);
	void lineRasterizer(Vec4 &vec1, Vec4 &vec2);

	Color interpolateColor(const Color &c1, const Color &c2, double t);

	bool isTriangleBackFacing(const Triangle& triangle, Camera *camera);
};

#endif
