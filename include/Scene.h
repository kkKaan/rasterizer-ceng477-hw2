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
#include <iostream>

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

	Matrix4 createModelTransformationMatrix(Mesh *mesh);
	Matrix4 createCameraTransformationMatrix(Camera *camera);

	bool liangBarskyClip(Camera *camera, Vec4 &p0, Vec4 &p1);
	void rasterizeLine(Vec4 &v1, Vec4 &v2, Color c1, Color c2, std::vector<std::vector<Color>> &image, int horRes, int verRes);
	void rasterizeTriangle(Vec4 &v0, Vec4 &v1, Vec4 &v2, Color &c0, Color &c1, Color &c2, Camera *camera, std::vector<std::vector<double>> &depthBuffer);

	// void midpoint1(Vec4 &vec1, Vec4 &vec2);
	// void midpoint2(Vec4 &vec1, Vec4 &vec2);
	// void lineRasterizer(Vec4 &vec1, Vec4 &vec2);

	Color interpolateColor(const Color &c1, const Color &c2, double t);

	bool isTriangleBackFacing(Vec4 &v0, Vec4 &v1, Vec4 &v2);
};

#endif
