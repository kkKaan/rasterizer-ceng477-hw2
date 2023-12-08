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

#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Helpers.h"
#include "Triangle.h"
#include "Matrix4.h"
#include "tinyxml2.h"

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

	void clipTriangle(Camera *camera, Triangle& triangle);
	bool isTriangleBackFacing(const Triangle& triangle, Camera *camera);
};

#endif
