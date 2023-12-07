#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "../include/tinyxml2.h"
#include "../include/Triangle.h"
#include "../include/Helpers.h"
#include "../include/Scene.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command, directory;
	
	// TODO: Change implementation if necessary.
	command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}

/*
	Creates orthographic projection matrix.
*/
Matrix4 createOrthographicProjectionMatrix(Camera *camera)
{
    double orthoMatrixValues[4][4] = {
        {2 / (camera->right - camera->left), 0, 0, -(camera->right + camera->left) / (camera->right - camera->left)},
        {0, 2 / (camera->top - camera->bottom), 0, -(camera->top + camera->bottom) / (camera->top - camera->bottom)},
        {0, 0, -2 / (camera->far - camera->near), -(camera->far + camera->near) / (camera->far - camera->near)},
        {0, 0, 0, 1}
    };

    return Matrix4(orthoMatrixValues);
}

/*
	Creates perspective projection matrix.
*/
Matrix4 createPerspectiveProjectionMatrix(Camera *camera) // ??????
{
    double aspectRatio = static_cast<double>(camera->horRes) / camera->verRes;
    double fovyRadians = 2 * atan((camera->top - camera->bottom) / (2 * camera->near)); // Field of view in radians
    double f = 1.0 / tan(fovyRadians / 2); // Focal length

    double perspectiveMatrixValues[4][4] = {
        {f / aspectRatio, 0, 0, 0},
        {0, f, 0, 0},
        {0, 0, (camera->far + camera->near) / (camera->near - camera->far), 2 * camera->far * camera->near / (camera->near - camera->far)},
        {0, 0, -1, 0}
    };

    return Matrix4(perspectiveMatrixValues);
}

/*
	Applies model transformations (translation, rotation, scaling) to the mesh.
*/
void Scene::applyModelTransformations(Mesh *mesh)
{
	for (int i = 0; i < mesh->numberOfTriangles; i++)
	{
		Triangle *triangle = &mesh->triangles[i];

		for (int j = 0; j < 3; j++)
		{
			Vec3 *vertex = this->vertices[triangle->vertexIds[j] - 1];

			for (int k = 0; k < mesh->numberOfTransformations; k++)
			{
				char transformationType = mesh->transformationTypes[k];
				int transformationId = mesh->transformationIds[k];

				switch (transformationType)
				{
					case 't':
					{
						Translation *translation = this->translations[transformationId - 1];
						*vertex = vertex->translateVec3(*vertex, *translation);
						break;
					}
					case 's':
					{
						Scaling *scaling = this->scalings[transformationId - 1];
						*vertex = vertex->scaleVec3(*vertex, *scaling);
						break;
					}
					case 'r':
					{
						Rotation *rotation = this->rotations[transformationId - 1];
						*vertex = vertex->rotateVec3(*vertex, *rotation);
						break;
					}
				}
			}
		}
	}
}

/*
	Applies camera transformations (translation, rotation) to the mesh.
*/
void Scene::applyCameraTransformations(Mesh *mesh, Camera *camera)
{
    // Constructing the rotation part of the view matrix
    double rotationMatrixValues[4][4] = {
		{camera->u.x, camera->u.y, camera->u.z, 0},
		{camera->v.x, camera->v.y, camera->v.z, 0},
		{camera->w.x, camera->w.y, camera->w.z, 0},
		{0, 0, 0, 1}
	};
	Matrix4 rotationMatrix = Matrix4(rotationMatrixValues);

    // Constructing the translation part of the view matrix
    double translationMatrixValues[4][4] = {
		{1, 0, 0, -camera->position.x},
		{0, 1, 0, -camera->position.y},
		{0, 0, 1, -camera->position.z},
		{0, 0, 0, 1}
	};
	Matrix4 translationMatrix = Matrix4(translationMatrixValues);

	// Constructing the projection matrix
	Matrix4 projectionMatrix = (camera->projectionType == 0) ? createOrthographicProjectionMatrix(camera) : createPerspectiveProjectionMatrix(camera);

    // Combining rotation and translation into the view matrix
    Matrix4 viewMatrix = rotationMatrix * translationMatrix;

    for (int i = 0; i < mesh->numberOfTriangles; i++)
    {
        Triangle *triangle = &mesh->triangles[i];

        for (int j = 0; j < 3; j++)
        {
            Vec3 *vertex = this->vertices[triangle->vertexIds[j] - 1];

            // Apply view transformation
            Vec4 vertexHomogeneous(vertex->x, vertex->y, vertex->z, 1); // Convert to homogeneous coordinates
            Vec4 transformedVertex = vertexHomogeneous.multiplyMatrixVec4(viewMatrix, vertexHomogeneous); // Apply view transformation

			// Apply projection transformation
			Vec4 projectedVertex = transformedVertex.multiplyMatrixVec4(projectionMatrix, transformedVertex);

			// Perspective divide for perspective projection
			if (camera->projectionType == 1)
			{
                projectedVertex.x /= projectedVertex.t;
                projectedVertex.y /= projectedVertex.t;
                projectedVertex.z /= projectedVertex.t;
				projectedVertex.t /= projectedVertex.t;
            }

			// Convert back to 3D coordinates
			*vertex = Vec3(projectedVertex.x, projectedVertex.y, projectedVertex.z);
        }
    }
}

/*
	Applies viewport transformation to the mesh.
*/
void Scene::applyViewportTransformation(Mesh *mesh, Camera *camera)
{
	double viewportMatrixValues[4][4] = { // viewport in the origin. xmin, ymin 0
		{camera->horRes / 2.0, 0, 0, (camera->horRes - 1) / 2.0},
		{0, camera->verRes / 2.0, 0, (camera->verRes - 1) / 2.0},
		{0, 0, 0.5, 0.5},
		{0, 0, 0, 1}
	};
	Matrix4 viewportMatrix = Matrix4(viewportMatrixValues);

	for (int i = 0; i < mesh->numberOfTriangles; i++)
	{
		Triangle *triangle = &mesh->triangles[i];

		for (int j = 0; j < 3; j++)
		{
			if (triangle->vertexIds[j] == -1)
			{
				continue;
			}

			Vec3 *vertex = this->vertices[triangle->vertexIds[j] - 1];

			// Apply viewport transformation
			Vec4 vertexHomogeneous(vertex->x, vertex->y, vertex->z, 1); // Convert to homogeneous coordinates
			Vec4 transformedVertex = vertexHomogeneous.multiplyMatrixVec4(viewportMatrix, vertexHomogeneous); // Apply viewport transformation

			// Convert back to 3D coordinates
			*vertex = Vec3(transformedVertex.x + 0.5, transformedVertex.y + 0.5, transformedVertex.z);
		}
	}
}

/*
	Clips triangles that are outside of the view volume. ????????????????????????
*/
uint8_t computeOutcode(double x, double y, Camera *camera)
{
    uint8_t code = 0;
    
    if (x < camera->left) code |= 1;
    else if (x > camera->right) code |= 2;
    if (y < camera->bottom) code |= 4;
    else if (y > camera->top) code |= 8;

    return code;
}

void clipLine(Vec3 &p0, Vec3 &p1, Camera *camera)
{
    double x0 = p0.x, y0 = p0.y;
    double x1 = p1.x, y1 = p1.y;

    uint8_t outcode0 = computeOutcode(x0, y0, camera);
    uint8_t outcode1 = computeOutcode(x1, y1, camera);
    bool accept = false;

    while (true)
	{
        if (!(outcode0 | outcode1))
		{ // Both points inside
            accept = true;
            break;
        }
		else if (outcode0 & outcode1)
		{ // Both points share an outside zone
            break;
        }
		else
		{
            // At least one endpoint is outside the clipping window
            double x, y;
            uint8_t outcodeOut = outcode0 ? outcode0 : outcode1;

            // Find intersection point using formulas y = y0 + slope * (x - x0), x = x0 + (1/slope) * (y - y0)
			switch (outcodeOut)
			{
				case 1:
					y = y0 + (y1 - y0) * (camera->left - x0) / (x1 - x0);
					x = camera->left;
					break;
				case 2:
					y = y0 + (y1 - y0) * (camera->right - x0) / (x1 - x0);
					x = camera->right;
					break;
				case 4:
					x = x0 + (x1 - x0) * (camera->bottom - y0) / (y1 - y0);
					y = camera->bottom;
					break;
				case 8:
					x = x0 + (x1 - x0) * (camera->top - y0) / (y1 - y0);
					y = camera->top;
					break;
			}

            // Replace point outside clipping area with intersection point
            if (outcodeOut == outcode0)
			{
                x0 = x; y0 = y;
                outcode0 = computeOutcode(x0, y0, camera);
            }
			else
			{
                x1 = x; y1 = y;
                outcode1 = computeOutcode(x1, y1, camera);
            }
        }
    }
    if (accept)
	{
        // Update the points of the line
        p0.x = x0; p0.y = y0;
        p1.x = x1; p1.y = y1;
    }
    // Else, the line is outside the clipping area and should not be drawn ????????????????????????????????
}

void Scene::clipTriangles(Mesh *mesh, Camera *camera)
{
    for (Triangle &triangle : mesh->triangles)
	{
        // Apply clipping to each edge of the triangle
        clipLine(*vertices[triangle.vertexIds[0] - 1], *vertices[triangle.vertexIds[1] - 1], camera);
        clipLine(*vertices[triangle.vertexIds[1] - 1], *vertices[triangle.vertexIds[2] - 1], camera);
        clipLine(*vertices[triangle.vertexIds[2] - 1], *vertices[triangle.vertexIds[0] - 1], camera);
    }
}

/*
	Checks if triangle is back facing.
*/
bool Scene::isTriangleBackFacing(const Triangle& triangle, Camera *camera)
{
    // Calculate the normal of the triangle
    Vec3 v0 = *vertices[triangle.vertexIds[0] - 1];
    Vec3 v1 = *vertices[triangle.vertexIds[1] - 1];
    Vec3 v2 = *vertices[triangle.vertexIds[2] - 1];

    Vec3 edge1 = v1 - v0;
    Vec3 edge2 = v2 - v0;
    Vec3 normal = crossProductVec3(edge1, edge2);

    // Calculate view direction (from triangle to camera)
    Vec3 viewDir = camera->position - v0;

    // Check if the dot product of the normal and view direction is positive
    return dotProductVec3(normal, viewDir) > 0;
}

/*
	Rasterizes triangle.
*/


/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
    for (Mesh *mesh : this->meshes)
    {
        applyModelTransformations(mesh);
        applyCameraTransformations(mesh, camera);
		
		if (mesh->type == 0) // Wireframe
        {
            clipTriangles(mesh, camera); // ??????
        }

		applyViewportTransformation(mesh, camera);

        for (Triangle& triangle : mesh->triangles)
        {
            if (!this->cullingEnabled || !isTriangleBackFacing(triangle, camera))
            {
                rasterize(triangle, camera);
            }
        }
    }
}

