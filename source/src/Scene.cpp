#include "../include/Scene.h"

#define PI 3.14159265

using namespace tinyxml2;
using namespace std;

/////////////////////////////////////////////////////////////////////////////////
bool isVisible(double den, double num, double &t_E, double &t_L) {
    double t = num / den;
    if (den > 0) {
        if (t > t_L) {
            return false;
        }
        if (t > t_E) {
            t_E = t;
        }
    } else if (den < 0) {
        if (t < t_E) {
            return false;
        }
        if (t < t_L) {
            t_L = t;
        }
    } else if (num > 0) {
        return false;
    }
    return true;
}


bool clippedLine(Scene *scene,Vec4 & v1_t, Vec4 & v2_t, Color & c1, Color & c2) {
    bool result = false;
    Vec4 v0 = v1_t, v1 = v2_t;
	int bflag = 0;
    double tE = 0, tL = 1;
	double d[]= {v1.x-v0.x,v1.y-v0.y,v1.z-v0.z};
	double  min[] = {-1,-1,-1};
	double max[] = {1,1,1};
	double * v[] = {&v0.x,&v0.y,&v0.z};
	Color dcolor;
	dcolor.r = c2.r - c1.r;
	dcolor.g = c2.g - c1.g;
	dcolor.b = c2.b - c1.b;
	Color *c1_t =new Color(c1);
	Color *c2_t = new Color(c2);
	int newColorId1 = v0.colorId;
	int newColorId2 = v1.colorId;

	
	
	for (int i = 0; i < 3; i++) {
		if (isVisible(d[i], min[i]-(*v[i]), tE, tL) && isVisible(-d[i], (*v[i])-max[i], tE, tL)) {
			bflag ++;
		}
	}
		if(bflag == 3){
			result = true;
			/* Line is visible */
			if (tL < 1.0) {
				v1.x = v0.x + (d[0] * tL);
				v1.y = v0.y + (d[1] * tL);
				v1.z = v0.z + (d[2] * tL);
				c2_t->r = c1_t->r + (dcolor.r * tL);
				c2_t->g = c1_t->g + (dcolor.g * tL);
				c2_t->b = c1_t->b + (dcolor.b * tL);
				scene->colorsOfVertices.push_back(c2_t);
				newColorId2= scene->colorsOfVertices.size();

			}
			if (tE > 0.0) {
				v0.x = v0.x + (d[0] * tE);
				v0.y = v0.y + (d[1] * tE);
				v0.z = v0.z + (d[2] * tE);
				c1_t->r = c1_t->r + (dcolor.r * tE);
				c1_t->g = c1_t->g + (dcolor.g * tE);
				c1_t->b = c1_t->b + (dcolor.b * tE);
				scene->colorsOfVertices.push_back(c1_t);
				newColorId1= scene->colorsOfVertices.size();
			}
		}

	v1_t = v0;
	v1_t.colorId = newColorId1;
	v2_t = v1;
	v2_t.colorId = newColorId2;
	return result;
}

void drawLine(Scene *scene, Vec4&v1, Vec4 &v2, Camera *camera){
	int x,y;
	double d;
	double y0,y1;
	double x0,x1;
	Color c;
	Color dc;
	double slope = 0;
	if(x1!=x0) slope = abs(v1.y-v2.y)/abs(v1.x-v2.x);
	if(slope<=1){
		if( v1.x <= v2.x ){
		x0 = v1.x;
		x1 = v2.x;
		y0 = v1.y;
		y1 = v2.y;
		c = *scene->colorsOfVertices[v1.colorId-1];
		dc.r = (scene->colorsOfVertices[v2.colorId-1]->r - scene->colorsOfVertices[v1.colorId-1]->r)/(x1-x0);
		dc.g = (scene->colorsOfVertices[v2.colorId-1]->g - scene->colorsOfVertices[v1.colorId-1]->g)/(x1-x0);
		dc.b = (scene->colorsOfVertices[v2.colorId-1]->b - scene->colorsOfVertices[v1.colorId-1]->b)/(x1-x0);
		}
		else{
		x0 = v2.x;
		x1 = v1.x;
		y0 = v2.y;
		y1 = v1.y;
		c = *scene->colorsOfVertices[v2.colorId-1];
		dc.r = (scene->colorsOfVertices[v1.colorId-1]->r - scene->colorsOfVertices[v2.colorId-1]->r)/(x1-x0);
		dc.g = (scene->colorsOfVertices[v1.colorId-1]->g - scene->colorsOfVertices[v2.colorId-1]->g)/(x1-x0);
		dc.b = (scene->colorsOfVertices[v1.colorId-1]->b - scene->colorsOfVertices[v2.colorId-1]->b)/(x1-x0);
		}
		y = y0;
		int dy = y0<y1? (y0-y1) : (y1-y0);
		d = 2*dy + 0.5*(x1-x0);

		for(x = x0; x <= x1; x+=1){
			if(x<0) continue;
			if(y<0){
				if(y0<y1) {y+=1; continue;}
				else break;
			}
			if(x >= camera->horRes) break;
			if(y >= (camera->verRes)){
				if(y0<y1) break;
				else y-=1;
				continue;
			}
			scene->image[x][y] = c;
			if(d<0){
				if(y0<y1) y+=1;
				else y-=1;
				d+= 2*(dy + (x1-x0));
			}
			else{
				d+= 2*dy;
			}
			c.r += dc.r;
			c.g += dc.g;
			c.b += dc.b;
		}
	}
	else{
		if( v1.y <= v2.y ){
			y0 = v1.y;
			y1 = v2.y;
			x0 = v1.x;
			x1 = v2.x;
			c = *scene->colorsOfVertices[v1.colorId-1];
			dc.r = (scene->colorsOfVertices[v2.colorId-1]->r - scene->colorsOfVertices[v1.colorId-1]->r)/(y1-y0);
			dc.g = (scene->colorsOfVertices[v2.colorId-1]->g - scene->colorsOfVertices[v1.colorId-1]->g)/(y1-y0);
			dc.b = (scene->colorsOfVertices[v2.colorId-1]->b - scene->colorsOfVertices[v1.colorId-1]->b)/(y1-y0);
		}
		else{
			y0 = v2.y;
			y1 = v1.y;
			x0 = v2.x;
			x1 = v1.x;
			c = *scene->colorsOfVertices[v2.colorId-1];
			dc.r = (scene->colorsOfVertices[v1.colorId-1]->r - scene->colorsOfVertices[v2.colorId-1]->r)/(y1-y0);
			dc.g = (scene->colorsOfVertices[v1.colorId-1]->g - scene->colorsOfVertices[v2.colorId-1]->g)/(y1-y0);
			dc.b = (scene->colorsOfVertices[v1.colorId-1]->b - scene->colorsOfVertices[v2.colorId-1]->b)/(y1-y0);
		}
		x = x0;
		int dy = x0<x1? (x0-x1) : (x1-x0);
		d = 2*dy + 0.5*(y1-y0);

		for(y = y0; y <= y1; y+=1){
			if(y<0) continue;
			if(x<0){
				if(x0<x1) {x+=1; continue;}
				else break;
			}
			if(y >= camera->verRes) break;
			if(x >= (camera->horRes)){
				if(x0<x1) break;
				else x-=1;
				continue;
			}
			
			scene->image[x][y] = c;
			if(d<0){
				if(x0<x1) x+=1;
				else x-=1;
				d+= 2*(dy + (y1-y0));
			}
			else{
				d+= 2*dy;
			}
			c.r += dc.r;
			c.g += dc.g;
			c.b += dc.b;
		}
	}
	
}

bool clipping(Scene& scene, Vec4 &vertex1, Vec4 &vertex2) { //Liang-Barsky Algorithm is implemented
    Color *color1 = scene.colorsOfVertices[vertex1.colorId - 1];
    Color *color2 = scene.colorsOfVertices[vertex2.colorId - 1];

    double t_E = 0;
    double t_L = 1;

    double dx = abs(vertex2.x - vertex1.x);
    double dy = abs(vertex2.y - vertex1.y);
    double dz = abs(vertex2.z - vertex1.z);

    Color dc;
    dc.r = color2->r - color1->r;
    dc.g = color2->g - color1->g;
    dc.b = color2->b - color1->b;

    double x_min = -1;
    double x_max = 1;
    double y_min = -1;
    double y_max = 1;
    double z_min = -1;
    double z_max = 1;

    bool lineVisible = true;
    if (isVisible(dx, x_min - vertex1.x, t_E, t_L)) {//left
        if (isVisible(-dx, vertex1.x - x_max, t_E, t_L)) {//right
            if (isVisible(dy, y_min - vertex1.y, t_E, t_L)) {//bottom
                if (isVisible(-dy, vertex1.y - y_max, t_E, t_L)) {//top
                    if (isVisible(dz, z_min - vertex1.z, t_E, t_L)) {//front
                        if (isVisible(-dz, vertex1.z - z_max, t_E, t_L)) {//back
                            lineVisible = true;
                            if (t_L < 1) {
                                vertex2.x = vertex1.x + t_L * dx;
                                vertex2.y = vertex1.y + t_L * dy;
                                vertex2.z = vertex1.z + t_L * dz;
                                color2->r = color1->r + t_L * dc.r;
                                color2->g = color1->g + t_L * dc.g;
                                color2->b = color1->b + t_L * dc.b;
                            }
                            if (t_E > 0) {
                                vertex1.x = vertex1.x + t_E * dx;
                                vertex1.y = vertex1.y + t_E * dy;
                                vertex1.z = vertex1.z + t_E * dz;
                                color1->r = color1->r + t_E * dc.r;
                                color1->g = color1->g + t_E * dc.g;
                                color1->b = color1->b + t_E * dc.b;
                            }
                        }
                    }
                }
            }
        }
    }

    return lineVisible;
}

/////////////////////////////////////////////////////////////////////////////////

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
Matrix4 createPerspectiveProjectionMatrix(Camera *camera)
{
    double perspectiveMatrixValues[4][4] = {
        {(2 * camera->near) / (camera->right - camera-> left), 0, (camera->right + camera-> left) / (camera->right - camera-> left), 0},
        {0, (2 * camera->near) / (camera->top - camera-> bottom), (camera->top + camera-> bottom) / (camera->top - camera-> bottom), 0},
        {0, 0, -(camera->far + camera->near) / (camera->far - camera->near), -(2 * camera->far * camera->near) / (camera->far - camera->near)},
        {0, 0, -1, 0}
    };

    return Matrix4(perspectiveMatrixValues);
}

Matrix4 createViewportMatrix(Camera *camera)
{
	double viewportMatrixValues[4][4] = {
		{camera->horRes / 2.0, 0, 0, (camera->horRes - 1) / 2.0},
		{0, camera->verRes / 2.0, 0, (camera->verRes - 1) / 2.0},
		{0, 0, 0.5, 0.5},
		{0, 0, 0, 1.0}
	};

	return Matrix4(viewportMatrixValues);
}

Matrix4 Scene::createModelTransformationMatrix(Mesh *mesh)
{
	Matrix4 returnMatrix = getIdentityMatrix();

	for (int i = 0; i < mesh->numberOfTransformations; i++)
	{
		char transformationType = mesh->transformationTypes[i];
		int transformationId = mesh->transformationIds[i] - 1;

		switch (transformationType)
		{
			case 't':
			{
				Translation *translation = this->translations[transformationId];
				double translationMatrixValues[4][4] = 
				{
					{1, 0, 0, translation->tx},
					{0, 1, 0, translation->ty},
					{0, 0, 1, translation->tz},
					{0, 0, 0, 1}
				};
				Matrix4 translationMatrix(translationMatrixValues);
				returnMatrix = multiplyMatrixWithMatrix(translationMatrix, returnMatrix);
				break;
			}
			case 's':
			{
				Scaling *scaling = this->scalings[transformationId];
				double scalingMatrixValues[4][4] = 
				{
					{scaling->sx, 0, 0, 0},
					{0, scaling->sy, 0, 0},
					{0, 0, scaling->sz, 0},
					{0, 0, 0, 1}
				};
				Matrix4 scalingMatrix(scalingMatrixValues);
				returnMatrix = multiplyMatrixWithMatrix(scalingMatrix, returnMatrix);
				break;
			}
			case 'r':
			{
				Rotation *rotation = this->rotations[transformationId];
				double angle = rotation->angle * PI / 180.0;
				double rotationMatrixValues[4][4] = 
				{
					{rotation->ux * rotation->ux * (1 - cos(angle)) + cos(angle), rotation->ux * rotation->uy * (1 - cos(angle)) - rotation->uz * sin(angle), rotation->ux * rotation->uz * (1 - cos(angle)) + rotation->uy * sin(angle), 0},
					{rotation->uy * rotation->ux * (1 - cos(angle)) + rotation->uz * sin(angle), rotation->uy * rotation->uy * (1 - cos(angle)) + cos(angle), rotation->uy * rotation->uz * (1 - cos(angle)) - rotation->ux * sin(angle), 0},
					{rotation->uz * rotation->ux * (1 - cos(angle)) - rotation->uy * sin(angle), rotation->uz * rotation->uy * (1 - cos(angle)) + rotation->ux * sin(angle), rotation->uz * rotation->uz * (1 - cos(angle)) + cos(angle), 0},
					{0, 0, 0, 1}
				};
				Matrix4 rotationMatrix(rotationMatrixValues);
				returnMatrix = multiplyMatrixWithMatrix(rotationMatrix, returnMatrix);
				break;
			}
		}
	}

	return returnMatrix;
}

Matrix4 Scene::createCameraTransformationMatrix(Camera* camera)
{
	double rotationMatrixValues[4][4] = 
	{
		{camera->u.x, camera->u.y, camera->u.z, 0},
		{camera->v.x, camera->v.y, camera->v.z, 0},
		{camera->w.x, camera->w.y, camera->w.z, 0},
		{0, 0, 0, 1}
	};
	Matrix4 rotationMatrix(rotationMatrixValues);

	double translationMatrixValues[4][4] = 
	{
		{1, 0, 0, -camera->position.x},
		{0, 1, 0, -camera->position.y},
		{0, 0, 1, -camera->position.z},
		{0, 0, 0, 1}
	};
	Matrix4 translationMatrix(translationMatrixValues);

	return multiplyMatrixWithMatrix(rotationMatrix, translationMatrix);
}

bool cull_triangle(Vec4 &vertex1, Vec4 &vertex2, Vec4 &vertex3){
	//Take the vertices and calculate the normal vector.
	Vec3 edge2_1 = Vec3(vertex2.x - vertex1.x, vertex2.y - vertex1.y, vertex2.z - vertex1.z, -1);
	Vec3 edge3_1 = Vec3(vertex3.x - vertex1.x, vertex3.y - vertex1.y, vertex3.z - vertex1.z, -1);

	Vec3 normal_vector = normalizeVec3(crossProductVec3(edge2_1, edge3_1));

	double dot_product = dotProductVec3(normal_vector, Vec3(vertex1.x, vertex1.y, vertex1.z, -1));

	return dot_product < 0;
}


/*
	Applies model transformations (translation, rotation, scaling) to the mesh.
*/
void Scene::applyModelTransformation(Mesh *mesh, Triangle& triangle)
{
	Vec3 *vertex1 = this->vertices[triangle.vertexIds[0] - 1];
	Vec3 *vertex2 = this->vertices[triangle.vertexIds[1] - 1];
	Vec3 *vertex3 = this->vertices[triangle.vertexIds[2] - 1];

	for (int k = 0; k < mesh->numberOfTransformations; k++)
	{
		char transformationType = mesh->transformationTypes[k];
		int transformationId = mesh->transformationIds[k];

		switch (transformationType)
		{
			case 't':
			{
				Translation *translation = this->translations[transformationId - 1];
				*vertex1 = vertex1->translateVec3(*translation);
				*vertex2 = vertex2->translateVec3(*translation);
				*vertex3 = vertex3->translateVec3(*translation);
				break;
			}
			case 's':
			{
				Scaling *scaling = this->scalings[transformationId - 1];
				*vertex1 = vertex1->scaleVec3(*scaling);
				*vertex2 = vertex2->scaleVec3(*scaling);
				*vertex3 = vertex3->scaleVec3(*scaling);
				break;
			}
			case 'r':
			{
				Rotation *rotation = this->rotations[transformationId - 1];
				*vertex1 = vertex1->rotateVec3(*rotation);
				*vertex2 = vertex2->rotateVec3(*rotation);
				*vertex3 = vertex3->rotateVec3(*rotation);
				break;
			}
		}
	}
}

/*
	Applies camera transformations (translation, rotation) to the triangle.
*/
void Scene::applyCameraTransformation(Camera *camera, Triangle& triangle)
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
    Matrix4 viewMatrix = multiplyMatrixWithMatrix(rotationMatrix, translationMatrix);

	Vec3 *vertex1 = this->vertices[triangle.vertexIds[0] - 1];
	Vec3 *vertex2 = this->vertices[triangle.vertexIds[1] - 1];
	Vec3 *vertex3 = this->vertices[triangle.vertexIds[2] - 1];

	// Convert to homogeneous coordinates
	Vec4 vertex1Homogeneous(vertex1->x, vertex1->y, vertex1->z, 1, vertex1->colorId); 
	Vec4 vertex2Homogeneous(vertex2->x, vertex2->y, vertex2->z, 1, vertex2->colorId); 
	Vec4 vertex3Homogeneous(vertex3->x, vertex3->y, vertex3->z, 1, vertex3->colorId); 

	// Apply view transformation
	Vec4 transformedVertex1 = multiplyMatrixWithVec4(viewMatrix, vertex1Homogeneous);
	Vec4 transformedVertex2 = multiplyMatrixWithVec4(viewMatrix, vertex2Homogeneous);
	Vec4 transformedVertex3 = multiplyMatrixWithVec4(viewMatrix, vertex3Homogeneous);

	// Apply projection transformation
	Vec4 projectedVertex1 = multiplyMatrixWithVec4(projectionMatrix, transformedVertex1);
	Vec4 projectedVertex2 = multiplyMatrixWithVec4(projectionMatrix, transformedVertex2);
	Vec4 projectedVertex3 = multiplyMatrixWithVec4(projectionMatrix, transformedVertex3);

	Vec4 tempProjectedVertex1 = projectedVertex1;
	Vec4 tempProjectedVertex2 = projectedVertex2;
	Vec4 tempProjectedVertex3 = projectedVertex3;

	// Perspective division
	if(camera->projectionType == 1){
		projectedVertex1.x /= projectedVertex1.t;
		projectedVertex1.y /= projectedVertex1.t;
		projectedVertex1.z /= projectedVertex1.t;

		projectedVertex2.x /= projectedVertex2.t;
		projectedVertex2.y /= projectedVertex2.t;
		projectedVertex2.z /= projectedVertex2.t;

		projectedVertex3.x /= projectedVertex3.t;
		projectedVertex3.y /= projectedVertex3.t;
		projectedVertex3.z /= projectedVertex3.t;
	}

	// Convert back to 3D coordinates
	*vertex1 = Vec3(projectedVertex1.x, projectedVertex1.y, projectedVertex1.z, projectedVertex1.colorId);
	*vertex2 = Vec3(projectedVertex2.x, projectedVertex2.y, projectedVertex2.z, projectedVertex2.colorId);
	*vertex3 = Vec3(projectedVertex3.x, projectedVertex3.y, projectedVertex3.z, projectedVertex3.colorId);
}

/*
	Applies viewport transformation to the triangle.
*/
void Scene::applyViewportTransformation(Camera *camera, Triangle& triangle)
{
	double viewportMatrixValues[4][4] = { // viewport in the origin. xmin, ymin 0 
		{camera->horRes / 2.0, 0, 0, (camera->horRes - 1) / 2.0},
		{0, camera->verRes / 2.0, 0, (camera->verRes - 1) / 2.0},
		{0, 0, 0.5, 0.5},
		{0, 0, 0, 1}
	};
	Matrix4 viewportMatrix = Matrix4(viewportMatrixValues);

	Vec3 *vertex1 = this->vertices[triangle.vertexIds[0] - 1];
	Vec3 *vertex2 = this->vertices[triangle.vertexIds[1] - 1];
	Vec3 *vertex3 = this->vertices[triangle.vertexIds[2] - 1];

	// Convert to homogeneous coordinates
	Vec4 vertex1Homogeneous(vertex1->x, vertex1->y, vertex1->z, 1, vertex1->colorId);
	Vec4 vertex2Homogeneous(vertex2->x, vertex2->y, vertex2->z, 1, vertex2->colorId);
	Vec4 vertex3Homogeneous(vertex3->x, vertex3->y, vertex3->z, 1, vertex3->colorId);

	// Apply viewport transformation
	Vec4 transformedVertex1 = multiplyMatrixWithVec4(viewportMatrix, vertex1Homogeneous);
	Vec4 transformedVertex2 = multiplyMatrixWithVec4(viewportMatrix, vertex2Homogeneous);
	Vec4 transformedVertex3 = multiplyMatrixWithVec4(viewportMatrix, vertex3Homogeneous);

	// Convert back to 3D coordinates
	*vertex1 = Vec3(transformedVertex1.x, transformedVertex1.y, transformedVertex1.z, transformedVertex1.colorId);
	*vertex2 = Vec3(transformedVertex2.x, transformedVertex2.y, transformedVertex2.z, transformedVertex2.colorId);
	*vertex3 = Vec3(transformedVertex3.x, transformedVertex3.y, transformedVertex3.z, transformedVertex3.colorId);
}

/*
	Liang-Barsky Algorithm for clipping	
*/
bool Scene::liangBarskyClip(Camera *camera, Vec4 &p0, Vec4 &p1) /// ??????????????
{
    // double dx = abs(p1.x - p0.x);
    // double dy = abs(p1.y - p0.y);
	// double dz = abs(p1.z - p0.z);
    // double p[6], q[6];
    // double t0 = 0.0;
    // double t1 = 1.0;

	// p[0] = dx;  q[0] = 0 - 1 - p0.x; // Left
	// p[1] = -dx; q[1] = p0.x - 1; // Right
	// p[2] = dy;  q[2] = 0 - 1 - p0.y; // Bottom
	// p[3] = -dy; q[3] = p0.y - 1;  // Top
	// p[4] = dz;  q[4] = 0 - 1 - p0.z; // Front
	// p[5] = -dz; q[5] = p0.z - 1;  // Back

    // for (int i = 0; i < 6; i++)
    // {
    //     if (p[i] > 0)
    //     {
    //         double t = q[i] / p[i];
	// 		if (t > t1) return false; // Line is outside
	// 		if (t > t0) t0 = t; // Line enters clipping area
    //     }
    //     else if (p[i] < 0)
	// 	{
	// 		double t = q[i] / p[i];
	// 		if (t < t0) return false; // Line is outside
	// 		if (t < t1) t1 = t; // Line exits clipping area
	// 	}
	// 	else if (q[i] > 0) return false; // Line is outside
    // }

    // if (t1 < 1)
    // {
    //     p1.x = p0.x + t1 * dx;
    //     p1.y = p0.y + t1 * dy;
	// 	p1.z = p0.z + t1 * dz;
    // }
    // if (t0 > 0)
    // {
    //     p0.x += t0 * dx;
    //     p0.y += t0 * dy;
	// 	p0.z += t0 * dz;
    // }

    // return true; // Line is inside or intersects the clipping area

	double x_min = -1;
	double x_max = 1;
	double y_min = -1;
	double y_max = 1;
	double z_min = -1;
	double z_max = 1;

	double t_E = 0;
	double t_L = 1;

	double dx = p1.x - p0.x;
	double dy = p1.y - p0.y;
	double dz = p1.z - p0.z;

	double p[6] = {dx, -dx, dy, -dy, dz, -dz}; //Den

	double q[6] = {x_min - p0.x, p0.x - x_max, y_min - p0.y, p0.y - y_max, z_min - p0.z, p0.z - z_max}; //Num

	for(int i = 0; i < 6; i++){
		if(p[i] > 0){
			double t = q[i] / p[i];

			if(t > t_L){
				return false;
			}
			else if(t > t_E){
				t_E = t;
			}
		}
		else if(p[i] < 0){
			double t = q[i] / p[i];

			if(t < t_E){
				return false;
			}
			else if(t < t_L){
				t_L = t;
			}
		}
		else if(q[i] > 0){
			return false;
		}
	}

	if(t_L < 1){
		p1.x = p0.x + t_L * dx;
		p1.y = p0.y + t_L * dy;
		p1.z = p0.z + t_L * dz;
	}

	if(t_E > 0){
		p1.x = p0.x + t_E * dx;
		p1.y = p0.y + t_E * dy;
		p1.z = p0.z + t_E * dz;
	}
	

	return true;
}

Color Scene::interpolateColor(const Color &c1, const Color &c2, double t)
{
    double r = c1.r + (c2.r - c1.r) * t;
    double g = c1.g + (c2.g - c1.g) * t;
    double b = c1.b + (c2.b - c1.b) * t;
    return Color(makeBetweenZeroAnd255(r), makeBetweenZeroAnd255(g), makeBetweenZeroAnd255(b));
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
    Vec3 normal = normalizeVec3(crossProductVec3(edge1, edge2));

    // Calculate view direction (from triangle to camera)
    Vec3 viewDir = camera->position - v0;

    // Check if the dot product of the normal and view direction is positive
    return dotProductVec3(normal, viewDir) > 0; // ?????
}

/*
	Rasterizes triangle.
*/
void rasterize_line(Vec4 &v1, Vec4 &v2, Color c1, Color c2, vector< vector<Color> > &image, int horRes, int verRes){
	double dx = abs(v2.x - v1.x);
	double dy = abs(v2.y - v1.y);

	int step;
	Color color, dc;
	double d, x, y;

	// slope m = (0,1]
	if(dx >= dy){

		if(v2.x < v1.x){
			swap(v1, v2);
			swap(c1, c2);
		}

		if(v2.y >= v1.y){
			step = 1;
		}
		else{
			step = -1;
		}

		color = c1;
		d = (v1.y - v2.y) + 0.5 * step * (v2.x - v1.x);
		dc = (c2 - c1) / (v2.x - v1.x);
		y = v1.y;

		for(x = v1.x; x <= v2.x; x++){
			if(x >= 0 && x <= horRes && y >= 0 && y <= verRes)  {
				image[x][y] = color;
			}

			if(d * step < 0){
				y = y + step;
				d += (v1.y - v2.y) + step * (v2.x - v1.x);
			}
			
			else{
				d += (v1.y - v2.y);
			}
			color = color + dc;
		}

	}

	// slope m = (1, inf)
	else{
		if(v2.y < v1.y){
			swap(v1, v2);
			swap(c1, c2);
		}

		if(v2.x >= v1.x){
			step = 1;
		}
		else{
			step = -1;
		}

		color = c1;
		d = (v2.x - v1.x) + 0.5 * step * (v1.y - v2.y);
		dc = (c2 - c1) / (v2.y - v1.y);
		x = v1.x;

		for(y = v1.y; y <= v2.y; y++){
			if(x >= 0 && x <= horRes && y >= 0 && y <= verRes)  {
				image[x][y] = color;
			}

			if(d * step > 0){
				x = x + step;
				d += (v2.x - v1.x) + step * (v1.y - v2.y);
			}
			
			else{
				d += (v2.x - v1.x);
			}
			color = color + dc;
		}
	}

}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	Vec3 xAxis = crossProductVec3(camera->gaze, camera->v);
	camera->u = normalizeVec3(xAxis);
	Matrix4 cameraMatrix = createCameraTransformationMatrix(camera);
	Matrix4 projectionMatrix = (camera->projectionType == 0) ? createOrthographicProjectionMatrix(camera) : createPerspectiveProjectionMatrix(camera);
	Matrix4 viewportMatrix = createViewportMatrix(camera);

    for (Mesh *mesh : this->meshes)
    {
		// printf ("triangle count: %d\n", mesh->triangles.size());
		Matrix4 modelMatrix = createModelTransformationMatrix(mesh);
		Matrix4 modelViewMatrix = multiplyMatrixWithMatrix(cameraMatrix, modelMatrix);
		Matrix4 modelViewProjectionMatrix = multiplyMatrixWithMatrix(projectionMatrix, modelViewMatrix);

		for (int i = 0; i < mesh->triangles.size(); ++i)
        {
			Triangle triangle = mesh->triangles[i];
            // Apply transformations
            // applyModelTransformation(mesh, triangle);
            // applyCameraTransformation(camera, triangle);

            if (mesh->type == 0) // wireframe
            {
				// printf ("wireframe, camera %d\n", camera->cameraId);
                Vec3 &v0 = *vertices[triangle.vertexIds[0] - 1];
                Vec3 &v1 = *vertices[triangle.vertexIds[1] - 1];
                Vec3 &v2 = *vertices[triangle.vertexIds[2] - 1];

				Color c0 = *colorsOfVertices[v0.colorId - 1]; Color c0temp = *colorsOfVertices[v0.colorId - 1];
				Color c1 = *colorsOfVertices[v1.colorId - 1]; Color c1temp = *colorsOfVertices[v1.colorId - 1];
				Color c2 = *colorsOfVertices[v2.colorId - 1]; Color c2temp = *colorsOfVertices[v2.colorId - 1];

				Vec4 v0Homogeneous(v0.x, v0.y, v0.z, 1, v0.colorId);
				Vec4 v1Homogeneous(v1.x, v1.y, v1.z, 1, v1.colorId);
				Vec4 v2Homogeneous(v2.x, v2.y, v2.z, 1, v2.colorId);

				v0Homogeneous = multiplyMatrixWithVec4(modelViewProjectionMatrix, v0Homogeneous);
				v1Homogeneous = multiplyMatrixWithVec4(modelViewProjectionMatrix, v1Homogeneous);
				v2Homogeneous = multiplyMatrixWithVec4(modelViewProjectionMatrix, v2Homogeneous);

				if(cullingEnabled){
					if(cull_triangle(v0Homogeneous, v1Homogeneous, v2Homogeneous)){
						continue;
					}
				}

				if (camera->projectionType == 1) // perspective
				{
					v0Homogeneous = v0Homogeneous / v0Homogeneous.t;
					v1Homogeneous = v1Homogeneous / v1Homogeneous.t;
					v2Homogeneous = v2Homogeneous / v2Homogeneous.t;
				}

				Vec4 v0tempHomogeneous = v0Homogeneous;
				Vec4 v1tempHomogeneous = v1Homogeneous;
				Vec4 v2tempHomogeneous = v2Homogeneous;

				bool edgeOneVisible = liangBarskyClip(camera, v0Homogeneous, v1Homogeneous);
				bool edgeTwoVisible = liangBarskyClip(camera, v1tempHomogeneous, v2Homogeneous);
				bool edgeThreeVisible = liangBarskyClip(camera, v2tempHomogeneous, v0tempHomogeneous);

				v0Homogeneous = multiplyMatrixWithVec4(viewportMatrix, v0Homogeneous);
				v1Homogeneous = multiplyMatrixWithVec4(viewportMatrix, v1Homogeneous);
				v2Homogeneous = multiplyMatrixWithVec4(viewportMatrix, v2Homogeneous);
				v0tempHomogeneous = multiplyMatrixWithVec4(viewportMatrix, v0tempHomogeneous);
				v1tempHomogeneous = multiplyMatrixWithVec4(viewportMatrix, v1tempHomogeneous);
				v2tempHomogeneous = multiplyMatrixWithVec4(viewportMatrix, v2tempHomogeneous);

                if (edgeOneVisible)
				{
					rasterize_line(v0Homogeneous, v1Homogeneous, c0, c1, this->image, camera->horRes, camera->verRes);
				}
				if (edgeTwoVisible)
				{
					rasterize_line(v1tempHomogeneous, v2Homogeneous, c1temp, c2, this->image, camera->horRes, camera->verRes);
				}
				if (edgeThreeVisible)
				{
					rasterize_line(v2tempHomogeneous, v0tempHomogeneous, c2temp, c0temp, this->image, camera->horRes, camera->verRes);
				}

				// printf("clipping\n");
				// if (clipping(*this, v0Homogeneous, v1Homogeneous))
				// {
				// 	printf ("lineRasterizer v0-v1 \n");
				// 	drawLine(this, v0Homogeneous, v1Homogeneous, camera);
				// }
				// if (clipping(*this, v1tempHomogeneous, v2Homogeneous))
				// {
				// 	printf ("lineRasterizer v1-v2 \n");
				// 	drawLine(this, v1tempHomogeneous, v2Homogeneous, camera);
				// }
				// if (clipping(*this, v2tempHomogeneous, v0tempHomogeneous))
				// {
				// 	printf ("lineRasterizer v2-v0 \n");
				// 	drawLine(this, v2tempHomogeneous, v0tempHomogeneous, camera);
				// }

				// printf("--------------------\n");

            }
            else // solid
            {
                // Rasterize the triangle
            }
        }
    }
}

// /*
// 	Clips triangles that are outside of the view volume.
// */
// uint8_t Scene::computeOutcode(Camera *camera, double x, double y)
// {
//     uint8_t code = 0;
    
//     if (x < camera->left) code |= 1;
//     else if (x > camera->right) code |= 2;
//     if (y < camera->bottom) code |= 4;
//     else if (y > camera->top) code |= 8;

//     return code;
// }

// void Scene::clipLine(Camera *camera, Vec3 &p0, Vec3 &p1)
// {
//     double x0 = p0.x, y0 = p0.y;
//     double x1 = p1.x, y1 = p1.y;

//     uint8_t outcode0 = computeOutcode(camera, x0, y0);
//     uint8_t outcode1 = computeOutcode(camera, x1, y1);
//     bool accept = false;

//     while (true)
// 	{
//         if (!(outcode0 | outcode1))
// 		{ // Both points inside
//             accept = true;
//             break;
//         }
// 		else if (outcode0 & outcode1)
// 		{ // Both points share an outside zone
//             break;
//         }
// 		else
// 		{
//             // At least one endpoint is outside the clipping window
//             double x, y;
//             uint8_t outcodeOut = outcode0 ? outcode0 : outcode1;

//             // Find intersection point using formulas y = y0 + slope * (x - x0), x = x0 + (1/slope) * (y - y0)
// 			switch (outcodeOut)
// 			{
// 				case 1:
// 					y = y0 + (y1 - y0) * (camera->left - x0) / (x1 - x0);
// 					x = camera->left;
// 					break;
// 				case 2:
// 					y = y0 + (y1 - y0) * (camera->right - x0) / (x1 - x0);
// 					x = camera->right;
// 					break;
// 				case 4:
// 					x = x0 + (x1 - x0) * (camera->bottom - y0) / (y1 - y0);
// 					y = camera->bottom;
// 					break;
// 				case 8:
// 					x = x0 + (x1 - x0) * (camera->top - y0) / (y1 - y0);
// 					y = camera->top;
// 					break;
// 			}

//             // Replace point outside clipping area with intersection point
//             if (outcodeOut == outcode0)
// 			{
//                 x0 = x; y0 = y;
//                 outcode0 = computeOutcode(camera, x0, y0);
//             }
// 			else
// 			{
//                 x1 = x; y1 = y;
//                 outcode1 = computeOutcode(camera, x1, y1);
//             }
//         }
//     }
//     if (accept)
// 	{
//         // Update the points of the line
//         p0.x = x0; p0.y = y0;
//         p1.x = x1; p1.y = y1;

// 		// Draw the line
// 		drawLine(camera, &p0, &p1, this->colorsOfVertices[p0.colorId - 1], this->colorsOfVertices[p1.colorId - 1]); // ????
//     }
//     // Else, the line is outside the clipping area and should not be drawn 
// }

// void Scene::clipTriangle(Camera *camera, Triangle& triangle)
// {
// 	// Apply clipping to each edge of the triangle
// 	clipLine(camera, *vertices[triangle.vertexIds[0] - 1], *vertices[triangle.vertexIds[1] - 1]);
// 	clipLine(camera, *vertices[triangle.vertexIds[1] - 1], *vertices[triangle.vertexIds[2] - 1]);
// 	clipLine(camera, *vertices[triangle.vertexIds[2] - 1], *vertices[triangle.vertexIds[0] - 1]);
// }

