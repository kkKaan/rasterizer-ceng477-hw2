#include "../include/Scene.h"

#include <iostream>
#include <vector>

using namespace std;

Scene* scene;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cout << "Please run the rasterizer as:" << endl << "\t./rasterizer <input_file_name>" << endl;
        return 1;
    }
    else
    {
        const char* xmlPath = argv[1];

        scene = new Scene(xmlPath);

        for (int i = 0; i < scene->cameras.size(); i++)
        {
            // initialize image with basic values
            scene->initializeImage(scene->cameras[i]);

            // do forward rendering pipeline operations
            scene->forwardRenderingPipeline(scene->cameras[i]);

            // generate PPM file
            scene->writeImageToPPMFile(scene->cameras[i]);

            // Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
            // Change/remove implementation if necessary.
            // scene->convertPPMToPNG(scene->cameras[i]->outputFilename);
        }

        // test
        Vec3 v1(1, 0, 0);

        // cout << "v1: " << v1 << endl;

        // Rotation r1(-1, 90, 0, 1, 0); // rotate 90 degrees around y axis
        // Rotation r2(-1, 45, 0, 0, 1); // rotate 45 degrees around z axis
        // Rotation r3(-1, -60, 0, 1, 0); // rotate 60 degrees around y axis
        // v1 = v1.rotateVec3(v1, r3);

        // Translation t1(-1, 1, 2, 1);
        // v1 = v1.translateVec3(v1, t1);

        // Scaling s1(-1, 2, 2, -31);
        // v1 = v1.scaleVec3(v1, s1);

        // cout << "v1: " << v1 << endl;

        return 0;
    }
}