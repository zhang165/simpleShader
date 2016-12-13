#include "parser.hpp"
#include "raytracer.hpp"
#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>


int main(int argc, char **argv)
{
	string scenefile;	//If aren't specified, it will render the default scene with only a cube
	string outfilename("scenes/default_out.bmp");	//NOTE: here we only support .bmp
	string outfilename_depth("scenes/default_depth_out.bmp");
	bool creative = false;

    // Scene filename specified
	if (argc > 1)
	{
		scenefile = string(argv[1]);
		creative = (argc > 2 && string(argv[2]) == "c");
		int dot = scenefile.find_last_of('.');
		outfilename = scenefile.substr(0, dot) + "_out.bmp";
		outfilename_depth = scenefile.substr(0, dot) + "_depth_out.bmp";
	}

    std::cout << "Rendering " << (scenefile.empty() ? "default scene" : scenefile) << std::endl;
    std::cout << "Output to " << outfilename << std::endl;
   
	Raytracer raytracer;
	if (scenefile.empty())	{
		//render default scene
		raytracer.render(
			outfilename.c_str(),
			outfilename_depth.c_str(),
			Scene(),
			creative);
	}
	else
	{
		// Parse the scene file
		Parser parser(new std::ifstream(scenefile.c_str()));
		if (!parser.parse()) {
			puts("Scene file can't be parsed. Use default scene.");
			raytracer.render(
				outfilename.c_str(),
				outfilename_depth.c_str(),
				Scene(),
				creative);
		}
		else
		{
			// Render the input scene with our raytracer.
			raytracer.render(
				outfilename.c_str(),
				outfilename_depth.c_str(),
				parser.scene,
				creative);
		}
	}
	
	// Use if you're running visual studio:
	//system("pause");
	return 0;
}

