#include "InputManager.h"
// #include "../DisplayGLFW/display.h"
#include "game.h"
#include "../res/includes/glm/glm.hpp"
#include <stb_image.h>
#include <fstream>
#include <iostream>

# define M_PI           3.14159265

class ObjectValues {
public:
	int objectIndex;
	float objectColors[4];
	float objectValues[4];
	int objectType;
};

// fucntion for assignment
std::pair<std::string, std::vector<float>> Split(std::string line, char delimiter);
float IntersectSphere(std::vector<float> sphere, std::vector<float> rayCasted, std::vector<float> camera);
float IntersectPlane(std::vector<float> plane, std::vector<float> rayCasted, std::vector<float> camera);
std::vector<float> ConstructRayThroughPixel(int y, int x);
std::pair<float[4], ObjectValues> FindIntersection(std::vector<float> ray, int objectIndexToPass, std::vector<float> camera);
std::vector<float> GetColor(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int times, std::vector<float> camera);
std::vector<float> ObjectGetColor(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, std::vector<float> camera);
std::vector<float> GetObjectColor(ObjectValues obj, float* hit);
std::vector<float> DiffuseReflection(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int i, std::vector<float> camera);
std::vector<float> SpecularReflection(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int i, std::vector<float> camera);
std::vector<float> Shadows(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int i, std::vector<float> camera);
std::vector<float> ReflectiveGetColor(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int times, std::vector<float> camera);
std::vector<float> TranparentGetColor(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int times, std::vector<float> camera);
std::vector<std::vector<float>> ConstructRayThroughPixelBonus(int y, int x);

std::vector<float> camera; // [0] - x , [1] - y , [2] - z , [3] - anti aliasing
std::vector<std::vector<float>> lightDirctions; // [0] - x , [1] - y , [2] - z , [3] - 0.0 for directional light and 1.0 for spotlight.
std::vector<float> ambientLight; // [0] - R , [1] - G , [2] - B , [3] - A
std::vector<std::vector<float>> positions; // [0] - x , [1] - y , [2] - z , [3] - w
std::vector<std::vector<float>> lightIntensities; // [0] - R , [1] - G , [2] - B , [3] - A
//          values, type(object=0,reflictive=1,tranparent=2).     x,y,z center position, r- radius always > 0.   coefficients of the plane equation when 'd' gets is a non - positive value
std::vector<std::pair<std::vector<float>, int>> spheresAndPlanes; // (x,y,z,r) for sphers ,                         (a,b,c,d) for planes
std::vector<std::vector<float>> colors; // [0] - R, [1] - G, [2] - B, [3] - A represents the shininess value


const int DISPLAY_WIDTH = 800;
const int DISPLAY_HEIGHT = 800;



int main(int argc,char *argv[])
{	
	const float CAMERA_ANGLE = 0.0f;
	const float NEAR = 1.0f;
	const float FAR = 100.0f;

	Game *scn = new Game(CAMERA_ANGLE,(float)DISPLAY_WIDTH/DISPLAY_HEIGHT,NEAR,FAR);
	
	Display display(DISPLAY_WIDTH, DISPLAY_HEIGHT, "OpenGL");
	
	Init(display);
	
	scn->Init();

	display.SetScene(scn);

	// first need to parse the scene file as suggested from tirgul 5
	// step 1 in the assignment in the terminal please enter the desired scene text file to parse
	std::string fileName;

	// Get user input for scene number (0-5)
	int sceneNumber = -1;
	while (sceneNumber < 0 || sceneNumber > 5)
	{
		std::cout << "Enter scene number:\n0 - custom_scene.txt\n1 - scene1.txt\n2 - scene2.txt\n3 - scene3.txt\n4 - scene4.txt\n5 - scene5.txt\nYour Choice: ";
		std::cin >> sceneNumber;
		if (sceneNumber < 0 || sceneNumber > 5)
		{
			std::cout << "Invalid scene number, please try again..." << std::endl;
		}
	}
	
	// open the scene file that the user choosen
	if (sceneNumber == 0)
	{
		fileName = "../assignment2/custom_scene.txt";
	}
	else
	{
		fileName = "../assignment2/scene" + std::to_string(sceneNumber) + ".txt";
	}

	// Open the file for reading
	std::ifstream sceneTxt(fileName);

	// Check if file can be opened
	if (!sceneTxt.is_open()) {
		std::cerr << "Error opening file: " << fileName << std::endl;
		return 1;
	}

	// Read and save the contents of the file to sceneContent and after close the sceneTxt input file
	std::string line;
	while (std::getline(sceneTxt, line))
	{
		std::pair<std::string, std::vector<float>> lineValues = Split(line, ' ');
		
		// check pair first value if it is: e , a , o , t , d , p
		if (lineValues.first == "e") // in case of eye or camera
		{
			camera = lineValues.second;
		}
		else if (lineValues.first == "d") // in case of Light direction
		{
			lightDirctions.push_back(lineValues.second);
		}
		else if (lineValues.first == "a") // in case of ambient lights
		{
			ambientLight = lineValues.second;
		}
		else if (lineValues.first == "p") // for spotlight, the position
		{
			positions.push_back(lineValues.second);
		}
		else if (lineValues.first == "i") // for light intensity
		{
			lightIntensities.push_back(lineValues.second);
		}
		// o / r / t for spheres and planes
		else if (lineValues.first == "o")
		{
			spheresAndPlanes.push_back(std::make_pair(lineValues.second, 0));
		}
		else if (lineValues.first == "r")
		{
			spheresAndPlanes.push_back(std::make_pair(lineValues.second, 1));
		}
		else if (lineValues.first == "t")
		{
			spheresAndPlanes.push_back(std::make_pair(lineValues.second, 2));
		}
		else if (lineValues.first == "c")
		{
			colors.push_back(lineValues.second);
		}
		std::cout << line << std::endl;
	}
	sceneTxt.close();
	
	// end of parsing the scene text file and initializing the variables

	// step 2 initialize the image data and applying the ray casting
	//unsigned char* data = (unsigned char*) malloc(800 * 800 * 4);
	unsigned char* data = new unsigned char[800 * 2 * 800 * 2];
	int maxTimes = 5;


	// bonus implementation
	bool bonus = camera[3] > 0;
	// Get user input for bonus (anti aliasing) run (0-1)
	int anti_aliasing = -1;
	while (anti_aliasing < 0 || anti_aliasing > 1)
	{
		std::cout << "Enter Run Option:\n0 - Without Anti Aliasing\n1 - With Anti Aliasing (Bonus)\nYour Choice: ";
		std::cin >> anti_aliasing;
		if (anti_aliasing < 0 || anti_aliasing > 1)
		{
			std::cout << "Invalid number, please try again..." << std::endl;
		}
	}
	
	if (!(anti_aliasing == 1 && bonus == true))
	{
		for (int y = 0; y < DISPLAY_HEIGHT; y++)
		{
			for (int x = 0; x < DISPLAY_WIDTH; x++)
			{
				// stages as in lecture 3 page 9 and 10
				// 1. shoot rays Shoot Rays from center of projection
				// through pixels
				std::vector<float> ray = ConstructRayThroughPixel(y, x);
				// 2. find intersection just like in lecture 3 page 19, broute force Approach
				std::pair<float[4], ObjectValues> hit = FindIntersection(ray, -1000, camera);
				// 3. get color
				std::vector<float> color = GetColor(hit, ray, maxTimes, camera);
				// update image pixel color values
				data[(800 * y + x) * 4 + 0] = (unsigned char)(color[0] * 255.0);
				data[(800 * y + x) * 4 + 1] = (unsigned char)(color[1] * 255.0);
				data[(800 * y + x) * 4 + 2] = (unsigned char)(color[2] * 255.0);
				data[(800 * y + x) * 4 + 3] = (unsigned char)(0 * 255.0);

			}
		}
	}

	else
	{
		// bonus implementation 
		for (int y = 0; y < DISPLAY_HEIGHT; y++)
		{
			for (int x = 0; x < DISPLAY_WIDTH; x++)
			{
				// constructRayThroughtPixelBonus
				// stages as in lecture 3 page 9 and 10
				// 1. shoot rays Shoot Rays from center of projection
				// through pixels
				std::vector<std::vector<float>> rays = ConstructRayThroughPixelBonus(y, x);


				float colors[4] = { 0,0,0 };
				for (int j=0;j<4;j++)
				{
					// 2. find intersection just like in lecture 3 page 19, broute force Approach
					std::pair<float[4], ObjectValues> hit = FindIntersection(rays[j], -1, camera);
					// 3. get color
					std::vector<float> color = GetColor(hit, rays[j], maxTimes, camera);
					// update image pixel color values
					colors[0] += (color[0] * 255.0);
					colors[1] += (color[1] * 255.0);
					colors[2] += (color[2] * 255.0);

				}

				data[(800 * y + x) * 4 + 0] = (unsigned char)(colors[0] / 4.0);
				data[(800 * y + x) * 4 + 1] = (unsigned char)(colors[1]  / 4.0);
				data[(800 * y + x) * 4 + 2] = (unsigned char)(colors[2] / 4.0);
				data[(800 * y + x) * 4 + 3] = (unsigned char)(0);
			}
		}
	}
	

	scn->AddTexture(800, 800, data);
	scn->SetShapeTex(0, 0); // plane number, texture number
	scn->Draw(1, 0, scn->BACK, true, false);


	scn->Motion();
	display.SwapBuffers();

	while(!display.CloseWindow())
	{
		/*scn->Draw(1,0,scn->BACK,true,false);
		scn->Motion();
		display.SwapBuffers();*/
		display.PollEvents();
	}
	delete scn;
	return 0;
}


std::pair<std::string, std::vector<float>> Split(std::string line, char delimiter) {
	// spliting each line to a pair of <letter, vactor of the remaining floats>
	std::vector<float> values;
	std::string type = "";
	std::string curr = "";
	

	for (int i = 0; i < line.length(); i++) 
	{
		if (line[i] == delimiter)
		{
			if (type == "")
			{
				type = curr;
			}
			else
			{
				values.push_back(stof(curr));
			}
			curr = "";
		}
		else
		{
			curr += line[i];
		}
	}
	if (curr != "")
	{
		values.push_back(stof(curr));
	}

	return std::make_pair(type,values);
}


std::vector<float> ConstructRayThroughPixel(int y, int x)
{
	// see https://computergraphics.stackexchange.com/questions/8479/how-to-calculate-ray
	// also tirgul 4 page 21 and page 22
	/*
	בכדי לחשב את התמונה על המישור שלנו, נבצע את הצעדים הבאים:
	1. נגדיר את גודל המישור שלנו, ואת כמות הפיקסלים שלו באורך וברוחב.
	2. נחשב האורך והרוחב של פיקסל בודד במישור.
	3. נירה קרן אל המרכז של פיקסל ראשוני ממנו אנחנו רוצים להתחיל, ונחשב באיזה צורות הוא פוגע כדי לבחור לו צבע.
	4. נזוז למרכזים של כל שאר הפיקסלים במסך (באמצעות קפיצה באורך והרוחב שחישבנו עבור פיקסל בודד), ונירה מהם קרן כדי לבדוק באילו צורות פגענו ואיזה צבע לבחור להם. 	
	*/
	// P = P_c + ( x - | R_x/2 | ) R V_right - ( y - | R_y/2 | ) R V_up
	// note from assignment pdf: 
	// The screen is located on z=0 plane. The right up corner of the screen is located at (1,1,0) and the
	// left bottom corner of the screen is located at(-1, -1, 0).All in the scene coordinates.
		
	float centerPixelHitPos[3] = { 0, 0, 0 };
	float rayCasted[3] = { 0, 0, 0 }; // ray direction values

	centerPixelHitPos[0] = -1.0 + (1.0 / float(DISPLAY_WIDTH));
	centerPixelHitPos[1] = 1.0 - (1.0 / float(DISPLAY_HEIGHT));

	centerPixelHitPos[0] += 2.0 * x / float(DISPLAY_WIDTH);
	centerPixelHitPos[1] -= 2.0 * y / float(DISPLAY_HEIGHT);


	// noramlizing 
	centerPixelHitPos[0] -= camera[0];
	centerPixelHitPos[1] -= camera[1];
	centerPixelHitPos[2] -= camera[2];

	float sizeOfVector = sqrt(centerPixelHitPos[0] * centerPixelHitPos[0] + centerPixelHitPos[1] * centerPixelHitPos[1] + centerPixelHitPos[2] * centerPixelHitPos[2]);

	rayCasted[0] = centerPixelHitPos[0] / sizeOfVector;
	rayCasted[1] = centerPixelHitPos[1] / sizeOfVector;
	rayCasted[2] = centerPixelHitPos[2] / sizeOfVector;

	std::vector<float> result;
	result.push_back(rayCasted[0]);
	result.push_back(rayCasted[1]);
	result.push_back(rayCasted[2]);
	return result;
}


std::vector<std::vector<float>> ConstructRayThroughPixelBonus(int y, int x)
{
	// see https://computergraphics.stackexchange.com/questions/8479/how-to-calculate-ray
	// also see explanation in the ConstructRayThroughPixel explanation
	// shoting 4 rays: top left, top right, buttom left, buttom right
	std::vector<std::vector<float>> rays;


	float centerPixelHitPos[3] = { 0, 0, 0 };
	float rayCasted[3] = { 0, 0, 0 }; // ray direction values

	centerPixelHitPos[0] = -1.0 + ((1.0 / float(DISPLAY_WIDTH)/2.0));
	centerPixelHitPos[1] = 1.0 - ((1.0 / float(DISPLAY_HEIGHT)/2.0));

	std::vector<float> ray1;
	ray1.push_back(centerPixelHitPos[0] + 2.0 * x / float(DISPLAY_WIDTH));
	ray1.push_back(centerPixelHitPos[1] - 2.0 * y / float(DISPLAY_HEIGHT));
	ray1.push_back(0.0);
	rays.push_back(ray1);

	std::vector<float> ray2;
	ray2.push_back(centerPixelHitPos[0] + 2.0 * x / float(DISPLAY_WIDTH) + 1.0 / float(DISPLAY_WIDTH));
	ray2.push_back(centerPixelHitPos[1] - 2.0 * y / float(DISPLAY_HEIGHT));
	ray2.push_back(0.0);
	rays.push_back(ray2);

	std::vector<float> ray3;
	ray3.push_back(centerPixelHitPos[0] + 2.0 * x / float(DISPLAY_WIDTH));
	ray3.push_back(centerPixelHitPos[1] - 2.0 * y / float(DISPLAY_HEIGHT) + 1.0 / float(DISPLAY_HEIGHT));
	ray3.push_back(0.0);
	rays.push_back(ray3);

	std::vector<float> ray4;
	ray4.push_back(centerPixelHitPos[0] + 2.0 * x / float(DISPLAY_WIDTH) + 1.0 / float(DISPLAY_WIDTH));
	ray4.push_back(centerPixelHitPos[1] - 2.0 * y / float(DISPLAY_HEIGHT) + 1.0 / float(DISPLAY_HEIGHT));
	ray4.push_back(0.0);
	rays.push_back(ray4);
	
	
	// noramlizing
	for (int j=0;j<4;j++)
	{
		rays[j][0] -= camera[0];
		rays[j][1] -= camera[1];
		rays[j][2] -= camera[2];

		float sizeOfVector = sqrt(rays[j][0] * rays[j][0] + rays[j][1] * rays[j][1] + rays[j][2] * rays[j][2]);

		rays[j][0] = rays[j][0] / sizeOfVector;
		rays[j][1] = rays[j][1] / sizeOfVector;
		rays[j][2] = rays[j][2] / sizeOfVector;
	}
	 
	return rays;
}



std::pair<float[4], ObjectValues> FindIntersection(std::vector<float> ray, int objectIndexToPass, std::vector<float> camera)
{
	// page 19 in lecture 3 algorithm
	float min_t = INFINITY;
	bool isMinPrimitiveNull = true;
	float minPrimitive[4] = { 1.0, 1.0, 1.0, 0 };
	float minPrimitiveType = -1; // 0 - o , 1 - r , 2 - t
	float minPrimitiveColors[4] = { 0,0,0,0 };
	
	int index = 0;
	for (int i = 0; i < spheresAndPlanes.size(); i++)
	{
		if (i == objectIndexToPass)
			continue;
		// t = Intersect(ray, primitive );
		float t = 0.0;
		if (spheresAndPlanes[i].first[3] > 0) // then this is sphere
		{
			t = IntersectSphere(spheresAndPlanes[i].first, ray, camera);
		}
		else
		{
			t = IntersectPlane(spheresAndPlanes[i].first, ray, camera);
		}

		if (t < 0)
			continue;
		if (t < min_t)
		{
			isMinPrimitiveNull = false;
			minPrimitive[0] = spheresAndPlanes[i].first[0];
			minPrimitive[1] = spheresAndPlanes[i].first[1];
			minPrimitive[2] = spheresAndPlanes[i].first[2];
			minPrimitive[3] = spheresAndPlanes[i].first[3];

			minPrimitiveType = spheresAndPlanes[i].second;
			index = i;
			minPrimitiveColors[0] = colors[i][0];
			minPrimitiveColors[1] = colors[i][1];
			minPrimitiveColors[2] = colors[i][2];
			minPrimitiveColors[3] = colors[i][3];

			min_t = t;
		}
	}

	std::pair<float[4], ObjectValues> intersectionHit;
	intersectionHit.first[0] = camera[0] + ray[0];
	intersectionHit.first[1] = camera[1] + ray[1];
	intersectionHit.first[2] = camera[2] + ray[2];

	
	(intersectionHit.second).objectValues[0] = minPrimitive[0];
	(intersectionHit.second).objectValues[1] = minPrimitive[1];
	(intersectionHit.second).objectValues[2] = minPrimitive[2];
	(intersectionHit.second).objectValues[3] = minPrimitive[3];

	(intersectionHit.second).objectColors[0] = minPrimitiveColors[0];
	(intersectionHit.second).objectColors[1] = minPrimitiveColors[1];
	(intersectionHit.second).objectColors[2] = minPrimitiveColors[2];
	(intersectionHit.second).objectColors[3] = minPrimitiveColors[3];

	(intersectionHit.second).objectIndex = index;
	(intersectionHit.second).objectType = minPrimitiveType;


	if (!isMinPrimitiveNull)
	{
		intersectionHit.first[0] = camera[0] + ray[0] * min_t;
		intersectionHit.first[1] = camera[1] + ray[1] * min_t;
		intersectionHit.first[2] = camera[2] + ray[2] * min_t;
	}

	return intersectionHit;
}


float IntersectSphere(std::vector<float> sphere, std::vector<float> rayCasted, std::vector<float> camera)
{
	// tirgul 4 pages 4 and 5
	/*
	Ray-Sphere Intersection I
	על מנת לחשב חיתוך של קרן עם כדור, נבצע את החישובים הבאים:
	1. נייצג נקודה 𝑃 על הקרן באמצעות משוואת הקרן, ובעזרת 𝑡.
	2. נציב במשוואת הכדור את משוואת הנקודה 𝑃.
	3. נפתור את המשוואה, וניעזר בנוסחת השורשים.
	4. ניקח את הפתרונות עבורם 𝑡>0 (כיוון שעבור ערכים אלו הנקודות נמצאות עם כיוון הקרן).
	5. נשים לב כי ייתכנו 3 מקרים:
	5.1. לא נקבל אף פיתרון ואז אין חיתוך.
	5.2. נקבל פיתרון יחיד ואז יש חיתוך יחיד.
	5.3. נקבל 2 פיתרונות ואז יש 2 נקודות חיתוך.
	*/
	
	// sphere -> O
	// camera -> P_0
	// rayCasted -> V
	// here we calculate P_0 - O
	float rayCastedHitPointMinusSphereCenter[3] = { 0,0,0 };
	rayCastedHitPointMinusSphereCenter[0] = camera[0] - sphere[0];
	rayCastedHitPointMinusSphereCenter[1] = camera[1] - sphere[1];
	rayCastedHitPointMinusSphereCenter[2] = camera[2] - sphere[2];

	// quadratic equation : at^2 + bt + c = 0
	// a = |V| = 1
	float a = 1.0;
	float b = 0.0;
	float c = 0.0;

	// b = 2V * (P_0 - O)
	for (int i = 0; i < 3; i++)
	{
		b += rayCastedHitPointMinusSphereCenter[i] * rayCasted[i];
	}
	b = b * 2;

	// c = |P_0 - O|^2 - r^2
	for (int i = 0; i < 3; i++)
	{
		c += rayCastedHitPointMinusSphereCenter[i] * rayCastedHitPointMinusSphereCenter[i];
	}
	c -= sphere[3] * sphere[3]; // r^2
	// נוסחאת השורשים delta calculation
	float d = (b * b) - 4 * a * c;
	if (d < 0)
		return -1; // cannot be minus
	
	float x1 = (-b + sqrt(d)) / 2 * a;
	float x2 = (-b - sqrt(d)) / 2 * a;

	// take only t>0
	float min = x1 < x2 ? x1 : x2;
	if (min <= 0.001f) // threashold for transparent spheres
		return x1 < x2 ? x2 : x1;
	return min;
}



float IntersectPlane(std::vector<float> plane, std::vector<float> rayCasted, std::vector<float> camera)
{
	// tirgul 4 pages 10 and 11
	/*
	Ray-Plane Intersection
	על מנת לחשב חיתוך של קרן עם מישור, נבצע את החישובים הבאים:
	1. נייצג נקודה 𝑃 על הקרן באמצעות משוואת הקרן, ובעזרת 𝑡.
	2. נציב במשוואת המישור את משוואת הנקודה 𝑃.
	3. נפתור את המשוואה.
	4. נשים לב כי ייתכנו 2 מקרים:
	4.1. לא נקבל אף פיתרון ואז אין חיתוך.
	4.2. נקבל פיתרון יחיד ואז יש חיתוך יחיד.
	*/
	// Plane -> N
	// camera -> (Q_0 - P_0)
	// rayCasted -> V
	// t = N * (Q_0 - P_0) / N*V
	//                     N.a          x_0          N.b      y_0        N.c          z_0      N.d=distance
	float intersection = plane[0] * camera[0] + plane[1] * camera[1] + plane[2] * camera[2] + plane[3];
	//                 N.a          V.x           N.b          V.y           N.c         V.z
	intersection /= (plane[0] * rayCasted[0] + plane[1] * rayCasted[1] + plane[2] * rayCasted[2]);
	return -1 * intersection; // -1 * intersection because all values of plane in the text file are negative
}


std::vector<float> GetColor(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int times, std::vector<float> camera)
{
	// Get Color of Object / Reflective Object / Transparent Object
	std::vector<float> color;
	if (intersectionHit.second.objectType == 0) // object=0
	{
		std::vector<float> result = ObjectGetColor(intersectionHit, ray, camera);
		color.push_back(result[0] > 1.0 ? 1.0 : result[0]); // clip
		color.push_back(result[1] > 1.0 ? 1.0 : result[1]);
		color.push_back(result[2] > 1.0 ? 1.0 : result[2]);
		return color;
	}

	if (intersectionHit.second.objectType == 1) // reflective=1
	{
		std::vector<float> result = ReflectiveGetColor(intersectionHit, ray, times, camera);
		color.push_back(result[0] > 1.0 ? 1.0 : result[0]); // clip
		color.push_back(result[1] > 1.0 ? 1.0 : result[1]);
		color.push_back(result[2] > 1.0 ? 1.0 : result[2]);
		return color;
	}

	if (intersectionHit.second.objectType == 2) // tranparent=2
	{
		std::vector<float> result = TranparentGetColor(intersectionHit, ray, times, camera);
		color.push_back(result[0] > 1.0 ? 1.0 : result[0]); // clip
		color.push_back(result[1] > 1.0 ? 1.0 : result[1]);
		color.push_back(result[2] > 1.0 ? 1.0 : result[2]);
		return color;
	}
	return color;
}

std::vector<float> DiffuseReflection(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int i, std::vector<float> camera)
{
	// tirgul 4 page 24
	// calculate diffuse reflection
	// cosine law 
	// N * L = |N||L|cos a
	// N^ * L^ = cos a
	// I_D = K_D ( N^ * L^ ) I_L
	/*
	על מנת לחשב את ההחזרה הדיפיוזית, נבצע את החישובים הבאים:
	1. נחשב את המכפלה הסקאלרית של וקטור הנורמל למישור עם הוקטור מהמישור למקור האור.
	2. נכפיל את התוצאה במקדם הדיפיוזי 𝐾_𝐷 ובעוצמת מקור האור 𝐼_𝐿.
	3. את התוצאה נסמן כ- 𝐼_𝐷.
	*/
	bool diffuseHit = true;


	// diffuse
	std::vector<float> diffuse = { 0,0,0 };
	float factor = intersectionHit.second.objectValues[3] >= 0.0 ? 1 : -1;
	// L^ -> normalizing ray direction
	float lightDirection[3] = { factor * lightDirctions[i][0],factor * lightDirctions[i][1],factor * lightDirctions[i][2] };
	// noramlizing vector
	float vecSize = sqrt(lightDirection[0] * lightDirection[0] + lightDirection[1] * lightDirection[1] + lightDirection[2] * lightDirection[2]);
	lightDirection[0] /= vecSize;
	lightDirection[1] /= vecSize;
	lightDirection[2] /= vecSize;

	if (lightDirctions[i][3]) //[3] == 1.0 for spotlight 
	{
		float spolightRay[3] = { intersectionHit.first[0] - positions[i][0],
		intersectionHit.first[1] - positions[i][1],
		intersectionHit.first[2] - positions[i][2] };
		// noramlizing vector
		float vecSize2 = sqrt(spolightRay[0] * spolightRay[0] + spolightRay[1] * spolightRay[1] + spolightRay[2] * spolightRay[2]);
		spolightRay[0] /= vecSize2;
		spolightRay[1] /= vecSize2;
		spolightRay[2] /= vecSize2;
		// light cos value
		float cos = 0.0;
		for (int j = 0; j < 3; j++)
		{
			cos += factor * spolightRay[j] * lightDirection[j];
		}

		if (cos < positions[i][3])
		{
			// spotlight did not hit the object
			diffuse[0] = 0.0;
			diffuse[1] = 0.0;
			diffuse[2] = 0.0;
			diffuseHit = false;
			return diffuse;
		}
		else
		{
			lightDirection[0] = factor * spolightRay[0];
			lightDirection[1] = factor * spolightRay[1];
			lightDirection[2] = factor * spolightRay[2];
		}

	}
	if (diffuseHit)
	{
		// in case of sphere
		float noramlizeObjectVector[3] = { 0,0,0 }; // N^
		// get normal of object
		if (intersectionHit.second.objectValues[3] > 0) // then this is sphere
		{
			// get normal of sphere
			noramlizeObjectVector[0] = intersectionHit.first[0] - intersectionHit.second.objectValues[0];
			noramlizeObjectVector[1] = intersectionHit.first[1] - intersectionHit.second.objectValues[1];
			noramlizeObjectVector[2] = intersectionHit.first[2] - intersectionHit.second.objectValues[2];
		}
		else
		{
			// get normal of plane
			noramlizeObjectVector[0] = intersectionHit.second.objectValues[0];
			noramlizeObjectVector[1] = intersectionHit.second.objectValues[1];
			noramlizeObjectVector[2] = intersectionHit.second.objectValues[2];
		}

		float sizeOfVector = sqrt(noramlizeObjectVector[0] * noramlizeObjectVector[0] + noramlizeObjectVector[1] * noramlizeObjectVector[1] + noramlizeObjectVector[2] * noramlizeObjectVector[2]);
		noramlizeObjectVector[0] = noramlizeObjectVector[0] / sizeOfVector;
		noramlizeObjectVector[1] = noramlizeObjectVector[1] / sizeOfVector;
		noramlizeObjectVector[2] = noramlizeObjectVector[2] / sizeOfVector;

		// N^ * L^ = cos a
		float cos = 0.0;
		for (int j = 0; j < 3; j++)
		{
			cos += -1 * lightDirection[j] * noramlizeObjectVector[j];
		}

		// I_D = (N ^ *L^) K_D I_L
		std::vector<float> objectColor = GetObjectColor(intersectionHit.second, intersectionHit.first);

		diffuse[0] = cos * objectColor[0] * lightIntensities[i][0];
		diffuse[1] = cos * objectColor[1] * lightIntensities[i][1];
		diffuse[2] = cos * objectColor[2] * lightIntensities[i][2];

		// max of two vectors
		diffuse[0] = diffuse[0] > 0 ? diffuse[0] : 0;
		diffuse[1] = diffuse[1] > 0 ? diffuse[1] : 0;
		diffuse[2] = diffuse[2] > 0 ? diffuse[2] : 0;

	}
	return diffuse;
}

std::vector<float> SpecularReflection(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int i, std::vector<float> camera)
{
	// calculate Specular Reflection
	// tirgul 4 pages 26 + 27
	// phong model
	// cos(a)^n
	// V^ * R^ = max( 0 , V^ * R^ )
	// I_S = K_S ( V^ * R^ )^n I_L
	/*
	על מנת לחשב את החזר הראי, נבצע את החישובים הבאים:
	1. נחשבת את המכפלה הסקאלרית בין הוקטור מהמישור למקור האור עם הוקטור מהמישור אל הצופה (או במקרה שלנו אל המצלמה), וניקח את המקסימום בין ערך זה ל- 0.
	2. את התוצאה נכפיל בחזקת מקדם הברק 𝑛.
	3. נכפיל את התוצאה במקדם הראי 𝐾_𝑆 ובעוצמת מקור האור 𝐼_𝐿.
	4. את התוצאה נסמן כ- 𝐼_𝑆.
	*/

	// specular
	std::vector<float> specular = { 0,0,0 };
	bool specularHit = true;
	// normalizing ray direction
	float lightDirection[3] = { lightDirctions[i][0], lightDirctions[i][1], lightDirctions[i][2] };
	// noramlizing vector
	float vecSize = sqrt(lightDirection[0] * lightDirection[0] + lightDirection[1] * lightDirection[1] + lightDirection[2] * lightDirection[2]);
	lightDirection[0] /= vecSize;
	lightDirection[1] /= vecSize;
	lightDirection[2] /= vecSize;

	// same calculation as in diffuse
	if (lightDirctions[i][3]) //[3] == 1.0 for spotlight 
	{
		float spolightRay[3] = { intersectionHit.first[0] - positions[i][0],
		intersectionHit.first[1] - positions[i][1],
		intersectionHit.first[2] - positions[i][2] };
		// noramlizing vector
		float vecSize2 = sqrt(spolightRay[0] * spolightRay[0] + spolightRay[1] * spolightRay[1] + spolightRay[2] * spolightRay[2]);
		spolightRay[0] /= vecSize2;
		spolightRay[1] /= vecSize2;
		spolightRay[2] /= vecSize2;
		// light cos value
		float cos = 0.0;
		for (int j = 0; j < 3; j++)
		{
			cos += spolightRay[j] * lightDirection[j];
		}

		if (cos < positions[i][3])
		{
			// spotlight did not hit the object
			specular[0] = 0.0;
			specular[1] = 0.0;
			specular[2] = 0.0;
			specularHit = false;
			return specular;
		}
		else
		{
			lightDirection[0] = spolightRay[0];
			lightDirection[1] = spolightRay[1];
			lightDirection[2] = spolightRay[2];
		}

	}
	if (specularHit)
	{
		// in case of sphere
		float noramlizeObjectVector[3] = { 0,0,0 };
		// get normal of object
		if (intersectionHit.second.objectValues[3] > 0) // then this is sphere
		{
			// get normal of sphere
			noramlizeObjectVector[0] = intersectionHit.first[0] - intersectionHit.second.objectValues[0];
			noramlizeObjectVector[1] = intersectionHit.first[1] - intersectionHit.second.objectValues[1];
			noramlizeObjectVector[2] = intersectionHit.first[2] - intersectionHit.second.objectValues[2];
		}
		else
		{
			// get normal of plane
			noramlizeObjectVector[0] = intersectionHit.second.objectValues[0];
			noramlizeObjectVector[1] = intersectionHit.second.objectValues[1];
			noramlizeObjectVector[2] = intersectionHit.second.objectValues[2];
		}

		float sizeOfVector = sqrt(noramlizeObjectVector[0] * noramlizeObjectVector[0] + noramlizeObjectVector[1] * noramlizeObjectVector[1] + noramlizeObjectVector[2] * noramlizeObjectVector[2]);
		noramlizeObjectVector[0] = noramlizeObjectVector[0] / sizeOfVector;
		noramlizeObjectVector[1] = noramlizeObjectVector[1] / sizeOfVector;
		noramlizeObjectVector[2] = noramlizeObjectVector[2] / sizeOfVector;

		
		//
		float normalToPlane[3] = { 0,0,0 };
		normalToPlane[0] = camera[0] - intersectionHit.first[0];
		normalToPlane[1] = camera[1] - intersectionHit.first[1];
		normalToPlane[2] = camera[2] - intersectionHit.first[2];
		sizeOfVector = sqrt(normalToPlane[0] * normalToPlane[0] + normalToPlane[1] * normalToPlane[1] + normalToPlane[2] * normalToPlane[2]);
		normalToPlane[0] = normalToPlane[0] / sizeOfVector;
		normalToPlane[1] = normalToPlane[1] / sizeOfVector;
		normalToPlane[2] = normalToPlane[2] / sizeOfVector;

		float dotProduct = 0;
		for (int j = 0; j < 3; j++)
		{
			dotProduct += noramlizeObjectVector[j] * lightDirection[j];
		}

		float planeToLightSource[3] = { 0, 0, 0 };
		for (int j = 0; j < 3; j++)
		{
			planeToLightSource[j] = lightDirection[j] - 2.0f * noramlizeObjectVector[j] * dotProduct;
		}


		// V^ * R^ = max( 0 , V^ * R^ )
		float cos = 0.0;
		for (int j = 0; j < 3; j++)
		{
			cos += planeToLightSource[j] * normalToPlane[j];
		}
		cos = cos > 0 ? cos : 0;

		// cos(a)^n                            shiness value
		cos = pow(cos, intersectionHit.second.objectColors[3]);
		
		//       K_S
		float specular_value = 0.7; // from assignment pdf: The specular value of an object is always (0.7,0.7,0.7)
		// I_S      =      K_S     ( V^ * R^ )^n        I_L
		specular[0] = specular_value * cos * lightIntensities[i][0];
		specular[1] = specular_value * cos * lightIntensities[i][1];
		specular[2] = specular_value * cos * lightIntensities[i][2];

		// max of two vectors
		specular[0] = specular[0] > 0 ? specular[0] : 0;
		specular[1] = specular[1] > 0 ? specular[1] : 0;
		specular[2] = specular[2] > 0 ? specular[2] : 0;

		return specular;
	}





}


std::vector<float> Shadows(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int i, std::vector<float> camera)
{
	// calculate Shadow term
	// tirgul 4 pages 31 and 32

	/*
	מקדם הצל אומר לנו איזה מקורות אור חסומים עבור הצורה שלנו (כלומר קרן האור לא מצליחה לפגוע בהם).
	לכן בכדי לחשב את הצל על הצורה שלנו, כל סכום של ההחזרה הדיפיוזית והחזר הראי, נכפיל במקדם הצל 𝑆_𝑖.
	מקדם הצל יהיה שווה ל- 0 אם מקור האור חסום עבור הצורה שלנו.
	אחרת, מקדם הצל יהיה שווה ל- 1.
	*/

	std::vector<float> shadows = { 0, 0, 0 };
	bool shadowHit = true;
	float min_t = INFINITY;
	// normalizing ray direction
	float lightDirection[3] = { lightDirctions[i][0], lightDirctions[i][1], lightDirctions[i][2] };
	// noramlizing vector
	float vecSize = sqrt(lightDirection[0] * lightDirection[0] + lightDirection[1] * lightDirection[1] + lightDirection[2] * lightDirection[2]);
	lightDirection[0] /= vecSize;
	lightDirection[1] /= vecSize;
	lightDirection[2] /= vecSize;

	// same calculation as in diffuse
	if (lightDirctions[i][3]) //[3] == 1.0 for spotlight 
	{
		float spolightRay[3] = { intersectionHit.first[0] - positions[i][0],
		intersectionHit.first[1] - positions[i][1],
		intersectionHit.first[2] - positions[i][2] };
		// noramlizing vector
		float vecSize2 = sqrt(spolightRay[0] * spolightRay[0] + spolightRay[1] * spolightRay[1] + spolightRay[2] * spolightRay[2]);
		spolightRay[0] /= vecSize2;
		spolightRay[1] /= vecSize2;
		spolightRay[2] /= vecSize2;
		// light cos value
		float cos = 0.0;
		for (int j = 0; j < 3; j++)
		{
			cos += spolightRay[j] * lightDirection[j];
		}

		if (cos < positions[i][3])
		{
			// spotlight did not hit the object
			//shadows[0] = 0.0;
			//shadows[1] = 0.0;
			//shadows[2] = 0.0;
			shadowHit = false;
			return shadows;
		}
		else
		{
			lightDirection[0] = spolightRay[0];
			lightDirection[1] = spolightRay[1];
			lightDirection[2] = spolightRay[2];

			// calculate min_t and update its value
			float dotOfLightDirectionAndPosition = 0;
			for (int j = 0; j < 3; j++)
			{
				dotOfLightDirectionAndPosition += - lightDirection[j] * positions[i][j];
			}
			dotOfLightDirectionAndPosition = abs(dotOfLightDirectionAndPosition);

			float dotOfHitAndPosition = 0;
			for (int j = 0; j < 3; j++)
			{
				dotOfHitAndPosition += intersectionHit.first[j] * positions[i][j];
			}
			dotOfHitAndPosition *= -1;

			min_t = dotOfHitAndPosition / dotOfLightDirectionAndPosition;
		}

	}
	if (shadowHit)
	{
		// loop through other objects and check
		for (int y = 0; y < spheresAndPlanes.size(); y++)
		{
			if (y == intersectionHit.second.objectIndex)
				continue;


			std::vector<float> ray1;
			ray1.push_back(-lightDirection[0]);
			ray1.push_back(-lightDirection[1]);
			ray1.push_back(-lightDirection[2]);

			// t = Intersect(ray, primitive );
			std::vector<float> hitVector;
			hitVector.push_back(intersectionHit.first[0]);
			hitVector.push_back(intersectionHit.first[1]);
			hitVector.push_back(intersectionHit.first[2]);
			hitVector.push_back(intersectionHit.first[3]);

			float t = 0.0;
			if (spheresAndPlanes[y].first[3] > 0) // then this is sphere
			{
				t = IntersectSphere(spheresAndPlanes[y].first, ray1, hitVector);
			}
			else
			{
				t = IntersectPlane(spheresAndPlanes[y].first, ray1, hitVector);
			}

			if ((t > 0) && (t < min_t))
				return shadows;
		}

	}

	shadows[0] = 1.0;
	shadows[1] = 1.0;
	shadows[2] = 1.0;
	return shadows;
}

std::vector<float> ObjectGetColor(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, std::vector<float> camera)
{
	// Calculates Diffuse Reflection, Specular Reflection And Shadow Term and update the color
	std::vector<float> color = GetObjectColor(intersectionHit.second, intersectionHit.first);
	color[0] *= ambientLight[0];
	color[1] *= ambientLight[1];
	color[2] *= ambientLight[2];

	for (int i = 0; i < lightDirctions.size(); i++)
	{
		std::vector<float> diffuse = DiffuseReflection(intersectionHit, ray, i, camera);
		std::vector<float> specular = SpecularReflection(intersectionHit, ray, i, camera);
		std::vector<float> shadows = Shadows(intersectionHit, ray, i, camera);

		color[0] += shadows[0] * (diffuse[0] + specular[0]);
		color[1] += shadows[1] * (diffuse[1] + specular[1]);
		color[2] += shadows[2] * (diffuse[2] + specular[2]);
	}

	return color;
}


std::vector<float> GetObjectColor(ObjectValues obj, float* hit)
{
	// if it is sphere simply return the color
	// if it is plane calculate color same as tirgul 5 page 4 
	std::vector<float> color;
	if (obj.objectValues[3] > 0.0)
	{
		// then this is sphere
		color.push_back(obj.objectColors[0]);
		color.push_back(obj.objectColors[1]);
		color.push_back(obj.objectColors[2]);
		return color;
	}
	else { // then this is Plane
		// from tirgul 5 page 4 same code with few adjustments
		float scale_parameter = 0.5f;
		float chessboard = 0;

		if (hit[0] < 0)
		{
			chessboard += floor((0.5 - hit[0]) / scale_parameter);
		}
		else {
			chessboard += floor(hit[0] / scale_parameter);
		}

		if (hit[1] < 0)
		{
			chessboard += floor((0.5 - hit[1]) / scale_parameter);
		}
		else {
			chessboard += floor(hit[1] / scale_parameter);
		}

		chessboard = (chessboard * 0.5) - int(chessboard * 0.5);
		chessboard *= 2;
		if (chessboard > 0.5) {
			color.push_back(0.5 * obj.objectColors[0]);
			color.push_back(0.5 * obj.objectColors[1]);
			color.push_back(0.5 * obj.objectColors[2]);
			return color;
		}
		color.push_back(obj.objectColors[0]);
		color.push_back(obj.objectColors[1]);
		color.push_back(obj.objectColors[2]);
		return color;

	}
}


std::vector<float> ReflectiveGetColor(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int times, std::vector<float> camera)
{
	/*
	Trace secondary ray in direction of mirror reflection
	Evaluate radiance along secondary ray and include it into illumination model
	במידה ויש לנו חומר מראתי, נחשב את הצבע שלו על ידי כך שנוסיף לסכום את החישוב הרקורסיבי של פגיעת קרן האור בגופים נוספים לאחר שהיא נשברה והמשיכה מחוץ לצורה שלנו.

	*/
	std::vector<float> color = { 0,0,0 };
	if (times == 0)
		return color;

	// in case of sphere
	float noramlizeObjectVector[3] = { 0,0,0 };
	// get normal of object
	if (intersectionHit.second.objectValues[3] > 0) // then this is sphere
	{
		// get normal of sphere
		noramlizeObjectVector[0] = intersectionHit.first[0] - intersectionHit.second.objectValues[0];
		noramlizeObjectVector[1] = intersectionHit.first[1] - intersectionHit.second.objectValues[1];
		noramlizeObjectVector[2] = intersectionHit.first[2] - intersectionHit.second.objectValues[2];
	}
	else
	{
		// get normal of plane
		noramlizeObjectVector[0] = intersectionHit.second.objectValues[0];
		noramlizeObjectVector[1] = intersectionHit.second.objectValues[1];
		noramlizeObjectVector[2] = intersectionHit.second.objectValues[2];
	}

	float sizeOfVector = sqrt(noramlizeObjectVector[0] * noramlizeObjectVector[0] + noramlizeObjectVector[1] * noramlizeObjectVector[1] + noramlizeObjectVector[2] * noramlizeObjectVector[2]);
	noramlizeObjectVector[0] = noramlizeObjectVector[0] / sizeOfVector;
	noramlizeObjectVector[1] = noramlizeObjectVector[1] / sizeOfVector;
	noramlizeObjectVector[2] = noramlizeObjectVector[2] / sizeOfVector;


	float dotProduct = 0;
	for (int j = 0; j < 3; j++)
	{
		dotProduct += noramlizeObjectVector[j] * ray[j];
	}


	std::vector<float> reflectionRay;
	reflectionRay.push_back(ray[0] - 2.0f * noramlizeObjectVector[0] * dotProduct);
	reflectionRay.push_back(ray[1] - 2.0f * noramlizeObjectVector[1] * dotProduct);
	reflectionRay.push_back(ray[2] - 2.0f * noramlizeObjectVector[2] * dotProduct);

	// note: the new eye/ camera is the previous hit point and the new ray is the reflected ray
	std::vector<float> newCamera;
	newCamera.push_back(intersectionHit.first[0]);
	newCamera.push_back(intersectionHit.first[1]);
	newCamera.push_back(intersectionHit.first[2]);
	newCamera.push_back(intersectionHit.first[3]);
	std::pair<float[4], ObjectValues> hit = FindIntersection(reflectionRay, intersectionHit.second.objectIndex, newCamera);

	if (hit.second.objectType == -1) // == Space //delete this comment later
	{
		std::vector<float> resulted = { 0,0,0 };
		return resulted; // zeros
	}

	std::vector<float> result = GetColor(hit, reflectionRay, times - 1, newCamera);
	return result;
}




std::vector<float> TranparentGetColor(std::pair<float[4], ObjectValues> intersectionHit, std::vector<float> ray, int times, std::vector<float> camera)
{
	/*
	במידה ויש לנו חומר שקוף, נחשב את הצבע שלו על ידי כך שנוסיף לסכום את החישוב הרקורסיבי של פגיעת קרן האור בגופים נוספים לאחר שהיא נשברה ועברה דרך הצורה שלנו.
	*/

	std::vector<float> color = { 0,0,0 };
	if (times == 0)
		return color;

	if (intersectionHit.second.objectValues[3] < 0.0) // then this is plane
	{
		std::vector<float> transparencyRay;
		transparencyRay.push_back(ray[0]);
		transparencyRay.push_back(ray[1]);
		transparencyRay.push_back(ray[2]);
		// note: the new eye/ camera is the previous hit point and the new ray is the reflected ray
		std::vector<float> newCamera;
		newCamera.push_back(intersectionHit.first[0]);
		newCamera.push_back(intersectionHit.first[1]);
		newCamera.push_back(intersectionHit.first[2]);
		newCamera.push_back(intersectionHit.first[3]);
		std::pair<float[4], ObjectValues> hit = FindIntersection(transparencyRay, intersectionHit.second.objectIndex, newCamera);

		if (hit.second.objectType == -1)
		{
			std::vector<float> zeros = { 0,0,0 };
			return zeros;
		}
		return GetColor(hit, transparencyRay, times - 1, newCamera);
	}
	else
	{
		// then this is sphere
		// https://personal.math.ubc.ca/~cass/courses/m309-01a/chu/Fundamentals/snell.htm
		// https://www.physicsclassroom.com/class/refrn/Lesson-2/Snell-s-Law
		// from assignment pdf: Use Snell's law with refractive index of 1.5 to determine the direction of the light rays intersect
		// with the sphere(air have refractive index of 1)
		float refractive = (1.0f / 1.5f);
		// snell's law
		// n_r * sin a_r = n_i * sin a_i
		// T = ( ( n_i / n_r ) * cos a_i - cos a_r ) * N - ( n_i / n_r ) * L )

		// in case of sphere
		float noramlizeObjectVector[3] = { 0,0,0 };
		// get normal of object
		if (intersectionHit.second.objectValues[3] > 0) // then this is sphere
		{
			// get normal of sphere
			noramlizeObjectVector[0] = intersectionHit.first[0] - intersectionHit.second.objectValues[0];
			noramlizeObjectVector[1] = intersectionHit.first[1] - intersectionHit.second.objectValues[1];
			noramlizeObjectVector[2] = intersectionHit.first[2] - intersectionHit.second.objectValues[2];
		}
		else
		{
			// get normal of plane
			noramlizeObjectVector[0] = intersectionHit.second.objectValues[0];
			noramlizeObjectVector[1] = intersectionHit.second.objectValues[1];
			noramlizeObjectVector[2] = intersectionHit.second.objectValues[2];
		}

		float sizeOfVector = sqrt(noramlizeObjectVector[0] * noramlizeObjectVector[0] + noramlizeObjectVector[1] * noramlizeObjectVector[1] + noramlizeObjectVector[2] * noramlizeObjectVector[2]);
		noramlizeObjectVector[0] = noramlizeObjectVector[0] / sizeOfVector;
		noramlizeObjectVector[1] = noramlizeObjectVector[1] / sizeOfVector;
		noramlizeObjectVector[2] = noramlizeObjectVector[2] / sizeOfVector;

		
		float cos_i = 0;
		for (int j = 0; j < 3; j++)
		{
			cos_i += noramlizeObjectVector[j] * -1 * ray[j];
		}
		float a_i = acos(cos_i) * (180.0f / M_PI);
		float sin_i = sin(a_i * (M_PI / 180.0f));

		float sin_r = sin_i * refractive;
		float a_r = asin(sin_r) * (180.0f / M_PI);
		float cos_r = cos(a_r * (M_PI / 180.0f));

		std::vector<float> transparentRay;
		transparentRay.push_back((refractive * cos_i - cos_r) * noramlizeObjectVector[0] - refractive * -ray[0]);
		transparentRay.push_back((refractive * cos_i - cos_r) * noramlizeObjectVector[1] - refractive * -ray[1]);
		transparentRay.push_back((refractive * cos_i - cos_r) * noramlizeObjectVector[2] - refractive * -ray[2]);
		// normalizing
		sizeOfVector = sqrt(transparentRay[0] * transparentRay[0] + transparentRay[1] * transparentRay[1] + transparentRay[2] * transparentRay[2]);
		transparentRay[0] = transparentRay[0] / sizeOfVector;
		transparentRay[1] = transparentRay[1] / sizeOfVector;
		transparentRay[2] = transparentRay[2] / sizeOfVector;

		// note: the new eye/ camera is the previous hit point and the new ray is the reflected ray
		std::vector<float> newCamera;
		newCamera.push_back(intersectionHit.first[0]);
		newCamera.push_back(intersectionHit.first[1]);
		newCamera.push_back(intersectionHit.first[2]);
		newCamera.push_back(intersectionHit.first[3]);


		std::pair<float[4], ObjectValues> hit = FindIntersection(transparentRay, -1, newCamera);
		if (hit.second.objectIndex != intersectionHit.second.objectIndex)
		{
			return GetColor(hit, transparentRay, times - 1, newCamera);
		}
		else
		{
			// t = Intersect(ray, primitive );
			float t = 0.0;
			if (intersectionHit.second.objectValues[3] > 0) // then this is sphere
			{
				t = IntersectSphere(spheresAndPlanes[intersectionHit.second.objectIndex].first, transparentRay, newCamera);
			}
			else
			{
				t = IntersectPlane(spheresAndPlanes[intersectionHit.second.objectIndex].first, transparentRay, newCamera);
			}

			std::vector<float> secondHitCamera;
			secondHitCamera.push_back(newCamera[0] + transparentRay[0] * t);
			secondHitCamera.push_back(newCamera[1] + transparentRay[1] * t);
			secondHitCamera.push_back(newCamera[2] + transparentRay[2] * t);


			// normalizing secondHitCamera vector
			// in case of sphere
			float noramlizeObjectVector2[3] = { 0,0,0 };
			// get normal of object
			if (intersectionHit.second.objectValues[3] > 0) // then this is sphere
			{
				// get normal of sphere
				noramlizeObjectVector2[0] = secondHitCamera[0] - intersectionHit.second.objectValues[0];
				noramlizeObjectVector2[1] = secondHitCamera[1] - intersectionHit.second.objectValues[1];
				noramlizeObjectVector2[2] = secondHitCamera[2] - intersectionHit.second.objectValues[2];
			}
			else
			{
				// get normal of plane
				noramlizeObjectVector2[0] = intersectionHit.second.objectValues[0];
				noramlizeObjectVector2[1] = intersectionHit.second.objectValues[1];
				noramlizeObjectVector2[2] = intersectionHit.second.objectValues[2];
			}

			float sizeOfVector = sqrt(noramlizeObjectVector2[0] * noramlizeObjectVector2[0] + noramlizeObjectVector2[1] * noramlizeObjectVector2[1] + noramlizeObjectVector2[2] * noramlizeObjectVector2[2]);
			noramlizeObjectVector2[0] = noramlizeObjectVector2[0] / sizeOfVector;
			noramlizeObjectVector2[1] = noramlizeObjectVector2[1] / sizeOfVector;
			noramlizeObjectVector2[2] = noramlizeObjectVector2[2] / sizeOfVector;

			refractive = 1 / refractive; //

			cos_i = 0;
			for (int j = 0; j < 3; j++)
			{
				cos_i += noramlizeObjectVector2[j] * transparentRay[j];
			}
			a_i = acos(cos_i) * (180.0f / M_PI);
			sin_i = sin(a_i * (M_PI / 180.0f));

			sin_r = sin_i * refractive;
			a_r = asin(sin_r) * (180.0f / M_PI);
			cos_r = cos(a_r * (M_PI / 180.0f));


			//
			std::vector<float> transparentRay2;
			transparentRay2.push_back((refractive * cos_i - cos_r) * -noramlizeObjectVector[0] - refractive * -transparentRay[0]);
			transparentRay2.push_back((refractive * cos_i - cos_r) * -noramlizeObjectVector[1] - refractive * -transparentRay[1]);
			transparentRay2.push_back((refractive * cos_i - cos_r) * -noramlizeObjectVector[2] - refractive * -transparentRay[2]);
			// normalizing
			sizeOfVector = sqrt(transparentRay2[0] * transparentRay2[0] + transparentRay2[1] * transparentRay2[1] + transparentRay2[2] * transparentRay2[2]);
			transparentRay2[0] = transparentRay2[0] / sizeOfVector;
			transparentRay2[1] = transparentRay2[1] / sizeOfVector;
			transparentRay2[2] = transparentRay2[2] / sizeOfVector;

			// note: the new eye/ camera is the previous hit point and the new ray is the transparent ray
			/*std::vector<float> newCamera2;
			newCamera2.push_back(intersectionHit.first[0]);
			newCamera2.push_back(intersectionHit.first[1]);
			newCamera2.push_back(intersectionHit.first[2]);
			newCamera2.push_back(intersectionHit.first[3]);*/

			secondHitCamera.push_back(0);
			std::pair<float[4], ObjectValues> hit = FindIntersection(transparentRay2, intersectionHit.second.objectIndex, secondHitCamera);

			if (hit.second.objectType == -1)
			{
				std::vector<float> zeros = { 0,0,0 };
				return zeros;
			}
			return GetColor(hit, transparentRay2, times - 1, newCamera);
		}
	}
}