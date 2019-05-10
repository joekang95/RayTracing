/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Joe Yu-Ho Chang
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
// Include to Enhance Code
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
// Basic Includes
#include <vector>
#include <iostream>
#include <math.h> 
#include <algorithm>   

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0
#define PI 3.14159265
#define RADIAN (fov / 2.0) * (PI / 180.0)

using namespace std;

unsigned char buffer[HEIGHT][WIDTH][3];
bool antialias = false, softshadow = false, animation = false, motionblur = false;
bool complete = false, nextstep = false;
glm::dvec3 CAMERA = glm::dvec3(0, 0, 0);
static int MAX_RAY_DEPTH = 0, NUM_RAND_LIGHTS = 1;

struct Vertex {
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double normal[3];
	double shininess;
};

struct Triangle {
	Vertex v[3];
};

struct Sphere {
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double shininess;
	double radius;
};

glm::dvec3 sphereCenter(Sphere &sphere, double time);

struct Light {
	double position[3];
	double color[3];
};

class Color {
public:
	double r;
	double g;
	double b;

	Color(double R, double G, double B) : r(R), g(G), b(B) {}

	// Addition Operator to Clamp Values
	Color& operator += (const Color &c) {
		r += c.r;
		if (r > 1.0) { r = 1.0; }
		if (r < 0.0) { r = 0.0; }

		g += c.g;
		if (g > 1.0) { g = 1.0; }
		if (g < 0.0) { g = 0.0; }

		b += c.b;
		if (b > 1.0) { b = 1.0; }
		if (b < 0.0) { b = 0.0; }
		return *this;
	}
	Color& operator / (const double &d) {
		r /= d;
		g /= d;
		b /= d;
		return *this;
	}
	Color& operator - (const  Color &c) {
		r -= c.r;
		if (r > 1.0) { r = 1.0; }
		if (r < 0.0) { r = 0.0; }

		g -= c.g;
		if (g > 1.0) { g = 1.0; }
		if (g < 0.0) { g = 0.0; }

		b -= c.b;
		if (b > 1.0) { b = 1.0; }
		if (b < 0.0) { b = 0.0; }
		return *this;
	}
};

class Ray {
public:
	glm::dvec3 origin;
	glm::dvec3 direction;
	float time;

	//Ray(glm::dvec3 origin, glm::dvec3 direction) : origin(origin), direction(direction) {}
	Ray(glm::dvec3 origin, glm::dvec3 direction, float time) : origin(origin), direction(direction), time(time) {}

	bool intersectSphere(Sphere &sphere, glm::dvec3 &intersection) {

		glm::dvec3 center = sphereCenter(sphere, time);
		glm::dvec3 distance = origin - center;
		double b = 2 * glm::dot(direction, distance);
		double c = pow(glm::length(distance), 2) - pow(sphere.radius, 2);
		double discriminant = pow(b, 2) - 4 * c;

		double t0 = (-b + sqrt(discriminant)) / 2;
		double t1 = (-b - sqrt(discriminant)) / 2;

		if (discriminant < 0) { return false; }
		if (t0 < 0 && t1 < 0) { return false; }
		if (t0 >= 0 || t1 >= 0) {
			double min;
			if (t0 >= 0 && t1 >= 0) { min = std::min(t0, t1); }
			else if (t0 >= 0) { min = t0; }
			else if (t1 >= 0) { min = t1; }
			intersection = origin + (direction * min);
			return true;
		}
		return false;
	}

	bool intersectTriangle(Triangle &triangle, glm::dvec3 &intersection) {

		glm::dvec3 VERTEX_A = glm::dvec3(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
		glm::dvec3 VERTEX_B = glm::dvec3(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
		glm::dvec3 VERTEX_C = glm::dvec3(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);

		glm::dvec3 NORMAL = glm::normalize(glm::cross(VERTEX_B - VERTEX_A, VERTEX_C - VERTEX_A));
		// -(n.p0) + (d coefficient)
		double NUMERATOR = glm::dot(NORMAL, VERTEX_A - origin);
		// n.d (ray direction)
		double DENOMINATOR = glm::dot(NORMAL, direction);
		// DBL_EPSILON or 1e-16
		if (abs(DENOMINATOR) < DBL_EPSILON) { return false;	}

		double t = NUMERATOR / DENOMINATOR;
		// If t <= 0, the Intersection is Behind Ray Origin
		if (t <= 0) { return false; }

		intersection = origin + (direction * t);

		// Barycentric Coordinates
		// |u x v| = |u x (va + vb)| = |u x va + 0| = |u||va| = |u||v|sin(theta) = Area
		// Negative Area = Outside of Triangle
		glm::dvec3 EDGE_1 = VERTEX_B - VERTEX_A;
		glm::dvec3 EDGE_2 = intersection - VERTEX_A;
		if (glm::dot(glm::cross(EDGE_1, EDGE_2), NORMAL) < 0) { return false; }

		EDGE_1 = VERTEX_C - VERTEX_B;
		EDGE_2 = intersection - VERTEX_B;
		if (glm::dot(glm::cross(EDGE_1, EDGE_2), NORMAL) < 0) {	return false; }

		EDGE_1 = VERTEX_A - VERTEX_C;
		EDGE_2 = intersection - VERTEX_C;
		if (glm::dot(glm::cross(EDGE_1, EDGE_2), NORMAL) < 0) {	return false; }

		return true;
	}
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

vector<Light> LIGHT_SOURCE;
void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);
Color traceRay(Ray &ray, int depth);

// Calculate Sphere Center
glm::dvec3 sphereCenter(Sphere &sphere, double time) {
	glm::dvec3 center = glm::dvec3(sphere.position[0], sphere.position[1], sphere.position[2]);
	if (motionblur) {
		// If Motion Blur is Enabled, Return Calculated Center
		return center + glm::dvec3(0.0, time / 20.0 * 0.25, 0.0);
	}
	return center;
}

// Sphere Phong Shading
Color spherePhong(Sphere &sphere, Light &light, glm::dvec3 &intersection, int depth, double time) {

	// Normal of a Point on a Sphere is the Vector Pointing to that Point from Center of the Sphere
	/********************************************************************  
		 /-----\
		|   .-->|(intersection)--> (normal)               <---- . (light)
		 \-----/
	**********************************************************************/
	glm::dvec3 NORMAL = glm::normalize(intersection - sphereCenter(sphere, time));

	// Setup Light Ray
	glm::dvec3 LIGHT_ORIGIN = glm::dvec3(light.position[0], light.position[1], light.position[2]);
	glm::dvec3 LIGHT_DIRECTION = glm::normalize(LIGHT_ORIGIN - intersection);

	// Magnitude of a Light at Intersection is its Projection onto the Normal
	double LIGHT = glm::dot(LIGHT_DIRECTION, NORMAL);
	if (LIGHT > 1.0) { LIGHT = 1.0; }
	if (LIGHT < 0.0) { LIGHT = 0.0; }

	// r = 2(l . n) n ¡V l
	// Reflect to Camera/Origin
	glm::dvec3 REFLECT = glm::normalize(2.0 * LIGHT * NORMAL - LIGHT_DIRECTION);
	glm::dvec3 DIRECTION = glm::normalize(-intersection);
	double REFLECTION = glm::dot(REFLECT, DIRECTION);
	if (REFLECTION > 1.0) { REFLECTION = 1.0; }
	if (REFLECTION < 0.0) { REFLECTION = 0.0; }

	// Initialize r, g, b, Diffuse Color, Specular Color, and Shininess
	double r = 0, g = 0, b = 0;
	Color REFLECTION_COLOR = Color(0, 0, 0);
	Color DIFFUSE(sphere.color_diffuse[0], sphere.color_diffuse[1], sphere.color_diffuse[2]);
	Color SPECULAR(sphere.color_specular[0], sphere.color_specular[1], sphere.color_specular[2]);
	double SHININESS = sphere.shininess;

	// Local Phong Shading, I = L(kd * (l . n) + ks (r . v)^alpha)
	r = light.color[0] * (DIFFUSE.r * LIGHT + (SPECULAR.r * pow(REFLECTION, SHININESS)));
	g = light.color[1] * (DIFFUSE.g * LIGHT + (SPECULAR.g * pow(REFLECTION, SHININESS)));
	b = light.color[2] * (DIFFUSE.b * LIGHT + (SPECULAR.b * pow(REFLECTION, SHININESS)));

	// Recursive Ray Tracing on Sphere
	if (depth < MAX_RAY_DEPTH) {
		// Add Offset to Intersection to Avoid Z-buffer Fighting
		Ray REFLECTION_RAY = Ray(intersection + 0.000001 * NORMAL, REFLECT, time);
		// Trace Reflection Ray
		Color RETURN_COLOR = traceRay(REFLECTION_RAY, depth + 1);
		REFLECTION_COLOR += RETURN_COLOR;
		Color WEIGHT = Color(SPECULAR.r, SPECULAR.g, SPECULAR.b);
		//(1 - ks) * localPhongColor + ks * colorOfReflectedRay
		Color FINAL_COLOR = Color(
			(1 - WEIGHT.r) * r + WEIGHT.r * REFLECTION_COLOR.r,
			(1 - WEIGHT.g) * g + WEIGHT.g * REFLECTION_COLOR.g,
			(1 - WEIGHT.b) * b + WEIGHT.b * REFLECTION_COLOR.b);
		return FINAL_COLOR;
	}

	// Return Color
	return Color(r, g, b);
}

// Triangle Phong Shading
Color trianglePhong(Triangle &triangle, Light &light, glm::dvec3 &intersection, int depth) {

	// Get Vertices Positions
	glm::dvec3 VERTEX_A = glm::dvec3(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
	glm::dvec3 VERTEX_B = glm::dvec3(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
	glm::dvec3 VERTEX_C = glm::dvec3(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);

	// Triangle Area = 1/2 * Noraml of Triangle
	// |u x v| = |u x (va + vb)| = |u x va + 0| = |u||va| = |u||v|sin(theta) = Area
	glm::dvec3 TRIANGLE_ABC = glm::cross(VERTEX_B - VERTEX_A, VERTEX_C - VERTEX_A);
	glm::dvec3 TRIANGLE_BCI = glm::cross(VERTEX_C - VERTEX_B, intersection - VERTEX_B);
	glm::dvec3 TRIANGLE_CAI = glm::cross(VERTEX_A - VERTEX_C, intersection - VERTEX_C);

	// p = alpha * p0 + beta * p1 + gamma * p2
	double alpha = glm::dot(TRIANGLE_ABC, TRIANGLE_BCI) / glm::dot(TRIANGLE_ABC, TRIANGLE_ABC);
	double beta = glm::dot(TRIANGLE_ABC, TRIANGLE_CAI) / glm::dot(TRIANGLE_ABC, TRIANGLE_ABC);
	double gamma = 1.0 - alpha - beta;

	// Get Vertices Normals
	glm::dvec3 NORMAL_A = glm::dvec3(triangle.v[0].normal[0], triangle.v[0].normal[1], triangle.v[0].normal[2]);
	glm::dvec3 NORMAL_B = glm::dvec3(triangle.v[1].normal[0], triangle.v[1].normal[1], triangle.v[1].normal[2]);
	glm::dvec3 NORMAL_C = glm::dvec3(triangle.v[2].normal[0], triangle.v[2].normal[1], triangle.v[2].normal[2]);

	// Calculate Intersection Normal
	glm::dvec3 INTERPOLATE_NORMAL = glm::normalize(alpha * NORMAL_A + beta * NORMAL_B + gamma * NORMAL_C);

	// Setup Light Ray
	glm::dvec3 LIGHT_ORIGIN = glm::dvec3(light.position[0], light.position[1], light.position[2]);
	glm::dvec3 LIGHT_DIRECTION = glm::normalize(LIGHT_ORIGIN - intersection);

	// Magnitude of a Light at Intersection is its Projection onto the Normal
	double LIGHT = glm::dot(LIGHT_DIRECTION, INTERPOLATE_NORMAL);
	if (LIGHT > 1.0) { LIGHT = 1.0; }
	else if (LIGHT < 0.0) { LIGHT = 0.0; }

	// r = 2(l . n) n ¡V l
	// Reflect to Camera/Origin
	glm::dvec3 REFLECT = glm::normalize(2.0 * LIGHT * INTERPOLATE_NORMAL - LIGHT_DIRECTION);
	glm::dvec3 DIRECTION = glm::normalize(-intersection);
	double REFLECTION = glm::dot(REFLECT, DIRECTION);
	if (REFLECTION > 1.0) { REFLECTION = 1.0; }
	else if (REFLECTION < 0.0) { REFLECTION = 0.0; }

	// Faster Without Defining Three Colors and a New Operator
	Color DIFFUSE(
		alpha * triangle.v[0].color_diffuse[0] + beta * triangle.v[1].color_diffuse[0] + gamma * triangle.v[2].color_diffuse[0],
		alpha * triangle.v[0].color_diffuse[1] + beta * triangle.v[1].color_diffuse[1] + gamma * triangle.v[2].color_diffuse[1],
		alpha * triangle.v[0].color_diffuse[2] + beta * triangle.v[1].color_diffuse[2] + gamma * triangle.v[2].color_diffuse[2]);

	Color SPECULAR(
		alpha * triangle.v[0].color_specular[0] + beta * triangle.v[1].color_specular[0] + gamma * triangle.v[2].color_specular[0],
		alpha * triangle.v[0].color_specular[1] + beta * triangle.v[1].color_specular[1] + gamma * triangle.v[2].color_specular[1],
		alpha * triangle.v[0].color_specular[2] + beta * triangle.v[1].color_specular[2] + gamma * triangle.v[2].color_specular[2]);

	double SHININESS = alpha * triangle.v[0].shininess + beta * triangle.v[1].shininess + gamma * triangle.v[2].shininess;

	// Local Phong Shading, I = L(kd * (l . n) + ks (r . v)^alpha)
	double r = light.color[0] * (DIFFUSE.r * LIGHT + (SPECULAR.r * pow(REFLECTION, SHININESS)));
	double g = light.color[1] * (DIFFUSE.g * LIGHT + (SPECULAR.g * pow(REFLECTION, SHININESS)));
	double b = light.color[2] * (DIFFUSE.b * LIGHT + (SPECULAR.b * pow(REFLECTION, SHININESS)));

	// Return Color
	return Color(r, g, b);
}

// Generate Random Unit Vector for Sphere
glm::dvec3 randomUnitVector() {
	// theta in range [0, 2PI)
	double theta = ((double)(rand() % 360) * PI / 180.0);
	// u in range [-1, 1]
	double u = ((double)rand() / (RAND_MAX)) * 2.0 - 1.0;
	double s = sqrt(1.0 - pow(u, 2));
	return glm::dvec3(u, s * cos(theta), s * sin(theta));
}

// Setup Light Source
vector<Light> createLightSource(Light lights[], double num) {
	vector<Light> LIGHT_SOURCE;
	NUM_RAND_LIGHTS = (int)num;
	// If Soft Shadow is Enabled, Create Spherical Light with num of Lights
	// Else Store Orginal Light(s)
	if (softshadow) {
		double num_spots = num;
		// For Each Original Light, Generate num of Lights
		for (int i = 0; i < num_lights; i++) {
			for (int j = 0; j < (int)num_spots; j++) {
				Light NEW_LIGHT;
				glm::dvec3 DIRECTION = glm::normalize(randomUnitVector());
				glm::dvec3 LIGHT_ORIGIN = glm::dvec3(lights[i].position[0], lights[i].position[1], lights[i].position[2]);
				glm::dvec3 LIGHT_POSITION = LIGHT_ORIGIN + DIRECTION / 10.0;
				Color LIGHT_COLOR = Color(lights[i].color[0] / num_spots, lights[i].color[1] / num_spots, lights[i].color[2] / num_spots);

				NEW_LIGHT.position[0] = LIGHT_POSITION.x;
				NEW_LIGHT.position[1] = LIGHT_POSITION.y;
				NEW_LIGHT.position[2] = LIGHT_POSITION.z;
				NEW_LIGHT.color[0] = LIGHT_COLOR.r;
				NEW_LIGHT.color[1] = LIGHT_COLOR.g;
				NEW_LIGHT.color[2] = LIGHT_COLOR.b;

				LIGHT_SOURCE.push_back(NEW_LIGHT);
			}
		}
	}
	else {
		for (int i = 0; i < num_lights; i++) {
			LIGHT_SOURCE.push_back(lights[i]);
		}
	}
	return LIGHT_SOURCE;
}

//Trace Spheres
void traceSphere(Ray &ray, Color &PIXEL_COLOR, double &MAX_DISTANCE, int depth) {
	// For Each Sphere
	for (int i = 0; i < num_spheres; i++) {
		// Initial Intersection
		glm::dvec3 INTERSECTION = glm::dvec3(0, 0, MAX_DISTANCE);
		// Detect Intersection With Sphere and Change Interscetion
		if (ray.intersectSphere(spheres[i], INTERSECTION) && INTERSECTION.z > MAX_DISTANCE) {

			// Change Pixel to White
			PIXEL_COLOR = Color(0.0, 0.0, 0.0);

			//Cast Shadow Ray from Surface Point to Each Light
			for (int j = 0; j < LIGHT_SOURCE.size(); j++) {
				glm::dvec3 LIGHT_POSITION = glm::dvec3(LIGHT_SOURCE[j].position[0], LIGHT_SOURCE[j].position[1], LIGHT_SOURCE[j].position[2]);
				glm::dvec3 RAY_ORIGIN = INTERSECTION;
				glm::dvec3 RAY_DIRECTION = glm::normalize(LIGHT_POSITION - RAY_ORIGIN);

				Ray SHADOW_RAY(RAY_ORIGIN, RAY_DIRECTION, ray.time);

				// If Shadow Ray Hits Opaque Object, No Contribution from that Light
				// If Shadow Ray Can Reach to the Light, Apply a Standard Phong Model

				bool within = false;

				// Check Every Shadow of Other Spheres
				for (int k = 0; k < num_spheres; k++) {
						glm::dvec3 SHADOW_INTERSECTION = glm::dvec3(0, 0, MAX_DISTANCE);
						if (SHADOW_RAY.intersectSphere(spheres[k], SHADOW_INTERSECTION) && i != k) {
							double SPHERE_DISTANCE = glm::length(SHADOW_INTERSECTION - RAY_ORIGIN);
							double LIGHT_DISTANCE = glm::length(LIGHT_POSITION - RAY_ORIGIN);
							if (SPHERE_DISTANCE < LIGHT_DISTANCE) {
								within = true;
								break;
							}
						}
					}

				// Check Every Shadow of Triangles
				if (!within) {
					for (int k = 0; k < num_triangles; k++) {
						glm::dvec3 SHADOW_INTERSECTION = glm::dvec3(0, 0, MAX_DISTANCE);
						if (SHADOW_RAY.intersectTriangle(triangles[k], SHADOW_INTERSECTION)) {
							double TRIANGLE_DISTANCE = glm::length(SHADOW_INTERSECTION - RAY_ORIGIN);
							double LIGHT_DISTANCE = glm::length(LIGHT_POSITION - RAY_ORIGIN);
							if (TRIANGLE_DISTANCE < LIGHT_DISTANCE) {
								within = true;
								break;
							}
						}
					}
				}
					
				// If Not in Shadow, Apply Phong Shading
				if (!within) {
					PIXEL_COLOR += spherePhong(spheres[i], LIGHT_SOURCE[j], INTERSECTION, depth, ray.time);
				}

				if (MAX_RAY_DEPTH != 0 && within) {
					Color tmp = spherePhong(spheres[i], LIGHT_SOURCE[j], INTERSECTION, depth, ray.time);
					PIXEL_COLOR.r = (PIXEL_COLOR.r + tmp.r) / 2.0;
					PIXEL_COLOR.g = (PIXEL_COLOR.g + tmp.g) / 2.0;
					PIXEL_COLOR.b = (PIXEL_COLOR.b + tmp.b) / 2.0;
				}
			}
			// Set Max Distance to Current Intersection
			MAX_DISTANCE = INTERSECTION.z;
		}
	}
}

//Trace Triangles
void traceTriangle(Ray &ray, Color &PIXEL_COLOR, double &MAX_DISTANCE, int depth) {
	// For Each Triangle
	for (int i = 0; i < num_triangles; i++) {
		// Initial Intersection
		glm::dvec3 INTERSECTION = glm::dvec3(0, 0, MAX_DISTANCE);
		// Detect Intersection With Triangle and Change Interscetion
		if (ray.intersectTriangle(triangles[i], INTERSECTION) && INTERSECTION.z > MAX_DISTANCE) {

			// Change Pixel to White
			PIXEL_COLOR = Color(0.0, 0.0, 0.0);

			//Cast Shadow Ray from Surface Point to Each Light
			for (int l = 0; l < LIGHT_SOURCE.size(); l++) {
				glm::dvec3 LIGHT_POSITION = glm::dvec3(LIGHT_SOURCE[l].position[0], LIGHT_SOURCE[l].position[1], LIGHT_SOURCE[l].position[2]);
				glm::dvec3 RAY_ORIGIN = INTERSECTION;
				glm::dvec3 RAY_DIRECTION = glm::normalize(LIGHT_POSITION - RAY_ORIGIN);

				Ray SHADOW_RAY(RAY_ORIGIN, RAY_DIRECTION, ray.time);

				// If Shadow Ray Hits Opaque Object, No Contribution from that Light
				// If Shadow Ray Can Reach to the Light, Apply a Standard Phong Model

				bool within = false;

				// Check Every Shadow of Spheres
				for (int k = 0; k < num_spheres; k++) {
						glm::dvec3 SHADOW_INTERSECTION = glm::dvec3(0, 0, MAX_DISTANCE);
						if (SHADOW_RAY.intersectSphere(spheres[k], SHADOW_INTERSECTION)) {
							double SPHERE_DISTANCE = glm::length(SHADOW_INTERSECTION - RAY_ORIGIN);
							double LIGHT_DISTANCE = glm::length(LIGHT_POSITION - RAY_ORIGIN);
							if (SPHERE_DISTANCE < LIGHT_DISTANCE) {
								within = true;
								break;
							}
						}
					}

				// Check Every Shadow of Other Triangles
				if (!within) {
					for (int k = 0; k < num_triangles; k++) {
						glm::dvec3 SHADOW_INTERSECTION = glm::dvec3(0, 0, MAX_DISTANCE);
						if (SHADOW_RAY.intersectTriangle(triangles[k], SHADOW_INTERSECTION) && i != k) {
							double TRIANGLE_DISTANCE = glm::length(SHADOW_INTERSECTION - RAY_ORIGIN);
							double LIGHT_DISTANCE = glm::length(LIGHT_POSITION - RAY_ORIGIN);
							if (TRIANGLE_DISTANCE < LIGHT_DISTANCE) {
								within = true;
								break;
							}
						}
					}
				}

				// If Not in Shadow, Apply Phong Shading
				if (!within) {
						PIXEL_COLOR += trianglePhong(triangles[i], LIGHT_SOURCE[l], INTERSECTION, depth);
					}
			}
			// Set Max Distance to Current Intersection
			MAX_DISTANCE = INTERSECTION.z;
		}
	}
}

// Trace Ray
Color traceRay(Ray &ray, int depth) {

	double MAX_DISTANCE = -INFINITY; // Negative Z Direction
	Color AMBIENT_COLOR(ambient_light[0], ambient_light[1], ambient_light[2]);
	Color PIXEL_COLOR(1.0, 1.0, 1.0);

	// Trace Spheres and Triagnles
	traceSphere(ray, PIXEL_COLOR, MAX_DISTANCE, depth);
	traceTriangle(ray, PIXEL_COLOR, MAX_DISTANCE, depth);

	PIXEL_COLOR += AMBIENT_COLOR;
	if (softshadow && depth != 0) {
		// Lower the Weight of Reflection Pixel Color
		return PIXEL_COLOR / (double)NUM_RAND_LIGHTS;
	}
	return PIXEL_COLOR;
}

// Setup Ray reference: Scratchapixel
Ray setupRay(double x, double y) {

	// Two Quadrant Size of Moving Half Pixel
	// Since Screen Space Coordinate (-1, -1) to (1, 1)
	double SCREEN_X = (2 * ((x + 0.5) / WIDTH)) - 1;
	double SCREEN_Y = (2 * ((y + 0.5) / HEIGHT)) - 1;

	double ASPECT_RATIO = (double) WIDTH / HEIGHT;

	double CAMERA_X = SCREEN_X * ASPECT_RATIO * tan(RADIAN);
	double CAMERA_Y = SCREEN_Y * tan(RADIAN);

	glm::dvec3 origin = CAMERA;
	glm::dvec3 direction = glm::normalize(glm::dvec3(CAMERA_X, CAMERA_Y, -1.0));

	//double time = ((double)rand() / (RAND_MAX));

	// Shoot Ray at Time 0 ~ 10
	double time = ((double)rand() / (RAND_MAX)) * 10.0;

	return Ray(origin, direction, time);
}

//MODIFY THIS FUNCTION
void draw_scene(double offsetX, double offsetY, double offsetZ) {
	// If Animation is Enabled, Move Sphere 1
	if (animation) {
		//CAMERA = glm::dvec3(0.0 + offsetX, 0.0 + offsetY, 0.0 + offsetZ);
		if (num_spheres != 0) {
			spheres[0].position[0] += offsetX;
			spheres[0].position[1] += offsetY;
			spheres[0].position[2] += offsetZ;
		}
	}

	//Offset to Split into Rays
	vector<double> offsetRow;
	vector<double> offsetCol;

	if (!motionblur) {
		double list[2] = { -0.25, 0.25 };
		//double list[4] = { -0.375, -0.125, 0.125, 0.375 };
		//double list[6] = { -0.4375, -0.1875, -0.0625, 0.0625, 0.1875, 0.4375 };
		//double list[8] = { -0.46875, -0.21875, -0.09375, -0.03125, 0.03125, 0.09375, 0.21875, 0.46875 };
		for (int i = 0; i < 2; i++) {
			offsetRow.push_back(list[i]);
			offsetCol.push_back(list[i]);
		}

	}

	// If Motion Blur is Enabled, Need More Rays to Look Smooth
	if (motionblur) {
		double list[9] = { -0.46875, -0.21875, -0.09375, -0.03125, 0.0, 0.03125, 0.09375, 0.21875, 0.46875 };
		for (int i = 0; i < 9; i++) {
			offsetRow.push_back(list[i]);
			offsetCol.push_back(list[i]);
		}
	}

	LIGHT_SOURCE = createLightSource(lights, 20.0);

	// For Every Coordinate
	for (unsigned int x = 0; x < WIDTH; x++) {
		glPointSize(2.0);
		glBegin(GL_POINTS);
		for (unsigned int y = 0; y < HEIGHT; y++) {
			Color color(0, 0, 0);
			// If Antialiasing is Enabled
			if (antialias) {
				// Antialiasing Using Multiple Rays
				for (int k = 0; k < offsetRow.size(); k++) {
					for (int l = 0; l < offsetCol.size(); l++) {
						// Fire Ray to Center of Subpixels
						Ray ray = setupRay(x + offsetRow[k], y + offsetCol[l]);
						Color tmpcolor = traceRay(ray, 0);

						// Add Seperately to Prevent Clamping
						color.r += tmpcolor.r;
						color.g += tmpcolor.g;
						color.b += tmpcolor.b;
					}
				}
				// Average the Colors
				color = color / pow((double)offsetRow.size(), 2);
			}
			else {
				// Normal Mode
				Ray ray = setupRay(x, y);
				color = traceRay(ray, 0);
			}
			plot_pixel(x, y, color.r * 255, color.g * 255, color.b * 255);
		}
		glEnd();
		glFlush();
	}
	printf("Done!\n"); std::fflush(stdout);
}

// Draw Scene Longer Method
/*void draw_scene(double offsetX, double offsetY, double offsetZ) {

	double CAMERA_WIDTH = (double) WIDTH / HEIGHT * tan(RADIAN);
	double CAMERA_HEIGHT = tan(RADIAN);

	// Width of Two Quadrant Size
	double PIXEL = (double)(2 * CAMERA_WIDTH) / WIDTH;
	// Counter for  Pixels
	int ROW_PIXEL = 0;
	int COL_PIXEL = 0;

	if (animation) {
		//CAMERA = glm::dvec3(0.0 + offsetX, 0.0 + offsetY, 0.0 + offsetZ);
		if (num_spheres != 0) {
			spheres[0].position[0] += offsetX;
			spheres[0].position[1] += offsetY;
			spheres[0].position[2] += offsetZ;
		}
	}

	//Offset to Split into Rays
	vector<double> offsetRow;
	vector<double> offsetCol;

	if (!motionblur) {
		double list[2] = { -0.25, 0.25 };
		for (int i = 0; i < 2; i++) {
			offsetRow.push_back(list[i]);
			offsetCol.push_back(list[i]);
		}
	}

	// If Motion Blur is Enabled, Need More Rays to Look Smooth
	if (motionblur) {
		double list[9] = { -0.46875, -0.21875, -0.09375, -0.03125, 0.0, 0.03125, 0.09375, 0.21875, 0.46875 };
		for (int i = 0; i < 9; i++) {
			offsetRow.push_back(list[i]);
			offsetCol.push_back(list[i]);
		}
	}

	LIGHT_SOURCE = createLightSource(lights, 20.0);

	//Screen Space Coordinate (-1, -1) to (1, 1)
	for(double x = -CAMERA_WIDTH; x < CAMERA_WIDTH - PIXEL; x += PIXEL) {
		glPointSize(2.0);  
		glBegin(GL_POINTS);
		COL_PIXEL = 0;
		for(double y = -CAMERA_HEIGHT; y < CAMERA_HEIGHT - PIXEL; y += PIXEL) {
			Color color(0, 0, 0);
			if (antialias) {
				// Antialiasing Using Multiple Rays
				for (int k = 0; k < offsetRow.size(); k++) {
					for (int l = 0; l < offsetCol.size(); l++) {
					// Camera Origin
					glm::dvec3 origin = glm::dvec3(0.0, 0.0, 0.0);
					// Fire Ray to Center of Subpixels
					glm::dvec3 direction = glm::normalize(glm::dvec3((x + PIXEL * offsetRow[k]), (y + PIXEL * offsetCol[l]), -1.0));
					Ray ray(origin, direction);
					Color tmpcolor = traceRay(ray, 0);

					// Add Seperately to Prevent Clamping
					color.r += tmpcolor.r;
					color.g += tmpcolor.g;
					color.b += tmpcolor.b;
					}
				}
				// Average the Colors
				color = color / pow((double)offsetRow.size(),2);
			}
			else {
				// Camera Origin
				glm::dvec3 origin = glm::dvec3(0.0, 0.0, 0.0);
				// Fire Ray to Center of Pixel
				glm::dvec3 direction = glm::normalize(glm::dvec3((x + PIXEL / 2), (y + PIXEL / 2), -1.0));
				Ray ray(origin, direction);
				color = traceRay(ray, 0);
			}

			plot_pixel(ROW_PIXEL, COL_PIXEL, color.r * 255, color.g * 255, color.b * 255);
			COL_PIXEL++;
		}
		glEnd();
		glFlush();
		ROW_PIXEL++;
	}
	printf("Done!\n"); fflush(stdout);
}*/

// Plot Pixel to Screen
void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
	glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
	glVertex2i(x, y);
}

// Set JPEG Buffer
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
	buffer[y][x][0] = r;
	buffer[y][x][1] = g;
	buffer[y][x][2] = b;
}

// Plot Pixel
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
	plot_pixel_display(x, y, r, g, b);
	if(mode == MODE_JPEG)
		plot_pixel_jpeg(x, y, r, g, b);
}

// Save JPEG with filename/path+filename
void save_jpg(const char * filename) {
	printf("Saving JPEG file: %s\n", filename);

	ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
	if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
		printf("Error in Saving\n");
	else 
	printf("File saved Successfully\n");
}

// Check Read Input
void parse_check(const char *expected, char *found) {
	if(strcasecmp(expected, found)) {
		printf("Expected '%s ' found '%s '\n", expected, found);
		printf("Parse error, abnormal abortion\n");
		exit(0);
	}
}

// Read Doubles
void parse_doubles(FILE* file, const char *check, double p[3]) {
	char str[100];
	fscanf(file, "%s", str);
	parse_check(check, str);
	fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
	printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

// Read Radius
void parse_rad(FILE *file, double *r) {
	char str[100];
	fscanf(file, "%s", str);
	parse_check("rad:", str);
	fscanf(file, "%lf", r);
	printf("rad: %f\n", *r);
}

// Read Shininess
void parse_shi(FILE *file, double *shi) {
	char s[100];
	fscanf(file, "%s", s);
	parse_check("shi:", s);
	fscanf(file, "%lf", shi);
	printf("shi: %f\n", *shi);
}

// Load Scene File
int loadScene(char *argv) {
	FILE * file = fopen(argv,"r");
	int number_of_objects;
	char type[50];
	Triangle t;
	Sphere s;
	Light l;
	fscanf(file,"%i", &number_of_objects);

	printf("number of objects: %i\n", number_of_objects);

	parse_doubles(file, "amb:", ambient_light);

	for(int i = 0 ; i < number_of_objects ; i++) {
		fscanf(file,"%s\n", type);
		printf("%s\n", type);
		if(strcasecmp(type, "triangle") == 0) {
			printf("found triangle\n");
			for(int j = 0 ; j < 3 ; j++) {
				parse_doubles(file, "pos:", t.v[j].position);
				parse_doubles(file, "nor:", t.v[j].normal);
				parse_doubles(file, "dif:", t.v[j].color_diffuse);
				parse_doubles(file, "spe:", t.v[j].color_specular);
				parse_shi(file, &t.v[j].shininess);
			}

			if(num_triangles == MAX_TRIANGLES) {
				printf("too many triangles, you should increase MAX_TRIANGLES!\n");
				exit(0);
			}
			triangles[num_triangles++] = t;
		}
		else if(strcasecmp(type, "sphere") == 0){
			printf("found sphere\n");

			parse_doubles(file, "pos:", s.position);
			parse_rad(file, &s.radius);
			parse_doubles(file, "dif:", s.color_diffuse);
			parse_doubles(file, "spe:", s.color_specular);
			parse_shi(file, &s.shininess);

			if(num_spheres == MAX_SPHERES) {
				printf("too many spheres, you should increase MAX_SPHERES!\n");
				exit(0);
			}
			spheres[num_spheres++] = s;
		}
		else if(strcasecmp(type,"light") == 0) {
			printf("found light\n");
			parse_doubles(file, "pos:", l.position);
			parse_doubles(file, "col:", l.color);

			if(num_lights == MAX_LIGHTS) {
				printf("too many lights, you should increase MAX_LIGHTS!\n");
				exit(0);
			}
			lights[num_lights++] = l;
		}
		else {
			printf("unknown type in scene description:\n%s\n", type);
			exit(0);
		}
	}
	return 0;
}

// Display Function
void displayFunc() {
}

// Keyboard Function
void keyboardFunc(unsigned char key, int x, int y) {
	switch (key) {
	case 27: // ESC key
		exit(0); // Exit the program
		break;
	}
}

// Used Function with Single Buffer
void idleFunc() {
	static int counter = 0;
	static int offset = 0;
	if (animation && !motionblur) {
		// If only Animation is Enabled
		if (counter < 60) {
			if (counter < 15) { draw_scene(0.0, offset++ / 60.0, 0.0); }
			else if(counter < 30){ draw_scene(0.0, -offset-- / 60.0, 0.0); }
			else if(counter < 45){ draw_scene(0.0, offset++ / 60.0, 0.0); }
			else if(counter < 60){ draw_scene(0.0, -offset-- / 60.0, 0.0); }

			if (mode == MODE_JPEG) {
				char anim_num[5];
				sprintf(anim_num, "%03d", ++counter);
				save_jpg(("./Animation/" + string(anim_num) + ".jpg").c_str());
			}
		}
		if (counter == 60) {
			cout << "Images Saved." << endl;
		}
	}
	else if (animation && motionblur) {
		// If Animation with Motion Blur
		if (counter < 15) {
			draw_scene(0.0, offset++ / 60.0, 0.0);

			if (mode == MODE_JPEG) {
				char anim_num[5];
				sprintf(anim_num, "%03d", ++counter);
				save_jpg(("./MotionBlur/" + string(anim_num) + ".jpg").c_str());
			}
		}
		if (counter == 15) {
			cout << "Images Saved." << endl;
		}
	}
	else {
		//hack to make it only draw once
		// Still Image
		if (!counter) {
			draw_scene(0.0, 0.0, 0.0);
			if (mode == MODE_JPEG) {
				save_jpg(filename);
			}
		}
		counter = 1;
	}
}

// Init Function
void init() {
	// Set Projection Matrix
	glMatrixMode(GL_PROJECTION);
	// Set Orthographic Projection
	glOrtho(0, WIDTH, 0, HEIGHT, 1, -1);
	// Initialize Modelview Matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// Reset Background Color
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);
}

int main(int argc, char ** argv) {
	if ((argc < 2) || (argc > 8)) {  
		printf ("Usage: %s <input scenefile> [output jpegname] [antialiasing: y/n] [softshadow: y/n] [recurse: y/n] [animation: y/n] [motionblur: y/n]\n", argv[0]);
		exit(0);
	}

	if(argc >= 3) {
		mode = MODE_JPEG;
		filename = argv[2];
		if (strlen(filename) > 4 && !strcmp(filename + strlen(filename) - 4, ".jpg")) {}
		else { exit(0); }
	}
	else if (argc == 2) {
		mode = MODE_DISPLAY;
	}

	// Input Arguments
	if (argc >= 4) { antialias = (tolower(*argv[3]) == 'y') ? true : false;	}
	if (argc >= 5) { softshadow = (tolower(*argv[4]) == 'y') ? true : false; }
	if (argc >= 6) { MAX_RAY_DEPTH = (tolower(*argv[5]) == 'y') ? 3 : 0; }
	if (argc >= 7) { animation = (tolower(*argv[6]) == 'y') ? true : false; }
	if (argc == 8) { motionblur = (tolower(*argv[7]) == 'y') ? true : false; }

	glutInit(&argc,argv);
	loadScene(argv[1]);

	glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(WIDTH, HEIGHT);
	int window = glutCreateWindow("Ray Tracer");
	glutDisplayFunc(displayFunc);
	glutIdleFunc(idleFunc);
	glutKeyboardFunc(keyboardFunc);

	init();
	glutMainLoop();
}

