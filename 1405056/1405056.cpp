#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <iostream>
#include <windows.h>
#include <GL/glut.h>
#include "bits/stdc++.h"
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))
#define deg_to_rad (0.0175)
using namespace std;
//bitmap_image b_img("texture.bmp");
double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
ifstream inFile("description.txt");
double near_p, far_p,fovy, aspectratio;
int level_of_recursion;
int number_of_pixels;
int checker_board_width;
double ambient_c, diffuse_c, reflection_c;
int number_of_objects;
double sphere_centre_x, sphere_centre_y, sphere_centre_z;
double sphere_radius;
double sphere_r, sphere_g, sphere_b;
double ambient_s, diffuse_s, specular_s, reflection_s, shininess_s;
double lowest_x, lowest_y, lowest_z;
double width_p, height_p;
double pyramid_r, pyramid_g, pyramid_b;
double ambient_p, diffuse_p, specular_p, reflection_p, shininess_p;
double screen_height, screen_width;

class point
{
public:
	double x, y, z;
	point(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	point() {

	}
	point operator*(double m) {
		point p(this->x*m, this->y*m, this->z*m);
		return p;
	}

	static double dot(point m, point n) {
		return n.x*m.x + n.y*m.y + n.z*m.z;
	}

	void normalize() {
		double sqrt_val = sqrt(pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2));
		this->x = this->x / sqrt_val;
		this->y = this->y / sqrt_val;
		this->z = this->z / sqrt_val;
	}
	point operator+(point p) {
		point p1(this->x + p.x, this->y + p.y, this->z + p.z);
		return p1;
	}

	void print() {
		cout << this->x << " " << this->y << " " << this->z << endl;
	}
};

class Square {
public:
	point a, b, c, d;
	Square(point a, point b, point c, point d) {
		this->a = a;
		this->b = b;
		this->c = c;
		this->d = d;
	}
	Square() {

	}
};

class Triangle {
public:
	point a, b, c;
	Triangle(point a, point b, point c) {
		this->a = a;
		this->b = b;
		this->c = c;
	}
	Triangle() {

	}
};

//class Pyramid {
//public:
//	Triangle t1, t2, t3, t4;
//	Square s;
//	Pyramid(Triangle t1, Triangle t2, Triangle t3, Triangle t4, Square s) {
//		this->t1 = t1;
//		this->t2 = t2;
//		this->t3 = t3;
//		this->t4 = t4;
//		this->s = s;
//	}
//};

//vector<Pyramid> pyramids;
class object {
public:
	string type;
	double sphere_centre_x, sphere_centre_y, sphere_centre_z;
	double sphere_radius;
	double r, g, b;
	double ambient_s, diffuse_s, specular_s, reflection_s, shininess_s;
	double lowest_x, lowest_y, lowest_z;
	double width_p, height_p;
	//double pyramid_r, pyramid_g, pyramid_b;
	double ambient_p, diffuse_p, specular_p, reflection_p, shininess_p;
	Triangle t1, t2, t3, t4;
	Square s;
	object(string type, double sphere_centre_x, double sphere_centre_y, double sphere_centre_z,
		double sphere_radius,
		double sphere_r, double sphere_g, double sphere_b,
		double ambient_s, double diffuse_s, double specular_s, double reflection_s, double shininess_s
		) {
		this->type = type;
		this->sphere_centre_x = sphere_centre_x;
		this->sphere_centre_y = sphere_centre_y;
		this->sphere_centre_z = sphere_centre_z;
		this->sphere_radius = sphere_radius;
		this->r = sphere_r;
		this->g = sphere_g;
		this->b = sphere_b;
		this->ambient_s = ambient_s;
		this->diffuse_s = diffuse_s;
		this->specular_s = specular_s;
		this->reflection_s = reflection_s;
		this->shininess_s = shininess_s;
	}
	object(string type, double lowest_x, double lowest_y, double lowest_z,
		double width_p, double height_p,
		double pyramid_r, double pyramid_g, double pyramid_b,
		double ambient_p, double diffuse_p, double specular_p, double reflection_p, double shininess_p,
		Triangle t1, Triangle t2, Triangle t3, Triangle t4, Square s) {
		this->type = type;
		this->lowest_x = lowest_x;
		this->lowest_y = lowest_y;
		this->lowest_z = lowest_z;
		this->width_p = width_p;
		this->height_p = height_p;
		this->r = pyramid_r;
		this->g = pyramid_g;
		this->b = pyramid_b;
		this->ambient_p = ambient_p;
		this->diffuse_p = diffuse_p;
		this->specular_p = specular_p;
		this->reflection_p = reflection_p;
		this->shininess_p = shininess_p;
		this->t1 = t1;
		this->t2 = t2;
		this->t3 = t3;
		this->t4 = t4;
		this->s = s;
	}
};
//vector<Sphere> spheres;
//vector<Pyramid> pyramids;
vector<object> objects;


point screen_mid_point;
vector<vector<point>> points;
double dr, du;

class Color {
public:
	double r, g, b;
	Color(double r, double g, double b) {
		this->r = r;
		this->g = g;
		this->b = b;
	}
	Color() {
	}


};


class matrix
{
public:
	double values[3][3];
	int num_rows, num_cols;

	// only set the number of rows and cols
	matrix(int rows, int cols)
	{
		assert(rows <= 4 && cols <= 4);
		num_rows = rows;
		num_cols = cols;
	}

	// prepare an nxn square matrix
	matrix(int n)
	{
		assert(n <= 4);
		num_rows = num_cols = n;
	}

	// prepare and return an identity matrix of size nxn
	static matrix make_identity(int n)
	{
		assert(n <= 4);
		matrix m(n);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j)
					m.values[i][j] = 1;
				else
					m.values[i][j] = 0;
			}
		}
		return m;
	}

	// print the matrix. exists for testing purposes
	void print()
	{
		cout << "Matrix:" << endl;
		for (int i = 0; i < num_rows; i++)
		{
			for (int j = 0; j < num_cols; j++)
			{
				cout << values[i][j] << "\t";
			}
			cout << endl;
		}
	}

};

//Color **textureBuffer;
//int height, width;
//
//height = b_img.height();
//width = b_img.width();
//textureBuffer = new Color*[width];
//for (int i = 0; i < width; i++) {
//	textureBuffer[i] = new Color[height];
//	for (int j = 0; j < height; j++) {
//		unsigned char r, g, b;
//		b_img.get_pixel(i, j, r, g, b);
//		Color c(r / 255.0, g / 255.0, b / 255.0);
//		textureBuffer[i][j] = c;
//	}
//}


point pos(100, 100, 20), u(0, 0, 1), r(-.707, .707, 0), l(-.707, -.707, 0);
point temp;
//double cameraAngle;


point getCrossProd(point a, point b) {
	double magnitude;
	point temp;
	temp.x = a.y*b.z - a.z*b.y;
	temp.y = a.z*b.x - a.x*b.z;
	temp.z = a.x*b.y - a.y*b.x;
	return temp;
}

point rotateVec(point a, point p) {
	a.x = a.x*cos(cameraAngle) + p.x*sin(cameraAngle);
	a.y = a.y*cos(cameraAngle) + p.y*sin(cameraAngle);
	a.z = a.z*cos(cameraAngle) + p.z*sin(cameraAngle);
	return a;
}

void drawAxes()
{
	if (drawaxes == 1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES); {
			glVertex3f(100, 0, 0);
			glVertex3f(-100, 0, 0);

			glVertex3f(0, -100, 0);
			glVertex3f(0, 100, 0);

			glVertex3f(0, 0, 100);
			glVertex3f(0, 0, -100);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if (drawgrid == 1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES); {
			for (i = -8; i <= 8; i++) {

				if (i == 0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i * 10, -90, 0);
				glVertex3f(i * 10, 90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i * 10, 0);
				glVertex3f(90, i * 10, 0);
			}
		}glEnd();
	}
}

void drawSquare(point a)
{
	//glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS); {
		glVertex3f(a.x, a.y, 0);
		glVertex3f(a.x, a.y + checker_board_width, 0);
		glVertex3f(a.x + checker_board_width, a.y + checker_board_width, 0);
		glVertex3f(a.x + checker_board_width, a.y, 0);
	}glEnd();
}

void drawPyramid(double width, double height, double lowest_x, double lowest_y, double lowest_z) {
	point apex = { lowest_x, lowest_y, lowest_z + height };
	double a = width / 2;
	glBegin(GL_TRIANGLES); {
		glVertex3f(apex.x, apex.y, apex.z);
		glVertex3f(lowest_x - a, lowest_y - a, lowest_z);
		glVertex3f(lowest_x + a, lowest_y - a, lowest_z);
	}glEnd();
	/*point b(lowest_x - a, lowest_y - a, lowest_z);
	point c(lowest_x + a, lowest_y - a, lowest_z);
	Triangle t1(apex, b, c);
*/
	glBegin(GL_TRIANGLES); {
		glVertex3f(apex.x, apex.y, apex.z);
		glVertex3f(lowest_x + a, lowest_y - a, lowest_z);
		glVertex3f(lowest_x + a, lowest_y + a, lowest_z);
	}glEnd();

	/*point b1(lowest_x + a, lowest_y - a, lowest_z);
	point c1(lowest_x + a, lowest_y + a, lowest_z);
	Triangle t2(apex, b1, c1);
	*/
	glBegin(GL_TRIANGLES); {
		glVertex3f(apex.x, apex.y, apex.z);
		glVertex3f(lowest_x + a, lowest_y + a, lowest_z);
		glVertex3f(lowest_x - a, lowest_y + a, lowest_z);
	}glEnd();

	/*point b2(lowest_x + a, lowest_y + a, lowest_z);
	point c2(lowest_x - a, lowest_y + a, lowest_z);
	Triangle t3(apex, b2, c2);*/
	
	glBegin(GL_TRIANGLES); {
		glVertex3f(apex.x, apex.y, apex.z);
		glVertex3f(lowest_x - a, lowest_y + a, lowest_z);
		glVertex3f(lowest_x - a, lowest_y - a, lowest_z);
	}glEnd();

	/*point b3(lowest_x - a, lowest_y + a, lowest_z);
	point c3(lowest_x - a, lowest_y - a, lowest_z);
	Triangle t4(apex, b3, c3);*/

	glBegin(GL_QUADS); {
		glVertex3f(lowest_x - a, lowest_y - a, lowest_z);
		glVertex3f(lowest_x + a, lowest_y - a, lowest_z);
		glVertex3f(lowest_x + a, lowest_y + a, lowest_z);
		glVertex3f(lowest_x - a, lowest_y + a, lowest_z);
	}glEnd();

	/*point a1(lowest_x - a, lowest_y - a, lowest_z);
	point b4(lowest_x + a, lowest_y - a, lowest_z);
	point c4(lowest_x + a, lowest_y + a, lowest_z);
	point d4(lowest_x - a, lowest_y + a, lowest_z);
	Square s(a1, b4, c4, d4);*/

	/*Pyramid p(t1, t2, t3, t4, s);
	pyramids.push_back(p);*/

}

//void drawCircle(double radius, int segments)
//{
//	int i;
//	struct point points[100];
//	glColor3f(0.7, 0.7, 0.7);
//	//generate points
//	for (i = 0; i <= segments; i++)
//	{
//		points[i].x = radius * cos(((double)i / (double)segments) * 2 * pi);
//		points[i].y = radius * sin(((double)i / (double)segments) * 2 * pi);
//	}
//	//draw segments using generated points
//	for (i = 0; i < segments; i++)
//	{
//		glBegin(GL_LINES);
//		{
//			glVertex3f(points[i].x, points[i].y, 0);
//			glVertex3f(points[i + 1].x, points[i + 1].y, 0);
//		}
//		glEnd();
//	}
//}

void drawCone(double radius, double height, int segments)
{
	int i;
	double shade;
	point points[100];
	//generate points
	for (i = 0; i <= segments; i++)
	{
		points[i].x = radius * cos(((double)i / (double)segments) * 2 * pi);
		points[i].y = radius * sin(((double)i / (double)segments) * 2 * pi);
	}
	//draw triangles using generated points
	for (i = 0; i < segments; i++)
	{
		//create shading effect
		if (i < segments / 2)shade = 2 * (double)i / (double)segments;
		else shade = 2 * (1.0 - (double)i / (double)segments);
		glColor3f(shade, shade, shade);

		glBegin(GL_TRIANGLES);
		{
			glVertex3f(0, 0, height);
			glVertex3f(points[i].x, points[i].y, 0);
			glVertex3f(points[i + 1].x, points[i + 1].y, 0);
		}
		glEnd();
	}
}


void drawSphere(double radius, int slices, int stacks)
{
	point points[100][100];
	int i, j;
	double h, r;
	//generate points
	for (i = 0; i <= stacks; i++)
	{
		h = radius * sin(((double)i / (double)stacks)*(pi / 2));
		r = radius * cos(((double)i / (double)stacks)*(pi / 2));
		for (j = 0; j <= slices; j++)
		{
			points[i][j].x = r * cos(((double)j / (double)slices) * 2 * pi);
			points[i][j].y = r * sin(((double)j / (double)slices) * 2 * pi);
			points[i][j].z = h;
		}
	}
	//draw quads using generated points
	for (i = 0; i < stacks; i++)
	{
		//glColor3f((double)i / (double)stacks, (double)i / (double)stacks, (double)i / (double)stacks);
		for (j = 0; j < slices; j++)
		{
			glBegin(GL_QUADS); {
				//upper hemisphere
				glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
				glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
				glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
				glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
				//lower hemisphere
				glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
				glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
				glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
				glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
			}glEnd();
		}
	}
}


//void drawSS()
//{
//	glColor3f(1, 0, 0);
//	drawSquare(20);
//
//	glRotatef(angle, 0, 0, 1);
//	glTranslatef(110, 0, 0);
//	glRotatef(2 * angle, 0, 0, 1);
//	glColor3f(0, 1, 0);
//	drawSquare(15);
//
//	glPushMatrix();
//	{
//		glRotatef(angle, 0, 0, 1);
//		glTranslatef(60, 0, 0);
//		glRotatef(2 * angle, 0, 0, 1);
//		glColor3f(0, 0, 1);
//		drawSquare(10);
//	}
//	glPopMatrix();
//
//	glRotatef(3 * angle, 0, 0, 1);
//	glTranslatef(40, 0, 0);
//	glRotatef(4 * angle, 0, 0, 1);
//	glColor3f(1, 1, 0);
//	drawSquare(5);
//}




double determinant(double matrix[3][3], int n) {
	double det = 0;
	double submatrix[3][3];
	if (n == 2)
		return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
	else {
		for (int x = 0; x < n; x++) {
			int subi = 0;
			for (int i = 1; i < n; i++) {
				int subj = 0;
				for (int j = 0; j < n; j++) {
					if (j == x)
						continue;
					submatrix[subi][subj] = matrix[i][j];
					subj++;
				}
				subi++;
			}
			det = det + (pow(-1, x) * matrix[0][x] * determinant(submatrix, n - 1));
		}
	}
	return det;
}

double triangle_t(Triangle t1,  point ray_dir, point ray_start) {
	matrix a = matrix::make_identity(3);
	a.values[0][0] = t1.a.x - t1.b.x;
	a.values[0][1] = t1.a.x - t1.c.x;
	a.values[0][2] = ray_dir.x;
	a.values[1][0] = t1.a.y - t1.b.y;
	a.values[1][1] = t1.a.y - t1.c.y;
	a.values[1][2] = ray_dir.y;
	a.values[2][0] = t1.a.z - t1.b.z;
	a.values[2][1] = t1.a.z - t1.c.z;
	a.values[2][2] = ray_dir.z;
	double det_a = determinant(a.values, 3);
	matrix b = matrix::make_identity(3);
	b.values[0][0] = t1.a.x - t1.b.x;
	b.values[0][1] = t1.a.x - t1.c.x;
	b.values[0][2] = t1.a.x - ray_start.x;
	b.values[1][0] = t1.a.y - t1.b.y;
	b.values[1][1] = t1.a.y - t1.c.y;
	b.values[1][2] = t1.a.y - ray_start.y;
	b.values[2][0] = t1.a.z - t1.b.z;
	b.values[2][1] = t1.a.z - t1.c.z;
	b.values[2][2] = t1.a.z - ray_start.z;
	double det_b = determinant(b.values, 3);
	double int_t = det_b / det_a;
	return int_t;
}

pair<int, int> intersect(point ray_dir, point ray_start) {
	double t = 10000;
	int color_indicator = -1;
	int obj_pos = -1;
	for (int i = 0; i < objects.size(); i++) {
		if (objects[i].type == "sphere") {
			ray_dir.normalize();
			double a = 1;
			double b = point:: dot(ray_dir, ray_start)* 2;
			double c = point::dot(ray_start, ray_start) - pow(objects[i].sphere_radius, 2);
			double d = sqrt(pow(b, 2) - 4 * a*c);
			double temp = (-b + d) / 2 * a;
			if (temp < 0) {
				temp = (-b - d) / 2 * a;
			}
			if (temp > 0 && temp < t ) {
				t = temp;
				obj_pos = i;
			}
		}
		else if (objects[i].type == "pyramid") {
			Triangle t1 = objects[i].t1;
			Triangle t2 = objects[i].t2;
			Triangle t3 = objects[i].t3;
			Triangle t4 = objects[i].t4;
			Square s = objects[i].s;
			double int_t = triangle_t(t1, ray_dir, ray_start);
			if (int_t > 0 && int_t < t) {
				t = int_t;
				obj_pos = i;
			}
			int_t = triangle_t(t2, ray_dir, ray_start);
			if (int_t > 0 && int_t < t) {
				t = int_t;
				obj_pos = i;
			}
			int_t = triangle_t(t3, ray_dir, ray_start);
			if (int_t > 0 && int_t < t) {
				t = int_t;
				obj_pos = i;
			}
			int_t = triangle_t(t4, ray_dir, ray_start);
			if (int_t > 0 && int_t < t) {
				t = int_t;
				obj_pos = i;
			}
			//pyramid square
			int_t = (s.a.z - ray_start.z) / ray_dir.z;
			if (int_t > 0 && int_t < t) {
				point intersection(ray_start.x + int_t * ray_dir.x, ray_start.y + t * ray_dir.y, ray_start.z + t * ray_dir.z);

				if (intersection.x >= s.a.x && intersection.x <= s.c.x && intersection.y >= s.a.y && intersection.y <= s.c.y) {
					t = int_t;
					obj_pos = i;
				}
			}
		}
	}
	//checker_board
	double int_t = -ray_start.z / ray_dir.z;
	if (int_t > 0 && int_t < t) {
		point bottom_left_c_point(-2475, -2475, 0);
		point intersection_c(ray_start.x+int_t*ray_dir.x, ray_start.y + t*ray_dir.y, ray_start.z+t*ray_dir.z);
		int i = floor(intersection_c.x / checker_board_width);
		int j = floor(intersection_c.y / checker_board_width);
		if ((abs(i) + abs(j)) % 2)color_indicator = 1;
		else color_indicator = 0;
	}
	return make_pair(color_indicator, obj_pos);
}

void grid_make() {
	double half_width = number_of_pixels / 2;
	//double incr = 0;
	if (!number_of_pixels % 2) {
		half_width -= 0.5;
	}
	for (int i = 0; i < number_of_pixels; i++) {
		for (int j = 0; j < number_of_pixels; j++) {
			point p(screen_mid_point.x + (i - half_width)*du*u.x + (j - half_width)*dr*r.x,
				screen_mid_point.y + (i - half_width)*du*u.y + (j - half_width)*dr*r.y,
				screen_mid_point.z + (i - half_width)*du*u.z + (j - half_width)*dr*r.z);
			points[number_of_pixels - 1 - i][j] = p;
		}
	}

	/*for (auto v: points) {
		for (auto p: v) {
			p.print();
		}
		cout << endl;
	}*/
}
void setColor() {
	Color background(0, 0, 0);
	Color** pixels = new Color*[number_of_pixels];
	for (int i = 0; i < number_of_pixels; i++) {
		pixels[i] = new Color[number_of_pixels];
		for (int j = 0; j < number_of_pixels; j++) {
			pixels[i][j] = background;
		}
	}

	for (int i = 0; i <(int) points.size(); i++) {
		for (int j = 0; j < (int)points.size(); j++) {
			point ray_dir(points[i][j].x - pos.x, points[i][j].y - pos.y, points[i][j].z - pos.z);
			pair<double, int> inters = intersect(ray_dir, points[i][j]);
			if (inters.first == -1) {
				//t[i][j] = inters.first;
				Color c(objects[inters.second].r, objects[inters.second].g, objects[inters.second].b);
				pixels[i][j] = c;
			}
			else {
				if (inters.first == 1) {
					Color c(1, 1, 1);
					pixels[i][j] = c;
				}
				else {
					Color c(0, 0, 0);
					pixels[i][j] = c;
				}
			}
		}

	}
	bitmap_image image(number_of_pixels, number_of_pixels);
	//bitmap_image image(number_of_pixels, number_of_pixels);
	for (int x = 0; x < number_of_pixels; x++) {
		for (int y = 0; y < number_of_pixels; y++) {
			image.set_pixel(x, y, pixels[x][y].r, pixels[x][y].g, pixels[x][y].b);
		}
	}
	image.save_image("out.bmp");

}

void keyboardListener(unsigned char key, int x, int y) {
	switch (key) {
	case '0':
		grid_make();
		setColor();
		break;
	case '1':
		temp = getCrossProd(u, l);
		l = rotateVec(l, temp);
		temp = getCrossProd(u, r);
		r = rotateVec(r, temp);
		break;

	case '2':
		temp = getCrossProd(l, u);
		l = rotateVec(l, temp);
		temp = getCrossProd(r, u);
		r = rotateVec(r, temp);
		break;

	case '3':
		temp = getCrossProd(r, l);
		l = rotateVec(l, temp);
		temp = getCrossProd(r, u);
		u = rotateVec(u, temp);
		break;

	case '4':
		temp = getCrossProd(l, r);
		l = rotateVec(l, temp);
		temp = getCrossProd(u, r);
		u = rotateVec(u, temp);
		break;

	case '5': //counterclock
		temp = getCrossProd(l, u);
		u = rotateVec(u, temp);
		temp = getCrossProd(l, r);
		r = rotateVec(r, temp);
		break;

	case '6':
		temp = getCrossProd(u, l);
		u = rotateVec(u, temp);
		temp = getCrossProd(r, l);
		r = rotateVec(r, temp);
		break;


	default:
		break;
	}
}


void specialKeyListener(int key, int x, int y) {
	switch (key) {
	case GLUT_KEY_UP:		// up arrow key
		pos.x += l.x;
		pos.y += l.y;
		pos.z += l.z;
		break;
	case GLUT_KEY_DOWN:		//down arrow key
		pos.x -= l.x;
		pos.y -= l.y;
		pos.z -= l.z;
		break;


	case GLUT_KEY_RIGHT:
		pos.x += r.x;
		pos.y += r.y;
		pos.z += r.z;
		break;
	case GLUT_KEY_LEFT:
		pos.x -= r.x;
		pos.y -= r.y;
		pos.z -= r.z;
		break;

	case GLUT_KEY_PAGE_UP:
		pos.x += u.x;
		pos.y += u.y;
		pos.z += u.z;
		break;
	case GLUT_KEY_PAGE_DOWN:
		pos.x -= u.x;
		pos.y -= u.y;
		pos.z -= u.z;
		break;

	case GLUT_KEY_INSERT:
		drawgrid = 1 - drawgrid;
		break;

	case GLUT_KEY_HOME:
		break;
	case GLUT_KEY_END:
		break;

	default:
		break;
	}
}


void mouseListener(int button, int state, int x, int y) {	//x, y is the x-y of the screen (2D)
	switch (button) {
	case GLUT_LEFT_BUTTON:
		if (state == GLUT_DOWN) {		// 2 times?? in ONE click? -- solution is checking DOWN or UP
			drawaxes = 1 - drawaxes;
		}
		break;

	case GLUT_RIGHT_BUTTON:
		//........
		break;

	case GLUT_MIDDLE_BUTTON:
		//........
		break;

	default:
		break;
	}
}



void display() {

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0, 0, 0, 0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);
	//gluLookAt(0, 0, 200, 0, 0, 0, 0, 1, 0);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	//drawAxes();
	//drawGrid();

	//glColor3f(1,0,0);

	//drawSquare(10);
	//int repeat = -10;
	
	
	
	//drawSS();

	//drawCircle(30,24);

	//drawCone(20,50,24);

	//drawSphere(30,24,20);

	/*for (int i = 0; i < (int)spheres.size(); i++) {
		glPushMatrix(); {
			glTranslatef(spheres[i].sphere_centre_x, spheres[i].sphere_centre_y, spheres[i].sphere_centre_z);
			glColor3f(spheres[i].sphere_r, spheres[i].sphere_g, spheres[i].sphere_b);
			drawSphere(spheres[i].sphere_radius, 24, 20);
		}
		glPopMatrix();
	}

	for (int i = 0; i < (int)pyramids.size(); i++) {
		glPushMatrix(); {
			glColor3f(pyramids[i].pyramid_r, pyramids[i].pyramid_g, pyramids[i].pyramid_b);
			drawPyramid(pyramids[i].width_p, pyramids[i].height_p, pyramids[i].lowest_x, pyramids[i].lowest_y, pyramids[i].lowest_z);
		}
		glPopMatrix();
	}*/

	for (int i = 0; i < (int)objects.size(); i++) {
		glPushMatrix(); {
			if (objects[i].type == "sphere") {
				glTranslatef(objects[i].sphere_centre_x, objects[i].sphere_centre_y, objects[i].sphere_centre_z);
				glColor3f(objects[i].r, objects[i].g, objects[i].b);
				drawSphere(objects[i].sphere_radius, 24, 20);
			}
			else if (objects[i].type == "pyramid") {
				glColor3f(objects[i].r, objects[i].g, objects[i].b);
				drawPyramid(objects[i].width_p, objects[i].height_p, objects[i].lowest_x, objects[i].lowest_y, objects[i].lowest_z);
			}
		}glPopMatrix();
		
	}

	//glPushMatrix(); {
		Color c(0, 0, 0);
		for (int i = -50; i <= 50; i++) {
			for (int j = -50; j <= 50; j++) {
				glPushMatrix(); {
					glColor3f(c.r, c.g, c.b);
					//glTranslatef(j*checker_board_width, i*checker_board_width, 0);
					point p(-checker_board_width*(i + 0.5), -checker_board_width*(j + 0.5), 0);
					//repeat -= 1;
					drawSquare(p);
					c.r = 1 - c.r;
					c.g = 1 - c.g;
					c.b = 1 - c.b;
				}glPopMatrix();
			}
		}
	//}
	//glPopMatrix();
	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}



void animate() {
	//angle += 0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init() {
	//codes for initialization
	drawgrid = 0;
	drawaxes = 1;
	cameraHeight = 150.0;
	cameraAngle = 0.05;
	angle = 0;

	inFile >> near_p >> far_p >> fovy >> aspectratio;
	inFile >> level_of_recursion >> number_of_pixels;
	inFile >> checker_board_width;
	inFile >> ambient_c >> diffuse_c >> reflection_c;
	inFile >> number_of_objects;
	screen_height = 2 * near_p * tan((fovy*deg_to_rad)/2);
	double fovx = fovy * aspectratio;
	screen_width = 2 * near_p * tan((fovx*deg_to_rad) / 2);
	screen_mid_point = pos + l * near_p;
	dr = screen_width / number_of_pixels;
	du = screen_height / number_of_pixels;
	points.resize(number_of_pixels, vector<point>(number_of_pixels));
	string command;
	for (int i = 0; i < number_of_objects; i++) {
		inFile >> command;
		if (command == "sphere") {
			inFile >> sphere_centre_x >> sphere_centre_y >> sphere_centre_z;
			inFile >> sphere_radius;
			inFile >> sphere_r >> sphere_g >> sphere_b;
			inFile >> ambient_s >> diffuse_s >> specular_s >> reflection_s >> shininess_s;
			object sphere("sphere", sphere_centre_x, sphere_centre_y, sphere_centre_z,
			sphere_radius,
			sphere_r, sphere_g, sphere_b,
			ambient_s, diffuse_s, specular_s, reflection_s, shininess_s);
			//spheres.push_back(sphere);
			objects.push_back(sphere);
		}
		else if (command == "pyramid") {
			inFile >> lowest_x >> lowest_y >> lowest_z;
			inFile >> width_p >> height_p;
			inFile >> pyramid_r >> pyramid_g >> pyramid_b;
			inFile >> ambient_p >> diffuse_p >> specular_p >> reflection_p >> shininess_p;

			point apex = { lowest_x, lowest_y, lowest_z + height_p };
			double a = width_p / 2;
			point b(lowest_x - a, lowest_y - a, lowest_z);
			point c(lowest_x + a, lowest_y - a, lowest_z);
			Triangle t1(apex, b, c);
			point b1(lowest_x + a, lowest_y - a, lowest_z);
			point c1(lowest_x + a, lowest_y + a, lowest_z);
			Triangle t2(apex, b1, c1);
			point b2(lowest_x + a, lowest_y + a, lowest_z);
			point c2(lowest_x - a, lowest_y + a, lowest_z);
			Triangle t3(apex, b2, c2);
			point b3(lowest_x - a, lowest_y + a, lowest_z);
			point c3(lowest_x - a, lowest_y - a, lowest_z);
			Triangle t4(apex, b3, c3);
			point a1(lowest_x - a, lowest_y - a, lowest_z);
			point b4(lowest_x + a, lowest_y - a, lowest_z);
			point c4(lowest_x + a, lowest_y + a, lowest_z);
			point d4(lowest_x - a, lowest_y + a, lowest_z);
			Square s(a1, b4, c4, d4);
			object pyramid("pyramid", lowest_x, lowest_y, lowest_z,
			width_p, height_p,
			pyramid_r, pyramid_g, pyramid_b,
			ambient_p, diffuse_p, specular_p, reflection_p, shininess_p, t1, t2, t3, t4, s);
			//pyramids.push_back(pyramid);
			objects.push_back(pyramid);
		}
	}

	

	//clear the screen
	glClearColor(0, 0, 0, 0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80, 1, 1, 1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv) {
	glutInit(&argc, argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();
	//grid_make();
	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
