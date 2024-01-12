#include <cmath>

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

static inline double sqr(double x) { return x * x; }

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}
	double& operator[](int i) { return coord[i]; }
	double operator[](int i) const { return coord[i]; }

	Vector& operator+=(const Vector& v) {
		coord[0] += v[0];
		coord[1] += v[1];
		coord[2] += v[2];
		return *this;
	}

	//already squared
	double norm2() const {
		return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
	}

	double coord[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector& a, double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator/(const Vector& a, double b) {
	return Vector(a[0]/b, a[1]/b, a[2]/b);
}

double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


class ray{
public:
	explicit ray(const Vector& center, const Vector& direction){
		O = center;
		// Error if norm(u) = 0
		u = direction / sqrt(direction.norm2());
	}
	Vector O;
	Vector u;

};

class sphere{
public:
	explicit sphere(double x,double y, double z, double r){
		C = Vector(x,y,z);
		R = r;
	}

	sphere(const Vector& v, double r){
		C = v;
		R = r;
	}

	bool intersect(const ray& r){
		double b = 2*dot(r.u , r.O - C) ;
		double c = (r.O - C).norm2() - R*R;
		double delta = pow(b,2) - 4 * c;
		double t = 0;
		double t2 = 0;

		if (delta < 0){
			return false;
		}
		else if (delta == 0){
			//t = -b/2;
			//return t>= 0;
			return b < 0;
		}
		else {
			//t, t2 = (-b-sqr(delta))/(2*a),(-b+sqr(delta))/(2*a);
			t , t2 = -b-sqrt(delta), -b+sqrt(delta);
			return t >= 0 || t2 >= 0 ;
		}

	};

	Vector C;
	double R;
};




int main() {
	int W = 512;
	int H = 512;

	Vector center(0, 0, 55);

	sphere boule(Vector(0,0,0),10);

	std::vector<unsigned char> image(W*H * 3, 0);


	ray r = ray(Vector(0,0,0),Vector(1,0,0));
	double x;
	double y;
	double z;
	double alpha = 2 * M_PI / 360 * 60;
#pragma omp parallel for
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			x = j - W/2 + 0.5;
			y = -i + W/2 - 0.5;
			z = -W/(2*tan(alpha/2));
			
			r = ray(center, Vector(x,y,z));

			if (boule.intersect(r)){
				image[(i*W + j) * 3 + 0] = 255;   // RED
				image[(i*W + j) * 3 + 1] = 255;  // GREEN
				image[(i*W + j) * 3 + 2] = 255 ;  // BLUE
			}
			else {
				image[(i*W + j) * 3 + 0] = 0;   // RED
				image[(i*W + j) * 3 + 1] = 0;  // GREEN
				image[(i*W + j) * 3 + 2] = 0 ; // BLUE
			}
			
			//Vector v(j / (double)W - 0.5, i / (double)H - 0.5, 0.);
			//double gaussianVal = exp(-(v-center).norm2()/(2*sqr(0.2)));

			
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}