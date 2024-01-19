#include <cmath>
#include <iostream>

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

	Vector normalize() const {
		double norm = sqrt(this->norm2());
		return Vector(coord[0]/norm,coord[1]/norm,coord[2]/norm);
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
		u = direction.normalize();
	}
	Vector O;
	Vector u;

};

class sphere{
public:
	explicit sphere(double x,double y, double z, double r,bool m){
		C = Vector(x,y,z);
		R = r;
		is_mirror = m;
	}

	sphere(const Vector& v, double r,bool m){
		C = v;
		R = r;
		is_mirror = m;
	}

	sphere(const Vector& v,double r, const Vector& a,bool m = false){
		C = v;
		R = r;
		albedo = a;
		is_mirror = m;
	}

	bool intersect(const ray& r, Vector &P, Vector &N, double &t){
		// renvoyer la normale : CP normalisé
		// renvoyer P et N: on pourrait créer une classe intersection.
		double b = 2*dot(r.u , r.O - C) ;
		double c = (r.O - C).norm2() - R*R;
		double delta = pow(b,2) - 4 * c;
		double t1 = 0;
		double t2 = 0;

		if (delta < 0){
			return false;
		}
		else {
			t1 = (-b-sqrt(delta))/2;
			t2 =  (-b+sqrt(delta))/2;
			//std::cout<<"t1:"<<t1<<" t2: "<<t2<<" ";
			if (t1 < 0 && t2 < 0) {return false;}

			// if at least one solution, compute N and P
			if (t1 > 0 && t2 > 0) {t = std::min(t1,t2);}
			else {t = std::max(t1,t2);}
			P = r.O + t * r.u;
			N = (P - C).normalize();
			return true;
		}

	};

	Vector C;
	double R;
	Vector albedo;
	bool is_mirror;
};


class scene{
public:
	scene(std::vector<sphere> sl){
		for (int i =0;i<sl.size();i++){
			sphere_list.push_back(sl[i]);
		}
	}
	
	//intersection with the closest sphere
	bool intersect(const ray& r,  Vector &Pscene, Vector &Nscene,double &tscene, int &index){
		double min_t = -1.0;
		int min_index = -1;
		for (int i=0; i<sphere_list.size(); i++){
			Vector P;
			Vector N;
			double t;
			if (sphere_list[i].intersect(r,P,N,t)){
				//std::cout<<" match! ";
				//t should be = 0 since its false
				if ( (t < min_t && min_t != -1.0 ) || min_t == -1 ){
					//new contestant, we update all
					min_t = t;
					min_index = i;
				}
			}
			//std::cout<<i<<" "<<t<<" "<<min_t<<" "<<min_index<<"\n";
		}
		//if not found, false
		if (min_index == -1) {return false;}

		//else update return values
		tscene =  min_t;
		index = min_index;
		Pscene = r.O + tscene * r.u;
		Nscene =  (Pscene - sphere_list[index].C).normalize();
		return true;
		
		// s'assurer que P, N,t et index correspondent bien au minimum de distance trouvé
	}

	std::vector<sphere> sphere_list;
};

Vector get_color(const ray& r, int iteration){

}


int main() {
	int W = 512;
	int H = 512;
	std::vector<unsigned char> image(W*H * 3, 0);

	Vector center(0, 0, 55);

	sphere boule(Vector(0,0,0),10,Vector(0.3,0.4,0.9),true);
	sphere boule1(Vector(0,0,-1000),940,Vector(0.1,0.8,0.2),false);
	sphere boule2(Vector(0,1000,0),940,Vector(0.85,0.2,0.2),false);
	sphere boule3(Vector(0,0,1000),940,Vector(0.1,0.8,0.2),false);
	sphere boule4(Vector(0,-1000,0),990,Vector(0.1,0.1,0.95),false);


	std::vector<sphere> sphere_list = {boule,boule1,boule2,boule3,boule4};
	scene scene_1(sphere_list);
	//update l'algo pour que la scène marche (boucle for)


	Vector light(-10,20,40);
	//Vector albedo(0.3,0.4,0.9);
	double I = 2E10;
	double gamma = 0.454;



	ray r = ray(Vector(0,0,0),Vector(1,0,0));
	double x;
	double y;
	double z;
	double alpha = 2 * M_PI / 360 * 60;

#pragma omp parallel for
	// pour chaque pixel, on crée un rayon 
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector P;
			Vector N;
			double t;
			int sphere_index;
			Vector albedo;

			x = j - W/2 + 0.5;
			y = -i + W/2 - 0.5;
			z = -W/(2*tan(alpha/2));
			
			r = ray(center, Vector(x,y,z));
			// si le rayon intersecte la scène
			if (scene_1.intersect(r,P,N,t,sphere_index)){
				// On regarde d'abord si le point d'intersection voit la lumière
				// Création d'un rayon de P vers L
				ray r2 = ray(P + 0.0001 * N,light-P);
				Vector P2;
				Vector N2;
				double t2;
				int sphere_index2;
				if (scene_1.intersect(r2,P2,N2,t2,sphere_index2) and sqr(t2) < (light- P).norm2()){
					image[(i*W + j) * 3 + 0] = 0;
					image[(i*W + j) * 3 + 1] = 0;
					image[(i*W + j) * 3 + 2] = 0;
				}
				else {
					// contient aussi le cas où pas d'inersection entre le 2ème rayon et la scène
					// intersection, il faut vérifier si l'intersection est avant la lumière
					Vector albedo = scene_1.sphere_list[sphere_index].albedo;
					Vector vect_lum = I * albedo / M_PI * std::max(dot((light-P).normalize(),N),0.) / ( 4*M_PI* ((P-light).norm2()) );

					image[(i*W + j) * 3 + 0] = std::min(pow(vect_lum[0],gamma),255.0) ; // RED
					image[(i*W + j) * 3 + 1] = std::min(pow(vect_lum[1],gamma),255.0) ; // GREEN
					image[(i*W + j) * 3 + 2] = std::min(pow(vect_lum[2],gamma),255.0) ; // BLUE
				}


			}
			else {
				image[(i*W + j) * 3 + 0] = 0;  // RED
				image[(i*W + j) * 3 + 1] = 0;  // GREEN
				image[(i*W + j) * 3 + 2] = 0 ; // BLUE
			}
						
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}