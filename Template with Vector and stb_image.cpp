#include <cmath>
#include <iostream>
#include <random>

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
Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1], -a[2]);
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

	sphere(const Vector& v,double r, const Vector& a,bool m = false,bool t = false,double refraction_index = 1.5){
		// vérifier qu'on ne puisse pas avoir m ET t en meme temps
		C = v;
		R = r;
		albedo = a;
		is_mirror = m;
		is_transparent = t;
		n = refraction_index;
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
	bool is_transparent;
	double n;
};


class scene{
public:
	scene(std::vector<sphere> sl,Vector l,double i){
		for (int i =0;i<sl.size();i++){
			sphere_list.push_back(sl[i]);
		}
		L = l;
		I = i;
	}
	
	//intersection with the closest sphere
	bool intersect(const ray& r,  Vector &Pscene, Vector &Nscene,double &tscene, int &index, bool count_transp){
		double min_t = -1.0;
		int min_index = -1;
		for (int i=0; i<sphere_list.size(); i++){
			Vector P;
			Vector N;
			double t;
			// On ne compte PAS les intersections avec les objets transparents!!
			// !count & !t & i   ||   count & i
			// i & ( !count & !t || count)
			if ( sphere_list[i].intersect(r,P,N,t) && ( (!sphere_list[i].is_transparent && !count_transp) || count_transp )  ){
				//std::cout<<" match! ";
				//t should be = 0 since its false
				if ( (t < min_t && min_t != -1.0 ) || min_t == -1 ){
					//new contestant, we update all
					min_t = t;
					min_index = i;
				}
			}
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

	Vector get_color(const ray& r, int ray_depth){
		//Calcul de l'intersection
		Vector P;
		Vector N;
		double t;
		int sphere_index;
		Vector albedo;


		// Si trop de récursion : Noir
		if (ray_depth < 0){
			return Vector(0.,0.,0.);
		}

		// Si intersection : disjonction selon la nature de la surface
		// true car on compte aussi les réflexions
		if (this->intersect(r,P,N,t,sphere_index,true)){
			
			// Si intersection miroir : Récursion
			if (sphere_list[sphere_index].is_mirror){
				// rayon réfléchi
				// On devrait faire gaffe à la normale : Si on vient de l'intérieur d'un cercle on devrait utiliser -N
				ray reflected(P + 0.0001 * N, r.u - 2 * dot(r.u,N) * N); //Centre décalé, et on inverse la normale
				return this->get_color(reflected,ray_depth -1);
			}

			// Si intersection transparents
			if (sphere_list[sphere_index].is_transparent){
				double n1;
				double n2;
				double orientation = dot(r.u,N); //rayon incident vs normale sortante
				Vector N_reflexion;
				if (orientation < 0){ // au pif pour le < ou <=
					//on entre dans la boule
					n1 = 1;
					n2 = sphere_list[sphere_index].n;
					N_reflexion = N;
				}
				else{
					//on sort de la boule, pour etre dans le schéma de snell-descartes on utilise -N
					n1 = sphere_list[sphere_index].n;
					n2 = 1;
					N_reflexion = -N;
				}
				
				if (sqr(n1/n2) * (1-sqr(dot(r.u,N_reflexion))) > 1){
					//si réflexion totale, on renvoie le rayon réfléchi. Jamais pu observer ça.
					//On utilise la normale "signée" pour
					ray reflected(P + 0.0001 * N_reflexion, r.u - 2 * dot(r.u,N_reflexion) * N_reflexion);
					return this->get_color(reflected,ray_depth - 1);
				}
				else
				{
					//sinon on renvoie le rayon transmis
					//t_t = n1/n2 (i - <i,n>n)
					//t_n = - sqrt(1 – (n1/n2)² * ( 1 – <i,n>²  ) . n
					Vector transmitted_tang = n1/n2 * ( r.u - dot(r.u, N_reflexion) * N_reflexion);
					Vector transmitted_norm = -sqrt(1 - sqr(n1/n2) * (1-sqr(dot(r.u,N_reflexion)))) * N_reflexion;

					ray transmitted(P - 0.0001 * N_reflexion,transmitted_norm + transmitted_tang);
					return this->get_color(transmitted,ray_depth - 1);
					
				}
				

			}

			// Si intersection diffuse : Calcul classique

			// On regarde d'abord si le point d'intersection voit la lumière
			// Création d'un rayon de P vers L
			ray r2 = ray(P + 0.0001 * N,L-P);
			Vector P2;
			Vector N2;
			double t2;
			int sphere_index2;
			//ici on ne compte pas les transparentes
			if (this->intersect(r2,P2,N2,t2,sphere_index2,false) and sqr(t2) < (L- P).norm2()){
				return Vector(0.,0.,0.);
			}
			else {
				// contient aussi le cas où pas d'inersection entre le 2ème rayon et la scène
				// intersection, il faut vérifier si l'intersection est avant la lumière
				Vector albedo = sphere_list[sphere_index].albedo;
				Vector vect_lum = I * albedo / M_PI * std::max(dot((L-P).normalize(),N),0.) / ( 4*M_PI* ((P-L).norm2()) );

				return vect_lum;
			}


		}
		// Si pas d'intersection, rien (mais ne devrait pas arriver, c'est une scène fermée)
		else {
			//std::cout<<"Warning: No intersection for ray";
			return Vector(0.,0.,0.);
		}
			
	}


	std::vector<sphere> sphere_list;
	Vector L;
	double I;
};



int main() {
	// Générateur de nbres aléatoires
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);


	int W = 2048;
	int H = 2048;
	std::vector<unsigned char> image(W*H * 3, 0);

	Vector center(0, 0, 55);

	sphere boule(Vector(0,0,0),10,Vector(0.3,0.4,0.9),true,false,1.333);
	
	sphere background(Vector(0,0,-1000),940,Vector(0.1,0.8,0.2),false,false,1.333);
	sphere topwall(Vector(0,1000,0),940,Vector(0.85,0.2,0.2),false,false,1.333);
	sphere foreground(Vector(0,0,1000),940,Vector(0.1,0.8,0.2),false,false,1.333);
	sphere bottomwall(Vector(0,-1000,0),990,Vector(0.1,0.1,0.95),false,false,1.333);
	sphere rightwall(Vector(1000,0,0),940,Vector(0.9,0.9,0.1),false,false,1.333);
	sphere leftwall(Vector(-1000,0,0),940,Vector(0.1,0.9,0.9),false,false,1.333);

	sphere boule7(Vector(10,10,10),5,Vector(0.3,0.4,0.9),false,true,1.333);


	std::vector<sphere> sphere_list;
	sphere_list.push_back(boule);
	sphere_list.push_back(background);
	sphere_list.push_back(topwall);
	sphere_list.push_back(foreground);
	sphere_list.push_back(bottomwall);
	sphere_list.push_back(rightwall);
	sphere_list.push_back(leftwall);
	//sphere_list.push_back(boule7);

	Vector light(-10,20,40);
	double I = 2E9;

	scene scene_1(sphere_list,light,I);
	double gamma = 0.454;
	//&& !sphere_list[i].is_transparent
	//condition dans intersect pour compter les sphères transp : pour savoir si on est à l'ombre on en a pas besoin, sinon oui



	ray r = ray(Vector(0,0,0),Vector(1,0,0));
	double x;
	double y;
	double z;
	double alpha = 2 * M_PI / 360 * 60;

#pragma omp parallel for
	// pour chaque pixel, on crée un rayon 
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			
			// antialiasing : On remplace + 0.5 par + [0;1]
			x = j - W/2 + distribution(generator);
			y = -i + W/2 - distribution(generator);
			z = -W/(2*tan(alpha/2));
			
			r = ray(center, Vector(x,y,z));
			Vector color(scene_1.get_color(r,15));

			image[(i*W + j) * 3 + 0] = std::min(pow(color[0],gamma),255.0);  // RED
			image[(i*W + j) * 3 + 1] = std::min(pow(color[1],gamma),255.0);  // GREEN
			image[(i*W + j) * 3 + 2] = std::min(pow(color[2],gamma),255.0);  // BLUE
						
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}