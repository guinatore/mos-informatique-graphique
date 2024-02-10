#include <cmath>
#include <ctime>
#include <iostream>
#include <random>

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Générateur de nbres aléatoires
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);

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
Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, double b) {
	return Vector(a[0]/b, a[1]/b, a[2]/b);
}

double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector& a, const Vector& b){
	return Vector(a[1]*b[2]-b[1]*a[2],a[2]*b[0]-b[2]*a[0],a[0]*b[1]-b[0]*a[1]);
}

Vector random_cos(const Vector& N){
	// N doit etre normalisé !!
	double r1 = distribution(generator);
	double r2 = distribution(generator);

	double x = cos(2*M_PI*r1);
	double y = sin(2*M_PI*r2);
	double z = sqrt(r2);

	// Selon la plus petite valeur de N on construit T1 orthogonal à N
	Vector T1;
	if(std::abs(N[0])<= std::abs(N[1]) && std::abs(N[0])<= std::abs(N[2]) ){
		T1 = Vector(0,-N[2],N[1]); 
	}
	else if(std::abs(N[1])<= std::abs(N[2]) && std::abs(N[1])<= std::abs(N[0]) ){
		T1 = Vector(-N[2],0,N[0]); 
	}
	else {
		T1 = Vector(-N[1],N[0],0); 
	}

	T1 = T1.normalize();

	Vector T2 = cross(N,T1); // Pas besoin de le normaliser
	
	return z*N + x*T1 + y*T2;
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
		//centre, rayon, albedo, miroir,transparent, refraction
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
	scene(std::vector<sphere> sl,double i){
		for (int i =0;i<sl.size();i++){
			sphere_list.push_back(sl[i]);
		}
		I = i;
	}
	
	//intersection with the closest sphere
	bool intersect(const ray& r,  Vector &Pscene, Vector &Nscene,double &tscene, int &index, bool count_transp){
		// count_transp à true si une intersection transparente compte comme une intersection
		// count_transp à false si une intersection transparente est ignorée
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
		
	}

	Vector get_color(const ray& r, int ray_depth,bool was_diffuse_interaction){
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
				return this->get_color(reflected,ray_depth -1,false);
			}

			// Si intersection transparente
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
					return this->get_color(reflected,ray_depth - 1,false);
				}
				else
				{
					//sinon : on renvoie le rayon transmis ou réfléchi, selon la pondération des coefficients de transmission
					double k0 = sqr(n1-n2)/sqr(n1+n2);
					double R = k0 + (1-k0)*pow(1-abs(dot(N,r.u)),5); //proba de réflexion

					if (distribution(generator) < R){
						// Réflexion //check la formule
						ray reflected(P + 0.0001 * N_reflexion, r.u - 2 * dot(r.u,N_reflexion) * N_reflexion);
						return this->get_color(reflected,ray_depth - 1,false);

					}
					else{
						// Transmission
						//t_t = n1/n2 (i - <i,n>n)
						//t_n = - sqrt(1 – (n1/n2)² * ( 1 – <i,n>²  ) . n
						Vector transmitted_tang = n1/n2 * ( r.u - dot(r.u, N_reflexion) * N_reflexion);
						Vector transmitted_norm = -sqrt(1 - sqr(n1/n2) * (1-sqr(dot(r.u,N_reflexion)))) * N_reflexion;

						ray transmitted(P - 0.0001 * N_reflexion,transmitted_norm + transmitted_tang);
						return this->get_color(transmitted,ray_depth - 1,false);

					}

					
					
				}
				

			}

			// Si intersection diffuse : Calcul classique

			// On regarde d'abord si le point d'intersection voit la lumière
			/*
			// Création d'un rayon de P vers L
			ray r2 = ray(P + 0.0001 * N,L-P);
			Vector P2;
			Vector N2;
			double t2;
			int sphere_index2;
			//ici on ne compte pas les transparentes

			Vector I_direct(0.,0.,0.);
			if (this->intersect(r2,P2,N2,t2,sphere_index2,false) and sqr(t2) < (L- P).norm2()){
				// La lumière est bloquée par une surface non transparente, éclairage direct nul
				I_direct = Vector(0.,0.,0.);
			}
			else {
				// contient aussi le cas où pas d'intersection entre le 2ème rayon et la scène
				// intersection, il faut vérifier si l'intersection est avant la lumière
				Vector albedo = sphere_list[sphere_index].albedo;
				Vector vect_lum = I * albedo / M_PI * std::max(dot((L-P).normalize(),N),0.) / ( 4*M_PI* ((P-L).norm2()) );

				I_direct =  vect_lum;
			}
			//  If lumière = 0: return I/4pi²r²
			*/

			////// Cas de l'intersection avec une matière diffuse
			if (sphere_index == 0){
				// Si interaction diffuse : On renvoie 0 pour éviter de compter 2 fois la source (I_direct et I_indirect)
				if (was_diffuse_interaction){return Vector(0.,0.,0.);}
				// Si intersection avec la lumière : On renvoit la lumière ! (Pondérée par l'élément de surface)
				else {return I*sphere_list[0].albedo /(4 * M_PI * M_PI * sphere_list[0].R * sphere_list[0].R ) ;}
				
			}
			else{
				//// Calcul de l'éclairage direct : On intègre l'éclairage de la sphère source
				Vector I_direct(0.,0.,0.);

				Vector vect_light = (sphere_list[0].C - P).normalize(); // vecteur de P vers L
				// Point sur la lumière: on utilise random_cos qu'on scale.
				Vector N_prime = random_cos(-vect_light);
				Vector P_prime = sphere_list[0].C + sphere_list[0].R * N_prime; 
				// Rayon incident
				Vector wi_direct = P_prime - P;
				wi_direct = wi_direct.normalize();

				// On renvoie un rayon de P vers P'. pour le calcul du terme de visibilité
				bool visibility = false;
				ray r2(P + 0.0001 * N , wi_direct); // Rayon de P vers P'
				double d2 = (P_prime - P).norm2(); // distance PP'²
				Vector P2;
				Vector N2;
				double t2;
				int sphere_index2;
				// Il y aura toujours une intersection(avec la lumière), mais on veut savoir si elle se fait avant la lumière
				if (this->intersect(r2,P2,N2,t2,sphere_index2,true)){
					visibility = sqr(t2 + 0.01) > d2; // P voit la lumière ssi l'intersection est plus loin que la lumière
				}
				else{
					// Au cas où on n'intersecte rien (Ne devrait pas arriver)
					visibility = false;
				}


				// Contribution ou non de l'éclairage direct selon la visibilité. Pas besoin de faire les calculs si pas
				// de visibilité, on laisse I_direct à 0.
				if (visibility){
					double light_intensity = I / (4 * M_PI * M_PI * sphere_list[0].R * sphere_list[0].R); // intensité de la lumière, normalisée par l'élément de surface
					Vector albedo = sphere_list[sphere_index].albedo/M_PI; // albedo de la sphere considérée
					double form_factor = std::max(0.,dot(N,wi_direct))*std::max(0.,dot(N_prime,-wi_direct))/(d2);
					double pdf = std::max(0.,dot(N_prime,-vect_light)) / ( M_PI * sphere_list[0].R * sphere_list[0].R ) ;

					I_direct = light_intensity * albedo * form_factor / pdf;
				}

				//// Calcul de l'éclairage indirect
				Vector I_indirect(0.,0.,0.);
				ray wi_indirect(Vector(P + 0.0001*N),random_cos(N)); 
				// albedo* getcolor(wi) et wi suit cos(theta)/pi
				I_indirect += sphere_list[sphere_index].albedo * get_color(wi_indirect,ray_depth-1, true); // Seul cas où on signale qu'on ne compte pas la lumière
				return I_direct + I_indirect;
			}


		}
		// Si pas d'intersection, rien (mais ne devrait pas arriver, c'est une scène fermée)
		else {
			//std::cout<<"Warning: No intersection for ray";
			return Vector(0.,0.,0.);
		}
			
	}


	std::vector<sphere> sphere_list;
	double I;
};



int main() {
	time_t time_start = time(NULL);

	int W = 2048;
	int H = 2048;
	int ray_count = 1024;
	std::vector<unsigned char> image(W*H * 3, 0);

	Vector center(0, 0, 55);

	// La lumière est la première boule
	sphere light_sphere(Vector(-10,20,40),10,Vector(1.,1.,1.),false,false,1.333);
	sphere boule1(Vector(0,0,0),10,Vector(0.3,0.4,0.9),true,false,1.333); // miroir
	sphere boule2(Vector(-20,0,0),10,Vector(0.3,0.4,0.9),false,true,1.333); //transparente 
	sphere boule3(Vector(20,0,0),10,Vector(0.3,0.4,0.9),false,false,1.333); // solide
	sphere background(Vector(0,0,-1000),940,Vector(0.85,0.1,0.1),false,false,1.333);
	sphere topwall(Vector(0,1000,0),940,Vector(0.1,0.9,0.2),false,false,1.333);
	sphere foreground(Vector(0,0,1000),940,Vector(0.9,0.3,0.1),false,false,1.333);
	sphere bottomwall(Vector(0,-1000,0),990,Vector(0.1,0.1,0.95),false,false,1.333);
	sphere rightwall(Vector(1000,0,0),940,Vector(0.95,0.95,.07),false,false,1.333);
	sphere leftwall(Vector(-1000,0,0),940,Vector(0.07,0.95,0.95),false,false,1.333);


	std::vector<sphere> sphere_list;
	sphere_list.push_back(light_sphere);
	sphere_list.push_back(boule1);
	sphere_list.push_back(boule2);
	sphere_list.push_back(boule3);
	sphere_list.push_back(background);
	sphere_list.push_back(topwall);
	sphere_list.push_back(foreground);
	sphere_list.push_back(bottomwall);
	sphere_list.push_back(rightwall);
	sphere_list.push_back(leftwall);

	double I = 3E9;

	scene scene_1(sphere_list,I);
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
			
			Vector color(0.,0.,0.);
			for (int ray_index = 0;ray_index<ray_count;ray_index++){

				// antialiasing : On remplace + 0.5 par + [0;1]
				double r1 = distribution(generator);
				double r2 = distribution(generator);

				double g1 = sqrt(-2 * log( r1 )) * cos( 2 * M_PI*r2 ) * 1 ;
				double g2 = sqrt(-2 * log( r1 )) * sin( 2 * M_PI*r2 ) * 1 ;


				x = j - W/2 + 0.5 + g1;
				y = -i + W/2 - 0.5 + g2;
				z = -W/(2*tan(alpha/2));
				
				r = ray(center, Vector(x,y,z));
				
				// On additionne la contribution du chemin à l'éclairage.
				color += scene_1.get_color(r,15,false)/ray_count;
			}

			image[(i*W + j) * 3 + 0] = std::min(pow(color[0],gamma),255.0);  // RED
			image[(i*W + j) * 3 + 1] = std::min(pow(color[1],gamma),255.0);  // GREEN
			image[(i*W + j) * 3 + 2] = std::min(pow(color[2],gamma),255.0);  // BLUE
		
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	time_t time_end = time(NULL);
	double time_diff_min = difftime(time_end, time_start)/60;
	std::cout << "Time elapsed : "<< time_diff_min << " min.";

	return 0;
}