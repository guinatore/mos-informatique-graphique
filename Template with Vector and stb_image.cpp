#include <cmath>
#include <ctime>
#include <iostream>
#include <random>
#include <stdio.h>
#include <string>
#include <algorithm>

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Générateur de nbres aléatoires
std::default_random_engine generator(2078369);
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

	double x = cos(2*M_PI*r1)*sqrt(1 - r2);
	double y = sin(2*M_PI*r1)*sqrt(1 - r2);
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

class Geometry{
// Classe abstraite qui permet de factoriser les propriétés de texture des éléments
// mais laisse l'implémentation de intersect aux classes filles (mesh & sphere)
public:
	explicit Geometry(const Vector& a = Vector(0.9,0.8,0.03),bool m = false,bool t = false,double refraction_index = 1.5){
		albedo = a;
		is_mirror = m;
		is_transparent = t;
		n = refraction_index;
	}

	virtual bool intersect(const ray& r, Vector &P, Vector &N, double &t) = 0; 
	
	Vector albedo;
	bool is_mirror;
	bool is_transparent;
	double n;

};


class sphere : public Geometry {
public:
	sphere(const Vector& v, double r, const Vector& a, bool m = false, bool t = false, double refraction_index = 1.5): Geometry(a,m,t,refraction_index){
		//centre, rayon, albedo, miroir,transparent, refraction
		// vérifier qu'on ne puisse pas avoir m ET t en meme temps
		C = v;
		R = r;
	}

	virtual bool intersect(const ray& r, Vector &P, Vector &N, double &t){
		// renvoyer la normale : CP normalisé
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
			if (t1 < 0 && t2 < 0) {return false;}

			// if at least one solution, compute N and P
			if (t1 >= 0 && t2 >= 0) {t = std::min(t1,t2);}
			else {t = std::max(t1,t2);}

			P = r.O + t * r.u;
			N = (P - C).normalize();
			return true;
		}

	};

	Vector C;
	double R;
};


class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};

class BoundingBox{
	public:
	BoundingBox(){}
	BoundingBox(double max_x,double  max_y,double  max_z,double  min_x,double  min_y,double  min_z){
		this->max_x = max_x;
		this->max_y = max_y;
		this->max_z = max_z;
		this->min_x = min_x;
		this->min_y = min_y;
		this->min_z = min_z;
	}

	// Intersection boite/rayon , en prenant en compte le fait que les plans sont finis
	bool intersect(const ray& r) const {
		//Crédits à Baptiste Perreyon pour cette méthode

		Vector invU = Vector(1 / r.u[0], 1 / r.u[1], 1 / r.u[2]);

		// X axis, normal is (1,0,0) so we only keep first coordinate
		double tX1 = (min_x - r.O[0]) * invU[0];
		double tX2 = (max_x - r.O[0]) * invU[0];
		double tXMin = std::min(tX1, tX2);
		double tXMax = std::max(tX1, tX2);

		// Y axis
		double tY1 = (min_y - r.O[1]) * invU[1];
		double tY2 = (max_y - r.O[1]) * invU[1];
		double tYMin = std::min(tY1, tY2);
		double tYMax = std::max(tY1, tY2);

		// Z axis
		double tZ1 = (min_z - r.O[2]) * invU[2];
		double tZ2 = (max_z - r.O[2]) * invU[2];
		double tZMin = std::min(tZ1, tZ2);
		double tZMax = std::max(tZ1, tZ2);

		// intersection des intervalles non vide

		//il faut aussi que t soit positif!!!
		return std::min({ tXMax, tYMax, tZMax }) > std::max({ tXMin, tYMin, tZMin }) && std::max({ tXMin, tYMin, tZMin }) >= 0;
		
	}
	
	double max_x, max_y, max_z, min_x, min_y, min_z;
};

class TriangleMesh : public Geometry  {
public:
    ~TriangleMesh() {}
	TriangleMesh(): Geometry() {}
	//TriangleMesh(const Vector& a,bool m = false,bool t = false,double refraction_index = 1.5): Geometry(a,m,t,refraction_index) {};
    
	virtual bool intersect(const ray& r, Vector &P, Vector &N, double &t){
		if (!bbox.intersect(r)){return false;}

		// On ne renvoie pas l'indice du triangle encore. (Comment fera-t-on pour la couleur?)
		t = -1.0; // Min tracker
		for (int i = 0; i< indices.size() ; i++){ // pour tous les triangles
			Vector A = vertices[indices[i].vtxi];
			Vector e1 = vertices[indices[i].vtxj] - A; // B - A : j-i
			Vector e2 = vertices[indices[i].vtxk] - A;
			Vector temp_N = cross(e1,e2);

			double beta = dot(e2,cross(A-r.O , r.u))/dot(r.u,temp_N);
			double gamma = - dot(e1,cross(A-r.O , r.u))/dot(r.u,temp_N);
			double alpha = 1 - gamma - beta;

			if (0 <= beta && beta <= 1 && 0 <= alpha && alpha <= 1 && 0 <= gamma && gamma <= 1 ){
				// Triangle intercepted. But is it closer?
				double temp_t = dot(A - r.O, temp_N)/dot(r.u, temp_N);

				if (temp_t > 0){
					if ((temp_t < t && t != -1.0 ) || t == -1.0 ){
						//Update return values
						t = temp_t;
						N = temp_N;
						P = r.O + t * r.u;
					}
				}
			}
		}

		if (t == -1.0){return false;} // Si t n'a pas été modifié

		return true; // Sinon 
	}

	void affine(const Vector& translation, double scale){
		for (int i = 0;i<vertices.size();i++){
			vertices[i] = vertices[i] * scale + translation;
		}
	}

	void set_bounding_box(){
		double min_x = (*std::min_element(vertices.begin(), vertices.end(), [](const Vector& v1, const Vector& v2) { return v1[0] < v2[0]; }))[0];
        double min_y = (*std::min_element(vertices.begin(), vertices.end(), [](const Vector& v1, const Vector& v2) { return v1[1] < v2[1]; }))[1];
        double min_z = (*std::min_element(vertices.begin(), vertices.end(), [](const Vector& v1, const Vector& v2) { return v1[2] < v2[2]; }))[2];
        double max_x = (*std::max_element(vertices.begin(), vertices.end(), [](const Vector& v1, const Vector& v2) { return v1[0] < v2[0]; }))[0];
        double max_y = (*std::max_element(vertices.begin(), vertices.end(), [](const Vector& v1, const Vector& v2) { return v1[1] < v2[1]; }))[1];
        double max_z = (*std::max_element(vertices.begin(), vertices.end(), [](const Vector& v1, const Vector& v2) { return v1[2] < v2[2]; }))[2];

        bbox = BoundingBox(max_x,max_y,max_z,min_x,min_y,min_z);
	}

    void readOBJ(const char* obj) {
 
        char matfile[255];
        char grp[255];
 
        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;
 
            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());
 
            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }
 
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
 
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));
 
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
 
                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;
 
                char* consumedline = line + 1;
                int offset;
 
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }
 
                consumedline = consumedline + offset;
 
                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
 
            }
 
        }
        fclose(f);
 
    }
 
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
	BoundingBox bbox;
    
};


class scene{
public:
	scene(std::vector<Geometry*> sl,double i){
		for (int i =0;i<sl.size();i++){
			object_list.push_back(sl[i]);
		}
		I = i;
	}
	
	//intersection with the closest sphere/mesh
	bool intersect(const ray& r,  Vector &Pscene, Vector &Nscene,double &tscene, int &index){
		double min_t = -1.0;
		int min_index = -1;
		Vector min_N;
		for (int i=0; i<object_list.size(); i++){
			Vector P;
			Vector N;
			double t;
			if ( object_list[i]->intersect(r,P,N,t) ){
				if ( (t < min_t && min_t != -1.0 ) || min_t == -1.0 ){
					//new contestant, we update all
					min_t = t;
					min_index = i;
					min_N = N;
				}
			}
		}
		//if not found, false
		if (min_index == -1) {return false;}

		//else update return values
		tscene =  min_t;
		index = min_index;
		Pscene = r.O + tscene * r.u;

		// On récupère la normale 
		Nscene = min_N;
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
			//Fin de la récursion
			return Vector(0,0,0);
		}

		// Si intersection : disjonction selon la nature de la surface
		if (this->intersect(r,P,N,t,sphere_index)){
			
			// Si intersection miroir : Récursion
			if (object_list[sphere_index]->is_mirror){
				// rayon réfléchi
				// On devrait faire gaffe à la normale : Si on vient de l'intérieur d'un cercle on devrait utiliser -N
				ray reflected(P + 0.0001 * N, r.u - 2 * dot(r.u,N) * N); //Centre décalé, et on inverse la normale
				return this->get_color(reflected,ray_depth -1,false);
			}

			// Si intersection transparente
			if (object_list[sphere_index]->is_transparent){
				
				double n1;
				double n2;
				double orientation = dot(r.u,N); //rayon incident vs normale sortante
				Vector N_reflexion;
				if (orientation < 0){ // au pif pour le < ou <=
					//on entre dans la boule
					n1 = 1;
					n2 = object_list[sphere_index]->n;
					N_reflexion = N;
				}
				else{
					//on sort de la boule, pour etre dans le schéma de snell-descartes on utilise -N
					n1 = object_list[sphere_index]->n;
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

			////// Cas de l'intersection avec une matière diffuse

			// Interaction avec la lumière : On **sait** que l'objet 0 est une sphère
			sphere* light_sphere = dynamic_cast<sphere*>(object_list[0]);

			if (sphere_index == 0){
				// Si interaction diffuse : On renvoie 0 pour éviter de compter 2 fois la source (I_direct et I_indirect)
				if (was_diffuse_interaction){return Vector(0.,0.,0.);}
				// Si intersection avec la lumière : On renvoit la lumière ! (Pondérée par l'élément de surface)
				else {return I*light_sphere->albedo /(4 * sqr(M_PI * light_sphere->R ) ) ;}
				
			}
			else{
				//// Calcul de l'éclairage direct : On intègre l'éclairage de la sphère source
				Vector I_direct(0.,0.,0.);

				Vector vect_light = (light_sphere->C - P).normalize(); // vecteur de P vers L
				// Point sur la lumière: on utilise random_cos qu'on scale.
				Vector N_prime = random_cos(-vect_light); 
				Vector P_prime = light_sphere->C + light_sphere->R * N_prime; 
				// Rayon incident
				Vector P_epsilon = P + 0.0001 * N;
				Vector wi_direct = P_prime - P_epsilon;
				wi_direct = wi_direct.normalize();

				// On renvoie un rayon de P(epsilon) vers P'. pour le calcul du terme de visibilité
				bool visibility = false;
				ray r2(P_epsilon , wi_direct); // Rayon de P vers P'
				double d2 = (P_prime - P_epsilon).norm2(); // distance PP'² 
				Vector P2;
				Vector N2;
				double t2;
				int sphere_index2;
				// Il y aura toujours une intersection(avec la lumière), mais on veut savoir si elle se fait avant la lumière
				if (this->intersect(r2,P2,N2,t2,sphere_index2)){
					//visibility = sqr(t2 + 0.01) > d2; // P voit la lumière ssi l'intersection est plus loin que la lumière
					//plutot tester si l'objet le plus proche est la lumière
					//le test sur t2 ne semble pas etre assez fiable :  différences de 1 ou + avec d2 
					visibility = (sphere_index2 == 0);
					
				}
				else{
					// Au cas où on n'intersecte rien (Ne devrait pas arriver)
					std::cout<<"light not found on diffuse bounce"<<std::endl;
					visibility = false;
				}


				// Contribution ou non de l'éclairage direct selon la visibilité. Pas besoin de faire les calculs si pas
				// de visibilité, on laisse I_direct à 0.
				if (visibility){
					double light_intensity = I / (4 * sqr(M_PI * light_sphere->R) ); // intensité de la lumière, normalisée par l'élément de surface
					Vector albedo = object_list[sphere_index]->albedo/M_PI; // albedo de la sphere considérée
					double form_factor = std::max(0.,dot(N,wi_direct))*std::max(0.,dot(N_prime,-wi_direct))/(d2);
					double pdf = std::max(0.,dot(N_prime,-vect_light)) / ( M_PI * light_sphere->R * light_sphere->R ) ;

					I_direct = light_intensity * albedo * form_factor / pdf;
				}

				//// Calcul de l'éclairage indirect
				Vector I_indirect(0.,0.,0.);
				ray wi_indirect(Vector(P + 0.0001*N),random_cos(N)); 
				// albedo* getcolor(wi) et wi suit cos(theta)/pi
				I_indirect += object_list[sphere_index]->albedo * get_color(wi_indirect,ray_depth-1, true); // Seul cas où on signale qu'on ne compte pas la lumière
				return I_direct + I_indirect;
			}


		}

		// Si pas d'intersection, rien (mais ne devrait pas arriver, c'est une scène fermée)
		// En pratique, arrive 3-5 fois pour 512x512p, 128 rpp
		else {
			std::cout<<"Warning: No intersection for ray "<<ray_depth<<std::endl;
			return Vector(0.,0.,0.);
		}
			
	}

	std::vector<Geometry*> object_list;
	double I;
};



int main() {
	time_t time_start = time(NULL);

	int W = 256;
	int H = 256;
	int ray_count = 16;
	int ray_depth = 0;
	std::vector<unsigned char> image(W*H * 3, 0);

	Vector center(0, 0, 55);

	// La lumière est la première boule
	sphere* light_sphere = new sphere(Vector(-10,20,40),10,Vector(1.,1.,1.),false,false,1.333);
	sphere* miroir = new sphere(Vector(0,0,0),10,Vector(0.3,0.4,0.9),true,false,1.333); // miroir
	sphere* transparente = new sphere(Vector(-20,0,0),10,Vector(0.3,0.4,0.9),false,true,1.333); //transparente 
	sphere* solide = new sphere(Vector(20,0,0),10,Vector(0.3,0.4,0.9),false,false,1.333); // solide
	sphere* background = new sphere(Vector(0,0,-1000),940,Vector(0.85,0.1,0.1),false,false,1.333);
	sphere* topwall = new sphere(Vector(0,1000,0),940,Vector(0.1,0.9,0.2),false,false,1.333);
	sphere* foreground = new sphere(Vector(0,0,1000),940,Vector(0.9,0.3,0.1),false,false,1.333);
	sphere* bottomwall = new sphere(Vector(0,-1000,0),990,Vector(0.1,0.1,0.95),false,false,1.333);
	sphere* rightwall = new sphere(Vector(1000,0,0),940,Vector(0.95,0.95,.07),false,false,1.333);
	sphere* leftwall = new sphere(Vector(-1000,0,0),940,Vector(0.07,0.95,0.95),false,false,1.333);

	TriangleMesh* cat = new TriangleMesh();
	cat->readOBJ("cat.obj");
	cat->affine(Vector(0,-10,0),0.3);
	cat->set_bounding_box();



	std::vector<Geometry*> object_list;
	object_list.push_back(light_sphere);
	//object_list.push_back(miroir);
	//object_list.push_back(transparente);
	//object_list.push_back(solide);
	object_list.push_back(background);
	object_list.push_back(topwall);
	object_list.push_back(foreground);
	object_list.push_back(bottomwall);
	object_list.push_back(rightwall);
	object_list.push_back(leftwall);

	object_list.push_back(cat);

	double I = 3E9;

	scene scene_1(object_list,I);
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
				color += scene_1.get_color(r,ray_depth,false)/ray_count;
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