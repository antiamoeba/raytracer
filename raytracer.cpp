#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iterator> 
#include <math.h>
#include <fstream>
#include <string>
#include <time.h> 
#include <algorithm>
#include "Eigen/Dense"
#include "lodepng.h"

using Eigen::Matrix4f;
using Eigen::Vector4f;
using namespace std;

#define minVal 0.001
#define maxVal 9999999
#define PI 3.14159265
void split(string s, char delim, vector<string> *elems) {
    stringstream stream;
    stream.str(s);
    string str;
    while (std::getline(stream, str, delim)) {
        elems->push_back(str);
    }
}
void split(string s, vector<string> *elems) {
    istringstream buffer(s);
    copy(istream_iterator<string>(buffer), istream_iterator<string>(), back_inserter(*elems));
}
class Vector3d {
    public: 
        float x;
        float y;
        float z;
        Vector3d(float x1, float y1, float z1) : x(x1), y(y1), z(z1) {}
        Vector3d() : x(0), y(0), z(0) {}
        static Vector3d mult(Vector3d vec1, Vector3d vec2) {
            Vector3d multed(vec1.x * vec2.x, vec1.y * vec2.y, vec1.z * vec2.z);
            return multed;
        }
        static Vector3d mult_num(Vector3d vec1, float num) {
            Vector3d multed(vec1.x * num, vec1.y * num, vec1.z * num);
            return multed;
        }
        static Vector3d dot(Vector3d vec1, Vector3d vec2) {
            float val = vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
            Vector3d dotted(val, val, val);
            return dotted;
        }
        static float dot_num(Vector3d vec1, Vector3d vec2) {
            float val = vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
            return val;
        }
        static Vector3d add(Vector3d vec1, Vector3d vec2) {
            Vector3d sum(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
            return sum;
        }
        static Vector3d sub(Vector3d vec1, Vector3d vec2) {
            Vector3d diff(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
            return diff;
        }
        static Vector3d cross(Vector3d vec1, Vector3d vec2) {
            Vector3d crossed(vec1.y * vec2.z - vec1.z * vec2.y, vec1.z * vec2.x - vec1.x * vec2.z, vec1.x * vec2.y - vec1.y * vec2.x);
            return crossed;
        }
        float magnitude() {
            return sqrt(x * x + y * y + z * z);
        }
        float getAxis(int num) {
            if (num == 0) {
                return x;
            } else if(num == 1) {
                return y;
            } else {
                return z;
            }
        }
        void normalize() {
            float len = sqrt(x * x + y * y + z * z);
            x = x / len;
            y = y / len;
            z = z / len;
        }
        Vector3d copy() {
            return Vector3d(x, y, z);
        }
};
class Ray {
    public:
        Vector3d start;
        Vector3d vec;
        Ray(Vector3d s, Vector3d v): start(s), vec(v) {}
        Vector3d getPoint(float t) {
            return Vector3d::add(start, Vector3d::mult_num(vec, t));
        }
};
class Material {
    public:
        Vector3d ks;
        Vector3d kd;
        Vector3d ka;
        Vector3d km;
        Vector3d krfr;
        float n;
        float sp;
        Material(Vector3d ks1, Vector3d kd1, Vector3d ka1, Vector3d km1, float p1, Vector3d krfr1, float n1) : ks(ks1), kd(kd1), ka(ka1), km(km1), sp(p1), krfr(krfr1), n(n1) {}
        Material() {}
};
class Surface {
    public: 
        Matrix4f minv;
        float xmin;
        float xmax;
        float ymin;
        float ymax;
        float zmin;
        float zmax;
        Vector3d getLoc() {
            Vector3d loc((xmin + xmax)/2, (ymin + ymax)/2, (zmin + zmax)/2);
            return loc;
        }
        Surface* intersects(Ray r, float* thit, Vector3d* normal, Material* m) {

            float txmin;
            float txmax;
            float ax = 1 / r.vec.x;
            if (ax >= 0) {
                txmin = ax * (xmin - r.start.x);
                txmax = ax * (xmax - r.start.x);
            } else {
                txmax = ax * (xmin - r.start.x);
                txmin = ax * (xmax - r.start.x);
            }

            float tymin;
            float tymax;
            float ay = 1 / r.vec.y;
            if (ay >= 0) {
                tymin = ay * (ymin - r.start.y);
                tymax = ay * (ymax - r.start.y);
            } else {
                tymax = ay * (ymin - r.start.y);
                tymin = ay * (ymax - r.start.y);
            }

            float tzmin;
            float tzmax;
            float az = 1 / r.vec.z;
            if (az >= 0) {
                tzmin = az * (zmin - r.start.z);
                tzmax = az * (zmax - r.start.z);
            } else {
                tzmax = az * (zmin - r.start.z);
                tzmin = az * (zmax - r.start.z);
            }

            if (tzmin > txmax || tzmin > tymax || txmin > tymax || txmin > tzmax || tymin > tzmax || tymin > txmax) {
                return NULL;
            }

            Vector4f a;
            Vector4f b;
            a(0) = r.start.x;
            a(1) = r.start.y;
            a(2) = r.start.z;
            a(3) = 1;
            b(0) = r.vec.x;
            b(1) = r.vec.y;
            b(2) = r.vec.z;
            b(3) = 0;
            Vector4f na = minv * a;
            Vector4f nb = minv * b;

            Vector3d nstart(na(0), na(1), na(2));
            Vector3d nvec(nb(0), nb(1), nb(2));
            Ray nr(nstart, nvec);
            //bounding box
            Surface* sub = intersectsSurface(nr, thit, normal, m);
            if (sub != NULL) {
                Vector4f n;
                n(0) = normal->x;
                n(1) = normal->y;
                n(2) = normal->z;
                n(3) = 1;
                Vector4f on = minv.transpose() * n;
                normal->x = on(0);
                normal->y = on(1);
                normal->z = on(2);
                normal->normalize();
                return sub;
            }
            return NULL;
        }
        virtual Surface* intersectsSurface(Ray r, float* thit, Vector3d* normal, Material* m) = 0;
};
class ParentSurface: public Surface {
    public:
        vector<Surface* > surfaces;
        ParentSurface(Matrix4f matrix) {
            xmin = maxVal;
            xmax = -maxVal;

            ymin = maxVal;
            ymax = -maxVal;

            zmin = maxVal;
            zmax = -maxVal;
            minv = matrix.inverse();
        }
        void addSurface(Surface* surface) {
            surfaces.push_back(surface);
            Vector4f aei;
            aei << surface->xmin, surface->ymin, surface->zmin, 1;
            aei = minv.inverse() * aei;
            Vector4f bei;
            bei << surface->xmin, surface->ymax, surface->zmin, 1;
            bei = minv.inverse() * bei;
            Vector4f cei;
            cei << surface->xmax, surface->ymin, surface->zmin, 1;
            cei = minv.inverse() * cei;
            Vector4f dei;
            dei << surface->xmax, surface->ymax, surface->zmin, 1;
            dei = minv.inverse() * dei;
            Vector4f eei;
            eei << surface->xmin, surface->ymin, surface->zmax, 1;
            eei = minv.inverse() * eei;
            Vector4f fei;
            fei << surface->xmin, surface->ymax, surface->zmax, 1;
            fei = minv.inverse() * fei;
            Vector4f gei;
            gei << surface->xmax, surface->ymin, surface->zmax, 1;
            gei = minv.inverse() * gei;
            Vector4f hei;
            hei << surface->xmax, surface->ymax, surface->zmax, 1;
            hei = minv.inverse() * hei;

            xmin = min(xmin, min(aei(0), min(bei(0), min(cei(0), min(dei(0), min(eei(0), min(fei(0), min(gei(0), hei(0)))))))));
            xmax = max(xmax, max(aei(0), max(bei(0), max(cei(0), max(dei(0), max(eei(0), max(fei(0), max(gei(0), hei(0)))))))));
            ymin = min(ymin, min(aei(1), min(bei(1), min(cei(1), min(dei(1), min(eei(1), min(fei(1), min(gei(1), hei(1)))))))));
            ymax = max(ymax, max(aei(1), max(bei(1), max(cei(1), max(dei(1), max(eei(1), max(fei(1), max(gei(1), hei(1)))))))));
            zmin = min(zmin, min(aei(2), min(bei(2), min(cei(2), min(dei(2), min(eei(2), min(fei(2), min(gei(2), hei(2)))))))));
            zmax = max(zmax, max(aei(2), max(bei(2), max(cei(2), max(dei(2), max(eei(2), max(fei(2), max(gei(2), hei(2)))))))));
        }
        Surface* intersectsSurface(Ray r, float* thit, Vector3d* normal, Material* m) {
            bool flag = false;
            *thit = maxVal;
            Surface* surf = NULL;
            for (int i = 0; i < surfaces.size(); i++) {
                Surface* surface = surfaces[i];
                Material mat;
                Vector3d n;
                float tm;
                Surface* sub = surface->intersects(r, &tm, &n, &mat);
                if (sub != NULL && tm < *thit && tm > minVal) {
                    *thit = tm;
                    surf = sub;
                    normal->x = n.x;
                    normal->y = n.y;
                    normal->z = n.z;
                    m->kd = mat.kd;
                    m->ks = mat.ks;
                    m->ka = mat.ka;
                    m->km = mat.km;
                    m->sp = mat.sp;
                    m->krfr = mat.krfr;
                    m->n = mat.n;
                    flag = true;
                }
            }
            return surf;
        }
};
class Triangle: public Surface {
    public:
        Vector3d pa;
        Vector3d paN;
        Vector3d pb;
        Vector3d pbN;
        Vector3d pc;
        Vector3d pcN;
        Material material;
        Triangle(Vector3d a1, Vector3d b1, Vector3d c1, Material m, Matrix4f matrix) : material(m) {
            Vector3d one = Vector3d::sub(c1, b1);
            Vector3d two = Vector3d::sub(a1, b1);
            Vector3d normal = Vector3d::cross(one, two);

            Vector3d loc((a1.x + b1.x + c1.x)/3, (a1.y + b1.y + c1.y)/3, (a1.z + b1.z + c1.z)/3 );
            pa = Vector3d(a1.x - loc.x, a1.y - loc.y, a1.z - loc.z);
            pb = Vector3d(b1.x - loc.x, b1.y - loc.y, b1.z - loc.z);
            pc = Vector3d(c1.x - loc.x, c1.y - loc.y, c1.z - loc.z);

            Matrix4f backMatrix;
            backMatrix << 1, 0, 0, loc.x,
                          0, 1, 0, loc.y,
                          0, 0, 1, loc.z,
                          0, 0, 0, 1;

            minv = (backMatrix * matrix).inverse();

            normal.normalize();
            paN = normal;
            pbN = normal;
            pcN = normal;

            Vector4f aei;
            aei << pa.x, pa.y, pa.z, 1;
            Vector4f bei;
            bei << pb.x, pb.y, pb.z, 1;
            Vector4f cei;
            cei << pc.x, pc.y, pc.z, 1;

            aei = minv.inverse() * aei;
            bei = minv.inverse() * bei;
            cei = minv.inverse() * cei;

            xmin = min(aei(0), min(bei(0), cei(0)));
            ymin = min(aei(1), min(bei(1), cei(1)));
            zmin = min(aei(2), min(bei(2), cei(2)));

            xmax = max(aei(0), max(bei(0), cei(0)));
            ymax = max(aei(1), max(bei(1), cei(1)));
            zmax = max(aei(2), max(bei(2), cei(2)));
        }
        Surface* intersectsSurface(Ray r, float* thit, Vector3d* normal, Material *mat) {
            float a = pa.x - pb.x;
            float b = pa.y - pb.y;
            float c = pa.z - pb.z;
            float d = pa.x - pc.x;
            float e = pa.y - pc.y;
            float f = pa.z - pc.z;
            float g = r.vec.x;
            float h = r.vec.y;
            float i = r.vec.z;
            float j = pa.x - r.start.x;
            float k = pa.y - r.start.y;
            float l = pa.z - r.start.z;

            float m = a * (e * i - h * f) + b * (g * f - d * i) + c * (d * h - e * g);

            float t = -(f*(a*k - j*b) + e*(j*c - a*l) + d*(b*l - k*c)) / m;
            if (t < minVal || t > maxVal) {
                return NULL;
            }
            float gamma = (i*(a*k - j*b) + h*(j*c - a*l) + g*(b*l - k*c))/m;
            if (gamma < minVal || gamma > maxVal) {
                return NULL;
            }
            float beta = (j*(e*i - h*f) + k*(g*f - d*i) + l*(d*h - e*g))/m;
            if (beta < minVal || beta > 1 - gamma) {
                return NULL;
            }
            *thit = t;
            float lambda = 1 - beta - gamma;
            Vector3d ac = Vector3d::mult_num(paN, lambda);
            Vector3d bc = Vector3d::mult_num(pbN, beta);
            Vector3d cc = Vector3d::mult_num(pcN, gamma);
            Vector3d inN = Vector3d::add(ac, Vector3d::add(bc, cc));
            inN.normalize();
            normal->x = inN.x;
            normal->y = inN.y;
            normal->z = inN.z;

            mat->ks = material.ks;
            mat->kd = material.kd;
            mat->ka = material.ka;
            mat->km = material.km;
            mat->sp = material.sp;
            mat->krfr = material.krfr;
            mat->n = material.n;
            return this;
        }
};
class Polygon: public Surface {
    public:
        vector<Ray* > edges;
        Vector3d point;
        Vector3d normal;
        Material material;
        float oldxmin;
        float oldymin;
        float oldzmin;
        float oldxmax;
        float oldymax;
        float oldzmax;
        Polygon(vector<Vector3d *> points, Material m, Matrix4f matrix) : material(m) {

            float avgx = 0;
            float avgy = 0;
            float avgz = 0;
            for (int i = 0; i < points.size(); i++) {
                Vector3d* point = points[i];
                avgx += point->x;
                avgy += point->y;
                avgz += point->z;
            }
            avgx /= points.size();
            avgy /= points.size();
            avgz /= points.size();
            Vector3d avg(avgx, avgy, avgz);

            Matrix4f backMatrix;
            backMatrix << 1, 0, 0, avg.x,
                          0, 1, 0, avg.y,
                          0, 0, 1, avg.z,
                          0, 0, 0, 1;
            minv = (backMatrix * matrix).inverse();

            point = Vector3d::sub(*points[0], avg);
            Vector4f aei;
            aei << point.x, point.y, point.z, 1;
            aei = minv.inverse() * aei;
            xmin = aei(0);
            xmax = aei(0);
            ymin = aei(1);
            ymax = aei(1);
            zmin = aei(2);
            zmax = aei(2);

            oldxmin = point.x;
            oldymin = point.y;
            oldzmin = point.z;
            oldxmax = point.x;
            oldymax = point.y;
            oldzmax = point.z;
            Vector3d p = point;
            for (int i = 1; i < points.size(); i++) {
                Vector3d npoint = Vector3d::sub(*points[i], avg);
                Vector3d dir = Vector3d::sub(npoint, p);
                Ray* r = new Ray(p, dir);
                p = npoint;
                edges.push_back(r);


                Vector4f bei;
                bei << p.x, p.y, p.z, 1;
                bei = minv.inverse() * bei;
                xmin = min(xmin, bei(0));
                xmax = max(xmax, bei(0));
                ymin = min(ymin, bei(1));
                ymax = max(ymax, bei(1));
                zmin = min(zmin, bei(2));
                zmax = max(zmax, bei(2));

                oldxmin = min(oldxmin, p.x);
                oldxmax = max(oldxmax, p.x);
                oldymin = min(oldymin, p.y);
                oldymax = max(oldymax, p.y);
                oldzmin = min(oldzmin, p.z);
                oldzmax = max(oldzmax, p.z);
            }
            Vector3d dir = Vector3d::sub(point, p);
            Ray* r = new Ray(p, dir);
            edges.push_back(r);

            Vector3d n = Vector3d::cross(Vector3d::sub(edges[2]->start, edges[1]->start), Vector3d::sub(edges[0]->start, edges[1]->start));
            normal = n;
            normal.normalize();

        }
        Surface* intersectsSurface(Ray r, float* thit, Vector3d* n, Material *m) {
            float top = Vector3d::dot_num(Vector3d::sub(point, r.start), normal);
            float bot = Vector3d::dot_num(r.vec, normal);
            if (bot == 0) {
                return NULL;
            }

            float t = top/bot;
            if (t < minVal) {
                return NULL;
            }

            Vector3d p = r.getPoint(t);
            int count = 0;
            if (oldxmin == oldxmax) {
                Vector3d dir(0, 1, 0);
                Ray nr(p, dir);

                for(int i = 0; i < edges.size(); i++) {
                    Ray edge = *edges[i];
                    float u = nr.start.z - edge.start.z;
                    u = u / edge.vec.z;
                    float v = u * edge.vec.y - nr.start.y + edge.start.y;
                    if (u >= 0 && u <= 1 && v > minVal) {
                        count++;
                    }
                }
            } else if (oldymin == oldymax) {
                Vector3d dir(1, 0, 0);
                Ray nr(p, dir);

                for(int i = 0; i < edges.size(); i++) {
                    Ray edge = *edges[i];
                    float u = nr.start.z - edge.start.z;
                    u = u / edge.vec.z;
                    float v = u * edge.vec.x - nr.start.x + edge.start.x;
                    if (u >= 0 && u <= 1 && v > minVal) {
                        count++;
                    }
                }
            } else {
                Vector3d dir(1, 0, 0);
                Ray nr(p, dir);

                for(int i = 0; i < edges.size(); i++) {
                    Ray edge = *edges[i];
                    float u = nr.start.y - edge.start.y;
                    u = u / edge.vec.y;
                    float v = u * edge.vec.x - nr.start.x + edge.start.x;
                    if (u >= 0 && u <= 1 && v > minVal) {
                        count++;
                    }
                }
            }
            if (count%2 == 0) {
                return NULL;
            }
            n->x = normal.x;
            n->y = normal.y;
            n->z = normal.z;
            *thit = t;

            m->ks = material.ks;
            m->kd = material.kd;
            m->ka = material.ka;
            m->km = material.km;
            m->sp = material.sp;
            m->krfr = material.krfr;
            m->n = material.n;
            return this;
        }
};
class Sphere: public Surface {
    public:
        Vector3d center;
        float radius;
        Material material;
        Sphere(Vector3d c, float r, Material m, Matrix4f matrix): radius(r), material(m) {
            Matrix4f backMatrix;
            backMatrix << 1, 0, 0, c.x,
                          0, 1, 0, c.y,
                          0, 0, 1, c.z,
                          0, 0, 0, 1;
            minv = (backMatrix * matrix).inverse();
            Vector4f aei;
            aei << -r, -r, -r, 1;
            aei = minv.inverse() * aei;
            Vector4f bei;
            bei << -r, r, -r, 1;
            bei = minv.inverse() * bei;
            Vector4f cei;
            cei << r, -r, -r, 1;
            cei = minv.inverse() * cei;
            Vector4f dei;
            dei << r, r, -r, 1;
            dei = minv.inverse() * dei;
            Vector4f eei;
            eei << -r, -r, r, 1;
            eei = minv.inverse() * eei;
            Vector4f fei;
            fei << -r, r, r, 1;
            fei = minv.inverse() * fei;
            Vector4f gei;
            gei << r, -r, r, 1;
            gei = minv.inverse() * gei;
            Vector4f hei;
            hei << r, r, r, 1;
            hei = minv.inverse() * hei;

            xmin = min(aei(0), min(bei(0), min(cei(0), min(dei(0), min(eei(0), min(fei(0), min(gei(0), hei(0))))))));
            xmax = max(aei(0), max(bei(0), max(cei(0), max(dei(0), max(eei(0), max(fei(0), max(gei(0), hei(0))))))));
            ymin = min(aei(1), min(bei(1), min(cei(1), min(dei(1), min(eei(1), min(fei(1), min(gei(1), hei(1))))))));
            ymax = max(aei(1), max(bei(1), max(cei(1), max(dei(1), max(eei(1), max(fei(1), max(gei(1), hei(1))))))));
            zmin = min(aei(2), min(bei(2), min(cei(2), min(dei(2), min(eei(2), min(fei(2), min(gei(2), hei(2))))))));
            zmax = max(aei(2), max(bei(2), max(cei(2), max(dei(2), max(eei(2), max(fei(2), max(gei(2), hei(2))))))));
        }
        Surface* intersectsSurface(Ray r, float* thit, Vector3d* normal, Material *m) {
            float a = Vector3d::dot_num(r.vec, r.vec);
            Vector3d emc = Vector3d::sub(r.start, center);
            float b = 2 * Vector3d::dot_num(r.vec, emc);
            float c = Vector3d::dot_num(emc, emc) - radius * radius;

            float determinant = b * b - 4 * a * c;
            if (determinant < 0) {
                return NULL;
            }
            float t1 = (-b + sqrt(determinant)) / (2 * a);
            float t2 = (-b - sqrt(determinant)) / (2 * a);
            if (abs(t1) > abs(t2)) {
                *thit = t2;
                Vector3d point = Vector3d::sub(r.getPoint(t2), center);
                normal->x = point.x;
                normal->y = point.y;
                normal->z = point.z;
                normal->normalize();
            }
            if (abs(t2) > abs(t1)) {
                *thit = t1;
                Vector3d point = Vector3d::sub(r.getPoint(t1), center);
                normal->x = point.x;
                normal->y = point.y;
                normal->z = point.z;
                normal->normalize();
            }

            m->ks = material.ks;
            m->kd = material.kd;
            m->ka = material.ka;
            m->km = material.km;
            m->sp = material.sp;
            m->krfr = material.krfr;
            m->n = material.n;
            return this;
        }
};
class Light {
    public:
        Vector3d intensity;
        int num_samples;
        virtual float getDirection(Vector3d p, Vector3d* dir) = 0;
        Light(Vector3d i): intensity(i) {}
        virtual Vector3d getIntensity(Vector3d p) = 0;
};
class AmbientLight: public Light {
    public:
        AmbientLight(Vector3d i): Light(i) {
            num_samples = 0;
        }
        float getDirection(Vector3d p, Vector3d* dir) {
            return 0;
        }
        Vector3d getIntensity(Vector3d p) {
            return intensity;
        }
};
class PointLight: public Light {
    public:
        Vector3d point;
        int falloffType;
        PointLight(Vector3d p, Vector3d i, int fall): Light(i) {
            point.x = p.x;
            point.y = p.y;
            point.z = p.z;
            falloffType = fall;
            num_samples = 1;
        }
        float getDirection(Vector3d p, Vector3d* dir) {
            Vector3d vec = Vector3d::sub(point, p);
            *dir = vec;
            return 1;
        }
        Vector3d getIntensity(Vector3d p) {
            if (falloffType == 0) {
                return intensity;
            }
            
            Vector3d dist = Vector3d::sub(p, point);
            float magnitude = dist.magnitude();
            if (falloffType == 1) {
                return Vector3d::mult_num(intensity, 1/(magnitude));
            }
            if (falloffType == 2) {
                return Vector3d::mult_num(intensity, 1/(magnitude*magnitude));
            }
            return intensity;
        }

};
class AreaLight: public Light {
    public:
        Vector3d point;
        int falloffType;
        float radius;
        AreaLight(Vector3d p, Vector3d i, int fall, float size, int ns): Light(i) {
            point.x = p.x;
            point.y = p.y;
            point.z = p.z;
            falloffType = fall;
            radius = size;
            num_samples = ns;
        }
        Vector3d getPoint() {
            float dx = ((float) rand() / (RAND_MAX)) * 2 * radius - radius;
            float dy = ((float) rand() / (RAND_MAX)) * 2 * radius - radius;
            float dz = ((float) rand() / (RAND_MAX)) * 2 * radius - radius;
            Vector3d np(point.x + dx, point.y + dy, point.z + dz);
            return np;
        }
        float getDirection(Vector3d p, Vector3d* dir) {
            Vector3d vec = Vector3d::sub(getPoint(), p);
            *dir = vec;
            return 1;
        }
        Vector3d getIntensity(Vector3d p) {
            if (falloffType == 0) {
                return intensity;
            }
            
            Vector3d dist = Vector3d::sub(p, getPoint());
            float magnitude = dist.magnitude();
            if (falloffType == 1) {
                return Vector3d::mult_num(intensity, 1/(magnitude));
            }
            if (falloffType == 2) {
                return Vector3d::mult_num(intensity, 1/(magnitude*magnitude));
            }
            return intensity;
        }
};
class DirectionalLight: public Light {
    public:
        Vector3d direction;
        DirectionalLight(Vector3d d, Vector3d i): Light(i) {
            direction.x = d.x;
            direction.y = d.y;
            direction.z = d.z;
            direction.normalize();
            num_samples = 1;
        }
        float getDirection(Vector3d p, Vector3d* dir) {
            *dir = direction;
            return maxVal;
        }
        Vector3d getIntensity(Vector3d p) {
            return intensity;
        }
};
class Scene {
    public:
        Vector3d eye;
        Vector3d tl;
        Vector3d tr;
        Vector3d bl;
        Vector3d br;
        int h;
        int w;
        vector<Surface* > surfaces;
        vector<Light* > lights;
        vector<Light* > ambientLights;
        Scene(Vector3d e, Vector3d tl, Vector3d tr, Vector3d bl, Vector3d br, int height, int width) {
            setup(e, tl, tr, bl, br, height, width);
        }
        Scene() {}
        void setup(Vector3d e, Vector3d tl, Vector3d tr, Vector3d bl, Vector3d br, int height, int width) {
            eye.x = e.x;
            eye.y = e.y;
            eye.z = e.z;

            this->tl = tl;
            this->tr = tr;
            this->bl = bl;
            this->br = br;
            h = height;
            w = width;
        }
        void addSphere(Vector3d center, float radius, Material m, Matrix4f matrix = Matrix4f::Identity(), ParentSurface* surface = NULL) {
            Sphere* sphere = new Sphere(center, radius, m, matrix); 
            if (surface == NULL) {
                surfaces.push_back(sphere);
            } else {
                surface->addSurface(sphere);
            }
        }
        void addTriangle(Vector3d a, Vector3d b, Vector3d c, Material m, Matrix4f matrix = Matrix4f::Identity(), ParentSurface* surface = NULL) {
            Triangle* triangle = new Triangle(a, b, c, m, matrix);
            if (surface == NULL) {
                surfaces.push_back(triangle);
            } else {
                surface->addSurface(triangle);
            }
        }
        void addPolygon(vector<Vector3d* > points, Material m, Matrix4f matrix = Matrix4f::Identity(), ParentSurface* surface = NULL) {
            Polygon* polygon = new Polygon(points, m, matrix);
            if (surface == NULL) {
                surfaces.push_back(polygon);
            } else {
                surface->addSurface(polygon);
            }
        }
        void addPointLight(Vector3d point, Vector3d intensity, int falloffType = 0) {
            PointLight* light = new PointLight(point, intensity, falloffType);
            lights.push_back(light);
        }
        void addAreaLight(Vector3d point, Vector3d intensity, int falloffType = 0, float radius = 1, int num_samples = 10) {
            AreaLight* light = new AreaLight(point, intensity, falloffType, radius, num_samples);
            lights.push_back(light);
        }
        void addDirectionalLight(Vector3d direction, Vector3d intensity) {
            DirectionalLight* light = new DirectionalLight(direction, intensity);
            lights.push_back(light);
        }
        void addAmbientLight(Vector3d intensity) {
            AmbientLight* light = new AmbientLight(intensity);
            ambientLights.push_back(light);
        }
        void addObj(string filename, Material m, Matrix4f matrix = Matrix4f::Identity()) {
            ifstream file(filename);
            string str;
            vector<Vector3d> vertices;
            vector<Vector3d> normals;
            vector<Surface* > sList;
            while (getline(file, str)) {
                vector<string> elems;
                printf("%s\n", str.c_str());
                split(str, &elems);
                if (elems.size() >= 4) {
                    string command = elems[0];
                    if (command == "v") {
                        Vector3d vertex(stof(elems[1]), stof(elems[2]), stof(elems[3]));
                        vertices.push_back(vertex);
                    } else if (command == "vn") {
                        Vector3d normal(stof(elems[1]), stof(elems[2]), stof(elems[3]));
                        normals.push_back(normal);
                    } else if (command == "f") {
                        if (elems.size() - 1 == 3) {
                            vector<string> params;
                            split(elems[1], '/', &params);
                            int index = stoi(params[0]) - 1;
                            params.clear();
                            split(elems[2], '/', &params);
                            int index1 = stoi(params[0]) - 1;
                            params.clear();
                            split(elems[3], '/', &params);
                            int index2 = stoi(params[0]) - 1;
                            Triangle* triangle = new Triangle(vertices[index], vertices[index1], vertices[index2], m, Matrix4f::Identity());
                            sList.push_back(triangle);
                        } else {
                            vector<Vector3d* > pverts;
                            for (int i = 1; i < elems.size(); i++) {
                                int index = stoi(elems[i]) - 1;
                                pverts.push_back(&vertices[index]);
                            }
                            Polygon* polygon = new Polygon(pverts, m, Matrix4f::Identity());
                            sList.push_back(polygon);
                        }
                    }
                }
            }
            ParentSurface* ps = bvhHelper(0, sList, 0, -1, matrix);
            surfaces.push_back(ps);
        }
        void generateBVHTree(int maxDepth = -1) {
            if (maxDepth == 0) {
                return;
            }
            ParentSurface* ps = bvhHelper(0, surfaces, 0, maxDepth);
            surfaces.clear();
            surfaces.push_back(ps);
        }
        class SortNode {
            public:
                float val;
                Surface* surface;
                SortNode(float v, Surface* s): val(v), surface(s) {}
        };
        ParentSurface* bvhHelper(int axis, vector<Surface* > sList, int depth, int maxDepth, Matrix4f matrix = Matrix4f::Identity()) {
            ParentSurface* ps = new ParentSurface(matrix);
            if (sList.size() <= 2 || (maxDepth > -1 && depth >= maxDepth)) {
                for (int i = 0; i < sList.size(); i++) {
                    ps->addSurface(sList[i]);
                }
                return ps;
            }
            vector<SortNode* > toSort;
            for (int i = 0; i < sList.size(); i++) {
                float val;
                if (axis == 0) {
                    val = (sList[i]->xmin + sList[i]->xmax) / 2;
                } else if (axis == 1) {
                    val = (sList[i]->ymin + sList[i]->ymax) / 2;
                } else {
                    val = (sList[i]->zmin + sList[i]->zmax) / 2;
                }
                SortNode* sn = new SortNode(val, sList[i]);
                toSort.push_back(sn);
            }
            sort(toSort.begin(), toSort.end(), compareSurfaces);

            vector<Surface* > left;
            for (int i = 0; i < toSort.size()/2; i++) {
                left.push_back(toSort[i]->surface);
            }
            vector<Surface* > right;
            for (int i = toSort.size()/2; i < toSort.size(); i++) {
                right.push_back(toSort[i]->surface);
            }
            printf("left!\n");
            ps->addSurface(bvhHelper((axis + 1)%3, left, depth+1, maxDepth));
            printf("right!\n");
            ps->addSurface(bvhHelper((axis + 1)%3, right, depth+1, maxDepth));
            return ps;
        }
    private:
        static bool compareSurfaces(SortNode* i, SortNode* j) {
            return i->val < j->val;
        }
};
class Sampler {
    public:
        virtual void sample(Vector3d bottomLeft, Vector3d topRight, vector<Vector3d> *vec) = 0;
};
class CenterSampler: public Sampler {
    public:
        void sample(Vector3d bottomLeft, Vector3d topRight, vector<Vector3d> *vec) {
            float xmid = (topRight.x - bottomLeft.x) / 2;
            float ymid = (topRight.y - bottomLeft.y) / 2;
            Vector3d point(bottomLeft.x + xmid, bottomLeft.y + ymid, 0);
            vec->push_back(point);
        }
};
class JitterSampler: public Sampler {
    public:
        int rows;
        int cols;
        JitterSampler(int r, int c): rows(r), cols(c) {}
        void sample(Vector3d bottomLeft, Vector3d topRight, vector<Vector3d> *vec) {
            float xrange = (topRight.x - bottomLeft.x)/cols;
            float yrange = (topRight.y - bottomLeft.y)/rows;
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    float dx = xrange * i + ((float) rand() / (RAND_MAX)) * xrange;
                    float dy = yrange * j + ((float) rand() / (RAND_MAX)) * yrange;
                    Vector3d point(bottomLeft.x + dx, bottomLeft.y + dy, 0);
                    vec->push_back(point);
                }
            }
        }
};
class MonteCarloSampler: public Sampler {
    public:
        int num;
        MonteCarloSampler(int count): num(count) {}
        void sample(Vector3d bottomLeft, Vector3d topRight, vector<Vector3d> *vec) {
            float xrange = topRight.x - bottomLeft.x;
            float yrange = topRight.y - bottomLeft.y;
            for (int i = 0; i < num; i++) {
                float dx = ((float) rand() / (RAND_MAX)) * xrange;
                float dy = ((float) rand() / (RAND_MAX)) * yrange;
                Vector3d point(bottomLeft.x + dx, bottomLeft.y + dy, 0);
                vec->push_back(point);
            }
        }
};
class RayTracer {
    public:
        Sampler* sampler;
        RayTracer(Sampler* s) : sampler(s) {}
        void rayTrace(Scene* scene, string filename, int maxdepth, Vector3d sky, bool flag) {
            vector<vector<Vector3d> > image;
            for (int i = 0; i < scene->w; i++) {
                vector<Vector3d> col(scene->h);
                image.push_back(col);
            }
            int counter = 0;
            for (int i = 0; i < scene->w; i++) {
                for (int j = 0; j < scene->h ; j++) {
                    //compute ray
                    Vector3d bls(i, j, 0);
                    Vector3d trs(i+1, j+1, 0);
                    vector<Vector3d> vec;
                    sampler->sample(bls, trs, &vec);
                    Vector3d color;
                    for (int k = 0; k < vec.size(); k++) {
                        float x = vec[k].x;
                        float y = vec[k].y;
                        float u = x*1.0f/(scene->w);
                        float v = (scene->h - y)*1.0f/(scene->h);
                        Vector3d uvll = Vector3d::mult_num(Vector3d::mult_num(scene->bl,1-v), 1-u);
                        Vector3d uvul = Vector3d::mult_num(Vector3d::mult_num(scene->tl, v), 1-u);
                        Vector3d uvlr = Vector3d::mult_num(Vector3d::mult_num(scene->br,1-v), u);
                        Vector3d uvur = Vector3d::mult_num(Vector3d::mult_num(scene->tr, v), u);

                        Vector3d point = Vector3d::add(Vector3d::add(uvll, uvul), Vector3d::add(uvlr, uvur));
                        Vector3d vec = Vector3d::sub(point, scene->eye);
                        Ray r(scene->eye, vec);

                        Vector3d nc = traceHelper(r, scene, 0, maxdepth, sky, flag);
                        color.x += nc.x;
                        color.y += nc.y;
                        color.z += nc.z;
                        counter++;
                    }
                    image[i][j] = Vector3d::mult_num(color, 1.0/vec.size());
                }
            }
            printf("Saving file!\n");
            saveFile(filename, scene->w, scene->h, image);
        }
        void saveFile(string filename, int w, int h, vector<vector<Vector3d> > image) {
            vector<unsigned char> imagedata;
            imagedata.resize(w * h * 4);
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    imagedata[4 * w * y + 4 * x + 0] = image[x][y].x * 255;
                    imagedata[4 * w * y + 4 * x + 1] = image[x][y].y * 255;
                    imagedata[4 * w * y + 4 * x + 2] = image[x][y].z * 255;
                    imagedata[4 * w * y + 4 * x + 3] = 255;
                }
            }
            unsigned error = lodepng::encode(filename.c_str(), imagedata, w, h);
            if (error) {
                printf("There was an error saving to %s.\n", filename.c_str());
            } else {
                printf("Wrote image to %s!\n", filename.c_str());
            }
        }
    private:
        Surface* intersector(Ray r, Scene* scene, float *t, Vector3d *normal, Material *m) {
            *t = maxVal;
            bool flag = false;
            Surface* surf = NULL;
            for (int i = 0; i < scene->surfaces.size(); i++) {
                Surface* s = scene->surfaces[i];
                float lt = 0;
                Vector3d lnormal;
                Material material;
                Surface* sfc = s->intersects(r, &lt, &lnormal, &material);
                if (sfc != NULL && lt < *t && lt > minVal*10) {
                    *t = lt;
                    surf = sfc;
                    m->kd = material.kd;
                    m->ks = material.ks;
                    m->ka = material.ka;
                    m->km = material.km;
                    m->sp = material.sp;
                    m->krfr = material.krfr;
                    m->n = material.n;
                    normal->x = lnormal.x;
                    normal->y = lnormal.y;
                    normal->z = lnormal.z;
                    flag = true;
                }
            }
            return surf;
        }
        bool intersects(Ray r, Scene* scene, float maxT, Surface* curr) {
            for (int i = 0; i < scene->surfaces.size(); i++) {
                Surface* s = scene->surfaces[i];
                if (s != curr) {
                    float lt = 0;
                    Vector3d lnormal;
                    Material material;
                    Surface* surf = s->intersects(r, &lt, &lnormal, &material);
                    if (curr != surf && surf != NULL && lt < maxT && lt > minVal) {
                        return true;
                    }
                }
            }
            return false;
        }
        Vector3d traceHelper(Ray r, Scene* scene, int depth, int maxdepth, Vector3d sky, bool ambientFlag = false) {
            float t = 0;
            Vector3d normal;
            Material material;
            Surface* surface = intersector(r, scene, &t, &normal, &material);
            if (surface == NULL) {
                return sky;
            }
            Vector3d total;
            Vector3d point = r.getPoint(t);
            normal.normalize();
            if (depth < 1) {
                for (int i = 0; i < scene->ambientLights.size(); i++) {
                    Vector3d ambient = Vector3d::mult(material.ka, scene->ambientLights[i]->getIntensity(point));
                    total.x += ambient.x;
                    total.y += ambient.y;
                    total.z += ambient.z;
                }
            }
            Vector3d refractV(0, 0, 0);
            if (material.krfr.x != 0 || material.krfr.y != 0 || material.krfr.z != 0) {
                Vector3d d = r.vec.copy();
                d.normalize();
                Vector3d reflr = Vector3d::mult(Vector3d::mult_num(Vector3d::dot(d, normal), 2), normal);
                Vector3d refl = Vector3d::sub(d, reflr);
                Ray reflectionray(point, refl);
                float dn = Vector3d::dot_num(d, normal);
                float c = 0;
                Vector3d krefract;
                Vector3d tVec;
                if (dn < 0) {
                    refract(d, normal, material.n, &tVec);
                    c = -dn;
                    krefract.x = 1;
                    krefract.y = 1;
                    krefract.z = 1;

                } else {
                    krefract.x = exp(-log(material.krfr.x) * t);
                    krefract.y = exp(-log(material.krfr.y) * t);
                    krefract.z = exp(-log(material.krfr.z) * t);
                    if (refract(d, Vector3d::mult_num(normal, -1), 1/material.n, &tVec)) {
                        c = Vector3d::dot_num(tVec, normal);
                    } else {
                        return Vector3d::mult(krefract, traceHelper(reflectionray, scene, depth+1, maxdepth, sky, ambientFlag));
                    }
                }
                float n = material.n;
                float R0 = ((n-1)*(n-1))/((n+1)*(n+1));
                float R = R0 + (1 - R0) * pow(1 - c, 5);
                Vector3d color1 = Vector3d::mult_num(traceHelper(reflectionray, scene, depth+1, maxdepth, sky, ambientFlag), R);
                Ray refractray(point, tVec);
                Vector3d color2 = Vector3d::mult_num(traceHelper(refractray, scene, depth+1, maxdepth, sky, ambientFlag), 1-R);
                Vector3d finalR = Vector3d::mult(krefract, Vector3d::add(color1, color2));
                refractV = finalR;
                return Vector3d::add(total, refractV);
            }
            for (int i = 0; i < scene->lights.size(); i++) {
                Light* light = scene->lights[i];
                Vector3d ambient;
                if (ambientFlag && depth < 1) {
                    ambient = Vector3d::mult(material.ka, light->getIntensity(point));
                }
                Vector3d diffuse(0, 0, 0);
                Vector3d spec(0, 0, 0);
                for (int j = 0; j < light->num_samples; j++) {
                    Vector3d lv;
                    float maxt = light->getDirection(point, &lv);
                    Ray shadowray(point, lv);
                    if (!intersects(shadowray, scene, maxt, surface)) {
                        lv.normalize();
                        Vector3d dot = Vector3d::dot(lv, normal);
                        if (dot.x > 0) {
                            diffuse = Vector3d::add(diffuse, Vector3d::mult(Vector3d::mult(material.kd, light->getIntensity(point)), dot));
                        }
                        Vector3d rr = Vector3d::mult(dot, Vector3d(2, 2, 2));
                        rr = Vector3d::mult(rr, normal);
                        Vector3d r1 = Vector3d::sub(rr, lv);
                        r1.normalize();
                        Vector3d vdv(-r.vec.x, -r.vec.y, -r.vec.z);
                        vdv.normalize();
                        Vector3d dot2 = Vector3d::dot(r1, vdv);
                        if (dot2.x > 0) {
                            float power = material.sp;
                            Vector3d coeff(pow(dot2.x, power), pow(dot2.y, power), pow(dot2.z, power));
                            spec = Vector3d::add(spec, Vector3d::mult(Vector3d::mult(material.ks, light->getIntensity(point)), coeff)); 
                        }
                    }
                }
                diffuse = Vector3d::mult_num(diffuse, 1.0/(light->num_samples));
                spec = Vector3d::mult_num(spec, 1.0/(light->num_samples));
                total.x += ambient.x + diffuse.x + spec.x;
                total.y += ambient.y + diffuse.y + spec.y;
                total.z += ambient.z + diffuse.z + spec.z;
            }
            Vector3d reflect(0, 0, 0);
            if (depth < maxdepth && (material.km.x != 0 || material.km.y != 0 || material.km.z != 0)) {
                Vector3d d = r.vec.copy();
                d.normalize();
                Vector3d reflr = Vector3d::mult(Vector3d::mult_num(Vector3d::dot(d, normal), 2), normal);
                Vector3d refl = Vector3d::sub(d, reflr);
                Ray reflectionray(point, refl);
                reflect = Vector3d::mult(material.km, traceHelper(reflectionray, scene, depth+1, maxdepth, sky, ambientFlag));
            }
            total.x += reflect.x;
            total.y += reflect.y;
            total.z += reflect.z;
            if (total.x > 1) {
                total.x = 1;
            }
            if (total.y > 1) {
                total.y = 1;
            }
            if (total.z > 1) {
                total.z = 1;
            }
            return total;
        }
        bool refract(Vector3d vec, Vector3d normal, float n, Vector3d* t) {
            float dn = Vector3d::dot_num(vec, normal);
            Vector3d left = Vector3d::mult_num(Vector3d::sub(vec, Vector3d::mult_num(normal, dn)), 1/n);
            float rco = 1 - (1 - dn*dn)/(n*n);
            if (rco < 0) {
                return false;
            }
            Vector3d right = Vector3d::mult_num(normal, sqrt(rco));
            Vector3d total = Vector3d::sub(left, right);
            t->x = total.x;
            t->y = total.y;
            t->z = total.z;
            return true;
        }
};

int main(int argc, char* argv[]) {
    if (argc == 1) {
        printf("Not enough inputs\n");
        return EXIT_FAILURE;
    }
    string outputFile = "output.ppm";
    string inputFile = "";
    srand (time(NULL));
    bool flag = false;
    int counter = 1;
    int depth = 1;
    int bvhDepth = 0;
    Sampler* cs = new CenterSampler();
    while (counter < argc) {
        string arg = string(argv[counter]);
        if (arg == "-i") {
            printf("argument: %s\n", argv[counter]);
            inputFile = string(argv[counter+1]);
            counter = counter + 2;
        } else if (arg == "-o") {
            printf("argument: %s\n", argv[counter]);
            outputFile = string(argv[counter+1]);
            counter = counter + 2;
        } else if (arg == "-ambient") {
            printf("argument: %s\n", argv[counter]);
            flag = true;
            counter = counter + 1;
        } else if (arg == "-d") {
            printf("argument: %s\n", argv[counter]);
            depth = stoi(argv[counter+1]);
            counter = counter + 2;
        } else if (arg == "-j") {
            printf("argument: %s\n", argv[counter]);
            int rows = stoi(argv[counter+1]);
            int cols = stoi(argv[counter+2]);
            cs = new JitterSampler(rows, cols);
            counter = counter + 3;
        } else if (arg == "-mc") {
            printf("argument: %s\n", argv[counter]);
            int count = stoi(argv[counter+1]);
            cs = new MonteCarloSampler(count);
            counter = counter + 2;
        } else if (arg == "-bvh") {
            printf("argument: %s\n", argv[counter]);
            bvhDepth = stoi(argv[counter+1]);
            counter = counter + 2;
        }
        else {
            counter++;
        }
    }
    ifstream file(inputFile);
    string str;
    Scene* scene = new Scene();

    Vector3d kaStart(0.1, 0.1, 0.1);
    Vector3d kdStart(0.5, 0.5, 0.5);
    Vector3d ksStart(1, 1, 1);
    Vector3d kmStart(0, 0, 0);
    Vector3d krfrStart(0, 0, 0);
    Material material(ksStart, kdStart, kaStart, kmStart, 50, krfrStart, 0);

    Matrix4f transform = Matrix4f::Identity();

    RayTracer rt(cs);
    Vector3d sky;
    while (getline(file, str)) {
        printf("%s\n", str.c_str());
        vector<string> elems;
        split(str, ' ', &elems);
        if (elems.size() > 0) {
            if (elems[0] == "cam") {
                Vector3d eye(stof(elems[1]), stof(elems[2]), stof(elems[3]));
                Vector3d bl(stof(elems[4]), stof(elems[5]), stof(elems[6]));
                Vector3d br(stof(elems[7]), stof(elems[8]), stof(elems[9]));
                Vector3d tl(stof(elems[10]), stof(elems[11]), stof(elems[12]));
                Vector3d tr(stof(elems[13]), stof(elems[14]), stof(elems[15]));
                scene->setup(eye, tl, tr, bl, br, 1000, 1000);
            } else if (elems[0] == "sph") {
                Vector3d sphere(stof(elems[1]), stof(elems[2]), stof(elems[3]));
                scene->addSphere(sphere, stof(elems[4]), material, transform);
            } else if (elems[0] == "tri") {
                Vector3d a(stof(elems[1]), stof(elems[2]), stof(elems[3]));
                Vector3d b(stof(elems[4]), stof(elems[5]), stof(elems[6]));
                Vector3d c(stof(elems[7]), stof(elems[8]), stof(elems[9]));
                scene->addTriangle(a, b, c, material, transform);
            } else if (elems[0] == "obj") {
                scene->addObj(elems[1], material, transform);
            } else if (elems[0] == "ltp") {
                Vector3d point(stof(elems[1]), stof(elems[2]), stof(elems[3]));
                Vector3d intensity(stof(elems[4]), stof(elems[5]), stof(elems[6]));
                int falloff = 0;
                if (elems.size() >= 8) {
                    falloff = stoi(elems[7]);
                }
                scene->addPointLight(point, intensity, falloff);
            } else if (elems[0] == "ltd") {
                Vector3d dir(-stof(elems[1]), -stof(elems[2]), -stof(elems[3]));
                Vector3d intensity(stof(elems[4]), stof(elems[5]), stof(elems[6]));
                scene->addDirectionalLight(dir, intensity);
            } else if (elems[0] == "lta") {
                Vector3d intensity(stof(elems[1]), stof(elems[2]), stof(elems[3]));
                scene->addAmbientLight(intensity);
            } else if (elems[0] == "ltr") {
                Vector3d point(stof(elems[1]), stof(elems[2]), stof(elems[3]));
                Vector3d intensity(stof(elems[4]), stof(elems[5]), stof(elems[6]));
                float radius = stof(elems[7]);
                int num_samples = stof(elems[8]);
                int falloff = 0;
                if (elems.size() >= 10) {
                    falloff = stoi(elems[9]);
                }
                scene->addAreaLight(point, intensity, falloff, radius, num_samples);
            } else if (elems[0] == "mat") {
                Vector3d ka(stof(elems[1]), stof(elems[2]), stof(elems[3]));
                Vector3d kd(stof(elems[4]), stof(elems[5]), stof(elems[6]));
                Vector3d ks(stof(elems[7]), stof(elems[8]), stof(elems[9]));
                Vector3d kr(stof(elems[11]), stof(elems[12]), stof(elems[13]));
                Vector3d krfr;
                float n = 0;
                if (elems.size() >= 18) {
                    krfr.x = stof(elems[14]);
                    krfr.y = stof(elems[15]);
                    krfr.z = stof(elems[16]);
                    n = stof(elems[17]);
                }
                material = Material(ks, kd, ka, kr, stof(elems[10]), krfr, n);
            } else if (elems[0] == "xft") {
                Vector3d sV(stof(elems[1]), stof(elems[2]), stof(elems[3]));
                Matrix4f matrix;
                matrix << 1, 0, 0, sV.x,
                          0, 1, 0, sV.y,
                          0, 0, 1, sV.z,
                          0, 0, 0, 1;
                transform = transform * matrix;
            } else if (elems[0] == "xfr") {
                Vector3d sV(stof(elems[1]), stof(elems[2]), stof(elems[3]));
                float theta = sV.magnitude() * PI / 180.0;
                sV.normalize();
                Matrix4f rcross;
                rcross << 0, -sV.z, sV.y, 0,
                          sV.z, 0, -sV.x, 0,
                          -sV.y, sV.x, 0, 0,
                          0, 0, 0, 1;
                
                Matrix4f matrix = rcross * sin(theta) + Matrix4f::Identity() + rcross * rcross * (1 - cos(theta));
                matrix(3,3) = 1;
                transform = transform * matrix;
            } else if (elems[0] == "xfs") {
                Vector3d sV(stof(elems[1]), stof(elems[2]), stof(elems[3]));
                Matrix4f matrix;
                matrix << sV.x, 0, 0, 0,
                          0, sV.y, 0, 0,
                          0, 0, sV.z, 0,
                          0, 0, 0, 1;
                transform = transform * matrix;
            } else if (elems[0] == "xfz") {
                transform = Matrix4f::Identity();
            } else if (elems[0] == "sky") {
                sky = Vector3d(stof(elems[1]), stof(elems[2]), stof(elems[3]));
            }
        }
    }
    scene->generateBVHTree(bvhDepth);
    rt.rayTrace(scene, outputFile, depth, sky, flag);
    /*Vector3d eye(0, 0, 0);
    Vector3d tl(-1, 1, -3);
    Vector3d tr(1, 1, -3);
    Vector3d bl(-1, -1, -3);
    Vector3d br(1, -1, -3);
    Scene* scene = new Scene(eye, tl, tr, bl, br, 1000, 1000);


    Vector3d ka(0.1, 0.1, 0.1);
    Vector3d kd(1, 0, 1);
    Vector3d ks(1, 1, 1);
    Vector3d km(0, 0, 0);
    Material m(ks, kd, ka, km, 50);

    Vector3d sphere(0, 0, -20);
    scene->addSphere(sphere, 3, m);


    Vector3d ka1(0.1, 0.1, 0.1);
    Vector3d kd1(1, 1, 0);
    Vector3d ks1(1, 1, 1);
    Vector3d km1(0, 0, 0);
    Material m1(ks1, kd1, ka1, km1, 50);
    Vector3d sphere2(-2, 2, -15);
    scene->addSphere(sphere2, 1, m1);


    Vector3d ka2(0.1, 0.1, 0.1);
    Vector3d kd2(0, 1, 1);
    Vector3d ks2(1, 1, 1);
    Vector3d km2(0, 0, 0);
    Material m2(ks2, kd2, ka2, km2, 50);
    Vector3d sphere3(-2, -2, -15);
    scene->addSphere(sphere3, 1, m2);


    Vector3d ka3(0.1, 0.1, 0.1);
    Vector3d kd3(0.1, 0.1, 0.1);
    Vector3d ks3(1, 1, 1);
    Vector3d km3(1, 1, 1);
    Material m3(ks3, kd3, ka3, km3, 50);
    //scene->addObj("triangle.obj", m3);
    Vector3d v1(5, 5, -17);
    Vector3d v2(1, 4, -20);
    Vector3d v3(6, -1, -20);

    scene->addTriangle(v1, v2, v3, m3);

    vector<Vector3d* > points;
    points.push_back(&v2);
    points.push_back(&v1);
    points.push_back(&v3);
    points.push_back(&p1);
    points.push_back(&p2);
    points.push_back(&p3);
    points.push_back(&p4);
    scene->addPolygon(points, m3);

    Vector3d intensity(1, 1, 1);
    Vector3d dir(-0.57735027,0.57735027,0.57735027);
    scene->addDirectionalLight(dir, intensity);

    Vector3d intensity2(0, 0, 1);
    Vector3d dir2(-0.57735027,-0.57735027,0.57735027);
    scene->addDirectionalLight(dir2, intensity2);

    Sampler* cs = new CenterSampler();
    RayTracer rt(cs);
    rt.rayTrace(scene, "test.ppm", 1);*/
    return 0;
}