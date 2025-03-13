#ifndef geometry
#define geometry
#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <ctime>
int depth = 5;
struct Vector {
    long double x;
    long double y;
    long double z;
    Vector() : x(0), y(0), z(0){};
    Vector(const std::vector<long double>& src) : x(src[0]), y(src[1]), z(src[2]) {};
    long double Module() const {
        return pow(x*x+y*y+z*z, 0.5);
    }
    long double ModuleSqr() const  {
        return x*x+y*y+z*z;
    }
    long double Scal(const Vector& other) const  {
        return x*other.x+y*other.y+z*other.z;
    }
    Vector Vect(const Vector& other) const  {
        return Vector({y*other.z-z*other.y, z*other.x-x*other.z, x*other.y-y*other.x});
    }
    Vector operator/ (long double other) const {
        return Vector({x / other, y / other, z / other});
    }
    Vector operator+ (const Vector& other) const {
        return Vector({x + other.x, y + other.y, z + other.z});
    }
    Vector operator- (const Vector& other) const {
        return Vector({x - other.x, y - other.y, z - other.z});
    }
    Vector operator- () const {
        return Vector({-x, -y, -z});
    }
};
const Vector operator*(const Vector& v, long double x) {
    return Vector({v.x * x, v.y * x, v.z * x});
}
const Vector operator*(long double x, const Vector& v) {
    return Vector({v.x * x, v.y * x, v.z * x});
}
std::istream& operator>> (std::istream& is, Vector& v) {
    is >> v.x >> v.y >> v.z;
    return is;
}
std::ostream& operator<< (std::ostream& os, const Vector& v) {
    os << '{' << v.x << ", " << v.y << ", " << v.z << '}';
    return os;
}
struct Point {
    long double x;
    long double y;
    long double z;
    Point() : x(0), y(0), z(0){};
    Point(const std::vector<long double>& src) : x(src[0]), y(src[1]), z(src[2]) {};
    Point operator+(const Vector& v) const {
        return Point({x+v.x, y+v.y, z+v.z});
    }
    Point operator-(const Vector& v) const {
        return Point({x-v.x, y-v.y, z-v.z});
    }
    Vector operator-(const Point& v) const {
        return Vector({x-v.x, y-v.y, z-v.z});
    }
};
std::istream& operator>> (std::istream& is, Point& v) {
    is >> v.x >> v.y >> v.z;
    return is;
}
std::ostream& operator<< (std::ostream& os, const Point& v) {
    os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
    return os;
}
struct Line {
    Point r_0;
    Vector tau;
    Line(const Point& src_1, const Vector& src_2) : r_0(src_1), tau(src_2) {};
    bool PointOn(const Point& p)  const {
        Vector check({p.x-r_0.x, p.y-r_0.y, p.z-r_0.z});
        return (tau.Vect(check).Module() <= pow(10, -11));
    }
};
std::istream& operator>> (std::istream& is, Line& v) {
    is >> v.r_0 >> v.tau;
    return is;
}
std::ostream& operator<< (std::ostream& os, const Line& v) {
    os << v.r_0 << " + t * " << v.tau;
    return os;
}
struct Surface {
    std::vector<long double> reflection;
    std::vector<long double> refraction;
    long double roughness;
    long double n;
    virtual std::vector<Point> CheckInters(const Line& line) const = 0;
    virtual Vector Normal(const Point& point) const = 0;
    virtual bool PointIn(const Point& point) const = 0;
    virtual ~Surface() {};
};
void ReflectLine(Line& line, const Surface* surf, const Vector& n) {
    std::vector<Point> points = surf->CheckInters(line);
    if (points.empty()) {
        return;
    }
    Point point = points[0];
    Vector a = line.tau - n * line.tau.Scal(n);
    line.tau = line.tau - 2 * a;
    line.tau = line.tau * (-1);
    line.r_0 = point + line.tau * 0.01;
}
bool RefractLine(Line& line, const Surface* surf, const Vector& n) {
    std::vector<Point> points = surf->CheckInters(line);
    if (points.empty()) {
        return false;
    }
    Point point = points[0];
    Vector a = line.tau - n * line.tau.Scal(n);
    if (surf->PointIn(line.r_0)) {
        if ((a * surf->n).ModuleSqr() < 1) {
            line.tau = a * surf->n + n * pow(1 - (a * surf->n).ModuleSqr(), 0.5);
        } else {
            return false;
        }
    } else {
        line.tau = a / surf->n - n * pow(1 - (a / surf->n).ModuleSqr(), 0.5);
    }
    line.tau = line.tau / line.tau.Module();
    line.r_0 = point + line.tau * 0.01;
    return true;
}
int Rand() {
    return rand();
}
struct Ray {
    Line l;
    long double red;
    long double green;
    long double blue;
    int is_copied;
    int max_copied;
    Ray(const Line& L, long double r, long double g, long double b) : l(L), red(r), green(g), blue(b), is_copied(0), max_copied(depth) {};
    bool Reflect(const Surface* s, int prec) {
        Line save = l;
        std::vector<Point> p = s->CheckInters(save);
        if (!p.empty()) {
            is_copied++;
            Vector n = s->Normal(p[0]);
            if (n.Module() == 0) {
                return false;
            }
            if (s->roughness != 0) {
                Vector add_f = n.Vect(Vector({1, 0, 0}));
                Vector add_s = Vector({1, 0, 0});
                if (add_f.Module() != 0) {
                    add_s = n.Vect(Vector({0, 1, 0}));
                    if (add_s.Module() == 0) {
                        add_s = n.Vect(Vector({0, 0, 1}));
                    }
                } else {
                    add_f = n.Vect(Vector({0, 1, 0}));
                    add_s = n.Vect(Vector({0, 0, 1}));
                }
                add_f = add_f / add_f.Module();
                add_s = add_s / add_s.Module();
                int rnd1 = Rand() % 80;
                int rnd2 = Rand() % prec;
                long double alpha = 2 * 3.1415926 * rnd1 / 80;
                Vector adding = cos(alpha) * add_f + sin(alpha) * add_s;
                n = n * pow(1 - pow(s->roughness * rnd2 / prec, 2), 0.5) + adding * s->roughness / prec;
                n = n / n.Module();
            }
            Line check = save;
            ReflectLine(check, s, n);
            if (check.tau.Scal(s->Normal(p[0])) * save.tau.Scal(s->Normal(p[0])) >= 0) {
                return false;
            }
            ReflectLine(l, s, n);
            red *= s->reflection[0] * save.tau.Scal(n) / save.tau.Module();
            green *= s->reflection[1] * save.tau.Scal(n) / save.tau.Module();
            blue *= s->reflection[2] * save.tau.Scal(n) / save.tau.Module();
            red = ((red > 0) ? red : -red);
            green = ((green > 0) ? green : -green);
            blue = ((blue > 0) ? blue : -blue);
            return true;
        }
        return false;
    }
    bool Refract(const Surface* s) {
        Line save = l;
        std::vector<Point> p = s->CheckInters(save);
        bool result = true;
        if (!p.empty()) {
            is_copied++;
            Vector n = s->Normal(p[0]);
            if (n.Module() == 0) {
                return false;
            }
            result = RefractLine(l, s, n);
            if (result) {
                red *= s->refraction[0] * save.tau.Scal(n) / save.tau.Module();
                green *= s->refraction[1] * save.tau.Scal(n) / save.tau.Module();
                blue *= s->refraction[2] * save.tau.Scal(n) / save.tau.Module();
                red = ((red > 0) ? red : -red);
                green = ((green > 0) ? green : -green);
                blue = ((blue > 0) ? blue : -blue);
                return true;
            }
            return false;
        }
        return false;
    }
};
Line Transform(const Point& start, long double phi, long double theta) {
    Vector a = Vector({cos(phi) * cos(theta), sin(theta), cos(theta)*sin(phi)});
    Line res = Line(start, a);
    return res;
}
struct Pixel {
    long double red;
    long double green;
    long double blue;
    Pixel() : red(0), green(0), blue(0){};
    Pixel(long double r, long double g, long double b) : red(r), green(g), blue(b){};
    Pixel operator+ (const Pixel& other) const {
        return Pixel(red + other.red, green + other.green, blue + other.blue);
    }
};
struct Lamp {
    Surface* pos;
    long double red;
    long double green;
    long double blue;
    std::vector<long double> intens;
    Lamp(Surface* p, long double r, long double g, long double b) : pos(p), red(r), green(g), blue(b) {
        intens.push_back(r);
        intens.push_back(g);
        intens.push_back(b);
    }
};
#endif