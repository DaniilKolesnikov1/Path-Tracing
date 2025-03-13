#ifndef surfaces
#define surfaces
#include "geometry.cpp"
struct FuncSurface : public Surface {
    long double C;
    long double (*f) (long double, long double, long double);
    FuncSurface(long double c, long double (*func) (long double, long double, long double), const std::vector<long double>& a, long double r, long double m, const std::vector<long double>& b = {0, 0, 0}) : C(c), f(func) {
        roughness = r;
        n = m;
        reflection = a;
        refraction = b;
    }
    long double Function(const Point& p) const {
        return f(p.x, p.y, p.z);
    }
    bool PointOn(const Point& p) const {
        return (pow(Function(p) - C, 2) <= pow(10, -10));
    }
    bool PointIn(const Point& p) const {
        return (Function(p) < C);
    }
    std::pair<Vector, Vector> Tanget(const Point& p) const {
        long double fx = (Function(p + Vector({0.001, 0, 0})) - Function(p)) * 1000;
        long double fy = (Function(p + Vector({0, 0.001, 0})) - Function(p)) * 1000;
        long double fz = (Function(p + Vector({0, 0, 0.001})) - Function(p)) * 1000;
        Vector a({0, 0 ,0});
        Vector b({0, 0 ,0});
        if (fz != 0) {
            a = Vector({1, 0, -fx/fz});
            b = Vector({0, 1, -fy/fz});
        } else {
            if (fy != 0) {
                a = Vector({1, -fx/fy, 0});
                b = Vector({0, 0, 1});
            } else {
                a = Vector({0, 1, 0});
                b = Vector({0, 0, 1});
            }
        }
        return std::make_pair(a, b);
    }
    Vector Normal(const Point& p) const {
        std::pair<Vector, Vector> a = this->Tanget(p);
        Vector res = a.first.Vect(a.second);
        res = res * (1 / res.Module());
        return res;
    }
    std::vector<Point> CheckInters(const Line& line) const {
        int n = 1000;
        long double prec = 0.01;
        std::vector<long double> ff(n + 1);
        std::vector<int> res;
        for (int i = 0; i <= n && res.size() == 0; i++) {
            ff[i] = Function(line.r_0 + i * prec * line.tau);
            if (i != 0 && (ff[i - 1] - C) * (ff[i] - C) < 0 || ff[i] == C) {
                res.push_back(i);
                break;
            }
        }
        std::vector<Point> ans(res.size(), Point({0, 0, 0}));
        for (uint i = 0; i < res.size(); i++) {
            Point n_p = line.r_0+(res[i])*prec*line.tau;
            ans[i] = n_p;
        }
        return ans;
    }
    ~FuncSurface() {};
};
struct Plane : public Surface {
    long double C;
    Vector normal;
    Plane() : C(0), normal(){};
    Plane(long double c, const Vector& k, const std::vector<long double>& a, long double r, long double m, const std::vector<long double>& b = {0, 0, 0}) : C(c), normal(k / k.Module()) {
        roughness = r;
        n = m;
        reflection = a;
        refraction = b;
    }
    Plane(const Point& c, const Vector& k, const std::vector<long double>& a, long double r, long double m, const std::vector<long double>& b = {0, 0, 0}) : C(k.Scal(c-Point({0,0,0}))), normal(k / k.Module()) {
        roughness = r;
        n = m;
        reflection = a;
        refraction = b;
    }
    long double Function(const Point& p) const {
        return normal.Scal(Vector({p.x, p.y, p.z}));
    }
    Vector Normal(const Point& p) const {
        return normal;
    }
    bool PointIn(const Point& p) const {
        return Function(p) - C <= pow(10, -12);
    }
    std::vector<Point> CheckInters(const Line& line) const {
        long double der = Function(Point({0, 0, 0}) + line.tau);
        long double t = -1;
        if (der != 0) {
            t = (C - Function(line.r_0)) / der;
            if (t >= 0) {
                return {line.r_0 + t * line.tau};
            } else {
                return {};
            }
        }
        return {};
    }
    ~Plane() {};
};
struct Cube : public Surface {
    Point left_down_far;
    std::vector<Vector> i;
    long double a;
    std::vector<Plane> ps;
    Cube (const Point& p, const std::vector<Vector>& b, const std::vector<long double>& al, long double len, long double r, long double m, const std::vector<long double>& ra = {0, 0, 0}) : left_down_far(p) {
        i = std::vector<Vector>(3);
        i[0] = b[0] / b[0].Module();
        i[1] = b[1] / b[1].Module();
        i[2] = b[2] / b[2].Module();
        reflection = al;
        refraction = ra;
        a = len;
        roughness = r;
        n = m;
        ps = std::vector<Plane>(6);
        ps[0] = Plane(left_down_far, -i[0].Vect(i[1]), {0,0,0}, roughness, n);
        ps[1] = Plane(left_down_far, -i[1].Vect(i[2]), {0,0,0}, roughness, n);
        ps[2] = Plane(left_down_far, -i[2].Vect(i[0]), {0,0,0}, roughness, n);
        ps[3] = Plane(left_down_far + (i[0] + i[1] + i[2]) * a, i[0].Vect(i[1]), {0,0,0}, roughness, n);
        ps[4] = Plane(left_down_far + (i[0] + i[1] + i[2]) * a, i[1].Vect(i[2]), {0,0,0}, roughness, n);
        ps[5] = Plane(left_down_far + (i[0] + i[1] + i[2]) * a, i[2].Vect(i[0]), {0,0,0}, roughness, n);
    }
    Vector Normal(const Point& p) const {
        std::vector<long double> x(6);
        for (int q = 0; q < 6; q++) {
            x[q] = ps[q].Function(p) - ps[q].C;
            x[q] = ((x[q] > 0) ? x[q] : -x[q]);
        }
        std::vector<long double> y(3);
        for (int q = 0; q < 3; q++) {
            y[q] = std::min(x[q], x[q + 3]);
        }
        if (y[0] <= y[1] && y[0] <= y[2]) {
            if (x[0] < x[3]) {
                return ps[0].Normal(p);
            } else {
                return ps[3].Normal(p);
            }
        }
        if (y[1] <= y[2] && y[1] <= y[0]) {
            if (x[1] < x[4]) {
                return ps[1].Normal(p);
            } else {
                return ps[4].Normal(p);
            }
        }
        if (x[2] < x[5]) {
            return ps[2].Normal(p);
        }
        return ps[5].Normal(p);
    }
    std::vector<Point> CheckInters(const Line& line) const {
        std::vector<std::vector<Point>> points(6);
        for (int q = 0; q < 6; q++) {
            points[q] = ps[q].CheckInters(line);
        }
        Point ans = line.r_0 - line.tau;
        long double min = -1;
        for (int q = 0; q < 6; q++) {
            for (int w = 0; w < points[q].size(); w++) {
                long double t = (points[q][w] - line.r_0).Module() / line.tau.Module();
                if (t >= 0 && (min < 0 || t < min)) {
                    Vector check_f = points[q][w] - left_down_far;
                    Vector check_s = check_f - a * (i[0] + i[1] + i[2]);
                    if (q < 3 && check_f.Scal(i[(q) % 3]) <= a && check_f.Scal(i[(q) % 3]) >= 0 && check_f.Scal(i[(q + 1) % 3]) <= a && check_f.Scal(i[(q + 1) % 3]) >= 0) {
                        ans = points[q][w];
                        min = t;
                    }
                    if (q >= 3 && check_s.Scal(i[(q) % 3]) >= -a && check_s.Scal(i[(q) % 3]) <= 0 && check_s.Scal(i[(q + 1) % 3]) >= -a && check_s.Scal(i[(q + 1) % 3]) <= 0) {
                        ans = points[q][w];
                        min = t;
                    }
                }
            }
        }
        if (min >= 0) {
            return {ans};
        }
        return {};
    }
    bool PointIn(const Point& p) const {
        bool ans = true;
        for (int q = 0; q < 3; q++) {
            ans = ans && ps[q].PointIn(p) && ps[q + 3].PointIn(p);
        }
        return ans;
    }
    ~Cube() {};
};
struct Sphere : public Surface {
    long double radius;
    Point pos;
    Sphere() : radius(0), pos(){};
    Sphere(const Point& k, long double c, const std::vector<long double>& a, long double r, long double m, const std::vector<long double>& b = {0, 0, 0}) {
        radius = c;
        pos = k;
        roughness = r;
        n = m;
        reflection = a;
        refraction = b;
    }
    long double Function(const Point& p) const {
        return (p - pos).Module();
    }
    Vector Normal(const Point& p) const {
        return (p - pos) / (p - pos).Module();
    }
    bool PointIn(const Point& p) const {
        return Function(p) - radius <= pow(10, -12);
    }
    std::vector<Point> CheckInters(const Line& line) const {
        Vector add = pos-line.r_0;
        long double dis = -line.tau.Vect(add).ModuleSqr() + pow(radius, 2);
        if (dis < 0) { return {}; }
        long double t_f = -1;
        long double t_s = -1;
        t_f = line.tau.Scal(add) - pow(-line.tau.Vect(add).ModuleSqr() + pow(radius, 2), 0.5);
        t_s = line.tau.Scal(add) + pow(-line.tau.Vect(add).ModuleSqr() + pow(radius, 2), 0.5);
        if (t_f >= -pow(10, -3)) {
            return {line.r_0 + t_f * line.tau, line.r_0 + t_s * line.tau};
        }
        if (t_f < 0 && t_s >= -pow(10, -3)) {
            return {line.r_0 + t_s * line.tau};
        }
        if (PointIn(line.r_0)) { std::cout << radius << ' ' << Function(line.r_0) << ' ' << t_s << ' ' << line.tau.Scal(add) << ' ' << pow(dis, 0.5) << '\n'; }
        return {};
    }
    ~Sphere() {};
};
struct Cilindre : public Surface {
    Point down;
    long double radius;
    long double h;
    Vector tau;
    Cilindre(const Point& d, const Vector& t, long double rad, long double hh, const std::vector<long double>& a, long double r, long double m, const std::vector<long double>& b = {0, 0, 0}) {
        down = d;
        tau = t;
        radius = rad;
        h = hh;
        reflection = a;
        refraction = b;
        roughness = r;
        n = m;
    }
    std::vector<Point> CheckInters(const Line& line) const {
        Vector temp = (line.r_0 - down).Vect(tau);
        long double disc = pow(radius, 2) * line.tau.Vect(tau).ModuleSqr();
        disc -= (line.tau.Vect(tau)).Vect(temp).ModuleSqr();
        if (disc < 0) { return {}; }
        long double t0 = -line.tau.Vect(tau).Scal(temp);
        long double t1 = t0 + pow(disc, 0.5);
        long double t2 = t0 - pow(disc, 0.5);
        if (t1 < 0) { return {}; }
        if (t2 >= 0) {
            Point res = line.r_0 + t2 * line.tau;
            if ((res - down).Scal(tau) < h && (res - down).Scal(tau) > 0) {
                return {res};
            }
            return {};
        }
        Point res = line.r_0 + t1 * line.tau;
        if ((res - down).Scal(tau) < h && (res - down).Scal(tau) > 0) {
            return {res};
        }
        return {};
    }
    Vector Normal(const Point& point) const {
        Vector temp = point - down;
        return (temp - temp.Scal(tau) * tau) / (temp - temp.Scal(tau) * tau).Module();
    }
    bool PointIn(const Point& point) const {
        Vector temp = point - down;
        return (temp - temp.Scal(tau) * tau).Module() - radius <= pow(10, -12);;
    }
    ~Cilindre() {};
};
struct Lense : public Surface {
    Sphere* first;
    Sphere* second;
    Lense(Sphere* f, Sphere* s, const std::vector<long double>& a, long double r, long double m, const std::vector<long double>& b = {0, 0, 0}) {
        first = f;
        second = s;
        reflection = a;
        refraction = b;
        roughness = r;
        n = m;
        first->n = m;
        second->n = m;
        first->reflection = a;
        second->reflection = a;
        first->refraction = b;
        second->refraction = b;
        first->roughness = roughness;
        second->roughness = roughness;
    }
    std::vector<Point> CheckInters(const Line& line) const {
        std::vector<Point> f = first->CheckInters(line);
        std::vector<Point> s = second->CheckInters(line);
        bool FPI = first->PointIn(line.r_0);
        bool SPI = second->PointIn(line.r_0);
        if (FPI && !SPI) {
            if (!s.empty() && (s[0] - line.r_0).ModuleSqr() < (f[0] - line.r_0).ModuleSqr()) {
                return s;
            }
        }
        if (!FPI && SPI) {
            if (!f.empty() && (f[0] - line.r_0).ModuleSqr() < (s[0] - line.r_0).ModuleSqr()) {
                return f;
            }
        }
        if (FPI && SPI) {
            if ((s[0] - line.r_0).ModuleSqr() < (f[0] - line.r_0).ModuleSqr()) {
                return s;
            } else {
                return f;
            }
        }
        if (!FPI && !SPI) {
            if (!f.empty() && !s.empty()) {
                long double t_f = (f[0] - line.r_0).Module();
                long double t_s = (s[0] - line.r_0).Module();
                long double t = std::max(t_f, t_s);
                if (first->PointIn(line.r_0 + t * line.tau) && second->PointIn(line.r_0 + line.tau * t)) {
                    return {line.r_0 + t * line.tau};
                }
            }
        }
        return {};
    }
    Vector Normal(const Point& point) const {
        if (first->Function(point) / first->radius < second->Function(point) / second->radius) {
            return second->Normal(point);
        }
        return first->Normal(point);
    }
    bool PointIn(const Point& point) const {
        return first->PointIn(point) && second->PointIn(point);
    }
    ~Lense() {
        delete first;
        delete second;
    };
};
int foggy = 0;
struct Fog : public Surface {
    Surface* pos;
    long double conc;
    long double r;
    long double probability;
    Fog(Surface* k, long double h, long double c, long double p, const std::vector<long double>& a) {
        pos = k;
        conc = h;
        r = c;
        roughness = 0;
        n = 30;
        reflection = a;
        probability = p;
        refraction = {0, 0, 0};
    }
    Vector Normal(const Point& p) const {
        int sign = ((Rand() % 2 == 0) ? 1 : -1);
        int f = (Rand() % 1000) * sign;
        sign = ((Rand() % 2 == 0) ? 1 : -1);
        int s = (Rand() % 1000) * sign;
        sign = ((Rand() % 2 == 0) ? 1 : -1);
        int t = (Rand() % 1000) * sign;
        Vector res = Vector({static_cast<long double>(f), static_cast<long double>(s), static_cast<long double>(t)});
        res = res / res.Module();
        return res;
    }
    bool PointIn(const Point& p) const {
        return pos->PointIn(p);
    }
    std::vector<Point> CheckInters(const Line& line) const {
        std::vector<Point> res = pos->CheckInters(line);
        if (res.size() == 0) {
            return {};
        }
        if (!pos->PointIn(line.r_0)) {
            return res;
        }
        int q = Rand() % 10000000;
        long double harsh = conc * pow(r, 3) * 3.1415926;
        long double p = q * 0.0000001 * harsh;
        long double dist = r * log(p / harsh) / log(1 - harsh) / probability;
        if (pos->PointIn(line.r_0 + dist * line.tau / line.tau.Module())) {
            return {line.r_0 + dist * line.tau / line.tau.Module()};
        }
        return {};
    }
    ~Fog() {
        delete pos;
    }
};
struct Atmosphere : public Surface {
    Sphere* pos;
    long double conc;
    long double r;
    std::vector<long double> probability;
    Atmosphere(Sphere* k, long double h, long double c, const std::vector<long double>& p, long double m, const std::vector<long double>& a, const std::vector<long double>& b={1,1,1}) {
        pos = k;
        conc = h;
        r = c;
        roughness = 0;
        n = m;
        reflection = a;
        probability = p;
        refraction = b;
    }
    Vector Normal(const Point& p) const {
        return pos->Normal(p);
    }
    bool PointIn(const Point& p) const {
        return pos->PointIn(p);
    }
    std::vector<Point> CheckInters(const Line& line) const {
        return pos->CheckInters(line);
    }
    ~Atmosphere() {
        delete pos;
    }
};
#endif