#include "geometry.cpp"
#include "surfaces.cpp"
long double max_k(const Ray& ray, const std::vector<Surface*>& obj, const std::vector<Lamp>& lamps, long double r) {
    long double min_dist = pow(10, 100);
    for (int i = 0; i < obj.size(); i++) {
        std::vector<Point> cur = obj[i]->CheckInters(ray.l);
        if (!cur.empty()) {
            min_dist = std::min(min_dist, (cur[0] - ray.l.r_0).Module());
        }
    }
    for (int i = 0; i < lamps.size(); i++) {
        std::vector<Point> cur = lamps[i].pos->CheckInters(ray.l);
        if (!cur.empty()) {
            min_dist = std::min(min_dist, (cur[0] - ray.l.r_0).Module());
        }
    }
    if (min_dist == pow(10, 100)) { return 0; }
    long double res = min_dist / r;
    return res;
}
long double Power(long double k, const Ray& ray, Atmosphere* fog, Sphere* sun) {
    const Vector tau = ray.l.tau;
    long double r = fog->r;
    long double Rz = fog->pos->radius;
    long double H = (sun->pos - ray.l.r_0).Module();
    Vector l_ = sun->pos - ray.l.r_0;
    Vector l = l_ / l_.Module();
    Point start = ray.l.r_0 + ray.l.tau * k * r;
    Line final(start, l);
    std::vector<Point> inter = fog->CheckInters(final);
    if (inter.empty()) { return 0; }
    long double res = (inter[0] - start).Module() / r + k;
    return res;
}
long double Func(long double k, const Ray& ray, Atmosphere* fog, Sphere* sun, int q) {
    long double pp = fog->conc * pow(fog->r, 3) * 3.1415926 * fog->probability[q];
    long double r = fog->r;
    long double R = sun->radius;
    long double H = (sun->pos - ray.l.r_0).Module();
    long double mult = (R * R) / 2 / (H * H) * (0.5 + k * r / 2 / H);
    return pow((1 - pp), Power(k, ray, fog, sun)) * pp * mult;
}
Point RandomPoint() {
    Sphere sph(Point({0, 0, 0}), 1000, {0, 0, 0}, 0, 0);
    Point p = Point({-100000, 0, 0});
    while (!sph.PointIn(p) || (p - Point({0, 0, 0})).Module() == 0) {
        long double x = static_cast<long double>(rand() % 2000 - 1000);
        long double y = static_cast<long double>(rand() % 2000 - 1000);
        long double z = static_cast<long double>(rand() % 2000 - 1000);
        p = Point({x, y, z});
    }
    long double dist = (p - Point({0, 0, 0})).Module();
    Point res({p.x / dist, p.y / dist, p.z / dist});
    return res;
}
std::vector<std::vector<long double>> FindProb(const Ray& ray, Atmosphere* fog, Sphere* sun, const std::vector<Surface*>& obj, const std::vector<Lamp>& lamps, int q, int depth) {
    std::vector<std::vector<long double>> S(depth + 1, std::vector<long double>(3, 0));
    long double m = max_k(ray, obj, lamps, fog->r);
    if (depth >= 2) {
        int prec = 6;
        for (int i = 1; i < prec; i++) {
            int am_of_rays = 10;
            for (int qq = 0; qq < 3; qq++) {
                for (int j = 0; j < am_of_rays; j++) {
                    Point p = ray.l.r_0 + i * 1.0 / prec * m * fog->r * ray.l.tau;
                    Point next = RandomPoint();
                    Vector new_tau = next - Point({0, 0, 0});
                    Line new_line = Line(p, new_tau);
                    std::vector<std::vector<long double>> add(depth, std::vector<long double>(3, 0));
                    if (fog->PointIn(p)) { add = FindProb(Ray(new_line, 1, 1, 1), fog, sun, obj, lamps, q, depth - 1); }
                    long double pp = fog->conc * pow(fog->r, 3) * 3.1415926 * fog->probability[qq];
                    for (int w = 1; w <= depth; w++) {
                        long double mult = 1.0 / am_of_rays * pow(new_tau.Scal(ray.l.tau), 4);
                        if (w == 1) { mult = pow((sun->radius / 2 / (ray.l.r_0 - sun->pos).Module()), 2) * pow(new_tau.Scal(ray.l.tau), 2); }
                        long double modifier = add[w - 1][qq] * mult * m / (prec - 1) * pow((1 - pp), i * 1.0 / prec * m) * pp;
                        S[w][qq] += modifier;
                    }
                }
                Point p = ray.l.r_0 + i * 1.0 / prec * m * fog->r * ray.l.tau;
                Point next = sun->pos;
                Vector new_tau = next - Point({0, 0, 0});
                new_tau = new_tau / new_tau.Module();
                Line new_line = Line(p, new_tau);
                std::vector<std::vector<long double>> add(depth, std::vector<long double>(3, 0));
                if (fog->PointIn(p)) { add = FindProb(Ray(new_line, 1, 1, 1), fog, sun, obj, lamps, q, depth - 1); }
                long double pp = fog->conc * pow(fog->r, 3) * 3.1415926 * fog->probability[qq];
                for (int w = 1; w <= depth; w++) {
                    long double mult = pow((sun->radius / 2 / (ray.l.r_0 - sun->pos).Module()), 2) * pow(new_tau.Scal(ray.l.tau), 2);
                    S[w][qq] += add[w - 1][qq] * mult * m / (prec - 1) * pow((1 - pp), i * 1.0 / prec * m) * pp;
                }
            }
        }
        depth = 0;
    }
    if (depth == 1) {
        bool not_int_abs = true;
        int prec = 60;
        for (int i = 1; i < prec; i++) {
            Point p = ray.l.r_0 + i * 1.0 / prec * m * fog->r * ray.l.tau;
            Line new_line = Line(p, (sun->pos - p) / (sun->pos - p).Module());
            bool not_int = true;
            for (int j = 0; j < obj.size() && not_int; j++) {
                not_int = not_int && (dynamic_cast<Atmosphere*>(obj[j]) || obj[j]->CheckInters(new_line).empty());
                std::vector<Point> points = obj[j]->CheckInters(ray.l);
                not_int = not_int && (points.empty() || (points[0] - ray.l.r_0).Module() > i * 1.0 / prec * m * fog->r);
            }
            for (int j = 0; j < lamps.size() && not_int; j++) {
                not_int = not_int && (j == q || lamps[j].pos->CheckInters(new_line).empty());
            }
            if (not_int) {
                for (int qq = 0; qq < 3; qq++) {
                    long double f = Func(i * 1.0 / prec * m, ray, fog, sun, qq) * pow(new_line.tau.Scal(ray.l.tau), 2);
                    S[1][qq] += f * m / (prec - 1);
                    long double pp = fog->conc * pow(fog->r, 3) * 3.1415926 * fog->probability[qq];
                }
            }
            if (i == 0) { not_int_abs = not_int; }
        }
        std::vector<Point> fix = sun->CheckInters(ray.l);
        std::vector<Point> fix_ = fog->pos->CheckInters(ray.l);
        if (!fix.empty() && !fix_.empty() && not_int_abs) {
            for (int qq = 0; qq < 3; qq++) {
                S[0][qq] = 0;
                //S[0][q] = pow((1 - fog->conc * pow(fog->r, 3) * 3.1415926 * fog->probability[q]), (fix_[0] - ray.l.r_0).Module() / fog->r);
            }
        }
    }
    if (depth == 0) {
        Point p = ray.l.r_0;
        Line new_line = Line(p, (sun->pos - p) / (sun->pos - p).Module());
        bool not_int_abs = true;
        for (int j = 0; j < obj.size() && not_int_abs; j++) {
            not_int_abs = not_int_abs && (dynamic_cast<Atmosphere*>(obj[j]) || obj[j]->CheckInters(new_line).empty());
        }
        for (int j = 0; j < lamps.size() && not_int_abs; j++) {
            not_int_abs = not_int_abs && (j == q || lamps[j].pos->CheckInters(new_line).empty());
        }
        std::vector<Point> fix = sun->CheckInters(ray.l);
        std::vector<Point> fix_ = fog->pos->CheckInters(ray.l);
        if (!fix.empty() && !fix_.empty() && not_int_abs) {
            for (int qq = 0; qq < 3; qq++) {
                S[0][qq] += pow(10, -4) * pow((1 - fog->conc * pow(fog->r, 3) * 3.1415926 * fog->probability[qq]), (fix_[0] - ray.l.r_0).Module() / fog->r);
            }
        }
    }
    return S;
}