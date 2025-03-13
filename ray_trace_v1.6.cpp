#include <iostream>
#include <cmath>
#include "mpi.h"
#include "geometry.cpp"
#include "surfaces.cpp"
#include "atmo_func.cpp"
int precise_of_atm = 1;
bool Modify_Res_By_Atm(std::vector<long double>& res, const std::vector<int>& atmos, const std::vector<Surface*>& objects, const std::vector<Lamp>& lamps, const std::vector<int>& sun, const Ray& save) {
    bool inters_sun = false;
    for (int j = 0; j < atmos.size(); j++) {
        if (objects[atmos[j]]->PointIn(save.l.r_0)) {
            for (int i = 0; i < sun.size(); i++) {
                inters_sun = inters_sun || (!lamps[sun[i]].pos->CheckInters(save.l).empty());
                Atmosphere* cur = dynamic_cast<Atmosphere*>(objects[atmos[j]]);
                std::vector<std::vector<long double>> add = FindProb(save, cur, dynamic_cast<Sphere*>(lamps[sun[i]].pos), objects, lamps, i, precise_of_atm);
                for (int ii = 0; ii <= precise_of_atm; ii++) {
                    std::vector<long double> mult = {pow(cur->reflection[0], ii), pow(cur->reflection[1], ii), pow(cur->reflection[2], ii)};
                    if (ii == 0) { mult = {(1 - pow(cur->reflection[0], ii + 1)) / (1 - cur->reflection[0]),
                                           (1 - pow(cur->reflection[1], ii + 1)) / (1 - cur->reflection[1]),
                                           (1 - pow(cur->reflection[2], ii + 1)) / (1 - cur->reflection[2])};}
                    long double add_r = lamps[sun[i]].red * save.red * mult[0] * add[ii][0];
                    long double add_g = lamps[sun[i]].green * save.green * mult[1] * add[ii][1];
                    long double add_b = lamps[sun[i]].blue * save.blue * mult[2] * add[ii][2];
                    res[0] = pow(pow(add_r, 2) + pow(res[0], 2), 0.5);
                    res[1] = pow(pow(add_g, 2) + pow(res[1], 2), 0.5);
                    res[2] = pow(pow(add_b, 2) + pow(res[2], 2), 0.5);
                }
            }
        }
    }
    return inters_sun;
}
std::vector<long double> Trace(const Ray& ray, int n, int m, const Point& pov, const std::vector<Lamp>& lamps, const std::vector<Surface*>& objects, const std::vector<int>& sun, const std::vector<int>& atmos) {
    std::vector<long double> res = {0.4, 0.4, 0.4};
    if (ray.is_copied >= ray.max_copied) {
        return res;
    }
    if (ray.red == 0 && ray.green == 0 && ray.blue == 0) {
        return res;
    }
    std::vector<std::vector<Point>> points(objects.size() + lamps.size());
    long double nearest_dist = -1;
    int min_dist = -1;
    for (int q = 0; q < objects.size(); q++) {
        points[q] = objects[q]->CheckInters(ray.l);
        if (!points[q].empty() && (nearest_dist > (points[q][0] - ray.l.r_0).Module() || nearest_dist == -1)) {
            nearest_dist = (points[q][0] - ray.l.r_0).Module();
            min_dist = q;
        }
    }
    for (int q = 0; q < lamps.size(); q++) {
        points[q + objects.size()] = lamps[q].pos->CheckInters(ray.l);
        if (!points[q + objects.size()].empty() && (nearest_dist > (points[q + objects.size()][0] - ray.l.r_0).Module() || nearest_dist == -1)) {
            nearest_dist = (points[q + objects.size()][0] - ray.l.r_0).Module();
            min_dist = q + objects.size();
        }
    }
    if (min_dist != -1 && min_dist < objects.size()) {
        int q = min_dist;
        Ray save = ray;
        if (objects[q]->roughness == 0) {
            Ray first = save;
            bool inters_sun = Modify_Res_By_Atm(res, atmos, objects, lamps, sun, save);
            if (!dynamic_cast<Atmosphere*>(objects[q])) {
                first.Reflect(objects[q], 1);
                std::vector<long double> add = Trace(first, n, m, pov, lamps, objects, sun, atmos);
                res[0] = pow(pow(res[0], 2) + pow(add[0], 2), 0.5);
                res[1] = pow(pow(res[1], 2) + pow(add[1], 2), 0.5);
                res[2] = pow(pow(res[2], 2) + pow(add[2], 2), 0.5);
            } else {
                bool refract = first.Refract(objects[q]);
                if (refract && objects[q]->PointIn(save.l.r_0)) {
                    std::vector<long double> add = Trace(first, n, m, pov, lamps, objects, sun, atmos);
                    for (int qq = 0; qq < 3; qq++) {
                        Atmosphere* fog = dynamic_cast<Atmosphere*>(objects[q]);
                        long double pp = fog->conc * pow(fog->r, 3) * 3.1415926 * fog->probability[qq];
                        long double prob = pow(1 - pp, nearest_dist / fog->r);
                        add[qq] *= prob;
                        if (inters_sun) { add[qq] *= pow(10, -4); }
                    }
                    res[0] = pow(pow(res[0], 2) + pow(add[0], 2), 0.5);
                    res[1] = pow(pow(res[1], 2) + pow(add[1], 2), 0.5);
                    res[2] = pow(pow(res[2], 2) + pow(add[2], 2), 0.5);
                }
            }
        }
        if (objects[q]->roughness != 0) {
            int prec = abs(objects[q]->roughness / 0.01);
            Ray first = save;
            bool inters_sun = Modify_Res_By_Atm(res, atmos, objects, lamps, sun, save);
            bool reflected_ok = false;
            while (!reflected_ok) {
                bool refl = first.Reflect(objects[q], prec);
                reflected_ok = refl;
                for (int i = 0; i < sun.size() && reflected_ok; i++) {
                    reflected_ok = lamps[sun[i]].pos->CheckInters(first.l).empty(); //можем не отразится нужным образом
                }
                if (!reflected_ok) { first = save; }
            }
            std::vector<long double> add = Trace(first, n, m, pov, lamps, objects, sun, atmos);
            //std::cerr << min_dist << ' ' << ray.l << ' ' << n << ' ' << m << ' ' << sun.size() << ' ' << objects[min_dist]->roughness << " hey\n";
            for (int i = 0; i < sun.size(); i++) {
                Vector dist = dynamic_cast<Sphere*>(lamps[sun[i]].pos)->pos - first.l.r_0;
                Vector l = dist / dist.Module();
                Line l_0 = save.l;
                long double alpha = asin(objects[q]->Normal(first.l.r_0).Vect(-l_0.tau).Module());
                long double beta = asin(objects[q]->roughness);
                long double gamma = asin(objects[q]->Normal(first.l.r_0).Vect(l).Module());
                long double check = asin((-l_0.tau).Vect(l).Module());
                long double checker = alpha + gamma;
                if (check - checker < pow(10, -4) && check - checker > -pow(10, -4)) {
                    alpha = -alpha;
                }
                Line fixing(first.l.r_0, l);
                if (gamma + alpha >= -2 * beta && gamma + alpha <= 2 * beta) {
                    bool not_inters_others = true;
                    int atmo_n_in = -1;
                    for (int ii = 0; ii < objects.size() && not_inters_others; ii++) {
                        not_inters_others = dynamic_cast<Atmosphere*>(objects[ii]) || objects[ii]->CheckInters(fixing).empty();
                        if (dynamic_cast<Atmosphere*>(objects[ii]) && !objects[ii]->CheckInters(fixing).empty()) { atmo_n_in = ii; }
                    }
                    for (int ii = 0; ii < lamps.size() && not_inters_others; ii++) {
                        not_inters_others = (ii == sun[i] || lamps[ii].pos->CheckInters(fixing).empty());
                    }
                    if (not_inters_others) {
                        for (int qq = 0; qq < 3; qq++) {
                            long double prob = dist.Module();
                            prob = pow((dynamic_cast<Sphere*>(lamps[sun[i]].pos)->radius / 2 / prob), 2);
                            if (atmo_n_in != -1) {
                                Atmosphere* fog = dynamic_cast<Atmosphere*>(objects[atmo_n_in]);
                                long double pp = fog->conc * pow(fog->r, 3) * 3.1415926 * fog->probability[qq];
                                long double dist = 0;
                                std::vector<Point> where_int = fog->CheckInters(fixing);
                                if (where_int.size() == 2) { dist = (where_int[1] - where_int[0]).Module(); }
                                else { dist = (fixing.r_0 - where_int[0]).Module(); }
                                prob *= pow(1 - pp, dist / fog->r);
                            }
                            add[qq] = pow(pow(add[qq], 2) + pow(lamps[sun[i]].intens[qq] * prob * dist.Scal(save.l.tau) / dist.Module() * objects[q]->reflection[qq], 2), 0.5);
                        }
                    }
                }
            }
            res[0] = pow(pow(res[0], 2) + pow(add[0], 2), 0.5);
            res[1] = pow(pow(res[1], 2) + pow(add[1], 2), 0.5);
            res[2] = pow(pow(res[2], 2) + pow(add[2], 2), 0.5);
        }
        if (objects[q]->refraction[0] + objects[q]->refraction[1] + objects[q]->refraction[2] != 0) {
            if (!dynamic_cast<Atmosphere*>(objects[q]) || !objects[q]->PointIn(save.l.r_0 - save.l.tau)) {
                Ray first = save;
                bool refract = first.Refract(objects[q]);
                if (refract) {
                    std::vector<long double> add = Trace(first, n, m, pov, lamps, objects, sun, atmos);
                    res[0] = pow(pow(res[0], 2) + pow(add[0], 2), 0.5);
                    res[1] = pow(pow(res[1], 2) + pow(add[1], 2), 0.5);
                    res[2] = pow(pow(res[2], 2) + pow(add[2], 2), 0.5);
                }
            }
        }
    }
    if (min_dist != -1 && min_dist >= objects.size()) {
        int q = min_dist - objects.size();
        res[0] = ray.red * lamps[q].red;
        res[1] = ray.green * lamps[q].green;
        res[2] = ray.blue * lamps[q].blue;
    }
    return res;
}
std::vector<std::vector<Pixel>> Calculate(long double angle, const Vector& spin, int n, int m, const std::vector<Lamp>& lamps, const std::vector<Surface*>& objects, const Plane* screen, const Point& pov, int zerox, int zeroy, int dx, int dy, const std::vector<int>& sun, const std::vector<int>& atmos) {
    srand(std::time(nullptr));
    std::vector<std::vector<Pixel>> res(dx, std::vector<Pixel>(dy));
    Vector add = screen->Normal(pov).Vect(spin);
    for (int i = zerox; i < zerox + dx; i++) {
        for (int j = zeroy; j < zeroy + dy; j++) {
            Vector temp = screen->Normal(pov) / angle + 1.0 / n * (i - n / 2) * spin + 1.0 / m * (j - m / 2) * add;
            Line line = Line(pov, temp / temp.Module());
            Ray ray = Ray(line, 1, 1, 1);
            std::vector<long double> ans = Trace(ray, n, m, pov, lamps, objects, sun, atmos);
            res[i - zerox][j - zeroy].red = ans[0];
            res[i - zerox][j - zeroy].green = ans[1];
            res[i - zerox][j - zeroy].blue = ans[2];
        }
    }
    return res;
}
int proc = 0;
std::vector<std::vector<Pixel>> Unification(long double angle, const Vector& spin, int n, int m, const std::vector<Lamp>& lamps, const std::vector<Surface*>& objects, const Plane* screen, const Point& pov, const std::vector<int>& sun = {}) {
    int rank, size;
    std::vector<std::vector<Pixel>> res(n, std::vector<Pixel>(m));
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int dx = 100;
    int dy = 100;
    std::vector<int> atmos;
    for (int i = 0; i < objects.size(); i++) {
        if (dynamic_cast<Atmosphere*>(objects[i])) {
            atmos.push_back(i);
        }
    }
    if (rank == 0) {
        int i = 0;
        int j = 0;
        for (int k = 1; k < size; k++) {
            std::vector<int> quest = {i, j};
            if (i != n) {
                MPI_Send(&(quest[0]), 2, MPI_INT, k, 0, MPI_COMM_WORLD);
                j += dy;
                if (j == m) {
                    i += dx;
                    j = 0;
                    if (i == n) { j = m; }
                }
            } else {
                MPI_Send(&(quest[0]), 2, MPI_INT, k, 1, MPI_COMM_WORLD);
            }
        }
        for (int k = 0; k < n * m / dx / dy; k++) {
            MPI_Status status;
            std::vector<long double> result(3 * dx * dy + 2);
            MPI_Recv(&(result[0]), 3 * dx * dy + 2, MPI_LONG_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int l = static_cast<int>(result[0]);
            int up = static_cast<int>(result[1]);
            for (int ii = l; ii < l + dx; ii++) {
                for (int jj = up; jj < up + dy; jj++) {
                    res[ii][jj] = Pixel(result[2 + ((ii - l) * dy + jj - up) * 3], result[2 + ((ii - l) * dy + jj - up) * 3 + 1], result[2 + ((ii - l) * dy + jj - up) * 3 + 2]);
                }
            }
            if (i != n || j != m) {
                std::vector<int> quest = {i, j};
                MPI_Send(&(quest[0]), 2, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                j += dy;
                if (j == m) {
                    i += dx;
                    j = 0;
                    if (i == n) { j = m; }
                }
            } else {
                std::vector<int> quest = {i, j};
                MPI_Send(&(quest[0]), 2, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
            }
            long double pro = (k + 1) * 1.0 / n / m * dx * dy;
            if (pro * 100 > proc) {
                std::cout << proc << "% finished\n";
                proc += 10;
            }
        }
        std::cout << 100 << "% finished\n";
    } else {
        int k = 0;
        while (k == 0) {
            MPI_Status status;
            std::vector<long double> result(3 * dx * dy + 2);
            std::vector<int> quest(2);
            MPI_Recv(&(quest[0]), 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            k = status.MPI_TAG;
            if (k == 0) {
                result[0] = 1.0 * quest[0];
                result[1] = 1.0 * quest[1];
                std::vector<std::vector<Pixel>> pixels = Calculate(angle, spin, n, m, lamps, objects, screen, pov, quest[0], quest[1], dx, dy, sun, atmos);
                for (int i = 0; i < dx; i++) {
                    for (int j = 0; j < dy; j++) {
                        result[2 + (i * dy + j) * 3] = pixels[i][j].red;
                        result[2 + (i * dy + j) * 3 + 1] = pixels[i][j].green;
                        result[2 + (i * dy + j) * 3 + 2] = pixels[i][j].blue;
                    }
                }
                MPI_Send(&(result[0]), 3 * dx * dy + 2, MPI_LONG_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
        }
    }
    return res;
}
#include <fstream>
std::pair<std::vector<Lamp>, std::vector<Surface*>> CreateScene() {
    std::pair<std::vector<Lamp>, std::vector<Surface*>> res;
    Point centre = Point({2, 0.3, -0.3});
    Cube* cube = new Cube(centre, {Vector({1, 0, 0}), Vector({0, 1, 0}), Vector({0, 0, 1})}, {0, 1, 1}, 0.4, 0, 1.3, {0.2, 0.2, 0.7});
    Sphere* sphere1 = new Sphere(Point({2, -0.4, 0}), 0.3, {0.8, 0.8, 0.8}, 0, 5, {0.5, 0.5, 0.5});
    Sphere* sphere2 = new Sphere(Point({3.3, 0.2, -0.1}), 0.2, {0.3, 1, 1}, 0, 0.4);
    Plane* plane = new Plane(Point({0, 0, 10}), Vector({0, 0, 1}), {1, 0, 0}, 0.9, 1.3);
    Plane* plane1 = new Plane(Point({10, 0, -2}), Vector({-1, 0, 1}), {1, 0, 0}, 0, 1.3);
    Plane* plane2 = new Plane(Point({0, 0, -0.3}), Vector({0, 0, 1}), {0, 1, 0}, 0.8, 1.3);
    res.first = {Lamp(plane, 1, 1, 1)};
    res.second = {cube, sphere1, sphere2, plane2};
    return res;
}
//cd Documents/lab/ray_tracing
//mpic++ -std=c++11 ray_trace_v1.6.cpp
//mpirun -n 4 a.out
int main() {
    std::ofstream fout("res.txt");
    auto t1 = std::clock();
    MPI_Init(NULL, NULL);
    std::pair<std::vector<Lamp>, std::vector<Surface*>> scene = CreateScene();
    int n = 1700, m = 1700;
    long double angle = 1; //angle - d/F
    Point pov = Point({-0.1, 0, 0});
    Vector screen_n = Vector({1, 0, 0});
    Plane screen = Plane(pov + 0.2 * screen_n, screen_n, {1, 1, 1}, 0.3, 10);
    Vector spin({0, 1, 0});
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) { std::cout << "Scene created successfully\n"; }
    if (rank == 0) {
        std::ofstream fout_c("conf.txt");
        fout_c << n << ' ' << m;
    }
    std::vector<std::vector<Pixel>> result = Unification(angle, spin, n, m, scene.first, scene.second, &screen, pov, {});
    MPI_Finalize();
    auto t2 = std::clock();
    std::cout << "Whole calculation time is " << (t2 - t1) * 1.0 / CLOCKS_PER_SEC << " s\n";
    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                fout << result[i][j].red << ' ' << result[i][j].green << ' ' << result[i][j].blue << ' ';
            }
            fout << std::endl;
        }
    }
    for (int i = 0; i < scene.first.size(); i++) {
        delete scene.first[i].pos;
    }
    for (int i = 0; i < scene.second.size(); i++) {
        delete scene.second[i];
    }
}
//красный рассвет (несколько отражений от атмосферы)
//рендерить маленький прямоугольник