#include "geometry.h"
#include <cmath>
#include <cstdio>

bool inside_polygon(polygon_t &polygon, point_t point)
{
    if (polygon.size() == 0) return false;
    for (polygon_t::iterator p = polygon.begin(); p < polygon.end(); p++)
        if (*p == point) return true;

    if (polygon.size() == 1) return false;

    long double total_area = area(polygon);
    long double sub_area = 0;
    for (unsigned i = 0; i < polygon.size(); i++) {
        point_t a = point;
        point_t b = polygon[i];
        point_t c = polygon[(i + 1) % polygon.size()];

        sub_area += fabsl((a.first * b.second + b.first * c.second + c.first * a.second) - (a.second * b.first + b.second * c.first + c.second * a.first)) / 2.0;
    }
    return (fabsl(total_area - sub_area) < 1e-10L);
}

long double area(polygon_t &polygon)
{
    if (polygon.size() < 3) return 0;
    
    long positive_term = 0, negative_term = 0;
    for (unsigned i = 0; i < polygon.size() - 1; i++) {
        positive_term += polygon[i].first * polygon[i + 1].second;
        negative_term += polygon[i].second * polygon[i + 1].first;
    }
    positive_term += polygon.back().first * polygon[0].second;
    negative_term += polygon.back().second * polygon[0].first;

    return fabsl((positive_term - negative_term) / 2.0);
}

long double minimum_enclosing_circle(polygon_t &convex_hull)
{
    if (convex_hull.size() <= 1) return 0;
    if (convex_hull.size() == 2) return distance_square(convex_hull[0], convex_hull[1]) / 4.0;
    
    point_t s_a = convex_hull[0];
    point_t s_b = convex_hull[1];
    
    while (1) {
        long double alpha = 100;
        point_t v;
        for (polygon_t::iterator p = convex_hull.begin(); p != convex_hull.end(); p++) {
            if (*p == s_a || *p == s_b) continue;
            long double a = angle(*p, s_a, s_b);
            if (a < alpha) {
                alpha = a;
                v = *p;
            }
        }
        
        if (alpha >= M_PI / 2) return distance_square(s_a, s_b) / 4.0;

        // printf("s_a:(%ld, %ld) s_b:(%ld, %ld) v:(%ld, %ld)\n", s_a.first, s_a.second, s_b.first, s_b.second, v.first, v.second);
        // printf("angle v-s_a-s_b: %Lf\n", angle(s_a, v, s_b));

        if (angle(s_a, v, s_b) >= M_PI / 2) {
            s_a = v;
            continue;
        }

        // printf("angle v-s_b-s_a: %Lf\n", angle(s_b, v, s_a));        
        if (angle(s_b, v, s_a) >= M_PI / 2) {
            s_b = v;
            continue;
        }
        
        /* v, s_a, s_b */
        fpoint_t center = find_center(v, s_a, s_b);
        // printf("center: %Lf, %Lf\n", center.first, center.second);
        // printf("%Lf = %Lf, %Lf\n", distance_square(center, v), distance_square(center, s_a), distance_square(center, s_b));
        return distance_square(center, v);
    }
}

fpoint_t find_center(point_t a, point_t b, point_t c)
{
    fpoint_t center_ab((a.first + b.first) / 2.0, (a.second + b.second) / 2.0);
    fpoint_t center_bc((c.first + b.first) / 2.0, (c.second + b.second) / 2.0);
    
    point_t normal_ab(a.second - b.second, b.first - a.first);
    point_t normal_bc(c.second - b.second, b.first - c.first);
    
    // printf("c(ab): (%Lf, %Lf), c(bc): (%Lf, %Lf)\n", center_ab.first, center_ab.second, center_bc.first, center_bc.second);
    // printf("n(ab): (%ld, %ld), n(bc): (%ld, %ld)\n", normal_ab.first, normal_ab.second, normal_bc.first, normal_bc.second);
    
    // center_ab + r * normal_ab == center_bc + s * normal_bc
    long double r = (long double)(normal_bc.second * (center_bc.first - center_ab.first) - normal_bc.first * (center_bc.second - center_ab.second)) / (normal_bc.second * normal_ab.first - normal_bc.first * normal_ab.second);
    // printf("r: %Lf\n", r);
    return fpoint_t(center_ab.first + r * normal_ab.first, center_ab.second + r * normal_ab.second);
}

/* return angle a-o-b in radian. */
long double angle(point_t o, point_t a, point_t b) 
{
    point_t v1(a.first - o.first, a.second - o.second);
    point_t v2(b.first - o.first, b.second - o.second);
    
    return acosl(inner_product(v1, v2) / sqrtl(inner_product(v1, v1) * inner_product(v2, v2)));
}

long double inner_product(point_t v1, point_t v2)
{
    return v1.first * v2.first + v1.second * v2.second;
}

long cross_product(point_t a, point_t b, point_t c)
{
    long x1 = b.first - a.first, x2 = c.first - b.first;
    long y1 = b.second - a.second, y2 = c.second - b.second;
    return x1 * y2 - x2 * y1;
}