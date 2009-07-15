#ifndef YLLAN_GEOMETRY_PACKAGE
#define YLLAN_GEOMETRY_PACKAGE 

#include <vector>
#include <set>
#include <algorithm>

typedef std::pair<long, long> point_t;
typedef std::pair<long double, long double> fpoint_t;
typedef std::vector<point_t> polygon_t;


#define distance_square(a, b)   (((a).first - (b).first) * ((a).first - (b).first) + ((a).second - (b).second) * ((a).second - (b).second))

long double minimum_enclosing_circle(polygon_t &convex_hull);
long cross_product(point_t a, point_t b, point_t c);
long double inner_product(point_t v1, point_t v2);
fpoint_t find_center(point_t a, point_t b, point_t c);
long double angle(point_t o, point_t a, point_t b);

bool inside_polygon(polygon_t &polygon, point_t point);
long double area(polygon_t &polygon);
#endif