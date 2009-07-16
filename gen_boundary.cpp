#include <iostream>
#include <utility>
#include <vector>
#include <set>
#include <algorithm>
#include "geometry.h"
#include <pthread.h>
#include <unistd.h>
#include <deque>
#include "cmath"

using namespace std;

set<polygon_t> cache;

bool unique_distance(polygon_t &polygon, point_t p, set<long> &existing_distance)
{
    set<long> newly_introduced_distance;
    for (unsigned i = 0; i < polygon.size(); i ++) {
        long dist = distance_square(polygon[i], p);
        if (existing_distance.count(dist) > 0 || newly_introduced_distance.count(dist)) return false;
        newly_introduced_distance.insert(dist);
    }
    return true;
}

void add_distances_from_new_point(polygon_t &polygon, point_t p, set<long> &existing_distance)
{
    for (unsigned i = 0; i < polygon.size(); i ++) 
        existing_distance.insert(distance_square(polygon[i], p));
}

void remove_distances_from_new_point(polygon_t &polygon, point_t p, set<long> &existing_distance)
{
    for (unsigned i = 0; i < polygon.size(); i ++) 
        existing_distance.erase(distance_square(polygon[i], p));
}

polygon_t &touch_axis(polygon_t &polygon)
{
    if (polygon.size() == 0) return polygon;
    long min_x = polygon[0].first, min_y = polygon[0].second;
    for (polygon_t::iterator p = polygon.begin(); p < polygon.end(); p++) {
        if (p->first < min_x) min_x = p->first;
        if (p->second < min_y) min_y = p->second;
    }
    for (polygon_t::iterator p = polygon.begin(); p < polygon.end(); p++) {
        p->first -= min_x;
        p->second -= min_y;
    }
    return polygon;
}

polygon_t normalize(polygon_t &polygon)
{
    int sign[2] = {1, -1};
    vector<polygon_t> family;
    
    for (int xy = 0; xy < 4; xy++) {
        polygon_t mirror;
        polygon_t rotate;
        for (polygon_t::iterator p = polygon.begin(); p < polygon.end(); p++) {
            int x_sign = sign[xy & 0x01], y_sign = sign[xy >> 1];
            mirror.push_back(point_t(x_sign * p->first, y_sign * p->second));
            rotate.push_back(point_t(x_sign * p->second, y_sign * p->first));
        }
        touch_axis(mirror);
        touch_axis(rotate);
        sort(mirror.begin(), mirror.end());
        sort(rotate.begin(), rotate.end());
        family.push_back(mirror);
        family.push_back(rotate);
    }
    sort(family.begin(), family.end());
    return family[0];
}

void generate_boundary(polygon_t &boundary, unsigned length, int bound, long minimum_side_length, set<long> &existing_distance)
{
    long double circle_area = minimum_enclosing_circle(boundary);

    // if (circle_area >= max_circle_area) return;
    
    if (boundary.size() == length) { // finish
        polygon_t representation = normalize(boundary);
        if (cache.count(representation) > 0) return;
        cache.insert(representation);
        for (polygon_t::iterator p = representation.begin(); p < representation.end(); p++)
            cout << "(" << p->first << "," << p->second << ")";
        cout << " " << circle_area << endl;
        return;
    }

    vector<long> not_tested(bound + 1, 0);
    // fill(not_tested.begin(), not_tested.end(), 0);

    long x = 1, y = 0;
    not_tested[1] = 1;
    int full_pointer = 1;
    while (x < bound) {
        long dist = x * x + y * y;

        /* Find the appropriate vector */
        if (boundary.size() == 1) {
            boundary.push_back(point_t(x, y));
            existing_distance.insert(dist);
            generate_boundary(boundary, length, bound, dist, existing_distance);
            existing_distance.erase(dist);
            boundary.pop_back();
        } else if (dist > minimum_side_length) {
            long sign[2] = {-1, 1};
            deque<point_t> valid_expanded_point;
            
            point_t p_a = boundary[boundary.size() - 2];
            point_t p_b = boundary[boundary.size() - 1];
            
            for (int xy_sign = 0; xy_sign < 4; xy_sign++) {
                long dx = x * sign[xy_sign & 0x01], dy = y * sign[xy_sign >> 1];
                if (x == 0 && (xy_sign & 0x01)) continue;
                if (y == 0 && (xy_sign >> 1)) continue;
                point_t p_c = point_t(p_b.first + dx, p_b.second + dy);
                point_t p_d = point_t(p_b.first + dy, p_b.second + dx);
                
                if (cross_product(p_a, p_b, p_c) > 0 && cross_product(p_b, p_c, boundary[0]) > 0 && cross_product(p_c, boundary[0], boundary[1]) > 0 &&
                    unique_distance(boundary, p_c, existing_distance) && (boundary.size() + 1 != length || cross_product(p_b, p_c, boundary[0]) > 0))
                    valid_expanded_point.push_back(p_c);

                if (x == y) continue;
                if (cross_product(p_a, p_b, p_d) > 0 && cross_product(p_b, p_d, boundary[0]) > 0 && cross_product(p_d, boundary[0], boundary[1]) > 0 &&
                    unique_distance(boundary, p_d, existing_distance) && (boundary.size() + 1 != length || cross_product(p_b, p_c, boundary[0]) > 0))
                    valid_expanded_point.push_back(p_d);
            }
            
            while (valid_expanded_point.size() > 0) {
                long double min_angle = 2 * M_PI;
                for (unsigned i = 0; i < valid_expanded_point.size(); i++) {
                    long double roundness = fabsl(angle(p_b, p_a, valid_expanded_point[i]) - M_PI / 2);
                    if (roundness < min_angle) {
                        min_angle = roundness;
                        swap(valid_expanded_point[0], valid_expanded_point[i]);
                    }
                }
                
                point_t pt = valid_expanded_point.front();
                add_distances_from_new_point(boundary, pt, existing_distance);
                boundary.push_back(pt);
                generate_boundary(boundary, length, bound, minimum_side_length, existing_distance);
                boundary.pop_back();
                remove_distances_from_new_point(boundary, pt, existing_distance);
                
                valid_expanded_point.pop_front();
            }
        }

        
        /* Getting the next pair of distance */
        long min_x_pos = full_pointer, min_sq = min_x_pos * min_x_pos + not_tested[min_x_pos] * not_tested[min_x_pos];
        for (long x_pos = full_pointer + 1; x_pos < bound; x_pos++) {
            long sq = x_pos * x_pos + not_tested[x_pos] * not_tested[x_pos];
            if (sq < min_sq) 
                min_sq = sq, min_x_pos = x_pos;
            if (not_tested[x_pos] == 0) break;
        }
        x = min_x_pos, y = not_tested[min_x_pos];
        not_tested[min_x_pos]++;
        if ((min_x_pos == full_pointer) && (not_tested[full_pointer] > full_pointer))
            full_pointer++;
    }
}

void solve_for(int n, int max)
{
    polygon_t polygon;
    polygon.push_back(point_t(0, 0));
    set<long> distances;
    distances.insert(0);
    generate_boundary(polygon, n, max, 0, distances);
}

int main(int argc, char *argv[])
{
    int n, max;
    if (argc <= 2) {
        cout << "usage: gen_boundary n" << endl;
        return 0;
    }
    
    n = atoi(argv[1]);
    max = atoi(argv[2]);

    solve_for(n, max);
    return 0;
}