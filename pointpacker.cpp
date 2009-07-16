#include <iostream>
#include <utility>
#include <vector>
#include <set>
#include <algorithm>
#include "geometry.h"
#include <pthread.h>
#include <unistd.h>

#define BOUND 40

using namespace std;

long max_width = 5e18;
long max_height = 5e18;
long max_diagonal_square = 5e18;
long double max_circle_area = 5e18L;

pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;

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

polygon_t normalize(polygon_t &old_polygon, long min_x, long min_y)
{
    polygon_t new_polygon;
    for (polygon_t::iterator p = old_polygon.begin(); p < old_polygon.end(); p++)
        new_polygon.push_back(point_t(p->first - min_x, p->second - min_y));
    return new_polygon;
}

void print_solution(polygon_t &polygon)
{
    for (unsigned i = 0; i < polygon.size(); i++) {
        if (i > 0) cout << ", ";
        cout << "(" << polygon[i].first << ", " << polygon[i].second << ")";
    }
    cout << endl;
}

bool test_fill(const vector<point_t> &candidates, int index, int left, polygon_t &polygon, set<long> &existing_distance)
{
    if (left == 0) return true;
    
    for (unsigned i = index; i < candidates.size(); i++) {
        point_t p = candidates[i];
        if (!unique_distance(polygon, p, existing_distance)) continue;

        add_distances_from_new_point(polygon, p, existing_distance);
        polygon.push_back(p);
        
        if (test_fill(candidates, i + 1, left - 1, polygon, existing_distance)) 
            return true;
        
        polygon.pop_back();
        remove_distances_from_new_point(polygon, p, existing_distance);
    }
    
    return false;
}

polygon_t valid_interior(polygon_t &boundary, unsigned interiors, set<long> &existing_distance, long width, long height)
{
    vector<point_t> candidates;
    for (int w = 0; w < width; w++) {
        for (int h = 0; h < height; h++) {
            point_t p = point_t(w, h);
            if (!inside_polygon(boundary, p)) continue;
            if (!unique_distance(boundary, p, existing_distance)) continue;
            candidates.push_back(p);
        }
    }
    if (candidates.size() < interiors) return vector<point_t>();
    
    polygon_t solution = boundary;
    set<long> distances = existing_distance;
    for (unsigned i = 0; i <= candidates.size() - interiors; i++) {
        add_distances_from_new_point(boundary, candidates[i], distances);
        solution.push_back(candidates[i]);
        
        if (test_fill(candidates, i + 1, interiors - 1, solution, distances)) {
            return solution;
        }
        
        solution.pop_back();
        remove_distances_from_new_point(boundary, candidates[i], distances);
    }
    return vector<point_t>();
}

void generate_boundary(polygon_t &boundary, unsigned length, int interiors, long minimum_side_length, set<long> &existing_distance, long left_bound, long right_bound, long bottom_bound, long top_bound)
{
    long width = right_bound - left_bound, height = top_bound - bottom_bound;
    long double circle_area = minimum_enclosing_circle(boundary);

    if (max(height, width) * max(height, width) > max_diagonal_square) return;
    if (circle_area >= max_circle_area) return;
    
    if (boundary.size() == length) { // finish
        if (distance_square(boundary[0], boundary.back()) < distance_square(boundary[0], boundary[1])) return;

        polygon_t normalized = normalize(boundary, left_bound, bottom_bound);

        if (interiors > 0) {
            normalized = valid_interior(normalized, interiors, existing_distance, width, height);
            if (normalized.size() == 0) return;
        }
        
        pthread_mutex_lock(&mtx);
        if (circle_area < max_circle_area) {
            max_circle_area = circle_area;
            cout << length << "/" << normalized.size() << ": ";
            print_solution(normalized);
            cout << "circle area=" << circle_area << endl;
                
            if (width * width + height * height < max(max_width, max_height)) {
                max_height = height;
                max_width = width;
                max_diagonal_square = width * width + height * height;
            }
        }
        pthread_mutex_unlock(&mtx);
        
        return;
    }

    vector<long> not_tested;
    for (int i = 0; i < BOUND; i++)
        not_tested.push_back(0);

    long x = 1, y = 0;
    not_tested[1] = 1;
    int full_pointer = 1;
    while (x < BOUND) {
        // cout << x << ", " << y << " = " << x * x + y * y << endl;

        long dist = x * x + y * y;

        /* Find the appropriate vector */
        if (boundary.size() == 1) {
            boundary.push_back(point_t(x, y));
            existing_distance.insert(dist);
            generate_boundary(boundary, length, interiors, dist, existing_distance, 0, x, 0, y);
            existing_distance.erase(dist);
            boundary.pop_back();
        } else if (dist > minimum_side_length) {
            long sign[2] = {-1, 1};
            vector<point_t> valid_expanded_point;
            
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
            
            for (vector<point_t>::iterator pt = valid_expanded_point.begin(); pt < valid_expanded_point.end(); pt++) {
                long extend_left = min(left_bound, pt->first), extend_right = max(right_bound, pt->first);
                long extend_bottom = min(bottom_bound, pt->second), extend_top = max(top_bound, pt->second);
                
                add_distances_from_new_point(boundary, *pt, existing_distance);
                boundary.push_back(*pt);
                generate_boundary(boundary, length, interiors, minimum_side_length, existing_distance, extend_left, extend_right, extend_bottom, extend_top);
                boundary.pop_back();
                remove_distances_from_new_point(boundary, *pt, existing_distance);
            }
            
        }

        
        /* Getting the next pair of distance */
        long min_x_pos = full_pointer, min_sq = min_x_pos * min_x_pos + not_tested[min_x_pos] * not_tested[min_x_pos];
        for (long x_pos = full_pointer + 1; x_pos < BOUND; x_pos++) {
            long sq = x_pos * x_pos + not_tested[x_pos] * not_tested[x_pos];
            if (sq < min_sq) 
                min_sq = sq, min_x_pos = x_pos;
            if (not_tested[x_pos] == 0) break;
        }
        x = min_x_pos, y = not_tested[min_x_pos];
        not_tested[min_x_pos]++;
        if ((min_x_pos == full_pointer) && (not_tested[full_pointer] > full_pointer))
            full_pointer++;
        // cout << "next x:" << x << ", next y:" << y << endl;
    }
}

void *generate_wrapper(void *ptr)
{
    pair<int, int> args = *((pair<int, int> *)ptr);
    int n = args.first;
    int b = args.second;
    polygon_t polygon;
    polygon.push_back(point_t(0, 0));

    set<long> distances;
    distances.insert(0);
    
    generate_boundary(polygon, b, n - b, 0, distances, 0, 0, 0, 0);
    return NULL;
}

void solve(int n)
{
    vector<pair<int, int> > args(n - 2);
    pthread_t tid[n - 2];

    for (int b = 3; b <= n; b++) {
        int index = b - 3;
        args[index].first = n;
        args[index].second = b;
        pthread_create(tid + index, NULL, generate_wrapper, (void *)(&args[index]));
    }
    
    for (int t = 0; t < n - 2; t++)
        pthread_join(tid[t], NULL);
}

int main(int argc, char *argv[])
{
    int n;
    if (argc <= 1) {
        cout << "usage: bruteforce n [bound]" << endl;
        return 0;
    }
    n = atoi(argv[1]);
    if (argc > 2) {
        max_circle_area = atof(argv[2]);
        max_diagonal_square = max_circle_area * 8;
    }
    
    solve(n);
    // cout << minimum_enclosing_circle(polygon) << endl;
    return 0;
}
