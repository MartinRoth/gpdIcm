#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H

typedef double coord_t;
typedef double coord2_t;

struct Point {
  coord_t x, y;
  bool operator <(const Point &p) const {
    return x < p.x || (x == p.x && y < p.y);
  }
};

std::vector<Point> convex_hull(std::vector<Point> P);

#endif
