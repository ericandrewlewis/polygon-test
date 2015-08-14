// ===================================================================
// simple_polygon.h - Class for a polygon
// Written by Glenn Burkhardt (2014)
/*
 * simple_polygon.h
 *
 */

#ifndef SIMPLE_POLYGON_H_
#define SIMPLE_POLYGON_H_

#include <stdlib.h>

/**
 * A single x,y coordinate point.
 */
typedef struct {
    double x, y;
} Point;

/**
 * A Polygon data structure.
 */
struct Polygon {
    /**
     * @constructor
     * @param npts The number of distinct points on the polygon.
     */
    Polygon(int npts) {
        n = npts;
        V = (Point*)malloc(npts * sizeof(Point));
    }

    ~Polygon() {
        free(V);
    }

public:
    /**
     * The number of points of the polygon.
     */
    int n;
    /**
     * A pointer that should hold the points.
     *
     * The points should be a list of all distinct points.
     * The last point defined implicitly connects to the first point.
     */
    Point *V;
};

typedef struct Polygon Polygon;

/**
 * Test whether a polygon is simple or not (i.e. none of its lines intersect).
 *
 * @param Polygon
 * @return True if the polygon is simple, false if not.
 */
bool simple_Polygon( Polygon &Pn );

#endif /* SIMPLE_POLYGON_H_ */