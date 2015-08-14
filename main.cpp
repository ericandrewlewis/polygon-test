/**
 A program to test whether a Polygon is simple or not.
 */

#include "simple_Polygon.h"
#include "simple_Polygon.cpp"
#include <iostream>

int main(int argc, const char * argv[]) {
    // Create a pointer to a polygon.
    Polygon *a_polygon = new Polygon(4);

    // Define the points of the polygon.
    a_polygon->V[0].x = 0;
    a_polygon->V[0].y = 0;

    a_polygon->V[1].x = 0;
    a_polygon->V[1].y = 5;

    a_polygon->V[2].x = 3;
    a_polygon->V[2].y = 3;

    a_polygon->V[3].x = -3;
    a_polygon->V[3].y = 3;

    bool is_simple_polygon = simple_Polygon( *a_polygon );

    std::cout << is_simple_polygon;
    return 0;
}
