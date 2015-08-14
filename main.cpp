//
//  main.cpp
//  Polygon-test
//
//  Created by Eric Lewis on 8/4/15.
//  Copyright (c) 2015 Eric Lewis. All rights reserved.
//

#include "simple_Polygon.h"
#include <iostream>

int main(int argc, const char * argv[]) {
    Polygon *a_polygon = new Polygon(3);
    a_polygon->V[0].x = 4;
    a_polygon->V[0].y = 5;
    a_polygon->V[1].x = 10;
    a_polygon->V[1].y = 12;
    bool is_simple_polygon = simple_Polygon( *a_polygon );

    std::cout << is_simple_polygon;
    return 0;
}
