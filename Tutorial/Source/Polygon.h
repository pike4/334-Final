#pragma once
#include "Geometry.h"

class Polygon {

    std::vector<Intercept*> vertices;

    Polygon(std::vector<Intercept*> v) {
        vertices = v;
    }

    Intercept* operator[](int i) {
        return vertices[i];
    }
};