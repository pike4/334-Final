#pragma once
#include "Geometry.h"

class mPolygon {
public:
    std::vector<Intercept*> vertices;

    mPolygon (std::vector<Intercept*> v) {
        vertices = v;
    }

    mPolygon(const mPolygon& other)
    {
        vertices = other.vertices;
    }

    Intercept* operator[](int i) {
        return vertices[i];
    }
    
    operator std::vector<Intercept*>() {
        return vertices;
    }

    mPolygon topHat();
    mPolygon bottomHat();
    mPolygon leftHat();
    mPolygon rightHat();

    //Return this polygon's vertices a in clockwise or maybe counter-clockwise order
    mPolygon perimiterOrdered();

    //Return the area of the polygon
    double area();

    std::vector<mPolygon> addRoundabout();
	std::vector<Line>  addRoundabouts(int n);

    //Return a vector of recursively split polygons until the minimum size is reached for each
    std::vector<mPolygon> mPolygon::iceLatticeSplit();

    //Recursive step of iceSplit
    std::vector<mPolygon> mPolygon::iceRecurse(mPolygon cur);

    //Return the two polygons resulting from splitting this along a random line through the centroid
    std::vector<mPolygon> mPolygon::split();

	//Return a similar polygon such that the difference between each vertex and the centroid is linearly reduced by $offset
	mPolygon mPolygon::getBufferedBlock(double offset);

	//Return a similar polygon such that the difference between each vertex and the centroid is scaled by $ratio
	mPolygon mPolygon::shrinkBlock(double ratio);

	std::vector<Highway*> getHighways();

    Point centroid();
};