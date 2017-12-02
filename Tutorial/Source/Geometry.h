#pragma once

#include <vector>
#include <algorithm>

class Point
{
public:
    Point()
    {
        x = 0;
        y = 0;
    }
    Point(float a, float b)
    {
        x = a;
        y = b;
    }

    bool isInf()
    {
        return ( (x == 0xdeadbeef && y == 0xdeadbeef) );
    }

    float m()
    {
        return x;
    }

    float b()
    {
        return y;
    }

    float x;
    float y;
};

class Line
{
public:
    Line(Point a1, Point b1)
    {
        a = a1;
        b = b1;
    }

    Point formula()
    {
        if (a.x == b.x)
        {
            return Point(0xdeadbeef, 0xdeadbeef);
        }

        float rise = (b.y - a.y);
        float run =  (b.x - a.x);

        float slope = rise / run;
        float intercept = a.y - (slope * a.x);

        return Point(slope, intercept);
    }

    Point a;
    Point b;
};

class Highway;

class Intercept : public Point
{
public:
    
    Intercept(float a, float b) : Point(a, b)
    {

    }

    Intercept(float a, float b, Highway* H1, Highway* H2) : Point(a, b)
    {
        h1 = H1;
        h2 = H2;
    }

    Highway* h1;
    Highway* h2;

    //Return the hihgway pointer that is not the one given
    Highway* other(Highway* h)
    {
        if (h1 == h)
            return h2;
        return h1;
    }
};

bool sortInterceptX(Intercept* i, Intercept* j);
bool sortInterceptY(Intercept* i, Intercept* j);

class Highway : public Line
{
public:
    Highway(Point a1, Point b1) : Line(a1, b1) {}
    Highway(Line l): Line(l.a, l.b) {}
    std::vector<Intercept*> intercepts;

    void sortPoints()
    {
        std::sort(intercepts.begin(), intercepts.end(), sortInterceptX);
    }
};



class Solid
{
public:
    Solid(float W, float H, float L)
    {
        w = W;
        h = H;
        l = L;
    }

    float w, h, l;
};

class Building
{
    std::vector<Solid> solids;

    Solid boundingBox;
};

extern std::vector<Highway> lines;
extern std::vector< std::vector<Intercept*> > chunks;
extern std::vector<Point> intersections;