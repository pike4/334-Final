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

    Point intercept(Line other)
    {
        Point f1 = formula();
        Point f2 = other.formula();

        if (!f1.isInf() && !f2.isInf())
        {
            //Can't calculate if f1.m is - cause we divide by it
            if (f1.m() != 0)
            {
                // Now solve the system of equations to get intersection point
                float factor = f2.m() / f1.m();
                float newB = f2.b() - (f1.b() * factor);
                float a = 1 - factor;
                float newY = newB / a;
                float newX = (newY - f1.b()) / f1.m();

                return Point(newX, newY);
            }

            //l1 is horizontal
            else
            {
                //l1 and l2 are parallel
                if (f2.m() == 0)
                {
                    return Point(0xdeadbeef, 0xdeadbeef);
                }

                //l2 and l1 intersect somewhere
                else
                {
                    float newY = a.y;
                    float newX = (newY - f2.b()) / f2.m();

                    return Point(newX, newY);
                }
            }
        }

        //Both linesToDo are vertical, no intersection
        else if (f1.isInf() && f2.isInf())
        {
            return Point(0xdeadbeef, 0xdeadbeef);
        }

        //One line is vertical and the other is not
        else
        {
            Point a, b;
            Line vertLine = Line(a, b);
            Point otherFormula;

            if (f1.isInf())
            {
                vertLine = *this;
                otherFormula = f2;
            }

            else
            {
                vertLine = other;
                otherFormula = f1;
            }

            float newX = vertLine.a.x;
            float newY = (otherFormula.m() * newX) + otherFormula.b();

            return Point(newX, newY);
        }
    }

    bool contains(Point other) {
        bool b1 = (std::min(a.x, b.x) < (other.x + 0.001) );
        bool b2 = (std::max(a.x, b.x) > (other.x - 0.001) );

        return (b1 && b2 && !other.isInf());
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