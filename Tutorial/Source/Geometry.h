#pragma once

#include <vector>
#include <algorithm>

extern double HEIGHT;
extern double WIDTH;
extern double UNIT;

class Point
{
public:
    Point() {
        x = 0;
        y = 0;
    }

    Point(double a, double b) {
        x = a;
        y = b;
    }

    bool isInf() {
        return ( (x == 0xdeadbeef && y == 0xdeadbeef) );
    }

    double m() {
        return x;
    }

    double b() {
        return y;
    }

    double x;
    double y;
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

        double rise = (b.y - a.y);
        double run =  (b.x - a.x);

        double slope = rise / run;
        double intercept = a.y - (slope * a.x);

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
                double factor = f2.m() / f1.m();
                double newB = f2.b() - (f1.b() * factor);
                double a = 1 - factor;
                double newY = newB / a;
                double newX = (newY - f1.b()) / f1.m();

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
                    double newY = a.y;
                    double newX = (newY - f2.b()) / f2.m();

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

            double newX = vertLine.a.x;
            double newY = (otherFormula.m() * newX) + otherFormula.b();

            return Point(newX, newY);
        }
    }

    // Get the x value of this line at the given 
    double xIntercept(double otherY) {
        Point f = formula();

        return ( otherY - f.b() ) / f.m();
    }

    // Get the y value of line defined by this segment at the given x coordinate
    double yIntercept(double otherX) {
        Point f = formula();

        return (f.m() * otherX) + f.b();
    }

    // True if the given point is within the domain and range of this segment
    bool contains(Point other) {
        bool b1 = (std::min(a.x, b.x) < (other.x + 0.1 * UNIT) );
        bool b2 = (std::max(a.x, b.x) > (other.x - 0.1 * UNIT) );

        return (b1 && b2 && !other.isInf());
    }

	// True if the lines actually cross, false if they intersect at the endpoints
	bool crosses(Point other) {
		bool b1 = (minX() < (other.x -  UNIT));
		bool b2 = (maxX() > (other.x +  UNIT));

		return (b1 && b2 && !other.isInf());
	}

    //True if the given x coordinate is within the domain of this segment
    bool containsX(double otherX) {
        bool b1 = (std::min(a.x, b.x) < (otherX + 0.1 * UNIT));
        bool b2 = (std::max(a.x, b.x) > (otherX - 0.1 * UNIT));

        return (b1 && b2);
    }

    //True if the given y coordinate is within the range of this segment
    bool containsY(double otherY) {
        bool b1 = (std::min(a.y, b.y) < (otherY + 0.1 * UNIT));
        bool b2 = (std::max(a.y, b.y) > (otherY - 0.1 * UNIT));

        return (b1 && b2);
    }

    //Return minimum x value of this line segment
    double minX() {
        return std::min(a.x, b.x);
    }

    //Return maximum x value of this line segment
    double maxX() {
        return std::max(a.x, b.x);
    }

    //Return minimum y value of this line segment
    double minY() {
        return std::min(a.y, b.y);
    }

    //Return maximum y value of this line segment
    double maxY() {
        return std::max(a.y, b.y);
    }

    double length() {
        return sqrt(pow(maxX() - minX(), 2) + pow(maxY() - minY(), 2));
    }

    Point a;
    Point b;
};

class Highway;

class Intercept : public Point
{
public:
    
    Intercept(double a, double b) : Point(a, b)
    {

    }

    Intercept(double a, double b, Highway* H1, Highway* H2) : Point(a, b)
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
    Solid(double W, double H, double L)
    {
        w = W;
        h = H;
        l = L;
    }

    double w, h, l;
};

class Building
{
    std::vector<Solid> solids;

    Solid boundingBox;
};

class mPolygon;

double randRange(double min, double max);

extern std::vector<Highway> lines;
extern std::vector< mPolygon > chunks;
extern std::vector<Point> intersections;