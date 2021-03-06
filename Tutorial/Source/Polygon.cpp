#include "Polygon.h"
#include <ctime>
extern std::vector<Line> bullshitLines;
extern std::vector<Point> streetInts;

bool sortLineLength(Line i, Line j) {
	return i.length() < j.length();
}

std::vector<Line> getAllConnections(std::vector<Point> points) {
	std::vector<Line> ret;
	for (int i = 0; i < points.size(); i++) {
		for (int j = i + 1; j < points.size(); j++) {
			ret.push_back(Line(points[i], points[j]));
			//bullshitLines.push_back(Line(points[i], points[j]));
		}
	}

	return ret;
}

std::vector<Line> pairIntersections(std::vector<Point> a, std::vector<Point> b) {
	std::vector<Line> ret;
	for (int i = 0; i < a.size(); i++) {
		for (int j = 0; j < b.size(); j++) {
			ret.push_back(Line(a[i], b[j]));
		}
	}
	//bullshitLines.insert(bullshitLines.begin(), ret.begin(), ret.end());
	return ret;
}

//Return a list of lines from $connections such that the shortest lines are added back in order, where no longer lines cross any other lines
std::vector<Line> filterIntersections(std::vector<Line> connections, std::vector<Line> ret) {
	std::sort(connections.begin(), connections.end(), sortLineLength);
	for (int i = 0; i < connections.size(); i++) {
		bool good = true;

		//Check each line against all lines already verified, only add it if it doesn't connect with any other line
		for (int j = 0; j < ret.size(); j++) {
			Point intercept = connections[i].intercept(ret[j]);
			if (connections[i].crosses(intercept)) {
				good = false;
			}
		}
		if (good)
		{
			ret.push_back(connections[i]);
		}
	}

	return ret;
}

mPolygon mPolygon::perimiterOrdered() {

    //Top and bottom hat are already perimeter ordered, bottom points come after top points
    std::vector<Intercept*> top = topHat();
    std::vector<Intercept*> bottom = bottomHat();

    //The first and last points on the bottom are the same as the first and last points on top
    for (int i = bottom.size() - 1; i >= 0; i--) {
        if (bottom[i] != top[0] && bottom[i] != top[top.size() - 1]) {
            top.push_back(bottom[i]);
        }
        else {
            int a = 0;
        }
    }

    return top;
}

double mPolygon::area() {
    std::vector<Intercept*> per = perimiterOrdered();
    double ret = 0;
    for (int i = 0; i < per.size(); i++) {
        Intercept* j = per[i];
        Intercept* j2 = per[(i + 1) % per.size()];

        ret += (j->x * j2->y) - (j2->x * j->y);
    }
    if (ret < 0) {
        ret *= -1;
    }
    return ret / 2;
}

Point mPolygon::centroid() {
    std::vector<Intercept*> per = perimiterOrdered();

    double A = area();

    double retX = 0;
    double retY = 0;

    for (int i = 0; i < per.size(); i++) {
        Intercept* j = per[i];
        Intercept* j2 = per[(i + 1) % per.size()];

        retX += (j->x + j2->x) * (j->x * j2->y - j2->x * j->y);
        retY += (j->y + j2->y) * (j->x * j2->y - j2->x * j->y);
    }

    retX /= (-6 * A);
    retY /= (-6 * A);

    return Point(retX, retY);
}

std::vector<mPolygon> mPolygon::addRoundabout()
{
    std::vector<Line> ret;
    std::vector<Point> ints;
    Point center = centroid();
    //streetInts.push_back(center);
    std::vector<Intercept*> per = perimiterOrdered();

    std::vector<mPolygon> rett;

    for(int i = 0; i < per.size(); i++)
    { 
        Line cur = Line(*per[i], *per[(i + 1) % per.size()]);
        
        double min = cur.minX() + ((cur.maxX() - cur.minX()) / 3);
        double max = cur.minX() + ((cur.maxX() - cur.minX()) * 2 / 3);

        double newX = (min + max) / 2;// randRange(min, max);
        double newY = cur.yIntercept(newX);

        ret.push_back(Line(center, Point(newX, newY)));
        bullshitLines.push_back(Line(center, Point(newX, newY)));

        double len = cur.length();
        if (len > 20 * UNIT)
        {
            ret.push_back(Line(center, *per[i]));
			bullshitLines.push_back(Line(center, *per[i]));
            ints.push_back(*per[i]);
        }
        ints.push_back(Point(newX, newY));
    }

    for (int i = 0; i < ints.size(); i++) {
        std::vector<Intercept*> curr;
        curr.push_back(new Intercept(center.x, center.y));
        curr.push_back(new Intercept(ints[i].x, ints[i].y));
        curr.push_back(new Intercept(ints[(i + 1) % ints.size()].x, ints[(i + 1) % ints.size()].y));


        rett.push_back(curr);
    }
    

    return rett;
}

std::vector<mPolygon> mPolygon::addRoundabouts(double r) {
	Point center = centroid();
	int n = vertices.size(); 

   std::vector<Point> outer;
   std::vector<Point> inner;

   std::vector<Intercept*> per = perimiterOrdered();

   std::vector<mPolygon> ret;

	for (int i = 0; i < n; i++) {
       double ratio = r;

		double newX = per[i]->x - ((per[i]->x - center.x) * (1 - ratio));
		double newY = per[i]->y - ((per[i]->y - center.y) * (1 - ratio));

		//The offset for the current vertex
		Point roundabout = Point(newX, newY);

      bullshitLines.push_back(Line(roundabout, Point(per[i]->x, per[i]->y)));

		//Line connecting the current offset and ots vector
      inner.push_back(roundabout);
	}
	

	//std::vector<Line> connections = getAllConnections(firstScatter);
	
	// Get a list of all connections between lines ascending by length
	
   std::vector<Intercept*> innerPoly;
	for (int i = 0; i < inner.size(); i++) {
       innerPoly.push_back(new Intercept(inner[i].x, inner[i].y));
	}
   ret.push_back(innerPoly);

   for (int i = 0; i < inner.size(); i++) {
       std::vector<Intercept*> newTrapezoid;
       newTrapezoid.push_back(new Intercept(inner[i].x, inner[i].y));
       newTrapezoid.push_back(new Intercept(inner[(i + 1) % inner.size()].x, inner[(i + 1) % inner.size()].y));
       newTrapezoid.push_back(new Intercept(per[i]->x, per[i]->y));
       newTrapezoid.push_back(new Intercept(per[(i + 1) % per.size()]->x, per[(i + 1) % per.size()]->y));
       ret.push_back(newTrapezoid);
   }

	return ret;
}

#pragma region hats
mPolygon mPolygon::topHat()
{
    std::vector<Intercept*> topHat;
    std::vector<Intercept*> topHat2;

    std::sort(vertices.begin(), vertices.end(), sortInterceptX);

    topHat.push_back(vertices[0]);

    Line topLine = Line(Point(vertices[0]->x, vertices[0]->y), Point(vertices[vertices.size() - 1]->x, vertices[vertices.size() - 1]->y));
    Point topF = topLine.formula();
    for (int i = 1; i < vertices.size() - 1; i++)
    {
        if (vertices[i]->y >((vertices[i]->x * topF.m()) + topF.b()))
        {
            topHat.push_back(vertices[i]);
        }
    }
    topHat.push_back(vertices[vertices.size() - 1]);

    return topHat;
}

mPolygon mPolygon::bottomHat()
{
    std::vector<Intercept*> topHat;

    std::sort(vertices.begin(), vertices.end(), sortInterceptX);

    topHat.push_back(vertices[0]);
    Line topLine = Line(Point(vertices[0]->x, vertices[0]->y), Point(vertices[vertices.size() - 1]->x, vertices[vertices.size() - 1]->y));
    Point topF = topLine.formula();
    for (int i = 1; i < vertices.size() - 1; i++)
    {
        if (vertices[i]->y  < ((vertices[i]->x * topF.m()) + topF.b()))
        {
            topHat.push_back(vertices[i]);
        }
    }
    topHat.push_back(vertices[vertices.size() - 1]);
    return topHat;
}

mPolygon mPolygon::leftHat()
{
    std::vector<Intercept*> topHat;

    std::sort(vertices.begin(), vertices.end(), sortInterceptY);

    topHat.push_back(vertices[0]);
    Line topLine = Line(Point(vertices[0]->x, vertices[0]->y), Point(vertices[vertices.size() - 1]->x, vertices[vertices.size() - 1]->y));
    Point topF = topLine.formula();
    for (int i = 1; i < vertices.size() - 1; i++)
    {
        if (vertices[i]->x < ((vertices[i]->y - topF.b()) / topF.m()))
        {
            topHat.push_back(vertices[i]);
        }
    }
    topHat.push_back(vertices[vertices.size() - 1]);
    return topHat;
}

mPolygon mPolygon::rightHat()
{
    std::vector<Intercept*> topHat;

    std::sort(vertices.begin(), vertices.end(), sortInterceptY);

    topHat.push_back(vertices[0]);
    Line topLine = Line(Point(vertices[0]->x, vertices[0]->y), Point(vertices[vertices.size() - 1]->x, vertices[vertices.size() - 1]->y));
    Point topF = topLine.formula();
    for (int i = 1; i < vertices.size() - 1; i++)
    {
        if (vertices[i]->x >((vertices[i]->y - topF.b()) / topF.m()))
        {
            topHat.push_back(vertices[i]);
        }
    }
    topHat.push_back(vertices[vertices.size() - 1]);
    return topHat;
}
#pragma endregion

//Return a vector of recursively split polygons until the minimum size is reached for each
std::vector<mPolygon> mPolygon::iceLatticeSplit() {
    srand(time(NULL));
    return iceRecurse(*this);
}

//Recursive step of iceSplit
std::vector<mPolygon> mPolygon::iceRecurse(mPolygon cur) {
    std::vector<mPolygon> ret;
    
    // Split current in half
    std::vector<mPolygon> res = cur.split();

    for (int i = 0; i < res.size(); i++) {
        if (res[i].area() > 0.03) {
            //Push back of recursing on either half
            std::vector<mPolygon> splitted = iceRecurse(res[i]);
            ret.insert(ret.end(), splitted.begin(), splitted.end());
        }
		else {
			ret.insert(ret.end(), res.begin(), res.end());
		}
    }

    return ret;
}

//split in half along random line through centroid
std::vector<mPolygon> mPolygon::split()
{
    std::vector<mPolygon> ret;

    Point center = centroid();
    //streetInts.push_back(center);
    std::vector<Intercept*> per = perimiterOrdered();

    //int ind;
    //
    //for (int i = 0; i < 1000; i++) {
    //    ind = rand() % (per.size());
    //}
    Line intersectWithThis = Line(Point(0,0),Point(0,0));
    int indd = 0;
    for (int ind = 0; ind < per.size(); ind++) {
        //Get a random line on the hull to intersect with
        Point a = Point(per[ind]->x, per[ind]->y);
        Point b = Point(per[(ind + 1) % per.size()]->x, per[(ind + 1) % per.size()]->y);
        Line temp = Line(a, b);

        if (temp.length() > intersectWithThis.length()) {
            intersectWithThis = temp;
            indd = ind;
        }
    }

    //Pick a random point in the middle 3/5s of the line
    double min = intersectWithThis.minX() + ((intersectWithThis.maxX() - intersectWithThis.minX()) / 2);
    double max = intersectWithThis.minX() + ((intersectWithThis.maxX() - intersectWithThis.minX()) * 2 / 3);

    double midX = randRange(min, max);
    double midY = intersectWithThis.yIntercept(midX);

    Point firstPoint = Point(midX, midY);
    Point otherPoint;
   
    //Line between the centroid and midpoint of a random side
    Line splitLine = Line(center, firstPoint);
   
    bool found = false;
   
   
    //TODO: now we need to get an intercept with the other side
    for (int i = 0; i < per.size(); i++) {
   
        //the line being checked is the same one we're shooting off from
        if (i == indd) {
            continue;
        }
   
        //Manufacture a line from the vertices at this index
        Point p1 = Point(per[i]->x, per[i]->y);
        Point p2 = Point(per[(i + 1) % per.size()]->x, per[(i + 1) % per.size()]->y);
        Line side = Line(p1, p2);
   
        Point curIntercept = side.intercept(splitLine);
   
        if (side.contains(curIntercept)) {
            otherPoint = curIntercept;
            found = true;
            break;
        }
    }
   
    if (found) {
        Line realSplitLine = Line(firstPoint, otherPoint);
        bullshitLines.push_back(realSplitLine);
        std::vector<Intercept*> r1;
        std::vector<Intercept*> r2;
        
        //Both intercept points go on both new chunks
        r1.push_back(new Intercept(firstPoint.x, firstPoint.y));
        r1.push_back(new Intercept(otherPoint.x, otherPoint.y));
        r2.push_back(new Intercept(firstPoint.x, firstPoint.y));
        r2.push_back(new Intercept(otherPoint.x, otherPoint.y));

        for (int i = 0; i < per.size(); i++) {
            //We don't know much, but we do know that every point is either above, or below realSplitLine
            if (per[i]->y > realSplitLine.yIntercept(per[i]->x)) {
                r1.push_back(per[i]);
            }
            else {
                r2.push_back(per[i]);
            }
        }
        ret.push_back(r1);
        ret.push_back(r2);
    }
    else {
        ret.push_back(*this);
    }
    return ret;
}

mPolygon mPolygon::shrinkBlock(double ratio) {

	std::vector<Intercept*> per = perimiterOrdered();
	std::vector<Intercept*> ret;

	Point center = centroid();

	for (int i = 0; i < per.size(); i++) {
		Point cur = Point(per[i]->x - center.x, per[i]->y - center.y);

		cur.x *= (1 - ratio);
		cur.y *= (1 - ratio);

		ret.push_back(new Intercept(per[i]->x - cur.x, per[i]->y - cur.y));
		//streetInts.push_back(*ret[i]);
	}

	for (int i = 0; i < ret.size(); i++) {
		Point a = *ret[i];
		Point b = *ret[(i + 1) % ret.size()];

		bullshitLines.push_back(Line(a, b));
	}

	return ret;
}

mPolygon mPolygon::getBufferedBlock(double offset) {
	std::vector<Intercept*> per = perimiterOrdered();
	std::vector<Intercept*> ret;

	Point center = centroid();

	for (int i = 0; i < per.size(); i++) {
		double newX = 0;
		double newY = 0;
		if (per[i]->x > center.x) {
			newX = per[i]->x - offset;
			if (newX < center.x) {
				ret.clear();
				break;
			}
		}
		else if (per[i]->x < center.x) {
			newX = per[i]->x + offset;
			if (newX > center.x) {
				ret.clear();
				break;
			}
		}
		else {
			int ff = 9;
		}
		if (per[i]->y > center.y) {
			newY = per[i]->y - offset;
			if (newY < center.y) {
				ret.clear();
				break;
			}
		}
		else if (per[i]->y < center.y) {
			newY = per[i]->y + offset;
			if (newY > center.y) {
				ret.clear();
				break;
			}
		}
		else {
			int yyy = 0;
		}

		if (i == per.size() - 1) {
			int why = 0;
		}

		ret.push_back(new Intercept(newX, newY));
	}

	for (int i = 0; i < ret.size(); i++) {
		Point a = *ret[i];
		Point b = *ret[(i + 1) % ret.size()];

		bullshitLines.push_back(Line(a, b));
	}

	return ret;
}

std::vector<Highway*> mPolygon::getHighways()
{
	std::vector<Intercept*> per = perimiterOrdered();

	std::vector<Highway*> ret;

	for (int i = 0; i < per.size(); i++) {
		Point a = Point(per[i]->x, per[i]->y);
		Point b = Point(per[(i + 1) % per.size()]->x, per[(i + 1) % per.size()]->y);

		Highway* cur = new Highway(a, b);
	}

	return ret;
}