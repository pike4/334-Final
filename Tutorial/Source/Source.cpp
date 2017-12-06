#include "glew.h"
#include "glut.h"
#include "Geometry.h"
#include "Debug_Display.h"
#include "Polygon.h"
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <cmath>

#include <fstream>
#include <iostream>

std::vector<Highway> lines;
std::vector<Highway> streets;
std::vector<std::vector<Highway>> streetSets;

std::vector< mPolygon > chunks;

std::vector<Point> streetInts;
std::vector<Point> intersections;

char initialCondition;
int numSides;

std::vector<Line> bullshitLines;

static const GLfloat vertices[] =
{
    -1.0f, -1.0f, 0.0f,
    1.0f, -1.0f, 0.0f,
    0.0f, 1.0f, 0.0f
};

int NUM_VERTICAL_HIGHWAYS = 2;
int NUM_HORIZONTAL_HIGHWAYS = 2;
double MIN_BLOCK_AREA = 0.01;
const int NUM_MID_STREETS = 2;

double UNIT = 1.0f;
double WIDTH = UNIT * 1000.0f;
double HEIGHT = UNIT * 1000.0f;

double STREET_WIDTH;
double GRID_SIZE;

// Likelihood that a given highway will be parallel to an axis rather than at an angle
const double HIGHWAY_TAXI_FACTOR = 0.0;

struct Rule {
	double minSize;
	int depth;
	char rule;
};

std::vector<Rule> readRules() {
	std::ifstream inFile;
	std::vector<Rule> ret;

	inFile.open("grammar.txt");

	std::string line1, line2;

	double minSize;
	int h, v;
	char cur;
	while (1) {
		inFile >> cur;
		if (!inFile.good()) {
			break;
		}
		Rule newRule;
		
		if (cur == 'H') {
			inFile >> NUM_HORIZONTAL_HIGHWAYS;
		}
		else if (cur == 'V') {
			inFile >> NUM_VERTICAL_HIGHWAYS;
		}
		else if (cur == 'W') {
			inFile >> STREET_WIDTH;
			STREET_WIDTH *= UNIT;
		}
		else if( cur == 'R' || cur == 'I' || cur == 'M' || cur == 'B') {
			newRule.rule = cur;
			ret.push_back(newRule);
		}
		else if (cur == 'G') {
			newRule.rule = cur;
			inFile >> newRule.minSize;
			ret.push_back(newRule);
		}
		else if (cur == 'C' ) {
			initialCondition = cur;
			inFile >> numSides;
		}
		else if (cur == 'S') {
			initialCondition = cur;
		}
	}
	return ret;
}

#pragma region Highway generation
std::vector<Highway> genTicTacToe()
{
    std::vector<Highway> ret;
    ret.push_back(Highway(Point(-WIDTH, 0.2), Point(WIDTH, 0.2 * HEIGHT) ) );
    //ret.push_back(Highway(Point(-1.0f, 0.5), Point(1.0f, 0.2)));
    ret.push_back(Highway(Point(-1.0f, -0.5), Point(1.0f, -0.3)));
    

    ret.push_back(Highway(Point(0.2f, -1.0f), Point(0.5f, 1.0f)));
    ret.push_back(Highway(Point(-0.2, -1.0f), Point(-0.5, 1.0f)));
    return ret;
}

std::vector<Highway> genBoundary()
{
    std::vector<Highway> ret;

    ret.push_back(Highway(Point(-WIDTH, -999.0f), Point(WIDTH, -998.0f)));
    ret.push_back(Highway(Point(-WIDTH, 999.0f), Point(WIDTH, 999.0f)));
    //
    ret.push_back(Highway(Point(-0.98f * WIDTH, -HEIGHT), Point(-0.998f * WIDTH, HEIGHT)));
    ret.push_back(Highway(Point(0.98f * WIDTH, -HEIGHT), Point(0.99f * WIDTH, HEIGHT)));

    return ret;
}

std::vector<Highway> genLines()
{
    srand(time(NULL));
    std::vector<Highway> ret;
    for (int i = 0; i < NUM_VERTICAL_HIGHWAYS; i++)
    {
        Point p1;
        Point p2;

		//Limit placements of highways to the ith TOTAL_LENGTH / NUM_VERTICAL_HIGHWAYS length stretch of the road
		double min = -WIDTH + ( (WIDTH - (-WIDTH)) * i) / NUM_VERTICAL_HIGHWAYS;
		double max = -WIDTH + ( (WIDTH - (-WIDTH)) * (i + 1)) / NUM_VERTICAL_HIGHWAYS;

		p1.x = randRange(min, max);
        p1.y = -HEIGHT;

        if (((rand() % 1000) / 1000.0f) < HIGHWAY_TAXI_FACTOR)
        {
            p2.x = p1.x;
        }
        else
        {
			p2.x = randRange(min, max);
        }
        
        p2.y = HEIGHT;

        Line L = Line(p1, p2);

        ret.push_back(L);
    }

    for (int i = 0; i < NUM_HORIZONTAL_HIGHWAYS; i++)
    {
        Point p1;
        Point p2;

		//Limit placements of highways to the ith TOTAL_LENGTH / NUM_VERTICAL_HIGHWAYS length stretch of the road
		double min = ( (HEIGHT - (-HEIGHT)) * i) / NUM_HORIZONTAL_HIGHWAYS;
		double max = ( (HEIGHT - (-HEIGHT)) * (i + 1)) / NUM_HORIZONTAL_HIGHWAYS;

		p1.y = randRange(min, max);
        p1.x = -WIDTH;

        if (((rand() % 1000) / 1000.0f) < HIGHWAY_TAXI_FACTOR)
        {
            p2.y = p1.y;
        }
        else
        {
            p2.y = randRange(min, max);
        }

        p2.x = WIDTH;

        Line L = Line(p1, p2);

        ret.push_back(L);
    }

    return ret;
}

mPolygon genPerimeter(int count) {
	
	Point center = Point(0,0);
	std::vector<Point> hull;
	std::vector<Intercept*> ret;
	std::vector<double> rads;

	double maxRad = WIDTH * 0.99;
	double minRad = HEIGHT * 0.7;

	for (int i = 0; i < count; i++) {
		double curRad = randRange(minRad, maxRad);

		rads.push_back(curRad);
	}

	std::sort(rads.begin(), rads.end());
	rads[0]+= (std::max(rads[1], rads[rads.size() - 1]) - std::min(rads[1], rads[rads.size() - 1])) * 1 / 3;
	rads[rads.size() - 1] -= (std::max(rads[1], rads[rads.size() - 1]) - std::min(rads[1], rads[rads.size() - 1])) * 2 / 10;

	
	for (int i = 0; i < rads.size(); i++) {
		double curRad = rads[i];
		double curAngle = ((2 * 3.14159265) * i) / count;
		hull.push_back(Point(center.x + (curRad * cos(curAngle)), center.y + (curRad * sin(curAngle))));
	}


	for (int i = 0; i < hull.size(); i++) {
		Point a = hull[i];
		Point b = hull[(i + 1) % hull.size()];
		bullshitLines.push_back(Line(a, b));
		ret.push_back(new Intercept(a.x, a.y));
	}

	return ret;
}

std::vector<Highway> genOffshoots(std::vector<Highway> orig)
{
    for (int i = 0; i < NUM_MID_STREETS; i++)
    {
        int ind = rand() % orig.size();

        Point form = orig[ind].formula();

        double minX = std::min(orig[ind].a.x, orig[ind].b.x);
        double maxX = std::max(orig[ind].a.x, orig[ind].b.x);

        double offX = randRange(minX, maxX);
        double offY = offX * form.m() + form.b();

        double endX = 0;
        double endY = 0;

        double newX, newY;

        double slope = -1.0/form.m();
        double newB = offY - (slope * offX);

        newX = 1.0;
        newY = endX * slope + newB;

        Point A = Point(offX, offY);
        Point B = Point(newX, newY);
        Highway newH = Highway(A, B);

        orig.push_back(newH);
    }

    return orig;
}

#pragma endregion

// Store a vector of intersections with other highways in each of the highways in lines,
// return a list of intercept points
std::vector<Point> getIntersections(std::vector<Highway>* linesToDo)
{
    std::vector <Point> ret;
    for (int i = 0; i < linesToDo->size(); i++)
    {
        for (int j = i; j < linesToDo->size(); j++)
        {
            if (i == j)
                continue;

            Line l1 = (*linesToDo)[i];
            Line l2 = (*linesToDo)[j];

            //Get the point at which the lines defined by these line segments intersect
            Point f1 = l1.intercept(l2);
            
            // Add a new intersection if the returned point exists and lies on both line segments
            if(l1.contains(f1) && l2.contains(f1)) {
                Intercept* newIntercept = new Intercept(f1.x, f1.y, &(*linesToDo)[i], &(*linesToDo)[j]);

                (*linesToDo)[i].intercepts.push_back(newIntercept);
                (*linesToDo)[j].intercepts.push_back(newIntercept);

                ret.push_back(f1);
                printf("x: %f\ny: %f\n\n", f1.x, f1.y);
            }
        }
    }
    printf("total: %i\n", ret.size());
    return ret;
}

#pragma region Extract polygons
/*
    origin - origin point of the polygon, return when current == current
    polygon - the working collection of points included in the polygon
    current - the intercept to be checked
    whichLine - Each vertex corresponds to one side, so when we push back one vertex, switch to the other
                line participating in the intercept, we pass this because we don't know where we came from
    turn -      We choose the next point based on where we are turning, you can walk from adjacent vertex
                to adjacent vertex making all right or left turns
                1:  left turns
                -1: right turns
                0:  We haven't turned yet, branch to the right and left to get the polygons 
*/
bool polygonize(Intercept* origin, std::vector<Intercept*>* polygon, Intercept* current, Highway* whichLine, int turn)
{
    //These are the same actual intercept pointer
    if (current == origin)
    {
        return true;
    }

    //These both represent the same point but somehow it got copied?
    // I suspect that this is happening because some lines appear far too close, and the algorithm is somehow
    // "clipping" them together when they aren't actually a part of the same poly
    //else if (((current->x - origin->x) < 0.00000001) && ((current->y - origin->y) < 0.00000001))
    //{
    //    return false;
    //}
    //This means a stack overflow is coming
    else if (polygon->size() > 100) {
        return false;
    }

    Intercept* prev = (*polygon)[polygon->size() - 1];
    polygon->push_back(current);

    Highway* other = current->other(whichLine);

    //TODO get the next intercept
    Intercept* p1 = prev;
    Intercept* p2 = current;

    Intercept* q1 = NULL;
    Intercept* q2 = NULL;

    //Get next and previous intercepts on the line
    for (int i = 0; i < other->intercepts.size(); i++)
    {
        if (other->intercepts[i] == current)
        {
            if (i > 0)
            {
                q1 = other->intercepts[i - 1];
            }

            if (i < (other->intercepts.size() - 1))
            {
                q2 = other->intercepts[i + 1];
            }
        }
    }

    //Represents the vector from p1 to p2
    Point v_origin = Point(p2->x - p1->x, p2->y - p1->y);
    
    if (q1)
    {
        // The vector from p1 to q1
        Point o1 = Point(q1->x - p1->x, q1->y - p1->y);

        double dot = (v_origin.y * o1.x) - (v_origin.x * o1.y);

        //Turn is the same
        if ( (dot < 0 && turn < 0) || (dot > 0 && turn > 0) )
        {
            return polygonize(origin, polygon, q1, other, turn);
        }

        // q2 is not null, meaning the current intercept is not the last on this line
        else if (q2)
        {
            return polygonize(origin, polygon, q2, other, turn);
        }

        // 
        else
        {
            return false;
        }
    }

    else if (q2)
    {
        // The vector from p1 to q1
        Point o2 = Point(q2->x - p1->x, q2->y - p1->y);

        int dot = (v_origin.y * o2.x) - (v_origin.x * o2.y);

        //Turn is the same
        if ((dot < 0 && turn < 0) || (dot > 0 && turn > 0))
        {
            return polygonize(origin, polygon, q2, other, turn);
        }

        //Turn is wrong direction, but there is no other intercept to check
        else
        {
            return false;
        }
    }

    return false;
}

std::vector<mPolygon> getPolygons(std::vector<Highway>* mLines)
{
    std::vector<mPolygon> ret;
    //This won't work unless the points are ordered geometrically
    for (int i = 0; i < mLines->size(); i++)
    {
        (*mLines)[i].sortPoints();
    }

    for (int i = 0; i < mLines->size(); i++)
    {
        for (int j = 0; j < (*mLines)[i].intercepts.size() - 1; j++)
        {
            if ((*mLines)[i].intercepts.size() == 0) { break; }
            std::vector<Intercept*> vertices1;
            std::vector<Intercept*> vertices2;

            //We pass in the origin and the next point on the line, i.e. the first side
            vertices1.push_back((*mLines)[i].intercepts[j]);
            vertices2.push_back((*mLines)[i].intercepts[j]);

            //Do twice, once for right turn once for left
            if (polygonize((*mLines)[i].intercepts[j], &vertices1, (*mLines)[i].intercepts[j + 1], &((*mLines)[i]), -1))
            {
                std::sort(vertices1.begin(), vertices1.end(), sortInterceptX);
                bool add1 = true;
                for (int i = 0; i < ret.size(); i++)
                {
                    if (ret[i].vertices.size() == vertices1.size())
                    {
                        bool add = false;
                        for (int j = 0; j < vertices1.size(); j++)
                        {
                            if (vertices1[j]->x != ret[i][j]->x && vertices1[j]->y != ret[i][j]->y)
                            {
                                add = true;
                                break;
                            }
                        }
                        if (!add)
                        {
                            add1 = false;
                            break;
                        }
                    }
                }
                if (add1)
                {
                    ret.push_back(vertices1);
                }
            }

            if (polygonize((*mLines)[i].intercepts[j], &vertices2, (*mLines)[i].intercepts[j + 1], &((*mLines)[i]), 1))
            {
                std::sort(vertices2.begin(), vertices2.end(), sortInterceptX);
                bool add1 = true;
                for (int i = 0; i < ret.size(); i++)
                {
                    if (ret[i].vertices.size() == vertices2.size())
                    {
                        bool add = false;
                        for (int j = 0; j < vertices2.size(); j++)
                        {
                            if (vertices2[j]->x != ret[i][j]->x && vertices2[j]->y != ret[i][j]->y)
                            {
                                add = true;
                                break;
                            }
                        }
                        if (!add)
                        {
                            add1 = false;
                            break;
                        }
                    }
                }
                if (add1)
                {
                    ret.push_back(vertices2);
                }
            }
        }
    }

    printf("%d chunks\n", ret.size());
    return ret;
}


#pragma endregion

#pragma region hats

#pragma endregion

//Get from the given vector the line that the given coordinate lies on
Line getCorresponding(std::vector<Intercept*> hat, double value, char xy)
{
    for (int i = 0; i < hat.size() - 1; i++)
    {

        if (xy == 'x')
        {
            if (value > hat[i]->x && value < hat[i + 1]->x)
            {
                Point a = Point(hat[i]->x, hat[i]->y);
                Point b = Point(hat[i+1]->x, hat[i+1]->y);
                return Line(a, b);
            }
        }

        else
        {
            if (value > hat[i]->y && value < hat[i + 1]->y)
            {
                Point a = Point(hat[i]->x, hat[i]->y);
                Point b = Point(hat[i + 1]->x, hat[i + 1]->y);
                return Line(a, b);
            }
        }
    }
    Point a;
    Point b;
    return Line(a, b);
}

std::vector<Intercept*> rotateToRandom(std::vector<Intercept*> hull, double* returnAngle, double* oX, double* oY)
{
    std::vector<Intercept*> ret;

    for (int i = 0; i < hull.size(); i++)
    {
        ret.push_back(new Intercept(hull[i]->x, hull[i]->y));
    }    
    std::sort(hull.begin(), hull.end(), sortInterceptX);
    
    //Ret now contains copies of hull sorted by x

    int ind = rand() % ret.size();

    Intercept* p1 = ret[ind];
    Intercept* p2 = ret[ (ind+1) % ret.size() ];

    double ox = p1->x;
    double oy = p1->y;

    *oX = ox;
    *oY = oy;

    double angle = atan2((p2->y - p1->y), (p2->x - p1->x));
    *returnAngle = angle;
    double cosTheta = cos(-angle);
    double sinTheta = sin(-angle);

    for (int i = 0; i < ret.size(); i++)
    {
        ret[i]->x -= ox;
        ret[i]->y -= oy;


        double xx = (ret[i]->x * cosTheta) - (ret[i]->y * sinTheta);
        double yy = (ret[i]->x * sinTheta) + (ret[i]->y * cosTheta);

        ret[i]->x = xx;
        ret[i]->y = yy;

        ret[i]->x += ox;
        ret[i]->y += oy;
    }

    return ret;
}

// 
std::vector<Highway> getVerticalStreets(mPolygon hull)
{
    std::vector<Highway> ret;
    double angle, ox, oy;
    mPolygon rotated = rotateToRandom(hull, &angle, &ox, &oy);
    mPolygon top = rotated.topHat();
    mPolygon bottom = rotated.bottomHat();
    mPolygon right = rotated.rightHat();
    mPolygon left = rotated.leftHat();

    std::sort(rotated.vertices.begin(), rotated.vertices.end(), sortInterceptX);

    double minX = rotated[0]->x;
    double maxX = rotated[rotated.vertices.size() - 1]->x;

    double pen = minX + GRID_SIZE;

    std::vector<Point> roads;

    //Draw the "horizontal" roads
    while (pen < (maxX - GRID_SIZE) )
    {
        Line topLine = getCorresponding(top, pen, 'x');
        Line bottomLine = getCorresponding(bottom, pen, 'x');

        double topY = topLine.yIntercept(pen);
        double bottomY = bottomLine.yIntercept(pen);

        roads.push_back(Point(pen, topY));
        roads.push_back(Point(pen, bottomY));

        pen += GRID_SIZE;
    }

    // Sort the points of the hull by y value, get the top and bottom points, and set the pen
    std::sort(rotated.vertices.begin(), rotated.vertices.end(), sortInterceptY);
    double minY = rotated[0]->y;
    double maxY = rotated[rotated.vertices.size() - 1]->y;
    pen = minY + GRID_SIZE;

    while (pen < (maxY - GRID_SIZE))
    {
        Line leftLine = getCorresponding(left, pen, 'y');
        Line rightLine = getCorresponding(right, pen, 'y');

        double leftX = leftLine.xIntercept(pen);
        double rightX = rightLine.xIntercept(pen);

        roads.push_back(Point(leftX, pen));
        roads.push_back(Point(rightX, pen));

        pen += GRID_SIZE;
    }


    for (int i = 0; i < roads.size(); i++)
    {
        roads[i].x -= ox;
        roads[i].y -= oy;
        double xx = (roads[i].x * cos(angle) ) - (roads[i].y * sin(angle) );
        double yy = (roads[i].x * sin(angle) ) + (roads[i].y * cos(angle) );

        roads[i].x = xx;
        roads[i].y = yy;

        roads[i].x += ox;
        roads[i].y += oy;
    }
    if (roads.size() > 0)
    {
        for (int i = 0; i < (roads.size() - 1); i += 2)
        {
            bullshitLines.push_back(Line(roads[i], roads[i + 1]));
            //streets.push_back(Highway(roads[i], roads[i + 1]));
            ret.push_back(Highway(roads[i], roads[i + 1]));
        }
    }

    std::vector<Intercept*> topp = hull.topHat();
    std::vector<Intercept*> bott = hull.bottomHat();

    for (int i = 0; i < topp.size() - 1; i++)
    {
        ret.push_back(Highway(Point(topp[i]->x, topp[i]->y), Point(topp[i + 1]->x, topp[i + 1]->y)));
    }
    for (int i = 0; i < bott.size() - 1; i++)
    {
        ret.push_back(Highway(Point(bott[i]->x, bott[i]->y), Point(bott[i + 1]->x, bott[i + 1]->y)));
    }

    return ret;
}

std::vector<std::vector<Highway>> genSubStreets()
{
    std::vector<std::vector<Highway>> ret;
    for (int i = 0; i < chunks.size(); i++)
    {
        mPolygon hull = chunks[i];
        ////std::vector<Intercept*> topHat = getTopHat(hull);
        ////std::vector<Intercept*> topHat = getBottomHat(hull);
        //for (int j = 0; j < topHat.size() - 1; j++)
        //{
        //    Point a = Point(topHat[j]->x, topHat[j]->y);
        //    Point b = Point(topHat[j+1]->x, topHat[j+1]->y);

        //  //  bullshitLines.push_back(Line(a, b));
        //}
        std::vector<Highway> streetSet = getVerticalStreets(hull.vertices);

        int a = 0;
        ret.push_back(streetSet);
    }
    return ret;
}

std::vector<mPolygon> splitInHalf(mPolygon toSplit) {

    std::vector<mPolygon> ret;

    std::vector<Intercept*> top = toSplit.topHat();
    std::vector<Intercept*> bottom = toSplit.bottomHat();

    // Select a random line segment from the top and bottom hull
    int ind1 = (rand() % (top.size() - 1));
    int ind2 = (rand() % (bottom.size() - 1));

    // Select a random point on the top and bottom lines
    double offX = randRange(top[ind1]->x, top[ind1 + 1]->x);
    double offX2 = randRange(top[ind2]->x, top[ind2 + 1]->x);

    //
    Line topLine = getCorresponding(top, offX, 'x');
    Line bottomLine = getCorresponding(bottom, offX2, 'x');

    double offY = topLine.yIntercept(offX);
    double offY2 = bottomLine.yIntercept(offX2);

    std::vector<Intercept*> split1;
    std::vector<Intercept*> split2;

    Intercept* i1 = new Intercept(offX, offY);
    Intercept* i2 = new Intercept(offX2, offY2);

    Line divider = Line(Point(offX, offY), Point(offX2, offY2));
    Point pp = divider.formula();

    split1.push_back(i1);
    split2.push_back(i1);
    split1.push_back(i2);
    split2.push_back(i2);

    for (int j = 0; j < toSplit.vertices.size(); j++)
    {
        if (toSplit[j]->y > toSplit[j]->x * pp.m() + pp.b())
        {
            split1.push_back(toSplit[j]);
        }
        else
        {
            split2.push_back(toSplit[j]);
        }
    }
    ret.push_back(split1);
    ret.push_back(split2);

    lines.push_back(Highway(Point(offX, offY), Point(offX2, offY2)));

    return ret;
}

std::vector<mPolygon> splitPolygons(std::vector<mPolygon> chunkSet) {

    for (int i = 0; i < NUM_MID_STREETS; i++)
    {
        int ind = rand() % chunkSet.size();
        mPolygon toSplit = chunkSet[ind];
        chunkSet.erase(chunkSet.begin() + ind);

        std::vector<mPolygon> splitted = splitInHalf(toSplit);
        
        chunkSet.insert(chunkSet.end(), splitted.begin(), splitted.end());
    }

    return chunkSet;
}

std::vector<mPolygon> recurse(mPolygon cur, std::vector<Rule> rules) {
	
	std::vector<mPolygon> ret;
	if (rules.size() == 0) {
		return ret;
	}

	std::vector<Rule> next = std::vector<Rule>(rules.begin() + 1, rules.end());
	std::vector<mPolygon> result;

	// TODO: One roundabout
	if (rules[0].rule == 'R') {
		result = cur.addRoundabout();
	}

	// TODO: Ice ray split
	else if (rules[0].rule == 'I') {
		result = cur.split();
	}

	//TODO: convert to multiple roundabouts
	else if (rules[0].rule == 'M') {
		std::vector<Line> needs = cur.addRoundabouts(10000000);
		std::vector<Highway> stret;
		for (int i = 0; i < needs.size(); i++) {
			stret.push_back(Highway(needs[i]));
		}

		getIntersections(&stret);
		result = getPolygons(&stret);
	}

	//TODO: return current chunk as blocks
	else if (rules[0].rule == 'B') {
		GRID_SIZE = STREET_WIDTH;
		std::vector<Highway> newYork = getVerticalStreets(cur);
		getIntersections(&newYork);
		result = getPolygons(&newYork);
	}

	else if (rules[0].rule == 'G') {
		GRID_SIZE = rules[0].minSize;
		std::vector<Highway> newYork = getVerticalStreets(cur);
		getIntersections(&newYork);
		result = getPolygons(&newYork);
	}

	for (int i = 0; i < result.size(); i++) {
		if (result[i].area() > STREET_WIDTH * STREET_WIDTH) {
			std::vector<mPolygon> recurseResult = recurse(result[i], next);

			if (recurseResult.size() > 0) {
				ret.insert(ret.end(), recurseResult.begin(), recurseResult.end());
			}

			else {
				ret.insert(ret.end(), result.begin(), result.end());
			}
		}

		// If the current block is too small to continue splitting, push back the current block with recursing any further
		else {
			ret.push_back(result[i]);
		}
	}

	return ret;
}

int main(int argc, char** argv)
{
    GLuint vertexArrayID = 0;

    GLuint vertexBuffer;
	std::vector<Rule> rules = readRules();

	if (STREET_WIDTH < 1) {
		STREET_WIDTH = 10;
	}

    glutInit(&argc, argv);

    glutInitWindowSize(WIDTH, HEIGHT);
    glutInitWindowPosition(200, 0);

    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("title");

    glewExperimental = GL_TRUE;

    //Glew cannot be initialized until after the window...
    int a = glewInit();

    if (a != GLEW_OK) {
        printf("error: %s\n", glewGetErrorString(a));
        exit(-1);
    }

    lines = genLines();

	if (initialCondition == 'C') {
		chunks.push_back(genPerimeter(numSides));
	}

	else {
		std::vector<Highway> bounds = genBoundary();
		lines.insert(lines.end(), bounds.begin(), bounds.end());
		intersections = getIntersections(&lines);
		chunks = getPolygons(&lines);
	}
	
	printf("done!\n");
	std::vector<mPolygon> fin;
	
	for (int i = 0; i < chunks.size(); i++) {
		std::vector<mPolygon> res = recurse(chunks[i], rules);
	
		fin.insert(fin.end(), res.begin(), res.end());
	}
	
    std::vector<mPolygon> blocks;
    


	//for (int i = 0; i < rules.size(); i++) {
	//
	//}
	//
    //for (int i = 0; i < chunks.size(); i++) {
	////	std::vector<Line> rounds = chunks[i].addRoundabouts(20);
	//	streetInts.push_back(chunks[i].centroid());
	//	if (i % 2) {
	//		std::vector<mPolygon> newLines = chunks[i].addRoundabout();
	//		
	//		for (int j = 0; j < newLines.size(); j++) {
	//			std::vector<mPolygon> newBlock = (newLines[j].iceLatticeSplit());
	//			blocks.insert(blocks.end(), newBlock.begin(), newBlock.end());
	//
	//		}
	//		//blocks.insert(blocks.end(), newLines.begin(), newLines.end());
	//	}
	//
	//	else {
	//		//TODO: New-York-Ize the chunk instead
	//		std::vector<Highway> newYorkChunk = getVerticalStreets(chunks[i]);
	//		getIntersections(&newYorkChunk);
	//		std::vector<mPolygon> curNewYork = getPolygons(&newYorkChunk);
	//		blocks.insert(blocks.end(), curNewYork.begin(), curNewYork.end());
	//	}
	//
	//	
    //}


	//This shrinks all the blocks in the whole world
	//for (int i = 0; i < blocks.size(); i++) {
	////	mPolygon subs = blocks[i].shrinkBlock(0.9);
	//
	//	std::vector<Line> rounds = blocks[i].addRoundabouts(3);
	//	bullshitLines.insert(bullshitLines.end(), rounds.begin(), rounds.end());
	//	//if (subs.vertices.size() > 0) {
	//	//	fin.push_back(subs);
	//	//}
	//}

    //for (int i = 0; i < chunks.size(); i++) {
    //    streetInts.push_back(chunks[i].centroid());
    //
    //    std::vector<mPolygon> newBlock = (chunks[i].iceLatticeSplit());
    //    //blocks.insert(blocks.end(), newBlock.begin(), newBlock.end());
    //    //streetInts.push_back(chunks[i].split());
    //}

    //chunks = splitPolygons(chunks);
    //streetSets = genSubStreets();
    //
    //
    //
    //
    //
    //for (int i = 0; i < streetSets.size(); i++)
    //{
    //    std::vector<Point> cur = getIntersections(&streetSets[i]);
    //    std::vector<mPolygon> curChunks = getPolygons(&streetSets[i]);
    //
    //    blocks.insert(blocks.begin(), curChunks.begin(), curChunks.end());
    //    streetInts.insert(streetInts.begin(), cur.begin(), cur.end());
    //}
    //
    //streetInts.clear();
    //
    //for (int i = 0; i < blocks.size(); i++) {
    //    streetInts.push_back(blocks[i].centroid());
    //}

    while (1) 
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

        /*glUseProgram(program);

        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

        glDrawArrays(GL_TRIANGLES, 0, 3);
        glDisableVertexAttribArray(0);*/

        
        glBegin(GL_LINES);

        // Draw "highways"
        glColor3f(0, 0, 1);
        for (int i = 0; i < lines.size(); i++)
        {
            glVertex2f(lines[i].a.x, lines[i].a.y);
            glVertex2f(lines[i].b.x, lines[i].b.y);
        }

        glColor3f(1, 0, 0);

        glEnd();
        
        
        //drawChunks(chunks);
		//
        //drawChunks(blocks);
		drawChunks(fin);

        glColor3f(1, 0, 0);
        markIntersections(streetInts);
      //  markIntersections(intersections);
        
        glBegin(GL_LINES);
        
        glColor3f(0, 0, 0);
        for (int i = 0; i < bullshitLines.size(); i++)
        {
            glVertex2f(bullshitLines[i].a.x / WIDTH, bullshitLines[i].a.y / HEIGHT);
            glVertex2f(bullshitLines[i].b.x / WIDTH, bullshitLines[i].b.y / HEIGHT);
        }
        glEnd();

        glFlush();
        glutSwapBuffers();
    }
    
    return 0;
}