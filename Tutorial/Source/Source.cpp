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

std::vector<Highway> lines;
std::vector<Highway> streets;
std::vector<std::vector<Highway>> streetSets;

std::vector< std::vector<Intercept*> > chunks;

std::vector<Point> intersections;

std::vector<Line> bullshitLines;

static const GLfloat vertices[] =
{
    -1.0f, -1.0f, 0.0f,
    1.0f, -1.0f, 0.0f,
    0.0f, 1.0f, 0.0f
};

const int NUM_VERTICAL_HIGHWAYS = 1;
const int NUM_HORIZONTAL_HIGHWAYS = 1;
const int NUM_MID_STREETS = 2;

const float STREET_WIDTH = 0.15;

// Likelihood that a given highway will be parallel to an axis rather than at an angle
const float HIGHWAY_TAXI_FACTOR = 0.0;

GLuint loadShader(const char* path, const int type) {
    //Initialize a slot for a vertex shader
    GLuint shaderID = glCreateShader(type);

    //Read the file to a string
    std::string code;
    std::ifstream VertexShaderStream(path, std::ios::in);
    if (VertexShaderStream.is_open()) {
        std::string Line = "";
        while (getline(VertexShaderStream, Line))
            code += "\n" + Line;
        VertexShaderStream.close();
    }
    else {
        printf("Impossible to open %s. Are you in the right directory ? Don't forget to read the FAQ !\n", path);
        getchar();
        return 0;
    }
    GLint result;
    int logLength;

    const char* vertexSource = code.c_str();

    //Set the source code for the new shader to the newly read source
    glShaderSource(shaderID, 1, &vertexSource, NULL);

    //Compile the shader
    glCompileShader(shaderID);

    // get GL_COMPILE_STATUS (status code) from shaderID and store in result
    glGetShaderiv(shaderID, GL_COMPILE_STATUS, &result);
    // get GL_INFO_LOG_LENGTH (length of log for given shader
    glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &logLength);

    if (logLength > 0) {
        std::vector<char> FragmentShaderErrorMessage(logLength + 1);
        glGetShaderInfoLog(shaderID, logLength, NULL, &FragmentShaderErrorMessage[0]);
        printf("%s\n", &FragmentShaderErrorMessage[0]);
    }
    return shaderID;
}

GLuint linkShaders(std::vector<GLuint> shaders) {
    //Initialize a new program
    GLuint program = glCreateProgram();

    //Attach each of the given shaders to the program
    for (int i = 0; i < shaders.size(); i++) {
        glAttachShader(program, shaders[i]);
    }

    //Link the program with the attached shaders
    glLinkProgram(program);

    int result;
    int logLength;

    // Get link status code and info log length
    glGetProgramiv(program, GL_LINK_STATUS, &result);
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLength);

    //Print errors if any
    if (logLength > 0) {
        std::vector<char> msg(logLength + 1);
        glGetProgramInfoLog(program, logLength, NULL, &msg[0]);
        printf("%s\n", &msg[0]);
    }

    //Detach and delete shaders from the program object, effects will
    // not come into effect until program is relinked
    for (int i = 0; i < shaders.size(); i++) {
        glDetachShader(program, shaders[i]);
        glDeleteShader(shaders[i]);
    }
    

    return program;
}

#pragma region Highway generation
std::vector<Highway> genTicTacToe()
{
    std::vector<Highway> ret;
    ret.push_back(Highway(Point(-1.0f, 0.2), Point(1.0f, 0.2) ) );
    //ret.push_back(Highway(Point(-1.0f, 0.5), Point(1.0f, 0.2)));
    ret.push_back(Highway(Point(-1.0f, -0.5), Point(1.0f, -0.3)));
    

    ret.push_back(Highway(Point(0.2f, -1.0f), Point(0.5f, 1.0f)));
    ret.push_back(Highway(Point(-0.2, -1.0f), Point(-0.5, 1.0f)));
    return ret;
}

std::vector<Highway> genBoundary()
{
    std::vector<Highway> ret;

    ret.push_back(Highway(Point(-1.0f, -0.99f), Point(1.0f, -0.99f)));
    ret.push_back(Highway(Point(-1.0f, 0.99f), Point(1.0f, 0.99f)));
    //
    ret.push_back(Highway(Point(-0.98f, -1.0f), Point(-0.99f, 1.0f)));
    ret.push_back(Highway(Point(0.98, -1.0f), Point(0.99f, 1.0f)));

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

        p1.x = ((float)(rand() % 2000) / 1000.0f) - 1.0f;
        p1.y = -1.0f;

        if (((rand() % 1000) / 1000.0f) < HIGHWAY_TAXI_FACTOR)
        {
            p2.x = p1.x;
        }
        else
        {
            p2.x = ((float)(rand() % 2000) / 1000.0f) - 1.0f;
        }
        
        p2.y = 1.0f;

        Line L = Line(p1, p2);

        ret.push_back(L);
    }

    for (int i = 0; i < NUM_HORIZONTAL_HIGHWAYS; i++)
    {
        Point p1;
        Point p2;

        p1.y = ((float)(rand() % 2000) / 1000.0f) - 1.0f;
        p1.x = -1.0f;

        if (((rand() % 1000) / 1000.0f) < HIGHWAY_TAXI_FACTOR)
        {
            p2.y = p1.y;
        }
        else
        {
            p2.y = ((float)(rand() % 2000) / 1000.0f) - 1.0f;
        }

        p2.x = 1.0f;

        Line L = Line(p1, p2);

        ret.push_back(L);
    }

    return ret;
}

std::vector<Highway> genOffshoots(std::vector<Highway> orig)
{
    for (int i = 0; i < NUM_MID_STREETS; i++)
    {
        int ind = rand() % orig.size();

        Point form = orig[ind].formula();

        float minX = std::min(orig[ind].a.x, orig[ind].b.x);
        float maxX = std::max(orig[ind].a.x, orig[ind].b.x);

        float offX = minX + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (maxX - minX)));
        float offY = offX * form.m() + form.b();

        float endX = 0;
        float endY = 0;

        float newX, newY;

        float slope = -1.0/form.m();
        float newB = offY - (slope * offX);

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

        float dot = (v_origin.y * o1.x) - (v_origin.x * o1.y);

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

std::vector<std::vector<Intercept*>> getPolygons(std::vector<Highway>* mLines)
{
    std::vector<std::vector<Intercept*>> ret;
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
                    if (ret[i].size() == vertices1.size())
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
                    if (ret[i].size() == vertices2.size())
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
std::vector<Intercept*> getTopHat(std::vector<Intercept*> hull)
{
    std::vector<Intercept*> topHat;
    std::vector<Intercept*> topHat2;

    std::sort(hull.begin(), hull.end(), sortInterceptX);

    topHat.push_back(hull[0]);

    Line topLine = Line(Point(hull[0]->x, hull[0]->y), Point(hull[hull.size() - 1]->x, hull[hull.size() - 1]->y));
    Point topF = topLine.formula();
    for (int i = 1; i < hull.size() - 1; i++)
    {
        if (hull[i]->y > ( (hull[i]->x * topF.m()) + topF.b()) )
        {
            topHat.push_back(hull[i]);
        }
    }
    topHat.push_back(hull[hull.size() - 1]);

    return topHat;
}

std::vector<Intercept*> getBottomHat(std::vector<Intercept*> hull)

{
    std::vector<Intercept*> topHat;

    std::sort(hull.begin(), hull.end(), sortInterceptX);

    topHat.push_back(hull[0]);
    Line topLine = Line(Point(hull[0]->x, hull[0]->y), Point(hull[hull.size() - 1]->x, hull[hull.size() - 1]->y));
    Point topF = topLine.formula();
    for (int i = 1; i < hull.size() - 1; i++)
    {
        if (hull[i]->y  < ((hull[i]->x * topF.m()) + topF.b()))
        {
            topHat.push_back(hull[i]);
        }
    }
    topHat.push_back(hull[hull.size() - 1]);
    return topHat;
}

std::vector<Intercept*> getLeftHat(std::vector<Intercept*> hull)

{
    std::vector<Intercept*> topHat;

    std::sort(hull.begin(), hull.end(), sortInterceptY);

    topHat.push_back(hull[0]);
    Line topLine = Line(Point(hull[0]->x, hull[0]->y), Point(hull[hull.size() - 1]->x, hull[hull.size() - 1]->y));
    Point topF = topLine.formula();
    for (int i = 1; i < hull.size() - 1; i++)
    {
        if (hull[i]->x < ((hull[i]->y - topF.b()) / topF.m()))
        {
            topHat.push_back(hull[i]);
        }
    }
    topHat.push_back(hull[hull.size() - 1]);
    return topHat;
}

std::vector<Intercept*> getRightHat(std::vector<Intercept*> hull)

{
    std::vector<Intercept*> topHat;

    std::sort(hull.begin(), hull.end(), sortInterceptY);

    topHat.push_back(hull[0]);
    Line topLine = Line(Point(hull[0]->x, hull[0]->y), Point(hull[hull.size() - 1]->x, hull[hull.size() - 1]->y));
    Point topF = topLine.formula();
    for (int i = 1; i < hull.size() - 1; i++)
    {
        if (hull[i]->x > ((hull[i]->y - topF.b()) / topF.m()))
        {
            topHat.push_back(hull[i]);
        }
    }
    topHat.push_back(hull[hull.size() - 1]);
    return topHat;
}
#pragma endregion

//Get from the given vector the line that the given coordinate lies on
Line getCorresponding(std::vector<Intercept*> hat, float value, char xy)
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

std::vector<Intercept*> rotateToRandom(std::vector<Intercept*> hull, float* returnAngle, float* oX, float* oY)
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

    float ox = p1->x;
    float oy = p1->y;

    *oX = ox;
    *oY = oy;

    float angle = atan2((p2->y - p1->y), (p2->x - p1->x));
    *returnAngle = angle;
    float cosTheta = cos(-angle);
    float sinTheta = sin(-angle);

    for (int i = 0; i < ret.size(); i++)
    {
        ret[i]->x -= ox;
        ret[i]->y -= oy;


        float xx = (ret[i]->x * cosTheta) - (ret[i]->y * sinTheta);
        float yy = (ret[i]->x * sinTheta) + (ret[i]->y * cosTheta);

        ret[i]->x = xx;
        ret[i]->y = yy;

        ret[i]->x += ox;
        ret[i]->y += oy;
    }

    return ret;
}

// 
std::vector<Highway> getVerticalStreets(std::vector<Intercept*> hull)
{
    std::vector<Highway> ret;
    float angle, ox, oy;
    std::vector<Intercept*> rotated = rotateToRandom(hull, &angle, &ox, &oy);
    std::vector<Intercept*> top = getTopHat(rotated);
    std::vector<Intercept*> bottom = getBottomHat(rotated);
    std::vector<Intercept*> right = getRightHat(rotated);
    std::vector<Intercept*> left = getLeftHat(rotated);

    std::sort(rotated.begin(), rotated.end(), sortInterceptX);

    float minX = rotated[0]->x;
    float maxX = rotated[rotated.size() - 1]->x;

    float pen = minX + STREET_WIDTH;

    std::vector<Point> roads;

    //Draw the "horizontal" roads
    while (pen < (maxX - STREET_WIDTH) )
    {
        Line topLine = getCorresponding(top, pen, 'x');
        Line bottomLine = getCorresponding(bottom, pen, 'x');

        float topY = topLine.yIntercept(pen);
        float bottomY = bottomLine.yIntercept(pen);

        roads.push_back(Point(pen, topY));
        roads.push_back(Point(pen, bottomY));

        pen += STREET_WIDTH;
    }

    // Sort the points of the hull by y value, get the top and bottom points, and set the pen
    std::sort(rotated.begin(), rotated.end(), sortInterceptY);
    float minY = rotated[0]->y;
    float maxY = rotated[rotated.size() - 1]->y;
    pen = minY + STREET_WIDTH;

    while (pen < (maxY - STREET_WIDTH))
    {
        Line leftLine = getCorresponding(left, pen, 'y');
        Line rightLine = getCorresponding(right, pen, 'y');

        float leftX = leftLine.xIntercept(pen);
        float rightX = rightLine.xIntercept(pen);

        roads.push_back(Point(leftX, pen));
        roads.push_back(Point(rightX, pen));

        pen += STREET_WIDTH;
    }


    for (int i = 0; i < roads.size(); i++)
    {
        roads[i].x -= ox;
        roads[i].y -= oy;
        float xx = (roads[i].x * cos(angle) ) - (roads[i].y * sin(angle) );
        float yy = (roads[i].x * sin(angle) ) + (roads[i].y * cos(angle) );

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

    std::vector<Intercept*> topp = getTopHat(hull);
    std::vector<Intercept*> bott = getBottomHat(hull);

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
        std::vector<Intercept*> hull = chunks[i];
        ////std::vector<Intercept*> topHat = getTopHat(hull);
        ////std::vector<Intercept*> topHat = getBottomHat(hull);
        //for (int j = 0; j < topHat.size() - 1; j++)
        //{
        //    Point a = Point(topHat[j]->x, topHat[j]->y);
        //    Point b = Point(topHat[j+1]->x, topHat[j+1]->y);

        //  //  bullshitLines.push_back(Line(a, b));
        //}
        std::vector<Highway> streetSet = getVerticalStreets(hull);

        int a = 0;
        ret.push_back(streetSet);
    }
    return ret;
}

std::vector<std::vector<Intercept*>> splitInHalf(std::vector<Intercept*> toSplit) {

    std::vector<std::vector<Intercept*>> ret;

    std::vector<Intercept*> top = getTopHat(toSplit);
    std::vector<Intercept*> bottom = getBottomHat(toSplit);

    // Select a random line segment from the top and bottom hull
    int ind1 = (rand() % (top.size() - 1));
    int ind2 = (rand() % (bottom.size() - 1));

    // Select a random point on the top and bottom lines
    float offX = top[ind1]->x + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (top[ind1 + 1]->x - top[ind1]->x)));
    float offX2 = bottom[ind2]->x + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (bottom[ind2 + 1]->x - bottom[ind2]->x)));

    //
    Line topLine = getCorresponding(top, offX, 'x');
    Line bottomLine = getCorresponding(bottom, offX2, 'x');

    float offY = topLine.yIntercept(offX);
    float offY2 = bottomLine.yIntercept(offX2);

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

    for (int j = 0; j < toSplit.size(); j++)
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

std::vector<std::vector<Intercept*>> splitPolygons(std::vector<std::vector<Intercept*>> chunkSet) {

    for (int i = 0; i < NUM_MID_STREETS; i++)
    {
        int ind = rand() % chunkSet.size();
        std::vector<Intercept*> toSplit = chunkSet[ind];
        chunkSet.erase(chunkSet.begin() + ind);

        std::vector<std::vector<Intercept*>> splitted = splitInHalf(toSplit);
        
        chunkSet.insert(chunkSet.end(), splitted.begin(), splitted.end());
    }

    return chunkSet;
}

int main(int argc, char** argv)
{
    GLuint vertexArrayID = 0;

    GLuint vertexBuffer;

    glutInit(&argc, argv);

    glutInitWindowSize(800, 800);
    glutInitWindowPosition(200, 40);

    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("title");

    glewExperimental = GL_TRUE;

    //Glew cannot be initialized until after the window...
    int a = glewInit();

    if (a != GLEW_OK) {
        printf("error: %s\n", glewGetErrorString(a));
        exit(-1);
    }

    GLuint* aa = &vertexArrayID;
    glGenVertexArrays(1, aa);
    glBindVertexArray(vertexArrayID);

    glGenBuffers(1, &vertexBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);

    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    //Load shaders and link a program
//    std::vector<GLuint> shaders;
//    shaders.push_back(loadShader("vertex.glsl", GL_VERTEX_SHADER));
//    shaders.push_back(loadShader("fragment.glsl", GL_FRAGMENT_SHADER));
//    GLuint program = linkShaders(shaders);


    lines = genLines();
 //   lines = genOffshoots(lines);
    
    std::vector<Highway> bounds = genBoundary();

    lines.insert(lines.end(), bounds.begin(), bounds.end());

    intersections = getIntersections(&lines);

    chunks = getPolygons(&lines);
    chunks = splitPolygons(chunks);
    streetSets = genSubStreets();

    std::vector<Point> streetInts;
    std::vector<std::vector<Intercept*>> blocks;

    

    for (int i = 0; i < streetSets.size(); i++)
    {
        std::vector<Point> cur = getIntersections(&streetSets[i]);
        std::vector<std::vector<Intercept*>> curChunks = getPolygons(&streetSets[i]);

        blocks.insert(blocks.begin(), curChunks.begin(), curChunks.end());
        streetInts.insert(streetInts.begin(), cur.begin(), cur.end());
    }

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

        drawChunks(blocks);

        glColor3f(1, 0, 0);
      //  markIntersections(streetInts);
      //  markIntersections(intersections);
        
        glBegin(GL_LINES);
        
        glColor3f(0, 0, 0);
        for (int i = 0; i < bullshitLines.size(); i++)
        {
            glVertex2f(bullshitLines[i].a.x, bullshitLines[i].a.y);
            glVertex2f(bullshitLines[i].b.x, bullshitLines[i].b.y);
        }
        glEnd();

        glFlush();
        glutSwapBuffers();
    }
    
    return 0;
}