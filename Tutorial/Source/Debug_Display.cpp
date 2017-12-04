#include "glew.h"
#include "glut.h"
#include <random>
#include <time.h>
#include "Geometry.h"
#include "Polygon.h"

void markIntersections(std::vector<Point> it)
{
    // Mark intersection points
    glBegin(GL_TRIANGLES);

    for (int i = 0; i < it.size(); i++)
    {
        Point p = it[i];
        glVertex2f(p.x, p.y + 0.01);
        glVertex2f(p.x - 0.01, p.y - 0.005);
        glVertex2f(p.x + 0.01, p.y - 0.005);
    }
    glEnd();
}

void drawChunks(std::vector<mPolygon> cur)
{
    // Draw the chunks
    for (int i = 0; i < cur.size(); i++)
    {
        glBegin(GL_POLYGON);

        float c[3] = { 0, 0, 0 };
        int q = i % 3;
        c[q] = ((float)i / cur.size()) + 0.5;

        glColor3f(c[0], c[1], c[2]);

        for (int j = 0; j < cur[i].vertices.size(); j++)
        {
            glVertex2f(cur[i][j]->x, cur[i][j]->y);
        }
        glEnd();
    }
}

bool sortInterceptX(Intercept* i, Intercept* j) { return (i->x < j->x); }
bool sortInterceptY(Intercept* i, Intercept* j) { return (i->y < j->y); }