#include "glew.h"
#include "glut.h"
#include <random>
#include <time.h>
#include "Geometry.h"
#include "Polygon.h"

extern double WIDTH;
extern double HEIGHT;

#define RAND_COLOR_MODE 0

void markIntersections(std::vector<Point> it)
{
    // Mark intersection points
    glBegin(GL_TRIANGLES);

    for (int i = 0; i < it.size(); i++)
    {
        Point p = it[i];
        glVertex2f(p.x / WIDTH, p.y / HEIGHT+ 0.01);
        glVertex2f(p.x / WIDTH - 0.01, p.y / HEIGHT - 0.005);
        glVertex2f(p.x / WIDTH + 0.01, p.y / HEIGHT - 0.005);
    }
    glEnd();
}

void drawChunks(std::vector<mPolygon> cur)
{
    // Draw the chunks
	
    for (int i = 0; i < cur.size(); i++)
    {
        glBegin(GL_POLYGON);

        double c[3] = { 0, 0, 0 };
        
		if (RAND_COLOR_MODE) {
			c[0] = (double)(rand() % 1000) / 1000;
			c[1] = (double)(rand() % 1000) / 1000;
			c[2] = (double)(rand() % 1000) / 1000;
		}

		else {
			int q = i % 3;
			c[q] = ((double)i / cur.size()) + 0.5;
		}

        glColor3f(c[0], c[1], c[2]);

		std::vector<Intercept*> bl = cur[i].perimiterOrdered();

        for (int j = 0; j < bl.size(); j++)
        {
            glVertex2f(bl[j]->x / WIDTH, bl[j]->y / HEIGHT);
        }
        glEnd();
    }
}

bool sortInterceptX(Intercept* i, Intercept* j) { return (i->x < j->x); }
bool sortInterceptY(Intercept* i, Intercept* j) { return (i->y < j->y); }