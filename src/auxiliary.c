/* auxiliary.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "defines.h"
#include "struct.h"

extern grain_struct *grain;
extern point_struct *point;
extern bound_struct *bound;
extern contact_struct *contact;

extern float  t;
extern int    ng;
extern double xmin, xmax, ymin, ymax;


int inbound(double x, double y)
{
	float boundx, boundy;

	boundx = 0.0*xmax;
	boundy = 0.0*ymax;

	if (x>xmax+boundx) return(0);
	if (x<xmin-boundx) return(0);
	if (y>ymax+boundy) return(0);
	if (y<ymin-boundy) return(0);
	return(1);
}

double calc_grain_area(int grn)
{
	int 	i;
	double 	x1, y1, x2, y2;
	double 	area;

	area= 0;
	for (i=0; i<grain[grn].numbound; i++) {
		x1 = bound[grain[grn].bound[i]].x1;
		y1 = bound[grain[grn].bound[i]].y1;	
		x2 = bound[grain[grn].bound[i]].x2;
		y2 = bound[grain[grn].bound[i]].y2;
		area += x1*y2 - y1*x2;
	}
	area *= 0.5;

	return(area);
}

/*
This routine determines if the distance r is greater than
the separation between points x1,y1 and x2,y2.
*/

int meet(double x1, double y1, double x2, double y2, double r)
{
	if ( ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)) > r*r )
		return(0);

	return(1);
}
