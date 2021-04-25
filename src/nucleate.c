/* nucleate.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "defines.h"
#include "struct.h"

extern grain_struct   *grain;
extern point_struct   *point;
extern bound_struct   *bound;
extern contact_struct *contact;

extern float  t;
extern int    ng;
extern float  g0;
extern float  scaling;
extern double xmax, ymax;

/* nucleate */
void new_grain (double, double);

void nucleate_island (float *time_step)
{
	/* Random position on substrate for nucleation attempt */
	double xx, yy;
	/* Distance between existing island and position of nulceation attempt */
	double dx, dy, r;
	/* If random nucleation position not contained within existing island */
	int okay;
	int i;

	/* Random time between nucleation attempts */
	*time_step = 2.0*((1+rand())/(double)(1.0+RAND_MAX)) / (xmax*ymax);

	/* Increase total simulation by time step */
	t += *time_step;
	/* Random position for current nucleation attempt */
	xx = rand()/(double)(RAND_MAX) * xmax;
	yy = rand()/(double)(RAND_MAX) * ymax;

	/* Check that new nucleation site outside any existing islands. */
	okay = 1;
	for (i=0; i<ng && okay; i++) {
		r  = g0*(t - grain[i].t);
		dx = grain[i].x - xx;
		dy = grain[i].y - yy;
		if ( (dx*dx + dy*dy) < (r*r) )
			okay=0;
	}

	/*
	If okay then record the new nuclei. Duplicate the nuclei
	eight times (once for each direction from the center rectangle) .
	*/
	if (okay) {
		new_grain(xx,	   yy-ymax);
		new_grain(xx,	   yy);
		new_grain(xx,	   yy+ymax);
		new_grain(xx-xmax, yy-ymax);
		new_grain(xx-xmax, yy);
		new_grain(xx-xmax, yy+ymax);
		new_grain(xx+xmax, yy-ymax);
		new_grain(xx+xmax, yy);
		new_grain(xx+xmax, yy+ymax);
	}
}

void new_grain(double x, double y)
{
	grain[ng].t		= t;
	grain[ng].x		= x;
	grain[ng].y		= y;

	grain[ng].area		= 0;
	grain[ng].delete_flag = 0;
	grain[ng].numcontact	= 0;
	grain[ng].stress	= 0;

	ng++;
}
