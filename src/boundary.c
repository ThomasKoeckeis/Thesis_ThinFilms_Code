/* boundary.c */

/* contact.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "defines.h"
#include "struct.h"

extern grain_struct *grain;
extern point_struct *point;
extern bound_struct *bound;
extern contact_struct *contact;

extern float 	t;
extern int 	ng;
extern float 	scaling;
extern int long numbound;
extern double 	xmin, xmax, ymin, ymax;
extern float 	g0;

float xcursor, ycursor;

/* boundary.c */
int move_to (double, double);
int draw_to (double, double, int);

/*
Subroutine to draw the hyperbola that forms the boundary
between two grains which nucleate at different times.
This is done by a simplified method, developed Dec. 1984, which
uses a transformation of coordinates to align the local axes,
(u,v), with the hyperbola.
x1,y1,t1 specifies the first nucleus
xn1,yn1,tn1
x2,y2,t2
xn2,yn2,tn2
*/

int drawhyp(int g1, int g2, int p1, int p2)
{
	double x1, y1, t1, x2, y2, t2;
	double xn1, yn1, tn1, xn2, yn2, tn2;
	double dx,dy,dr,ct,xc,yc,st,un1,vn1,un2,vn2;
	double a,b,c,b2,v,u,dv,x,y,xx,yy;
	int i,points;

	x1 = grain[g1].x;	y1 = grain[g1].y;	t1 = grain[g1].t;
	x2 = grain[g2].x;	y2 = grain[g2].y;	t2 = grain[g2].t;

	xn1= point[p1].x;	yn1= point[p1].y;	tn1= point[p1].t;
	xn2= point[p2].x;	yn2= point[p2].y;	tn2= point[p2].t;

	/* Find the parameters for the coordinate transformation. */
	dx = x2 - x1;
	dy = y2 - y1;
	dr = sqrt(dx*dx + dy*dy);
	ct = dx/dr;
	st = dy/dr;
	xc = (x1 + x2)/2.;
	yc = (y1 + y2)/2.;

	/* Transform the coordinates of the triple point nodes. */
	un1 =  (xn1-xc)*ct + (yn1-yc)*st;
	vn1 = -(xn1-xc)*st + (yn1-yc)*ct;
	un2 =  (xn2-xc)*ct + (yn2-yc)*st;
	vn2 = -(xn2-xc)*st + (yn2-yc)*ct;

	/* Specify the hyperbola. */
	a = 0.5 * g0 * (t2-t1);
	c = 0.5 * dr;
	b2 = c*c - a*a;
	b = sqrt(b2);

	/*
	** Determine the number of segment points we need to keep
	** the segment point spacing close to the desired spacing. Don't
	** allow there to be zero segment points as GB_INIT will barf.
	** The "points" variable is actually equal the number of segment
	** points plus 1.
	*/
	if (TOOCLOSE > TOOSMALL3) {
		points = fabs(vn2-vn1)/TOOCLOSE;
		if (points < 2) points = 2;
	}
	else
	points = 8;

	/* Tell where to start */
	move_to(xn1,yn1);

	/* Increment over values of v to draw the segment. */
	dv = (vn2-vn1)/( (double) (points) );
	v = vn1;
	for (i=0;i<points-1;i++) {
		v = v + dv;
		u = a * sqrt(1. + v*v/b2);
		x = u*ct - v*st + xc;
		y = u*st + v*ct + yc;
	
		draw_to (x,y,g1);

		xx = x;	yy = y;

		if (x < xmin) xx = x + (xmax-xmin);
		if (x > xmax) xx = x - (xmax-xmin);
		if (y < ymin) yy = y + (ymax-ymin);
		if (y > ymax) yy = y - (ymax-ymin);
	}

	/* Tell where to stop */
	draw_to (xn2,yn2,g1);

	return(1) ;
}

int draw_circle_full(int g)
{
	double 	xn, yn, rn;
	int 	points;
	double 	theta, dt;
	double 	dx, dy;
	double 	x, y;
	int 	i;

	xn = grain[g].x;	yn = grain[g].y;	rn = grain[g].r;

	points= 2*PI*rn/TOOCLOSE;
	if (points < 2) points= 2;

	/* Tell where to stop */
	move_to (xn+rn,yn);

	dt    = 2*PI/( (double) (points) );
	theta = 0;
	for (i=0; i<points; i++) {
		dx = rn*cos(theta);
		dy = rn*sin(theta);
		x  = xn+dx;
		y  = yn+dy;

		draw_to(x,y,g);
		theta += dt;
	}

	/* Tell where to stop */
	draw_to(xn+rn,yn,g);
	return(1);
}

int draw_circle(int g, int p0, int p1)
{
	double xn, yn, rn;
	double x0, y0, theta0;
	double x1, y1, theta1;
	double dc[2][2];
	double dtheta;
	double arc_length;
	int points;
	double theta, dt;
	double dx_prime, dy_prime;
	double dx, dy;
	double x, y;
	int i;

	xn = grain[g].x;	yn = grain[g].y;	rn= grain[g].r;
	x0 = point[p0].x;	y0 = point[p0].y;	theta0 = point[p0].theta;
	x1 = point[p1].x;	y1 = point[p1].y;	theta1 = point[p1].theta;

	/* calculate direction cosine rotation matrix */
	dc[0][0] =  cos(theta0);
	dc[0][1] =  sin(theta0);
	dc[1][0] = -sin(theta0);
	dc[1][1] =  cos(theta0);

	dtheta= theta1 - theta0;
	if (dtheta < 0)
		dtheta += 2*PI;

	arc_length= dtheta*rn;

	points= arc_length/TOOCLOSE;
	if (points < 2) points= 2;

	/* Tell where to stop */
	move_to (x0,y0);

	dt 	= (dtheta)/( (double) (points) );
	theta	= 0;

	for (i=0; i<points-1; i++) {
		dx_prime = rn*cos(theta);
		dy_prime = rn*sin(theta);
		dx	 = dc[0][0]*dx_prime + dc[1][0]*dy_prime;
		dy	 = dc[0][1]*dx_prime + dc[1][1]*dy_prime;
		x	 = xn+dx;
		y	 = yn+dy;
		draw_to(x,y,g);
		theta += dt;
	}

	/* Tell where to stop */
	draw_to(x1,y1,g);

	return(1);
}

int move_to (double x, double y)
{
	xcursor = x;	ycursor = y;
	return(1) ;
}

int draw_to(double x, double y, int i)
{
	bound[numbound].x1 = xcursor;
	bound[numbound].y1 = ycursor;
	bound[numbound].x2 = x;
	bound[numbound].y2 = y;
	grain[i].bound[grain[i].numbound++] = numbound++;

	xcursor=x; ycursor=y;
	return (1) ;
}
