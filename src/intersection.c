/* intersection.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "defines.h"
#include "struct.h"

extern grain_struct   *grain;
extern point_struct   *point;
extern bound_struct   *bound;
extern contact_struct *contact;

extern float 	t;
extern float 	g0;
extern int 	ng;
extern double 	xmin, xmax, ymin, ymax;
extern int long	numbound;

/* auxiliary.c */
extern int meet (double, double, double, double, double);

/* intersection.c */
int  same_point(int, int);
void circle_intersect(int, int, double *, double *);
int  find_triple(int, int, int, int *, double, double, double, int *);
int  hyper(double, double, double, double, double, double, 
	   double, double, double, double *, double *, double *);
double tp_grain_angle(int, int, int);

/* boundary.c */
extern int drawhyp(int, int, int, int);
extern int draw_circle_full(int);
extern int draw_circle (int, int, int);

void island_intersection(int *numpoint) {

	double 	xi, yi, ri;
	double 	xj, yj, rj;
	double 	xk, yk, rk;
	int 	near[MAXSIDES*10];
	double 	xci[2], yci[2];
	int 	intersect;
	int 	i, j, k, p, nn, jj, kk;

	/* The "points" are those which results from the intersection of
	** two circles. Since the circles grow each timestep, they must be
	** recalculated each timestep, so initialize numpoint to zero.
	*/
	*numpoint = 0;

	/* Determine which grains intersect. */
	for (i=0; i<ng; i++) {
		/* Store data for current grain, i */
		xi = grain[i].x;	yi = grain[i].y;	ri = grain[i].r;

		/* Generate a list of possible neighbors. */
		nn = 0;
		for (j=0; j<ng; j++) {
			if (j != i) {
				xj = grain[j].x;	yj = grain[j].y;	rj = grain[j].r;
				if ( meet(xi, yi, xj, yj, ri+rj) ) {
					near[nn++]= j;
				}
			}
		}

		/* Find new points of intersection of grains */
		for (jj=0; jj<nn; jj++) {
			j = near[jj];
			/* Determine the x,y coords of the two points of intersection
			** of the grain, i, and the neighboring grain, j.
			*/
			circle_intersect(i,j,xci,yci);

			/* Determine if any other grains have intersected the x,y coord
			** of the intersection point. If no other intersections are
			** found, then the point is recorded. These type of points will
			** move each timestep and are therefore recalculated each timestep.
			** If the point is intersected, the point must be part of a triple
			point.
			*/
			for (p=0; p<2; p++) {
				intersect= 0;
	
				/* loop for the kth neighbor */
				for (kk=0; kk<nn; kk++) {
					k = near[kk];
					if (k != j) {
						xk = grain[k].x; yk = grain[k].y; rk = grain[k].r;
						if ( meet(xk,yk,xci[p],yci[p],rk) )
							intersect= 1;
					}
				}
				if (!intersect) {
					point[*numpoint].x	  = xci[p];
					point[*numpoint].y	  = yci[p];
					point[*numpoint].t	  = t;
					point[*numpoint].grain[0] = j;
					point[*numpoint].grain[1] = j;
					point[*numpoint].theta	  = atan2(yci[p]-yi,xci[p]-xi);
					point[*numpoint].cw	  = p % 2;
					point[*numpoint].ccw	  = (p+1) % 2;
					grain[i].point[grain[i].numpoint++] = (*numpoint)++;
				}
				else {
					find_triple(i,jj,nn,near,xci[p],yci[p],t,numpoint);
				}
			}
		}
	}/* loop for ith grain */
}

void organize_intersection_points()
{

	int 	pt;
	double 	theta0, theta1, thetac;
	int 	grn;
	int 	i, j, k, l, m;

	numbound= 0;

	/* Now organize the points around each grain */
	for (i=0; i<ng; i++) {
		/*
		** Determine ccw direction of grains around triple points.
		** If grain0 and grain1 equal, then only circle intersection point .
		*/
		for (j=0; j < grain[i].numpoint; j++) {
			pt = grain[i].point[j];
			if (point[pt].grain[0] != point[pt].grain[1]) {
				theta0= tp_grain_angle(i, pt, point[pt].grain[0]);
				theta1= tp_grain_angle(i, pt, point[pt].grain[1]);

				if (theta0 < theta1) {
					grn= point[pt].grain[0];
					point[pt].grain[0] = point[pt].grain[1];
					point[pt].grain[1] = grn;
				}
			}
		}

		/* Sort intersection points around grain */
		for (j=1; j<grain[i].numpoint; j++) {
			pt	= grain[i].point[j];
			thetac	= point[pt].theta;
			m	= 0;
			while ( (thetac>point[grain[i].point[m]].theta) && m<j) m++;
			for (l=j; l>m; l--)
				grain[i].point[l]= grain[i].point[l-1];

			grain[i].point[m]= pt;
		}

		/* Remove any duplicate points */
		for (j=0; j<grain[i].numpoint-1; j++)
			if ( same_point(grain[i].point[j], grain[i].point[j+1]) ) {
				for (k=j; k<grain[i].numpoint-1; k++)
					grain[i].point[k]= grain[i].point[k+1];
				grain[i].numpoint--;
				j--;
			}
	}
}

void draw_island_boundary()
{
	int pt_start, pt_stop;
	int i, j, k;

	for (i=0; i<ng; i++) {
		/* No intersection points then still a circular island */
		if (grain[i].numpoint == 0)
			draw_circle_full(i);
		/* Traverse intersection points around grain and calculate the length
		** of the grain boundary and the grain size.
		*/
		for (j=0; j<grain[i].numpoint; j++) {
			pt_start = grain[i].point[j];

			if (j<grain[i].numpoint-1)
				pt_stop = grain[i].point[j+1];
			else
				pt_stop = grain[i].point[0];
	
			k = point[pt_start].grain[0];

			if (same_point(pt_start, pt_stop) ) {
				/* printf("DUPLICATE POINTS IN STRUCTURE?\n"); */
				continue;
			}

			if ( (point[pt_start].ccw==1) && (point[pt_stop].cw==1) ) {
				drawhyp(i,k,pt_start,pt_stop) ;
			}
			if ( (point[pt_start].ccw==0) && (point[pt_stop].cw==0) ) {
				draw_circle(i,pt_start,pt_stop);
			}
		}
	}
}

int same_point(int p1, int p2)
{
	double x1, y1, x2, y2;
	double distance;

	x1 = point[p1].x;	y1 = point[p1].y;
	x2 = point[p2].x;	y2 = point[p2].y;
	distance = sqrt( (x2-x1)*(x2-x1) + (y2-y1)* (y2-y1) );
	if (distance < TOL1)
		return(1);
	else
		return(0);
}

void circle_intersect(int i, int j, double *xci, double *yci)
{
	double x1, y1, r1;
	double x2, y2, r2;
	double r3;
	double theta1, theta2, theta3;
	double theta;
	double dx_prime, dy_prime;
	double dx, dy;
	double dc[2][2];

	x1 = grain[i].x;	y1 = grain[i].y;	r1 = grain[i].r;
	x2 = grain[j].x;	y2 = grain[j].y;	r2 = grain[j].r;

	r3= sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );

	/* l:m:s r1:r2:r3 */
	if ( (r1>r2) && (r2>r3) && (r1>r3) ) {
		theta1= acos( (r2*r2 + r3*r3 - r1*r1)/(2.0*r2*r3) );
		theta2= asin( (r2/r1)*sin(theta1) );
		theta3= asin( (r3/r1)*sin(theta1) );
	}
	/* l:m:s r1:r3:r2 */
	if ( (r1>r2) && (r2<r3) && (r1>r3) ) {
		theta1= acos( (r2*r2 + r3*r3 - r1*r1)/(2.0*r2*r3) );
		theta2= asin( (r2/r1)*sin(theta1) );
		theta3= asin( (r3/r1)*sin(theta1) );
	}
	/* l:m:s r2:r1:r3 */
	if ( (r1<r2) && (r2>r3) && (r1>r3) ) {
		theta2= acos( (r1*r1 + r3*r3 - r2*r2)/(2.0*r1*r3) );
		theta1= asin( (r1/r2)*sin(theta2) );
		theta3= asin( (r3/r2)*sin(theta2) );
	}
	/* l:m:s r2:r3:r1 */
	if ( (r1<r2) && (r2>r3) && (r1<r3) ) {
		theta2= acos( (r1*r1 + r3*r3 - r2*r2)/(2.0*r1*r3) );
		theta1= asin( (r1/r2)*sin(theta2) );
		theta3= asin( (r3/r2)*sin(theta2) );
	}
	/* l:m:s r3:r1:r2 */
	if ( (r1>r2) && (r2<r3) && (r1<r3) ) {
		theta3= acos( (r1*r1 + r2*r2 - r3*r3)/(2.0*r1*r2) );
		theta1= asin( (r1/r3)*sin(theta3) );
		theta2= asin( (r2/r3)*sin(theta3) );
	}
	/* l:m:s r3:r2:r1 */
	if ( (r1<r2) && (r2<r3) && (r1<r3) ) {
		theta3= acos( (r1*r1 + r2*r2 - r3*r3)/(2*r1*r2) );
		theta1= asin( (r1/r3)*sin(theta3) );
		theta2= asin( (r2/r3)*sin(theta3) );
	}

	/* components of distance to intersection of circles */
	dx_prime= r1*cos(theta2);
	dy_prime= r1*sin(theta2);

	/* calculate direction cosine rotation matrix */
	theta= atan2(y2-y1,x2-x1);
	dc[0][0] =  cos(theta);
	dc[0][1] =  sin(theta);
	dc[1][0] = -sin(theta);
	dc[1][1] =  cos(theta);

	/* already in rotated coordinate system, so transform back */
	dx = dc[0][0]*dx_prime - dc[1][0]*dy_prime;
	dy = dc[0][1]*dx_prime - dc[1][1]*dy_prime;
	xci[0] = x1+dx;
	yci[0] = y1+dy;

	dx = dc[0][0]*dx_prime + dc[1][0]*dy_prime;
	dy = dc[0][1]*dx_prime + dc[1][1]*dy_prime;
	xci[1] = x1+dx;
	yci[1] = y1+dy;

/*	printf("theta: if %f if %f\n", theta1, theta2, theta3, theta);
	printf("%f %f %f %f %f %f %f\n", x1, y1, r1, x2, y2, r2, r3);
	printf("%f %f %f %f %f\n", theta, dx_prime, dy_prime, dx, dy);
	printf("%f %f %f %f\n", xci[0], yci[0], xci[1], yci[1]);
	
	}*/
}	
	
int find_triple(int i, int jj, int nn, int *near, double xci, double yci, double t, int *numpoint)
{
	int 	j, k, kk, nt, m, mm;	
	double 	xi, yi, ri, ti;
	double 	xj, yj, rj, tj;
	double 	xk, yk, rk, tk;
	double 	xm, ym, rm, tm;
	double 	xn, yn, tn;
	double	rch;
	int 	ntrip;
	double 	xp[2], yp[2], tp[2];
	double 	bound_x, bound_y;
	int 	okay;

	bound_x= DUPL_FRAC*xmax;
	bound_y= DUPL_FRAC*ymax;

	xi = grain[i].x;	yi = grain[i].y;	ri = grain[i].r;	ti = grain[i].t;

	j  = near[jj];
	xj = grain[j].x;	yj = grain[j].y;	rj = grain[j].r;	tj = grain[j].t;

	/* loop for the kth neighbor */
	for (kk=0;kk<nn;kk++) {
		k = near[kk];
		if (k != j) {
			xk= grain[k].x;	yk= grain[k].y;	rk= grain[k].r; tk= grain[k].t;
			if ( meet(xk,yk,xci,yci,rk) ) {
			
				/*
				** Determine the number, position, and formation time of triple
				** points associated with the three grains. Alert user if no
				** triple points were found.
				*/

				ntrip = hyper(xi,yi,ti,xj,yj,tj,xk,yk,tk,xp,yp,tp);

				/*
				** Process each triple point found. In particular check to see that
				** the triple point is a valid triple point and not erroneous due to
				** edge effects. Also check to make sure that no other grain could
				** have grown to the triple point before the three grains we have
				** found to form the triple point.
				*/

				for (nt=0;nt<ntrip;nt++) {
					xn = xp[nt]; yn = yp[nt]; tn = tp[nt];
					if ( (xn < xmin-bound_x) || (yn < ymin-bound_y) ||
					     (xn > xmax+bound_x) || (yn > ymax+bound_y) ) {
						/* printf("Triple out of bounds.\n"); */ 
					}
					else if (tn > t) {
						/* printf("Triple point time too large 
							   %d %d %f %f.\n",i,j,t); */
					}
					else {
						/*
						** Check to see that no other grain has reached the triple
						** point first.
						*/
						okay = 1;
						for (mm=0;mm<nn;mm++) {
							m = near [mm];
							if ( (m != j) && (m != k) ){
								xm = grain[m].x; ym = grain[m].y;
								rm = grain[m].r; tm = grain[m].t;
								rch = g0 * (tn - tm);
								if (rch>0.0 && meet(xn,yn,xm,ym,rch) ) {
									okay = 0;
								}
							}
						}
	
						/* If okay=1 we have found a triple point. */
						if (okay == 1) {
							/* Record triple point information */
							point[*numpoint].x= xn;
							point[*numpoint].y= yn;
							point[*numpoint].t= tn;
							point[*numpoint].grain[0] = j;
							point[*numpoint].grain[1] = k;
							point[*numpoint].theta = atan2(yn-yi,xn-xi);
							point[*numpoint].ccw = 1;
							point[*numpoint].cw = 1;
							grain[i].point[grain[i].numpoint++]= (*numpoint)++;
						}/* if valid triple point */
					}/* if triple in bounds */
				}/* loop over triples found */
			}/* if jth and kth could meet */
		}/* if (k != j) */
	}/* loop for kth neighbor */
	return (1) ;
}		
		
/*
Subroutine to calculate the intersection points of two hyperbolae.
These represent the boundaries between grains 1 and 2 and
grains 1 and 3, when they grow circularly, nucleated at different
times.
*/
int hyper(double x1, double y1, double t1, double  x2, double  y2, double  t2, 
	  double x3, double y3, double t3, double *xp, double *yp, double *tt)
{
	double aa1,aa2,bb1,bb2,w1,w2,den,aaa,bbb,ccc,ddd,at,bt;
	double ct,f2,root,cc1,cc2,g2;

	g2  = g0*g0;
	aa1 = 2.*(x2-x1);
	aa2 = 2.*(x3-x1);
	bb1 = 2.*(y2-y1);
	bb2 = 2.*(y3-y1);

	w1 = x2*x2-x1*x1+y2*y2-y1*y1+g2*(t1*t1-t2*t2);
	w2 = x3*x3-x1*x1+y3*y3-y1*y1+g2*(t1*t1-t3*t3);
	den = bb2*aa1 - bb1*aa2;

	if (den == 0.0) {
		printf("Special case: Three points in straight line.\n");
		return(0);
	}

	if (aa1 == 0.0) {
		/* Special case of xl = x2 */
		aaa = w1/bb1;
		bbb = -2.*g2*(t1-t2)/bb1;
		ccc = (w2-bb2*aaa)/aa2;
		ddd = (-2.*g2*(t1-t3) - bb2*bbb)/aa2;
	}
	else {
		aaa = (w2*aa1-w1*aa2)/den;
		bbb = 2.*g2*((t3-t1)*aa1 + (t1-t2)*aa2)/den;
		ccc = (w1-bb1*aaa)/aa1;
		ddd = (-2.*g2*(t1-t2) - bb1*bbb)/aa1;
	}

	at = g2 - ddd*ddd - bbb*bbb;
	bt = 2.*(-t1*g2-ccc*ddd+x1*ddd-aaa*bbb+y1*bbb);
	ct = t1*t1*g2-(ccc-x1)*(ccc-x1)-(aaa-y1)*(aaa-y1);

	if (at == 0.0) {
		/* special case of a=O; bt*T+c*T=O */
		printf ("Special case, t is not quadratic. \n" );
		tt[0] = -ct/bt;
		cc1   = w1 - 2.*g2*(t1-t2)*tt[0];
		cc2   = w2 - 2.*g2*(t1-t3)*tt[0];

		if (aa1 == 0.0) xp[0] = (cc2-bb2*yp[0])/aa2;
		else 		xp[0] = (cc1-bb1*yp[0])/aa1;

		return(1);
	}

	f2 = bt*bt - 4.*at*ct;
	if (f2 < 0.0) {
		/* Handle the special case of imaginary roots. */
		printf("Speeial case: imaginary roots.\n");

		/* Return to program if error was from roundoff. */
		if (fabs(f2/(bt*bt)) > 0.001) return(0);
		else {
			root = 0.0;
			printf("Continued with root O. \n");
		}
	}
	else root = sqrt(f2);

	tt[0] = (-bt+root)/(2.*at);
	if (tt[0] < t3) return(0);

	cc1   = w1 - 2.*g2*(t1-t2)*tt[0];
	cc2   = w2 - 2.*g2*(t1-t3)*tt[0];
	yp[0] = (cc2*aa1-cc1*aa2)/den;

	if (aa1 == 0.0)
		xp[0] = (cc2 - bb2*yp[0])/aa2;
	else
		xp[0] = (cc1 - bb1*yp[0])/aa1;

	tt[1] = (-bt-root)/(2.*at);
	if (tt[1] < t3) return(1);

	cc1   = w1 - 2.*g2*(t1-t2)*tt[1];
	cc2   = w2 - 2.*g2*(t1-t3)*tt[1];
	yp[1] =  (cc2*aa1-cc1*aa2)/den;

	if (aa1 == 0.0)
		xp[1] = (cc2 - bb2*yp[1])/aa2;
	else
		xp[1] = (cc1 - bb1*yp[1])/aa1;
	return(2);
}

double tp_grain_angle(int g1, int pt, int g2)
{
	double x1, y1;
	double xj, yj;
	double x2, y2, x2_new, y2_new;
	double theta;
	double dc[2][2];

	x1 	 = grain[g1].x;		y1 = grain[g1].y;
	xj 	 = point[pt].x;		yj = point[pt].y;
	theta	 = atan2(yj-y1,xj-x1);
	dc[0][0] =  cos(theta);
	dc[0][1] =  sin(theta);
	dc[1][0] = -sin(theta);
	dc[1][1] =  cos(theta);

	x2	= grain[g2].x;
	y2	= grain[g2].y;
	x2_new	= dc[0][0]*(x2-xj) + dc[0][1]*(y2-yj);
	y2_new	= dc[1][0]*(x2-xj) + dc[1][1]*(y2-yj);

	theta	= atan2(y2_new,x2_new);

	return (theta);
}
