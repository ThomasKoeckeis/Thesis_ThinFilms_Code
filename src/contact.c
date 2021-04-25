/* contact.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "defines.h"
#include "struct.h"
#include "length_height.h"

extern grain_struct *grain;
extern point_struct *point;
extern bound_struct *bound;
extern contact_struct *contact;

extern float t;
extern int ng;
extern float scaling;
extern int long numbound;

/* auxiliary.c */
extern int inbound (double, double);

/* contact.c */

/*float gb_length(int);
float zip_height(int, int, float, float, float);*/
float height(int, int);
float cluster_area(int);

void new_contacts (int *numcontact)
{
	int pt_start, pt_stop;
	int ct;
	int found;
	int i, j, k, m;

	/* Determine if new contact points */
	for (i=0; i<ng; i++) {

		/* Only process those grains that are inbounds */
		if (inbound(grain[i].x, grain[i].y)) {

			/* Search intersection points to see if new contact points */
			for (j=0; j<grain[i].numpoint; j++) {
				pt_start = grain[i].point[j];
				if (j<grain[i].numpoint-1)
					pt_stop= grain[i].point[j+1];
				else
					pt_stop= grain[i].point[0];
				k= point[pt_start].grain[0];

				if ( (point[pt_start].ccw==1) && (point[pt_stop].cw==1) ) {
					found= 0;
					for (m=0; m<grain[i].numcontact; m++) {
						if (contact[grain[i].contact[m]].grain_contact == k) {
							found= 1;
							ct= grain[i].contact[m];
						}
					}
					if (!found) {
						contact[*numcontact].t		   = t;
						contact[*numcontact].grain	   = i;
						contact[*numcontact].grain_contact = k;
						contact[*numcontact].point[0] 	   = pt_start;
						contact[*numcontact].point[1] 	   = pt_stop;

						grain[i].contact[grain[i].numcontact++] = (*numcontact)++;
					}
					if (found) {
						contact[ct].point[0]= pt_start;
						contact[ct].point[1]= pt_stop;
					}
				}
			}
		}
	}
}

float gb_length(int ct)
{
	int grn;
	double x1, y1;
	double x2, y2;
	double xb, yb;
	float distance;
	float length;
	int start, stop;
	int i;

	grn = contact[ct].grain;
	x1  = point[contact[ct].point[0]].x;
	y1  = point[contact[ct].point[0]].y;
	x2  = point[contact[ct].point[1]].x;
	y2  = point[contact[ct].point[1]].y;

	start	= stop	= -1;
	for (i=0; i<grain[grn].numbound; i++) {
		xb	 = bound[grain[grn].bound[i]].x1;
		yb	 = bound[grain[grn].bound[i]].y1;
		distance = sqrt( (xb-x1)*(xb-x1)+(yb-y1)*(yb-y1));
		if (distance < TOL1)
			start= i;
		distance = sqrt( (xb-x2)*(xb-x2)+(yb-y2)*(yb-y2));
		if (distance < TOL1)
			stop= i;
	}

	length= 0;
	if (start < stop)
		for (i=start; i<stop; i++) {
			x1 = bound[grain[grn].bound[i]].x1;
			y1 = bound[grain[grn].bound[i]].y1;
			x2 = bound[grain[grn].bound[i]].x2;
			y2 = bound[grain[grn].bound[i]].y2;
			length += sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		}
	else {
		for (i=start; i<grain[grn].numbound; i++) {
			x1 = bound[grain[grn].bound[i]].x1;
			y1 = bound[grain[grn].bound[i]].y1;
			x2 = bound[grain[grn].bound[i]].x2;
			y2 = bound[grain[grn].bound[i]].y2;
			length += sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		}
		for (i=0; i<stop; i++) {
			x1 = bound[grain[grn].bound[i]].x1;
			y1 = bound[grain[grn].bound[i]].y1;
			x2 = bound[grain[grn].bound[i]].x2;
			y2 = bound[grain[grn].bound[i]].y2;
			length += sqrt( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		}
	}
	return(length);
}

float zip_height(int g1, int g2, float ca, float E, float dg)
{
	float r0, r1, r2;
	float z0_r1, z0_r2;
	float y0_r1, y0_r2;
	float z0_r1_min;
	float min;
	float se1_unit, se2_unit, se;
	float egb;
	float etot;
	int i, steps;

	if ( (grain[g1].area < TOL1) || (grain[g2].area < TOL1) )
		return(0);

	r0= 100e-10;
	r1= scaling*sqrt(grain[g1].area/PI);
	r2= scaling*sqrt(grain[g2].area/PI);

	steps= 1000;
	min= 0;
	z0_r1_min= 1.0;

	for (i=1; i<steps; i++) {
		z0_r1	 = (1.0/steps)*i;
		y0_r1	 = 1.0 - sqrt(1.0/(sin(ca)*sin(ca)) - pow(z0_r1+cos(ca)/sin(ca), 2.0));
		z0_r2	 = z0_r1*r1/r2;
		y0_r2	 = 1.0 - sqrt(1.0/(sin(ca)*sin(ca)) - pow(z0_r2+cos(ca)/sin(ca), 2.0));
		se1_unit = E*0.2576*pow(y0_r1, 2.0);
		se2_unit = E*0.2576*pow(y0_r2, 2.0);
		se	 = se1_unit*r1*r1 + se2_unit*r2*r2;
		egb	 = dg*z0_r1*r1;
		etot	 = se-egb;
		/* printf("%f %e %e %e\n", zO_rl, se, egb, etot); */
		if (etot<min) {
			min	  = etot;
			z0_r1_min = z0_r1;
		}
	}
	if (z0_r1_min*r1/scaling>0.1) {
		/**/
		printf("\nz0_r1_min=%f, r1=%f, scaling=%f\n",z0_r1_min, r1, scaling);
		/**/
	}
	return(z0_r1_min*r1/scaling);
}

void update_contacts (int ct)
{
	int i, k;

	i = contact[ct].grain;
	k = contact[ct].grain_contact;

	contact[ct].length = gb_length(ct);
	contact[ct].height = height(i,k);
	contact[ct].radius = sqrt(grain[i].area/PI);
}

float height (int i, int j)
{
	double x1, y1, r1;
	double x2, y2, r2;
	double r3;
	double theta1, theta2, theta3;
	double h;

	x1 = grain[i].x;	y1 = grain[i].y;	r1 = grain[i].r;
	x2 = grain[j].x;	y2 = grain[j].y;	r2 = grain[j].r;

	r3= sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );

	/* l:m:s r1:r2:r3 */
	if ( (r1>r2) && (r2>r3) && (r1>r3) ) {
		theta1= acos( (r2*r2 + r3*r3 - r1*r1)/(2*r2*r3) );
		theta2= asin( (r2/r1)*sin(theta1) );
		theta3= asin( (r3/r1)*sin(theta1) );
	}
	/* l:m:s r1:r3:r2 */
	if ( (r1>r2) && (r2<r3) && (r1>r3) ) {
		theta1 = acos( (r2*r2 + r3*r3 - r1*r1)/(2*r2*r3) );
		theta2 = asin( (r2/r1)*sin(theta1) );
		theta3 = asin( (r3/r1)*sin(theta1) );
	}
	/* l:m:s r2:r1:r3 */
	if ( (r1<r2) && (r2>r3) && (r1>r3) ) {
		theta2 = acos( (r1*r1 + r3*r3 - r2*r2)/(2*r1*r3) );
		theta1 = asin( (r1/r2)*sin(theta2) );
		theta3 = asin( (r3/r2)*sin(theta2) );
	}
	/* l:m:s r2:r3:r1 */
	if ( (r1<r2) && (r2>r3) && (r1<r3) ) {
		theta2 = acos( (r1*r1 + r3*r3 - r2*r2)/(2*r1*r3) );
		theta1 = asin( (r1/r2)*sin(theta2) );
		theta3 = asin( (r3/r2)*sin(theta2) );
	}
	/* l:m:s r3:r1:r2 */
	if ( (r1>r2) && (r2<r3) && (r1<r3) ) {
		theta3 = acos( (r1*r1 + r2*r2 - r3*r3)/(2*r1*r2) );
		theta1 = asin( (r1/r3)*sin(theta3) );
		theta2 = asin( (r2/r3)*sin(theta3) );
	}
	/* l:m:s r3:r2:r1 */
	if ( (r1<r2) && (r2<r3) && (r1<r3) ) {
		theta3 = acos( (r1*r1 + r2*r2 - r3*r3)/(2*r1*r2) );
		theta1 = asin( (r1/r3)*sin(theta3) );
		theta2 = asin( (r2/r3)*sin(theta3) );
	}
	/* components of distance to intersection of circles */
	h = fabs(r1*cos(theta2));

	return(h);
}
