/* graphics.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "defines.h"
#include "struct.h"

#define MAXLINE 	400
#define UNIT 		144
#define TOOCLOSE2 	0.5
#define XM 		4
#define YM 		4
#define LINEWIDTH 	0.5
#define MARGIN 		0.5
#define HEADER 		"%!%!PS-Adobe"
#define FONT 		"Courier"
#define SIZEFONT 	16

extern grain_struct *grain;
extern int long numbound;

/* auxiliary.c */
extern int inbound (double, double);

extern grain_struct *grain;
extern bound_struct *bound;

extern double xmin, xmax, ymin, ymax;
extern double xsize, ysize, xcursor, ycursor;


void ps_write (	FILE *f, 
		float t, 
		float scaling,
		float thickness_unscaled,
		double coverage,
		double stress_c,
		double stress_t,
		double x_size_unscaled,
		int ng)
{
	float thickness;
	double stress_comp;
	double stress_tensile;
	double x_size;

	double x, y;
	int    xmpts, ympts;
	int    margx, margy;
	float  ax, ay;
	float  bx, by;
	float  tooclx, toocly;
	float  xpscurr1, ypscurr1;
	float  xsized2, ysized2;
	int    i, j, counter;

	/* Area plot */
	double maxarea, minarea, steparea;
	/* Stress plot */
	double maxstress, minstress, stepstress;
	double ylegstep = 0;
	int    imin = 0;
	int imax = 0;

	thickness      = 1e9*scaling*thickness_unscaled;
	stress_comp    = stress_c*scaling*thickness_unscaled;
	stress_tensile = stress_t*scaling*thickness_unscaled;
	x_size	       = 1e9*scaling*x_size_unscaled;

	printf ("...ps_write in function\n");

	/* initialization of the postscript coordinates of the page */
	xmpts 	= XM*UNIT;
	ympts 	= YM*UNIT;
	margx 	= (int) (MARGIN*UNIT);
	margy 	= (int) (MARGIN*UNIT);
	xsized2 = xsize/2.0;
	ysized2 = ysize/2.0;

	/* defines the coordinates transformation */
	ax 	= (xmpts-2*margx)/xmax;
	bx 	= margx;
	ay 	= (ympts-2*margy)/ymax;
	by 	= margy;
	tooclx 	= TOOCLOSE2/ax;
	toocly 	= TOOCLOSE2/ay;

	/* header of the postscript file */
	fprintf(f, "%s\n", HEADER);
	fprintf(f, "initgraphics\n");
	fprintf(f, "newpath \n");
	fprintf(f, "0 setlinecap\n");
	fprintf(f, "%.2f setlinewidth\n",LINEWIDTH);
/*	fprintf(f, "0.8 setgray\n");	*/

	maxarea   = -1000; minarea   = 1000;
	maxstress = -1000; minstress = 1000;

	for (i=0; i<ng; i++) {
		for (j=0; j<grain[i].numbound; j++) {
			x = bound[grain[i].bound[j]].x1;
			y = bound[grain[i].bound[j]].y1;
			if (inbound(x,y)) {
				/* Area plot */
				if (grain[i].area>maxarea) {
					maxarea=grain[i].area;
					imax = i;
				}
				if (grain[i].area<minarea) {
					minarea=grain[i].area;
					imin = i;
				}
				/* Stress plot */
				if (grain[i].stress*scaling>maxstress) {
					maxstress=grain[i].stress*scaling;
				}
				if (grain[i].stress*scaling<minstress) {
					minstress=grain[i].stress*scaling;
				}

			}
		}
	}
	steparea   = (maxarea  - minarea)   / 7.0;
	stepstress = (maxstress- minstress) / 7.0;

	for (i=0; i<ng; i++) {
		counter = MAXLINE;

		if (grain[i].stress*scaling > maxstress-1.0*stepstress) {
			fprintf (f, "1 0 1 setrgbcolor\n");
		} else if (grain[i].stress*scaling > maxstress-2.0*stepstress) {
			fprintf (f, "1 0 0 setrgbcolor\n");
		} else if (grain[i].stress*scaling > maxstress-3.0*stepstress) {
			fprintf (f, "1 0.5 0 setrgbcolor\n");
		} else if (grain[i].stress*scaling > maxstress-4.0*stepstress) {
			fprintf (f, "1 1 0 setrgbcolor\n");
		} else if (grain[i].stress*scaling > maxstress-5.0*stepstress) {
			fprintf (f, "0 1 0 setrgbcolor\n");
		} else if (grain[i].stress*scaling > maxstress-6.0*stepstress) {
			fprintf (f, "0 1 1 setrgbcolor\n");
		} else {
			fprintf (f, "0 0 1 setrgbcolor\n");
		}
		for (j=0; j<(grain[i].numbound); j++) {

			x = bound[grain[i].bound[j]].x1;
			y = bound[grain[i].bound[j]].y1;
			if (inbound(x,y)) {
				xpscurr1 = ax*x+bx;
				ypscurr1 = ay*y+by;
				counter++;

				if (counter > MAXLINE) {
					counter= 0;
/*					fprintf (f, "stroke\n");*/
					fprintf (f, "%.1f %.1f moveto\n", xpscurr1, ypscurr1);
				} else {
					fprintf (f, "%.1f %.1f lineto\n", xpscurr1, ypscurr1);
				}

				if (j==grain[i].numbound-1) {
					x = bound[grain[i].bound[j]].x2;
					y = bound[grain[i].bound[j]].y2;
					xpscurr1 = ax*x+bx;
					ypscurr1 = ay*y+by;
					fprintf (f , "%.1f %.1f lineto\n", xpscurr1, ypscurr1);
				}
			}
		}


/*		if (grain[i].area > maxarea-1.0*steparea) {
			fprintf (f, "1 0 1 setrgbcolor\n");
		} else if (grain[i].area > maxarea-2.0*steparea) {
			fprintf (f, "1 0 0 setrgbcolor\n");
		} else if (grain[i].area > maxarea-3.0*steparea) {
			fprintf (f, "1 0.5 0 setrgbcolor\n");
		} else if (grain[i].area > maxarea-4.0*steparea) {
			fprintf (f, "1 1 0 setrgbcolor\n");
		} else if (grain[i].area > maxarea-5.0*steparea) {
			fprintf (f, "0 1 0 setrgbcolor\n");
		} else if (grain[i].area > maxarea-6.0*steparea) {
			fprintf (f, "0 1 1 setrgbcolor\n");
		} else {
			fprintf (f, "0 0 1 setrgbcolor\n");
		} */

		fprintf (f, "fill\n");
		fprintf (f, "stroke\n");

	}
	printf ("Minarea = %f um^2, Minradius = %f um\n", minarea, grain[imin].r);
	printf ("Maxarea = %f um^2, Maxradius = %f um\n", maxarea, grain[imax].r);
	printf ("Minstress = %f, Maxstress = %f\n", minstress, maxstress);


	fprintf (f, "0.0 setgray\n");
	for (i=0; i<ng; i++) {
		counter = MAXLINE;

		for (j=0; j<grain[i].numbound; j++) {

			x = bound[grain[i].bound[j]].x1;
			y = bound[grain[i].bound[j]].y1;

			if (inbound(x,y)) {
				xpscurr1 = ax*x+bx;
				ypscurr1 = ay*y+by;
				counter++;

				if (counter > MAXLINE) {
					counter= 0;

					fprintf (f, "stroke\n");
					fprintf (f, "%.1f %.1f moveto\n", xpscurr1, ypscurr1);
				}
				else
					fprintf (f, "%.1f %.1f lineto\n", xpscurr1, ypscurr1);

				if (j==grain[i].numbound-1) {
					x = bound[grain[i].bound[j]].x2;
					y = bound[grain[i].bound[j]].y2;
					xpscurr1 = ax*x+bx;
					ypscurr1 = ay*y+by;
					fprintf (f , "%.1f %.1f lineto\n", xpscurr1, ypscurr1);
				}
			}
		}
		fprintf (f, "stroke\n");
	}
	/* draws the frames */
	fprintf(f, "0.0 setgray\n");
	fprintf(f, "1 setlinewidth\n");
	fprintf(f, "%d %d moveto\n", margx,margy);
	fprintf(f, "%d %d lineto\n", xmpts-margx,margy);
	fprintf(f, "%d %d lineto\n", xmpts-margx,ympts-margy);
	fprintf(f, "%d %d lineto\n", margx,ympts-margy);
	fprintf(f, "%d %d lineto\n", margx,margy);
	fprintf(f, "stroke\n");

	/* Prints the stress legend boxes */
	ylegstep = (ympts-2*margy-60)/7.0;
	fprintf(f, "0.5 setlinewidth\n");
	for (i=1; i<8; i++) {
		fprintf(f, "%d %f moveto\n", xmpts-margx+30, (i-1)* ylegstep+margy+20);
		fprintf(f, "%d %f lineto\n", xmpts-margx+40, (i-1)* ylegstep+margy+20);
		fprintf(f, "%d %f lineto\n", xmpts-margx+40, (i)  * ylegstep+margy+20);
		fprintf(f, "%d %f lineto\n", xmpts-margx+30, (i)  * ylegstep+margy+20);
		fprintf(f, "%d %f lineto\n", xmpts-margx+30, (i-1)* ylegstep+margy+20);
		switch(i) {
			case 1: fprintf(f, "0 0 1 setrgbcolor\n"); break;
			case 2: fprintf(f, "0 1 1 setrgbcolor\n"); break;
			case 3: fprintf(f, "0 1 0 setrgbcolor\n"); break;
			case 4: fprintf(f, "1 1 0 setrgbcolor\n"); break;
			case 5: fprintf(f, "1 0.5 0 setrgbcolor\n"); break;
			case 6: fprintf(f, "1 0 0 setrgbcolor\n"); break;
			case 7: fprintf(f, "1 0 1 setrgbcolor\n"); break;
			default: fprintf(f, "0 0 0 setrgbcolor\n"); break;
		}
		fprintf(f, "fill\n");
		fprintf(f, "stroke\n");
		fprintf(f, "0.0 setgray\n");
		fprintf(f, "%d %f moveto\n", xmpts-margx+30, (i-1)* ylegstep+margy+20);
		fprintf(f, "%d %f lineto\n", xmpts-margx+40, (i-1)* ylegstep+margy+20);
		fprintf(f, "%d %f lineto\n", xmpts-margx+40, (i)  * ylegstep+margy+20);
		fprintf(f, "%d %f lineto\n", xmpts-margx+30, (i)  * ylegstep+margy+20);
		fprintf(f, "%d %f lineto\n", xmpts-margx+30, (i-1)* ylegstep+margy+20);
		fprintf(f, "stroke\n");
	}

	/* Prints the stress legend text */
	fprintf(f, "/Times-Roman findfont\n");
	fprintf(f, "20 scalefont\n");
	fprintf(f, "setfont\n");
	fprintf(f, "newpath\n");
	fprintf(f, "%d %d moveto\n", xmpts-margx+3, ympts-margy-12);
	fprintf(f, "(Stress) show\n");/*"(Area (um  )) show\n");*/
	fprintf(f, "%d %d moveto\n", xmpts-margx+25, ympts-margy-32);
	fprintf(f, "(%.2f) show\n", maxstress);/*maxarea);*/
	fprintf(f, "%d %d moveto\n", xmpts-margx+30, margy);
	fprintf(f, "(%.2f) show\n", minstress);/*minarea);*/
/*	fprintf(f, "/Times-Roman findfont\n");
	fprintf(f, "15 scalefont\n");
	fprintf(f, "setfont\n");
	fprintf(f, "newpath\n");
	fprintf(f, "%d %d moveto\n", xmpts-margx+82, ympts-margy-8);
	fprintf(f, "(2) show\n");*/

	/* Prints the time, thickness, coverage statistics */
	fprintf(f,  "/Times-Roman findfont\n");
	fprintf(f,  "20 scalefont\n");
	fprintf(f,  "setfont\n");
	fprintf(f,  "newpath\n");

	fprintf(f,  "%d %d moveto\n", margx, ympts-margy+130);
	fprintf(f,  "(Time = %.2f s) show\n", t);
	fprintf(f,  "%d %d moveto\n", margx, ympts-margy+110);
	fprintf(f,  "(Thickness = %.2f nm) show\n", thickness);
	fprintf(f,  "%d %d moveto\n", margx, ympts-margy+90);
	fprintf(f,  "(Coverage = %.2f ) show\n", coverage);

	fprintf(f,  "%d %d moveto\n", margx, ympts-margy+50);
	fprintf(f,  "(Total area = %.0f nm2) show\n", x_size*x_size);
/*	fprintf(f,  "%d %d moveto\n", margx, ympts-margy+30);
	fprintf(f,  "(Compressive stress = %.2f N/m) show\n", stress_comp);*/
	fprintf(f,  "%d %d moveto\n", margx, ympts-margy+10);
	fprintf(f,  "(Tensile stress = %.2f N/m) show\n", stress_tensile);

	fprintf(f,  "%d %d moveto\n", 3*margx, margy-20);
	fprintf(f,  "(%.2f nm) show\n", x_size);

	/* prints the simulation parameters */
	/* 
	fprintf (f, "/%5 findfont %d scalefont setfont\n",FONT,SIZEFONT);
	fprintf (f, "%d %.1f moveto\n",margx, ympts+O.5*UNIT);
	fprintf (f, "(G_I_ratio: %f Delta: %.2f Seed: %d) show\n", g_i_ratio, delta, seed) ;
	fprintf (f, "%d %.1f moveto\n",margx, ympts+O.25*UNIT);
	fprintf (f, "(Average grain size: %f) show\n", d_avg);
	*/
	fprintf(f, "showpage\n");
	fclose(f);
}
