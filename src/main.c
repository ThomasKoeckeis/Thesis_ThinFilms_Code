/* main.c */

//#include <stdio.h>
//#include <math.h>
//#include <stdlib.h>
//#include "defines.h"
//#include "struct.h"
//#include "external.h"
//
//  int main_old (int argc, char *argv[])/*eigentlich main ohne 2*/
//{
//
//	int 	printed;
//
//	char	basename[50];		/* input file basename */
//	char	filename[100];		/* file name plus io directory */
//	char	text[80];		/* dummy text string */
//	FILE	*fin;			/* pointer to input file */
//	FILE	*fout;			/* pointer to output file */
//	FILE	*fps[10];		/* pointer to graphics file */
//	int	i_fps;
//	FILE	*fpsfinal;
//	FILE	*fpsthickness[2];
//
//	/* Parameter read from input file */
//	int	seed;			/* seed for random number generator, rand() */
//	int	number_grains;		/* desired number grains in final structure */
//	float	radius_real;		/* final grain size from input file */
//	float	dep_rate_real;		/* deposition rate in m/sec */
//	int	contact_angle;		/* island-substrate contact angle, degrees */
//	float	ca;			/* island-substrate contact angle, radians */
//	float	gamma;			/* surface stress of islands, J/m2 */
//	float	r_lock;			/* lock-down radius used in comp stress calc */
//	float	E;			/* Young's modulus of film materia.outl, Pa */
//	float	nu;			/* Poisson ratio of film material, unitless */
//	float	gamma_gb;		/* stress of grain boundaries J/m2 */
//	float	dg;			/* 2*gamma_surface - gamma_gb, J/m2 */
//	float	cluster_area_crit;	/* island sliding criteria */
//	float	temperature;		/* deposition temperature in K */
//	float	pre_exp;		/* Pre-exp factor to surface coble creep */
//	float	q_act;			/* Activation energy to surface coble creep */
//	float	graphics_coverage;	/* Substrate coverage at which to do graphics */
//	float	graphics_thickness;	/* Substrate thickness at which to do graphics */
//
//	/* Relate real dimensions are rate to simulation growth conditions */
//	float	geom;			/* geometric factor for spherical-cap island */
//	float	g0_i_ratio;		/* lateral growth rate to nucleation rate ratio */
//	float	dep_rate;		/* deposition rate in simulation */
//	float	size;			/* estimated average grain size in simulation */
//
//	/* film thickness during deposition */
//	float	thickness;
//
//	/* fractional substrate coverate */
//	double	coverage;
//	double  old_coverage; /*used to limit amount to print out */
//	double  stop_coverage;
//
//	/* Random time between nucleation attempts */
//	float	time_step;
//
//	/* Used to check if graphics dumped yet */
//	int	gdump;
//
//	/* Track number of intersections and contact points */
//	int	numpoint;
//	int	numcontact;
//	int	numcontact_old;
//
//	double	area;
//	double	stress_comp;
//	double	stress_tensile;
//	double	stress;
//	int	ng_inbound;
//
//	float	h_begin, h_final, h_step;
//
//	float	sim_radius;
//
//	int	i;
//	int	qi;
//
//	if (argc!=2) {
//		printf("USAGE:  ./grow results/BASENAME\n");
//		printf("where results/BASENAME.inp must have the correct format\n");
//		exit(0);
//	}
//
//	/* Read in name of input file in results directory */
//	sprintf(basename, "%s", argv[1]);
//	sprintf(filename, "%s.inp", basename);
//	
//	if ((fin=fopen(filename, "r"))==NULL) {
//		printf ("\nError opening input file %s\n", filename);
//		exit(0);
//	}
//	
//	/* Use in basename to create output file pointers */
//	sprintf(filename, "%s.out", basename);
//	
//	if ((fout=fopen(filename, "w"))==NULL) {
//		printf("\nError opening output file %s\n", filename);
//		exit(0);
//	}
//
//	for (i_fps=0;i_fps<10;i_fps++) {
//		sprintf(filename, "%s_%d.ps", basename, i_fps);
//	
//		if ((fps[i_fps]=fopen(filename, "w") )==NULL) {
//			printf("\nError opening graphics output file %s\n", filename);
//			exit(0);
//		}
//	}
//
//	sprintf(filename, "%s_final.ps", basename);
//	
//	if ((fpsfinal=fopen(filename, "w") )==NULL) {
//		printf("\nError opening graphics output file %s\n", filename);
//		exit(0);
//	}
//
//	for (i_fps=0;i_fps<2;i_fps++) {
//		sprintf(filename, "%s_thickness_%d.ps", basename, i_fps);
//	
//		if ((fpsthickness[i_fps]=fopen(filename, "w") )==NULL) {
//			printf("\nError opening graphics output file %s\n", filename);
//			exit(0);
//		}
//	}
//
//	/*
//	** Read simulation parameters from input file
//	*/
//
//	/* Random seed */
//	fgets(text,80,fin);
//	fscanf(fin, "%d", &seed);
//	fgets(text,80,fin);	fgets(text,80,fin);
//
//	/* Total number of grains once continuous film */
//	fscanf(fin, "%d", &number_grains);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//
//	/* Average grain radius once continuous (meters) */
//	fscanf(fin, "%f", &radius_real);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//
//	/* Deposition rate (meter/sec) */
//	fscanf(fin, "%f", &dep_rate_real);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//
//	/* Contact angle (degrees) */
//	fscanf(fin, "%d", &contact_angle);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//	/* Convert contact angle from degrees to radians */
//	ca = PI*contact_angle/180.0;
//
//	/* Surface stress (J/m2) */
//	fscanf(fin, "%f", &gamma);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//
//	/* Lock-down radius (meters) */
//	fscanf(fin, "%f", &r_lock);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//
//	/* Youngs' modulus of film material (Pa) */
//	fscanf(fin, "%f", &E);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//
//	/* Poisson ratio of film material (unitless) */
//	fscanf(fin, "%f", &nu);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//
//	/* Delta gamm = 2*gamma_surface - gamma_gb (J/m2) */
//	fscanf(fin, "%f", &gamma_gb);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//	dg = 2.0*gamma-gamma_gb;
///*	dg = gamma_gb;*/
//
//	/* Critical cluster area (m2) */
//	fscanf(fin, "%f", &cluster_area_crit);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//
//	/* Deposition temperature */
//	fscanf(fin, "%f", &temperature);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//
//	/* Pre-exponential to surface coble creep */
//	fscanf(fin, "%f", &pre_exp);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//
//	/* Activation energy for surface coble creep (eV) */
//	fscanf(fin, "%f", &q_act);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//
//	/* Graphics at specific substrate coverage
//	** If graphics are not desired, set above 1.0
//	*/
//	fscanf(fin, "%f", &graphics_coverage);
//	fgets(text, 80, fin);	fgets(text, 80, fin);
//
//	/* Graphics at specific substrate thickness
//	*/
//	fscanf(fin, "%f", &graphics_thickness);
//	fclose(fin);
//
//	printf("# INPUT FILE DATA:\n");
//	printf("# Seed:			%d\n", seed);
//	printf("# Number of grains:	%d\n", number_grains);
//	printf("# Real radius:		%e A\n", 1e10*radius_real);
//	printf("# Deposition rate:	%f nm/s\n", 1e9*dep_rate_real);
//	printf("# Contact angle:	%d degrees\n", contact_angle);
//	printf("# Surface stress:	%f J/m2\n", gamma);
//	printf("# Lock-down radius:	%f nm\n", 1e9*r_lock);
//	printf("# Young's modulus:	%f GPa\n", 1e-9*E);
//	printf("# Poisson ratio:	%f\n", nu);
//	printf("# Delta gamma:		%f\n", dg);
//	printf("# Crit. clust. area:	%e m2\n", cluster_area_crit);
//	printf("# Temperature:		%f K\n", temperature);
//	printf("# Pre exponential:	%e\n", pre_exp);
//	printf("# Activation energy:	%f\n", q_act);
//	printf("# Graphics coverage:	%f\n", graphics_coverage);
//	printf("# Graphics thickness:	%.2f nm\n\n", graphics_thickness);
//
//	/* Assign seed for pseudo-random number generator, rand(). */
//	srand((unsigned int) (seed));
//
//	/* Read C.V. Thompson, JMR 14, 3164 (1999) */
//	/* The final average grain radius can be estimated as:
//	** radius = O.602*g_i_ratio^(1/3)
//	** where g_i_ratio is the lateral growth rate/nucleation rate ratio.
//	**
//	** The scaling factor for the simulation relates simulation dimension
//	** to real dimensions:
//	** radius_real= scaling*radius_simulation
//	**
//	** From Thompson 1999, the lateral growth rate of an isolated island
//	** is the deposition rate time a geometric factor:
//	** dep_rate = g/geom;
//	**
//	** The real deposition rate supplied as an input to the simulation
//	** is therefore given by:
//	** dep_rate_real = (g_sim/geom)*scaling
//	**
//	** The choosen time step (shown later) gives a constant nucleation rate
//	** of one new nuclei per unit exposed area per unit time.
//	** I_sim = 1
//	**
//	** Combining the following three expressions from above:
//	** radius_real = scaling*radius_sim
//	** radius_sim  = O.602(g_sim/I_sim)^(1/3)
//	** g_sim       = dep_rate_real*geom/scaling
//	**
//	** radius_real = scaling*0.362*(dep_rate_real*geom/(I_sim*scaling))^(1/3)
//	** radius real = 0.602* (scaling^2*dep_rate_real*geom/(I_sim))^1/3
//	**
//	** Rearranging the above equation to solve for scaling;
//	** scaling = (radius_real/0.602)^(3/2)*(I_sim/dep_rate_real*geom)^(1/2)
//	*/
//
//	geom	   = 1.0;
//	geom 	  *= pow(sin(ca), 3.0);
//	geom 	  *= 1.0/pow(1.0 - cos(ca), 2.0);
//	geom 	  *= 1.0/(2.0 + cos(ca));
//	scaling	   = pow((radius_real/0.602), 1.5)*sqrt(1.0/(dep_rate_real*geom));
//	g0	   = dep_rate_real*geom/scaling;
//	g0_i_ratio = g0;
//
//	/* Calculate deposition rate for simulation */
//	dep_rate   = g0/geom;
//
//	/* Use average grain size to calculate size of substrate. */
//	size	= 2*0.602*pow(g0_i_ratio, (1.0/3.0));
//	xmin	= 0.0;	ymin= 0.0;
//
//	xmax	= sqrt(number_grains)*size;
//	ymax	= sqrt(number_grains)*size;
//	xsize	= xmax-xmin;
//	ysize	= ymax-ymin;
//
//	printf("# SIZE OF GRAINS: %f\n",      size);
//	printf("# xsize x ysize:  %f x %f\n", 1e9*scaling*xsize, 1e9*scaling*ysize);
//	printf("# dep rate        %f nm/s\n", dep_rate*scaling*1e9);
//	printf("# g0 		  %f\n",      g0);
//	printf("# geom 		  %f\n",      geom);
//
//	/* Allocate memory for structure arrays */
//	if (allocate_memory (number_grains, size) == 0)
//		exit(0);
//
//	/* Initialize to zero:
//	** number of grains (ng)
//	** time (t)
//	** film thickness
//	** fractional substrate coverage
//	** variable to check if graphics dumped yet
//	*/
//	ng	     = 0;
//	t	     = 0.0;
//	thickness    = 0.0;
//	coverage     = 0.0;
//	old_coverage = 0.0;
//	stop_coverage= 0.1;
//	numcontact   = 0;
//	gdump	     = 0;
//
//	/* Convert to scaling dimensions */
//	cluster_area_crit *= 1.0/(scaling*scaling);
//
//	fprintf(fout, "t\t");
//	fprintf(fout, "1e9*scaling*thickness\t");
//	fprintf(fout, "stress_comp*1e-6\t");
//	fprintf(fout, "stress_comp*scaling*thickness\t");
//	fprintf(fout, "0.5*stress_tensile*1e-6\t");
//	fprintf(fout, "stress_tensile*scaling*thickness\t");
//	fprintf(fout, "area\t");
//	fprintf(fout, "coverage\t");
//	fprintf(fout, "stress-thickness\t");
//	fprintf(fout, "stress\n");
//
//	printed=0; i_fps=0;
//	while ((stop_coverage > 0)&&(coverage!=stop_coverage)) {/*1.0001) {*/
//		if (((coverage - old_coverage) > 0.099999) && (coverage > 0.0)){
//			printf("coverage = %f, writing %s_%d.ps\n",coverage, basename, i_fps);
//			old_coverage = coverage;
//			ps_write(fps[i_fps++],	t, 
//						scaling,
//						thickness, 
//						coverage, 
//						stress_comp,
//						stress_tensile,
//						xsize,
//						ng);				
//		}
//
//		/* Do graphics dump at specified coverage */
//		/*if ((coverage > graphics_coverage) && !(gdump)) {
//			ps_write(fps, ng);
//			gdump= 1;
//		}*/
//		if ((1e9*scaling*thickness > graphics_thickness) && (gdump<1)) {
//			ps_write(fpsthickness[0], t, 
//						  scaling,
//						  thickness, 
//						  coverage, 
//						  stress_comp,
//						  stress_tensile,
//						  xsize,
//						  ng);
///*			ps_write(fpsthickness[0], t, 
//						  1e9*scaling*thickness, 
//						  coverage, 
//						  stress_comp*scaling*thickness,//*1e-6,
//						  stress_tensile*scaling*thickness,//*1e-6,
//						  1e9*scaling*xsize,
//						  ng);				*/
//			gdump = 1;
//		}
//		if ((1e9*scaling*thickness > 10.0*graphics_thickness) && (gdump<2)) {
//			ps_write(fpsthickness[1], t, 
//						  scaling,
//						  thickness, 
//						  coverage, 
//						  stress_comp,
//						  stress_tensile,
//						  xsize,
//						  ng);				
//			gdump = 2;
//		}
//
//		/*
//		** (1) Increment time step of simulation
//		** (2) Attempt to nucleate new island
//		*/
//
//		nucleate_island(&time_step) ;
//
//		/*
//		** At begining of new timestep:
//		** (1) Calculate radius of island (ignoring intersections)
//		** (2) Initialize number of intersection points to zero
//		** (3) Initialize number of boundary points to zero
//		*/
//		for (i=0; i<ng; i++) {
//			grain[i].r= g0*(t-grain[i].t);
//			grain[i].numpoint= 0;
//			grain[i].numbound= 0;
//		}
//
//		/* Determine intersection points between island */
//		island_intersection(&numpoint);
//
//		/* Organize intersection point in a CCW fashion around island */
//		organize_intersection_points();
//
//		/* Create boundary points between intersection points */
//		draw_island_boundary();
//
//		/* Calculate area of islands/grains */
//		for (i=0; i<ng; i++) {
//			grain[i].area= calc_grain_area(i);
//			if (grain[i].area < TOL1)
//				grain[i].area = TOL1;
//		}
//
//		/* Determine new contacts between island
//		** numcontact_added= numcontact_new - numcontact_old
//		*/
//		numcontact_old	= numcontact;
//		new_contacts(&numcontact);
//
//		/* Update old contact points */
//		for (i=0; i<numcontact_old; i++)
//			update_contacts(i) ;
//
//		/* calculate zipping stress for new contact points */
//		for (i=numcontact_old; i<numcontact; i++)
//			stress_contacts(i, ca, E, nu, dg, cluster_area_crit);
//
//		/* Calculate Laplace pressure for all island */
//		laplace_pressure(ca, gamma, r_lock);
//
//		/* Stres relaxation */
//		stress_relax_coble(time_step, temperature, pre_exp, q_act);
//
//		/* Area averaged stress */
//		stress_tensile	= 0;
//		stress_comp	= 0;
//		area		= 0;
//		ng_inbound	= 0;
//
//		for (i=0; i<ng; i++) {
//
//			if (inbound(grain[i].x, grain[i].y)) {
//				area 		+= grain[i].area;
//				stress_comp 	+= grain[i].area*grain[i].stress_comp;
//				stress_tensile	+= grain[i].area*grain[i].stress;
//				ng_inbound++;
//			}
//		}
//
//		stop_coverage	= (coverage>0)?coverage:stop_coverage;
//		coverage	= area/(xsize*ysize);
///*		if ((stop_coverage > 1) && (coverage == stop_coverage)) break; */
//		thickness 	+= dep_rate*time_step*coverage;
//		stress_comp 	*= 1.0/area;
//		stress_tensile 	*= 1.0/area;
//		fprintf(fout, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", t, 
//							1e9*scaling*thickness, 
//							stress_comp*1e-6, 
//							stress_comp*scaling*thickness, 
//							0.5*stress_tensile*1e-6, 
//							stress_tensile*scaling*thickness, 
//							1e18*scaling*scaling*xsize*ysize,
//							coverage,
//							(stress_comp+stress_tensile)*scaling*thickness,
//							(stress_comp+stress_tensile)*scaling//(stress_comp+stress_tensile)*1e-6
//							);
//
//		if ((printed==0) && (1e9*scaling*thickness>50)) {
//			printf("t\t");
//			printf("1e9*scaling*thickness\t");
//			printf("stress_comp*1e-6\t");
//			printf("stress_comp*scaling*thickness\t");
//			printf("0.5*stress_tensile*1e-6\t");
//			printf("stress_tensile*scaling*thickness\t");
//			printf("area\n");
//			printf("coverage\n");
//			printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", t, 
//						1e9*scaling*thickness,
//						stress_comp*1e-6,
//						stress_comp*scaling*thickness,
//						0.5*stress_tensile*1e-6,
//						stress_tensile*scaling*thickness,
//						1e18*scaling*scaling*area,
//						coverage);
//			printed = 1;
//		}
//			
//	}
//	ps_write(fps[9],	t, 
//				scaling,
//				thickness, 
//				coverage, 
//				stress_comp,
//				stress_tensile,
//				xsize,
//				ng);				
//
//
//	printf("FinalCoverage = %f\n",coverage);
//	h_step	  = thickness/10.0;//200.0;
//	h_begin	  = thickness + h_step;
//	h_final	  = 20.0*thickness;//200.0*thickness;
//	time_step = h_step/dep_rate;
//
//	printf("hstep=%f, hbegin=%f, hfinal=%f\n",h_step,h_begin,h_final);
//
//	for (thickness=h_begin; thickness<h_final; thickness+=h_step) {
//
//		if ((1e9*scaling*thickness > graphics_thickness) && (gdump<1)) {
//			ps_write(fpsthickness[0], t, 
//						  scaling,
//						  thickness, 
//						  coverage, 
//						  stress_comp,
//						  stress_tensile,
//						  xsize,
//						  ng);
//			gdump = 1;
//		}
//		if ((1e9*scaling*thickness > 10.0*graphics_thickness) && (gdump<2)) {
//			ps_write(fpsthickness[1], t, 
//						  scaling,
//						  thickness, 
//						  coverage, 
//						  stress_comp,
//						  stress_tensile,
//						  xsize,
//						  ng);
//			gdump = 2;
//		}						
//
//		t += time_step;
//		film_stress_relax_coble(time_step, temperature, pre_exp, q_act);
//
//		/* Area averaged stress */
//		stress_tensile	= 0;
//		stress_comp	= 0;
//		area		= 0;
//		ng_inbound	= 0;
//
//		for (i=0; i<ng; i++) {
//
//			if (inbound(grain[i].x, grain[i].y)) {
//				area 		+= grain[i].area;
//				stress_comp 	+= grain[i].area*grain[i].stress_comp;
//				stress_tensile 	+= grain[i].area*grain[i].stress;
//				ng_inbound++;
//			}
//		}
//
//		coverage	= area/(xsize*ysize);
//		stress_comp 	*= 1.0/area;
//		stress_tensile 	*= 1.0/area;
//
//		fprintf(fout, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", t, 
//							1e9*scaling*thickness,
//							stress_comp*1e-6, 
//							stress_comp*scaling*thickness, 
//							0.5*stress_tensile*1e-6, 
//							stress_tensile*scaling*thickness, 
//							1e18*scaling*scaling*area,
//							coverage,
//							(stress_comp+stress_tensile)*scaling*thickness,
//							(stress_comp+stress_tensile)*scaling//(stress_comp+stress_tensile)*1e-6
//							);
//	}
//	ps_write(fpsfinal,	t, 
//				scaling,
//				thickness, 
//				coverage, 
//				stress_comp,
//				stress_tensile,
//				xsize,
//				ng);				
//
//
//	return 0;
//}