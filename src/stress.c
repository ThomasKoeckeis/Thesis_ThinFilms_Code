/* stress.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "defines.h"
#include "struct.h"
#include "length_height.h"

extern grain_struct *grain;
extern contact_struct *contact;

extern float scaling;
extern float g0;
extern float t;
extern int   ng;

/* auxiliary.c */
extern int inbound(double, double);

/* contact.c */
extern float height(int, int);

/* stress.c */
float cluster_area(int);

void stress_contacts(int ct, float ca, float E, float nu, float dg, float cluster_area_crit)
{
	float z0_r, y0_r;
	float stress_tensile;
	int i, k;

	i= contact[ct].grain;
	k= contact[ct].grain_contact;

	contact[ct].z0	   = zip_height(i,k,ca,E,dg);
	contact[ct].length = gb_length(ct);
	contact[ct].height = height(i,k);
	contact[ct].radius = sqrt(grain[i].area/PI);

	z0_r 		= contact[ct].z0/contact[ct].radius;
	/**/
/*	if ((1.0/(sin(ca)*sin(ca)) - pow(z0_r + cos(ca)/sin(ca), 2.0))<0.0) {
		printf("\n\nz0=%f, z0_r=%f\n\n", contact[ct].z0, z0_r);
	} else {
		printf("\n\nz0=%f", contact[ct].z0);
	}
*/	/**/
	y0_r 		= 1.0 - sqrt(1.0/(sin(ca)*sin(ca)) - pow(z0_r + cos(ca)/sin(ca), 2.0));
	stress_tensile	= 0.5*(E/(1.0-nu*nu))*pow(y0_r, 1.3892);

	if ( (cluster_area(i) > cluster_area_crit) &&
	     (cluster_area(k) > cluster_area_crit) )
		grain[i].stress += stress_tensile;
	else
		grain[i].stress += 0.0;
}

void laplace_pressure(float ca, float gamma, float r_lock)
{
	float comp_prefactor;
	float radius;
	int i;

	/* Laplace Pressure
	** Since radius changing each timestep, must recalculate stress for
	** each grain.
	**
	** Two models for the origin of compressive stress resulting from
	** the Laplace pressure. The first says simply that the compressive
	** stress in an island is given by:
	** stress= -gamma/radius
	** The other model assumes that an island can easily slide when very
	** small, but above a critical size the island attaches to the
	** substrate with a certain lattice parameter. As the island grows
	** larger, the equilibrium lattice parameter dictated by the Lapalce
	** pressure increases. However the newly added material assumes the
	** lattice parameter of the existing material. Consequently, the
	** material is under a compressive stres since the lattice parameter
	** is smaller than the equilibrium lattice parameter:
	** stress= gamma(1/radius - 1/radius_lock)
	**
	** For gamma= 1.5 and radius lock= 50, both models give similar results.
	*/

	comp_prefactor	= 4.0*gamma*sin(ca);
	comp_prefactor *= 1.0/(1.0-cos(ca));
	comp_prefactor *= 1.0/(2.0+cos(ca));

	for (i=0; i<ng; i++) {

		radius= scaling*sqrt(grain[i].area/PI);

		if (radius < 5e-10)
			radius= 5e-10;
		grain[i].stress_comp = -0.5*comp_prefactor/radius;

/*		grain[i].stress_comp = 0.0;*/
		/*if (radius > r_lock)
			grain[i].stress_comp = comp_prefactor*(1/radius - 1/r_lock);
		   else
			grain[i].stress_comp = 0;
		*/
	}
}

void stress_relax_coble(float time_step, float temperature, float pre_exp, float q_act)
{
	float ht;
	float length;
	int   ct;
	float stress_relax;
	int   i,j;
	
	for (i=0; i<ng; i++) {
		if (inbound(grain[i].x, grain[i].y)) {

			ht= 0;
			length= 0;
			for (j=0; j<grain[i].numcontact; j++) {
				ct= grain[i].contact[j];
				if (contact[ct].height > contact[ct].z0)
					ht += contact[ct].height;
				else
					ht += contact[ct].z0;
	
				if (contact[ct].length > 2*contact[ct].z0)
					length += contact[ct].length;
				else
					length += 2*contact[ct].z0;
			}

			if (grain[i].numcontact > 0) {
				ht *= 1.0/grain[i].numcontact;

				stress_relax  = 1.0/pow(ht,3.0);
				stress_relax *= 1.0/pow(scaling,3.0);
				stress_relax *= pre_exp;
				stress_relax *= exp(-q_act*1.602e-19/(1.381e-23*temperature));
				stress_relax *= grain[i].stress;
				stress_relax *= time_step;

				if (stress_relax > grain[i].stress)
					grain[i].stress= 0;
				else
					grain[i].stress -= stress_relax;
			}
		}
	}
}

void film_stress_relax_coble(float time_step, float temperature, float pre_exp, float q_act)
{
	float ht;
	float length;
	int   ct;
	float radius;
	float stress_relax;
	int   i, j, k;

	for (i=0; i<ng; i++) {
		grain[i].r= g0*(t-grain[i].t);

		if (inbound(grain[i].x, grain[i].y)) {

			ht	= 0;
			length  = 0;
			for (j=0; j<grain[i].numcontact; j++) {

				ct = grain[i].contact[j];
				k  = contact[ct].grain_contact;

				contact[ct].height = height(i,k);
				contact[ct].radius = sqrt(grain[i].area/PI);

				radius= 0.5*(contact[ct].radius + sqrt(grain[k].area/PI));

				if (contact[ct].height > contact[ct].z0)
					ht += contact[ct].height;
				else
					ht += contact[ct].z0;
				length += contact[ct].length;
			}

			if (grain[i].numcontact > 0) {
				ht *= 1.0/grain[i].numcontact;

				stress_relax  = 1.0/pow(ht,3.0);
				stress_relax *= 1.0/pow(scaling,3.0);
				stress_relax *= pre_exp;
				stress_relax *= exp(-q_act*1.602e-19/(1.381e-23*temperature));
				stress_relax *= grain[i].stress;
				stress_relax *= time_step;

				if (stress_relax > grain[i].stress)
					grain[i].stress = 0;
				else
					grain[i].stress -= stress_relax;
			}
		}
	}
}

float cluster_area (int grn)
{
	int	*cluster_list;
	int	num_cl, num_cl_last;
	float	area;
	int	ct;
	int	found1, found2;
	int	i, j, k;

	if (grain[grn].numcontact==0)
		return (grain[grn].area);

	cluster_list = (int *) malloc(ng*sizeof(int));

	num_cl = 0;
	cluster_list[num_cl++] = grn;

	num_cl_last = 0;
	while (num_cl_last < num_cl) {
		num_cl_last = num_cl;
		for (i=0; i<num_cl; i++) {
			grn = cluster_list[i];
			for (j=0; j<grain[grn].numcontact; j++) {
				ct     = grain[grn].contact[j];

				found1 = found2 = 0;
				for (k=0; k<num_cl; k++) {
					if ( (contact[ct].grain == cluster_list[k]) )
						found1 = 1;
					if ( (contact[ct].grain_contact == cluster_list[k]) )
						found2= 1;
				}

				if (!found1)
					cluster_list[num_cl++] = contact[ct].grain;
				if (!found2)
					cluster_list[num_cl++] = contact[ct].grain_contact;
			}
		}
	}

	area = 0.0;
	for (i=0; i<num_cl; i++)
		area += grain[cluster_list[i]].area;

	free(cluster_list);
	return(area);
}
