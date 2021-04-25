
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "defines.h"
#include "struct.h"

extern grain_struct	*grain;
extern point_struct 	*point;
extern bound_struct 	*bound;
extern contact_struct 	*contact;

int allocate_memory(int number_grains, float size)
{
	/* Estimations of the the maximum size for array declaration */
	int maxgrains;
	int maxpoints;
	int maxbound;
	int maxbounds;
	int maxcontacts;

	/* Calculate total memory usage of allocated structures */
	int ram_mb;

	/* Estimate size of arrays */

	/* PEC mirror copies of the structure around central structure.
	** The 3X3 array of mirror copied structure gives 9*number grains
	** The factor of 2 in front is for a margain of safety
	*/
	maxgrains = 2*(9*number_grains);

	/* Each island impinges with on average 6 grains
	** The extra 100 if for a margin of safety
	*/
	maxpoints= 6*maxgrains + 100;

	/* Each island forms a grain boundary with on average 6 grains
	** The distance between boundary points is never less than TOOCLOSE
	** Each island must be alloted the max number of possible boundary points
	*/
	maxbound  = 6*size/TOOCLOSE;
	maxbounds = maxbound*maxgrains;

	/* Each island impinges with on average 6 grains
	** The extra 100 if for a margin of safety
	*/
	maxcontacts = 6*maxgrains + 100;

	printf("\n# maxgrain	%d\n", maxgrains);
	printf("# maxpoint 	%d\n", maxpoints);
	printf("# maxbound 	%d\n", maxbound);
	printf("# maxbounds 	%d\n", maxbounds);
	printf("# maxcontacts 	%d\n\n", maxcontacts);

	ram_mb	= 0;
	ram_mb += (int) maxgrains*sizeof(grain_struct);
	printf("# Using %.0f MB of RAM\n", ram_mb/1.0e6);
	ram_mb += (int) maxpoints*sizeof(point_struct);
	printf("# Using %.0f MB of RAM\n", ram_mb/1.0e6);
	ram_mb += (int) maxbounds*sizeof(bound_struct);
	printf("# Using %.0f MB of RAM\n", ram_mb/1.0e6);
	ram_mb += (int) maxcontacts*sizeof(contact_struct);
	printf("# Using %.0f MB of RAM\n", ram_mb/1.0e6);
	ram_mb *= 1.0e-6;

	/* Allocate memory */
	grain= (grain_struct *) malloc(maxgrains*sizeof(grain_struct));
	point= (point_struct *) malloc(maxpoints*sizeof(point_struct));
	bound= (bound_struct *) malloc(maxbounds*sizeof(bound_struct));
	contact= (contact_struct *) malloc(maxcontacts*sizeof(contact_struct));

	if (ram_mb < RAM_MAX)
		return(1) ;
	else {
		printf("%d MB of RAM allocated greater than %d MB limit.\n", ram_mb, RAM_MAX);
		return(0);
	}
}
