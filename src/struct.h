/* struct.h */
#include "defines.h"

/* Structure of island/grain information */
typedef struct
{ double x, y;			/* x,y nucleation coordinate */
  double r;			/* island radius */
  double t;			/* time of nucleation */
  int	 numpoint;		/* number of intersection points */
  int	 point[MAXSIDES];	/* intersection point indexed to point struct */
  int	 numbound;		/* number of boundary points */
  int	 bound[200*MAXSIDES];	/* boundary points indexed to bound struct */
  int	 numcontact;		/* number of contact points */
  int	 contact[2*MAXSIDES];	/* contact points indexed to contact struct */
  double area;			/* in-plane area of island/grain */
  float	 stress;		/* calculated stress in grain */
  float  stress_comp;		/* calculated comp. stress in grain */
  int	 delete_flag;		/* flag for island coarsening */
} grain_struct;

/* Structure of intersection point information
** Point of intersection between two circles (island/grains)
*/
typedef struct
{ double x, y;			/* x,y coordinate of intersection point */
  double t;			/* time intersection point formed */
  int	 cw;			/* next point in cw direction around island */
  int	 ccw;			/* next point in ccw direction around island */
  int  	 grain[2];		/* two grains involved with intersection point */
  double theta;			/* angle between int. point and nucl. point of island */
} point_struct;

/* Structure of boundary point information
** Boundary points are the evenly spaced points that form the grain
** boundary between intersection and/or triple points around the
** island/grain
*/
typedef struct
{
  double x1, y1;		/* x,y coordinate of boundary point */
  double x2, y2;		/* x,y coordinate of next boundary point in ccw direction */
} bound_struct;

/* Structure of contact point information
** Two contact points are formed the first time two islands contact
*/
typedef struct
{ float	t;		/* time contact formed */
  int	grain;		/* grain#1 associated with island contact */
  int	grain_contact;	/* grain#2 associated with island contact */
  int	point[2];	/* intersection points on either side of contact */
  float length;		/* length of grain boundary between ints. points */
  float radius;		/* radius of grain#1 */
  float z0;		/* zipping height due to island coalescence */
  float height;		/* height of grain boundary */
  float theta;		/* angle between contact point & nucl. point of grain#1 */
  float stress;		/* zipping stress associated with coalescence */
  float area;		/* grain boundary area */
} contact_struct;
