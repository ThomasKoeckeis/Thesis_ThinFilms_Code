/* external.h */
#include "struct.h"
#include <vector>
#include "Point.h"
#include "Polygon.h"

/* The main structures used in the simulation (see struct.h) */
grain_struct	*grain;
point_struct	*point;
int pointcount;
bound_struct	*bound;
contact_struct	*contact;

int numcontact_extern = 0;//count of contact points
std::vector<Point> grain_points = std::vector<Point>(); //points on floor
std::vector<Point> grain_points_zip = std::vector<Point>(); //zipping points
std::vector<Polygon> grain_polygons = std::vector<Polygon>(); //polygons of floor points
std::vector<Polygon> grain_polygons_zip = std::vector<Polygon>(); //polygons of zipping points
/*
** The scaling factor for the simulation relates simulation dimension
** to real dimensions;
** radius_real = scaling*radius_simulation
*/
float scaling;

/*
** The lateral growth rate of the growing island in the simulation
** g0 = dep_rate_real*geom/scaling
*/
float g0;

/* Dimensions of entire simulation field of grains */
double xmin, xmax, ymin, ymax;
double xsize, ysize;

/* Simulation time in seconds */
float t;

/* Current total number of grains including those due to PCB */
int ng;

int long numbound;

/* allocate.c */
extern int  allocate_memory(int, float);

/*nucleate.c*/
extern void nucleate_island(float *);

/* intersection.c */
extern void island_intersection(int *);
extern void organize_intersection_points();
extern void draw_island_boundary();

/* contact.c */
extern void new_contacts(int *);
extern void update_contacts(int);
extern void stress_contacts(int, float, float, float, float, float);

/*extern float gb_length(int);
extern float zip_height(int, int, float, float, float); */

/* auxiliary.c */
extern double calc_grain_area(int);
extern int    inbound(double, double);

/* coarsen.c */
extern void island_coalesce(int, int, double);
extern void island_coalesce_noshift(int, int, double);

/* graphics.c */
extern void ps_write(FILE *, float, float, float, double, double, double, double, int);

/* stress.c */
extern void laplace_pressure(float, float, float);
extern void stress_relax_coble(float, float, float, float);
extern void film_stress_relax_coble(float, float, float, float);

