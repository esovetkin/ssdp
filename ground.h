typedef struct topology {
	vec *points; // 3D coordinates
	int N;		 // number of points
	double d;	 // diameter of the points (all points have the same diameter)
				 // you'll have to create more than one topology to have different sizes
} topology;


void free_topo (topology *T);
void MakeHorizon(sky_grid *sky, topology T, double height);
void ClearHorizon(sky_grid *sky);
