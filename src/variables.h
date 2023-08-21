/*      POA_ABS     absolurte height and orientation
 *      POA_ZREL    height is relative to surface level
 *      POA_SURFREL height and orientation relative to surface level
 */
typedef struct simulation_config {
	/* a collection of data we need
	 * The only things not in here are 
	 * - time
	 * - GHI
	 * - DHI
	 */
	AOI_Model_Data M; 		// angle of incidence effects
	double lon, lat, E;		// longitude, latitude, elevation
	char sky_init, topo_init, grid_init; // flags indicating whether the structures have been initialized or not
	sky_grid S; 			// sky
	topology T; 			// topology
	topogrid Tx; 			// topogrid
	double albedo;			// ground albedo
	char loc_init; 			// flags indicating whether location arrays have been initialized or not
	double *x, *y, *z;		// locations of the module(s)
	sky_pos *o;				// orientation of the module(s)
	location *L;			// traced locations
	int Nl;					// number of locations
	struct horizoncache* hcache;	// horizon cache
	horizon **uH;	// list of pointers to unique horizons, always length Nl
	int *uHi;	// index uH -> index location, len(uHi) = Nl
} simulation_config;

typedef struct array {
	double *D; // array of double floats
	int N; // length of array
} array;

simulation_config InitConf();
void InitVars();
void ClearVars();
int lookupvar(char *name);
int LookupSimConf(char *name, simulation_config **C);
int LookupArray(char *name, array **a);
int AddSimConf(char *name, simulation_config C);
int AddArray(char *name, array a);
int RMVar(char *name);

void ListVars();
