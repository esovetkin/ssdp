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

// number of locations (and effective number Nl_eff that depends on
// chunks). Nl_o correspond to the current chunk shift
		int Nl, Nl_eff, Nl_o;

		// flags indicating whether the structures have been initialized
		char sky_init, topo_init, grid_init, loc_init;

// number of topogrids. Multiple topogrids can be setup and contribute
// to the horizon. Order matters.
		int nTx;

// number of ->S, equals to number of threads
		int nS;

// longitude, latitude, elevation
		double lon, lat, E;

// sky, one sky for each OMP thread, at least one
		sky_grid* S;

// topogrid
		topogrid* Tx;

// ground albedo
		double albedo;

// locations and orientation of the module(s)
		double *x, *y, *z;
		sky_pos *o;

// traced locations
		location *L;

// ordering locations according to the  pixels of the **first** topogrid pixels
		int *xyorder;

// horizon cache
		struct rtreecache* hcache;
		horizon **uH;	// list of pointers to unique horizons, always length Nl
		int *uHi;	// index uH -> index location, len(uHi) = Nl

// stcache: initial sky transfer cache for each location the initial
// sky is copied. uST is list of pointers to unique sky transfers,
// always length Nl. uSTi is the index uST -> index location,
// len(uSTi) = NluSTii is the inverse of uSTi, len(uSTii) = Nl
		struct rtreecache* stcache;
		sky_transfer **uST;
		int *uSTi;
		int *uSTii;

// number of points used in the approximate horizon
		int approx_n;

// approximate horizon sampling type
		enum SampleType approx_stype;

// number of locations do process in one chunk. If negative process
// everything in one chunk
		int chunked;

// angle of incidence effects
		AOI_Model_Data M;

// irregular topology
		topology T;

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
