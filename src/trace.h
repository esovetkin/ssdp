
//BEGIN_SSDP_EXPORT
typedef struct sky_transfer {
	double *t; // transfer sky-patch->POA
	double g;  // transfer ground to GHI->POA
	int N;
} sky_transfer;

typedef struct horizon {
		double *zen;
		float astep, minz, maxz;
		int N;
} horizon;
typedef struct skypoly {
	int N;
	double *x, *y, *z;
} skypoly;
struct rtreecache {
		struct rtree *rtree;
		double xydelta;
		double zdelta;
};
//END_SSDP_EXPORT


void FreeSkyTransfer(sky_transfer *T);
sky_transfer InitSkyTransfer(int Nz, sky_transfer *T);
horizon InitHorizon(int Nz);
void AtanHorizon(horizon *H); // also sets minz, maxz
void FreeHorizon(horizon *H);
int BelowHorizon(const horizon *H, sky_pos p);

void TransTrans(const sky_transfer *a, const sky_transfer *b, sky_transfer *c);
void HorizTrans(const sky_grid *sky, const horizon *a, const sky_transfer *b, sky_transfer *c);

/** initialise rtree

	@xydelta, @zdelta: coordinates within a box of size

        xydelta x xydelta x zdelta

	points to the same horizon
*/
struct rtreecache* rtreecache_init(double xydelta, double zdelta);


/** cleanup rtree **and** all the computed horizons/sky transfers
 */
void horizoncache_free(struct rtreecache* self);
void stcache_free(struct rtreecache* self);


/** get pointer to a horizon

	if the horizon at the given location is stored in the index, the
	horizon is returned. if the location is new, the new horizon is
	created (with h->zen == NULL).

	returns NULL if something went wrong
*/
horizon* horizoncache_get(struct rtreecache* self, double x, double y, double z);
sky_transfer* stcache_get(struct rtreecache* self, double x, double y, double z);
