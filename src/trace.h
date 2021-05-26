

//BEGIN_SSDP_EXPORT
typedef struct sky_transfer {
	double *t; // transfer sky-patch->POA
	double g;  // transfer ground to GHI->POA
	int N;
} sky_transfer;

typedef struct horizon {
	int N;
	double astep;
	double *zen;
} horizon;

typedef struct skypoly {
	int N;
	double *x, *y, *z;
} skypoly;
//END_SSDP_EXPORT


void FreeSkyTransfer(sky_transfer *T);
sky_transfer InitSkyTransfer(int Nz);
horizon InitHorizon(int Nz);
void AtanHorizon(horizon *H);
void FreeHorizon(horizon *H);
int BelowHorizon(const horizon *H, sky_pos p);

void TransTrans(const sky_transfer *a, const sky_transfer *b, sky_transfer *c);
void HorizTrans(const sky_grid *sky, const horizon *a, const sky_transfer *b, sky_transfer *c);
