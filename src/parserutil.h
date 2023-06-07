#ifndef _PARSEUTILS_H
#define _PARSEUTILS_H
#define RADPDEG 1.745329251994329e-02
#define DEGPRAD 5.729577951308232e+01
#define rad2deg(r) (DEGPRAD*(r))
#define deg2rad(d) (RADPDEG*(d))

#ifdef OPENMP
	extern double tic;
	double omp_get_wtime(void);
	#define TIC() (tic=omp_get_wtime())
	#define TOC() (omp_get_wtime()-tic)
#else // OPENMP
	extern clock_t tic;
	#define TIC() (tic=clock())
	#define TOC() ((double)(clock()-tic)/CLOCKS_PER_SEC)
#endif // OPENMP

struct cvec {
        char **s;
        int astep, a, n;
};

struct cvec* cvec_init(int a);
void cvec_free();
// returns -1 on realloc fails
int cvec_push(struct cvec*, char *word);
// returns -1 if failed to read line from fn
int read_filelist(char *fn, struct cvec *dst);

/* utility function to write parsers */
char * GetWord(const char *in, char *word);
int GetOption(const char *in, const char *opt, char *word);
int GetArg(const char *in, const char *opt, char *word);
int GetNumOption(const char *in, const char *opt, int i, char *word);

int FetchConfig(const char *in, const char *pat, char *str, simulation_config **a);
int FetchArray(const char *in, const char *pat, char *str, array **a);
int FetchOptArray(const char *in, const char *pat, char *str, array **a);

int FetchFloat(const char *in, const char *pat, char *str, double *a);
int FetchOptFloat(const char *in, const char *pat, char *str, double *a);
int FetchInt(const char *in, const char *pat, char *str, int *a);
int FetchOptInt(const char *in, const char *pat, char *str, int *a);

void InitConfigMask(simulation_config *C);
void InitConfigMaskNoH(simulation_config *C);

#define ProgressLen 40
#define ProgressTics 4
#endif //_PARSEUTILS_H

