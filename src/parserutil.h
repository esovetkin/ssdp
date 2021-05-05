#ifndef _PARSEUTILS_H
#define _PARSEUTILS_H
#define RADPDEG 1.745329251994329e-02
#define DEGPRAD 5.729577951308232e+01
#define rad2deg(r) (DEGPRAD*(r))
#define deg2rad(d) (RADPDEG*(d))
extern clock_t tic;
#define TOC() ((double)(clock()-tic)/CLOCKS_PER_SEC)
#define TIC() (tic=clock())

/* utility function to write parsers */
char * GetWord(const char *in, char *word);
int GetOption(const char *in, const char *opt, char *word);
int GetArg(const char *in, const char *opt, char *word);

int FetchConfig(const char *in, const char *pat, char *str, simulation_config **a);
int FetchArray(const char *in, const char *pat, char *str, array **a);

int FetchFloat(const char *in, const char *pat, char *str, double *a);
int FetchInt(const char *in, const char *pat, char *str, int *a);

void InitConfigMask(simulation_config *C);
void InitConfigMaskNoH(simulation_config *C);

#define ProgressLen 40
#define ProgressTics 4
#endif //_PARSEUTILS_H

