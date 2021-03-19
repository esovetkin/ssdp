#ifndef _UTIL_H
#define _UTIL_H
typedef enum {QUIET, NORMAL, VERBOSE} VERB;
extern VERB Verbosity;

void Print(VERB v, const char *format_str, ...);
void Fatal( const char *format_str, ...);
void Warning( const char *format_str, ...);
int ProgressBar(int pcn, int pco, int len, int tics);
#endif
