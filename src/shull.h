#ifndef shull_h
#define shull_h

#include <stdint.h>
#ifdef _WIN64
#define _UNISTD_H    1

/* This is intended as a drop-in replacement for unistd.h on Windows.
 * Please add functionality as neeeded.
 * https://stackoverflow.com/a/826027/1202830
 */

#include <stdlib.h>
#include <getopt.h> /* getopt at: https://gist.github.com/ashelly/7776712 */
#include <process.h> /* for getpid() and the exec..() family */
#include <direct.h> /* for _getcwd() and _chdir() */

#define srandom srand
#define random rand

/* Values for the second argument to access.
   These may be OR'd together.  */
#define R_OK    4       /* Test for read permission.  */
#define W_OK    2       /* Test for write permission.  */
//#define   X_OK    1       /* execute permission - unsupported in windows*/
#define F_OK    0       /* Test for existence.  */

#define access _access
#define dup2 _dup2
#define execve _execve
#define ftruncate _chsize
#define unlink _unlink
#define fileno _fileno
#define getcwd _getcwd
#define chdir _chdir
#define isatty _isatty
#define lseek _lseek
/* read, write, and close are NOT being #defined here, because while there are file handle specific versions for Windows, they probably don't work for sockets. You need to look at your app and consider whether to call e.g. closesocket(). */

#ifdef _WIN64
#define ssize_t __int64
#else
#define ssize_t long
#endif

#define STDIN_FILENO 0
#define STDOUT_FILENO 1
#define STDERR_FILENO 2

#else
#include <unistd.h>
#endif


#include "ll.h"


typedef struct sh_point sh_point;
typedef struct sh_edge sh_edge;
typedef struct sh_triangle sh_triangle;

struct sh_point {
	double x;
	double y;
	int index;
};

struct sh_edge {
	sh_point *p[2];
	sh_triangle *t[2];
	unsigned int flipcount;
};

struct sh_triangle {
	sh_point *p[3];
	sh_edge *e[3];
	sh_point cc;
	double ccr2;
};

typedef struct {
	ll_node *hull_edges;
	ll_node *internal_edges;
	ll_node *triangles;
} sh_triangulation_data;

int triangulate(sh_triangulation_data *td, sh_point *ps, size_t n);
int make_delaunay(sh_triangulation_data *td);
int delaunay(sh_triangulation_data *td, sh_point *ps, size_t n);

#endif
