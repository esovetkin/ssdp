# define ADDR_OF(X) ((int[]){X})
#include <stdlib.h>
#include <stdbool.h>

void print_array(const int *array, int length);

char* concat(const char *s1, const char *s2);

void make_dir(const char* path);

int recursive_delete(const char *dir);

char * make_tempfile(const char *dir, bool create);