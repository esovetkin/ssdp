#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ftw.h>
#include <fts.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"


void print_array(const int *array, int length){
    for (int i = 0; i < length; i++){
        printf("%i ",array[i]);
    }
    printf("\n");
}

void make_dir(const char* path){
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        mkdir(path, 0700);
    }
    else{
        fprintf(stderr, "Could not create directory under %s\n", path);
    }
}

/*
    create a temporary filename
    must be freed
    args:
        dir: name of a directory without a trailing /
*/
char * make_tempfile(const char *dir, bool create){
    int fd;
    char* filename;
    filename = concat(dir, "/fileXXXXXX.h5");
    fd = mkstemps(filename,3);
    if(!create){
        close(fd);
        remove(filename);
    }
    return filename;

}

char* concat(const char *s1, const char *s2)
{
    char *result = malloc(strlen(s1) + strlen(s2) + 1); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
    if(NULL == result){
        // IDK how to properly work here in case malloc does not work
        // if we only had a mechanism to raise Exceptions in C...
        return "";
    }
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

int recursive_delete(const char *dir)
{
    int ret = 0;
    FTS *ftsp = NULL;
    FTSENT *curr;

    // Cast needed (in C) because fts_open() takes a "char * const *", instead
    // of a "const char * const *", which is only allowed in C++. fts_open()
    // does not modify the argument.
    char *files[] = { (char *) dir, NULL };

    // FTS_NOCHDIR  - Avoid changing cwd, which could cause unexpected behavior
    //                in multithreaded programs
    // FTS_PHYSICAL - Don't follow symlinks. Prevents deletion of files outside
    //                of the specified directory
    // FTS_XDEV     - Don't cross filesystem boundaries
    ftsp = fts_open(files, FTS_NOCHDIR | FTS_PHYSICAL | FTS_XDEV, NULL);
    if (!ftsp) {
        fprintf(stderr, "%s: fts_open failed: %s\n", dir, strerror(errno));
        ret = -1;
        goto finish;
    }

    while ((curr = fts_read(ftsp))) {
        switch (curr->fts_info) {
        case FTS_NS:
        case FTS_DNR:
        case FTS_ERR:
            fprintf(stderr, "%s: fts_read error: %s\n",
                    curr->fts_accpath, strerror(curr->fts_errno));
            break;

        case FTS_DC:
        case FTS_DOT:
        case FTS_NSOK:
            // Not reached unless FTS_LOGICAL, FTS_SEEDOT, or FTS_NOSTAT were
            // passed to fts_open()
            break;

        case FTS_D:
            // Do nothing. Need depth-first search, so directories are deleted
            // in FTS_DP
            break;

        case FTS_DP:
        case FTS_F:
        case FTS_SL:
        case FTS_SLNONE:
        case FTS_DEFAULT:
            if (remove(curr->fts_accpath) < 0) {
                fprintf(stderr, "%s: Failed to remove: %s\n",
                        curr->fts_path, strerror(curr->fts_errno));
                ret = -1;
            }
            break;
        }
    }

finish:
    if (ftsp) {
        fts_close(ftsp);
    }

    return ret;
}