/*
  based on the code originally written by Michael Gordon:
  https://github.com/IEK-5/ssdp_hdf5_extention.git
*/

#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#ifdef WIN32
#include <io.h>
#define F_OK 0
#define access _access
#endif

#include "h5io.h"

static hid_t h5io_fopen(const char *);
static hid_t subgroups(void);
static hid_t float16(void);
static hid_t chunks(int, int, int);
static hid_t cache(int,int);
static hid_t str2dtype(const char*);
static int read_arr(hid_t, hid_t, double**, int, int);
static int write_arr(hid_t, hid_t, double*, int, int);


struct h5io *h5io_init(const char *fn)
{
        struct h5io *self = malloc(sizeof(*self));
        if (NULL == self) goto eself;

        if (H5I_INVALID_HID == (self->file=h5io_fopen(fn))) goto efile;

        h5io_setdataset(self, "data");
        h5io_setdtype(self, "float64");
        self->compression = 0;
        self->chunkarr = 1;
        self->cachemb = 64;
        self->cacheslots = 12421;

        return self;
efile:
        free(self);
eself:
        return NULL;
}


void h5io_free(struct h5io* self)
{
        H5Fclose(self->file);
        free(self);
        self = NULL;
}


void h5io_setdataset(struct h5io* self, const char* dataset)
{
        strncpy(self->dataset, dataset, 1024);
        self->dataset[1023] = '\0';
}


void h5io_setdtype(struct h5io* self, const char* dtype)
{
        strncpy(self->dtype, dtype, 1024);
        self->dtype[1023] = '\0';
}


int h5io_read(struct h5io* self, double ***data, int* arrlen, int narr)
{
        int i;
        hid_t dst, dsp;
        hsize_t dim[2], mdim[2];

        dst = H5Dopen(self->file, self->dataset, H5P_DEFAULT);
        if (H5I_INVALID_HID == dst) goto edst;
        if (H5I_INVALID_HID == (dsp=H5Dget_space(dst))) goto edsp;
        if (H5Sget_simple_extent_dims(dsp, dim, mdim) < 0) goto edim;
        if ((int)dim[0] < narr) {
                printf("Error: not enough arrays available in the chosen file!\n");
                goto enarr;
        }

        *arrlen = dim[1];
        if (NULL == (*data=malloc(narr*sizeof(**data)))) goto edata;

        for (i=0; i < narr; ++i)
                if (read_arr(dsp, dst, &((*data)[i]), *arrlen, i)) goto erow;

        H5Sclose(dsp);
        H5Dclose(dst);
        return 0;
erow:
        for (int j=0; j < i; ++j)
                free((*data)[j]);
        free(*data);
edata:
enarr:
edim:
        H5Sclose(dsp);
edsp:
        H5Dclose(dst);
edst:
        return -1;
}


int h5io_write(struct h5io* self, double **data, int arrlen, int narr)
{
        int i;
        hid_t dsp, gcp, dtype, dcp, dap, dst;

        dsp = H5Screate_simple(2, (hsize_t []){narr, arrlen}, NULL);
        if (H5I_INVALID_HID == dsp) goto edsp;
        if (H5I_INVALID_HID == (gcp=subgroups())) goto egcp;
        if (H5I_INVALID_HID == (dtype=str2dtype(self->dtype))) goto edtype;

        self->chunkarr = self->chunkarr > narr ? narr : self->chunkarr;
        dcp = chunks(self->chunkarr, arrlen, self->compression);
        if (H5I_INVALID_HID == dcp) goto edcp;
        if (H5I_INVALID_HID == (dap=cache(self->cacheslots,
                                          self->cachemb))) goto edap;

        dst = H5Dcreate(self->file, self->dataset, dtype, dsp, gcp, dcp, dap);
        if (H5I_INVALID_HID == dst) goto edst;

        for (i=0; i < narr; ++i)
                if (write_arr(dsp, dst, data[i], arrlen, i))
                        goto ewrite;

        H5Dclose(dst);
        H5Pclose(dap);
        H5Pclose(dcp);
        H5Tclose(dtype);
        H5Pclose(gcp);
        H5Sclose(dsp);
        return 0;
ewrite:
        H5Dclose(dst);
edst:
        H5Pclose(dap);
edap:
        H5Pclose(dcp);
edcp:
        H5Tclose(dtype);
edtype:
        H5Pclose(gcp);
egcp:
        H5Sclose(dsp);
edsp:
        return -1;
}


int h5io_isin(struct h5io* self)
{
        // I want to avoid any additional error messages from
        // HDF. That's why I test that every subgroup exist, starting
        // from the first one.  see:
        // https://docs.hdfgroup.org/hdf5/v1_12/group___h5_l.html#title11
        int res=1;
        // ensure appending delimeter last '/' is okay
        int n = strlen(self->dataset) + 1;
        char *token, *x, *p, delim[] = "/";
        if (NULL == (x=malloc((n+1)*sizeof(*x)))) goto ex;
        if (NULL == (p=malloc((n+1)*sizeof(*p)))) goto ep;
        strncpy(x, self->dataset, n+1);
        p[0] = '\0';

        token = strtok(x,"/");
        while (NULL != token) {
                strcat(p, token);

                res = H5Lexists(self->file, p, H5P_DEFAULT);
                if (res <= 0)
                        break;

                strcat(p, delim);
                token = strtok(NULL, "/");
        }

        free(p);
        free(x);
        return res;
ep:
        free(p);
ex:
        return -2;
}


static hid_t h5io_fopen(const char *fn)
{
        if (0 == access(fn, F_OK))
                return H5Fopen(fn, H5F_ACC_RDWR, H5P_DEFAULT);

        return H5Fcreate(fn, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
}


static hid_t subgroups(void)
{
        hid_t x;
        if (H5I_INVALID_HID == (x = H5Pcreate(H5P_LINK_CREATE))) goto ecreate;
        if (H5Pset_create_intermediate_group(x, 1) < 0) goto egroup;

        return x;
egroup:
        H5Pclose(x);
ecreate:
        return H5I_INVALID_HID;
}


static hid_t float16(void)
{
        hid_t x = H5Tcopy(H5T_NATIVE_FLOAT);
        if (H5Tset_fields(x, 15, 10, 5, 0, 10) < 0) goto err;
        if (H5Tset_size(x, 2) < 0) goto err;
        if (H5Tset_ebias(x, 15) < 0) goto err;

        return x;
err:
        H5Tclose(x);
        return H5I_INVALID_HID;
}


static hid_t chunks(int narr, int arrlen, int level)
{
        if (narr <= 0)
                return H5P_DEFAULT;

        hid_t x = H5Pcreate(H5P_DATASET_CREATE);
        if (H5Pset_chunk(x, 2, (hsize_t []){narr, arrlen}) < 0) goto err;

        level = level < 0 ? 0 : level > 9 ? 9 : level;
        if (H5Pset_deflate(x, level) < 0) goto err;

        return x;
err:
        H5Pclose(x);
        return H5I_INVALID_HID;
}


static hid_t cache(int cacheslots, int cachemb)
{
        // ? see: https://docs.hdfgroup.org/hdf5/v1_14/group___d_a_p_l.html#ga104d00442c31714ee073dee518f661f1
        hid_t x = H5Pcreate(H5P_DATASET_ACCESS);
        if (H5Pset_chunk_cache(x,cacheslots, cachemb*1024*1024,
                               H5D_CHUNK_CACHE_W0_DEFAULT) < 0) goto err;
        return x;
err:
        H5Pclose(x);
        return H5I_INVALID_HID;
}


static hid_t str2dtype(const char *name)
{
        if (0 == strcmp(name, "float64"))
                return H5Tcopy(H5T_NATIVE_DOUBLE);
        else if (0 == strcmp(name, "float32"))
                return H5Tcopy(H5T_NATIVE_FLOAT);
        else if (0 == strcmp(name, "float16"))
                return float16();
        else if (0 == strcmp(name, "int16"))
                return H5Tcopy(H5T_NATIVE_SHORT);
        else if (0 == strcmp(name, "uint16"))
                return H5Tcopy(H5T_NATIVE_USHORT);
        else if (0 == strcmp(name, "int32"))
                return H5Tcopy(H5T_NATIVE_INT);
        else if (0 == strcmp(name, "uint32"))
                return H5Tcopy(H5T_NATIVE_UINT);
        else if (0 == strcmp(name, "int64"))
                return H5Tcopy(H5T_NATIVE_LLONG);
        else if (0 == strcmp(name, "uint64"))
                return H5Tcopy(H5T_NATIVE_ULLONG);
        else
                printf("Error: unknown data type=%s\n", name);

        return H5I_INVALID_HID;
}


static int read_arr(hid_t dsp, hid_t dst, double** data, int nrow, int icol)
{
        herr_t status;
        hid_t hs = H5Screate_simple(2, (hsize_t []){1, nrow}, NULL);
        if (H5I_INVALID_HID == hs) goto ehs;

        status = H5Sselect_hyperslab(dsp, H5S_SELECT_SET,
                                     (hsize_t []){icol, 0}, NULL,
                                     (hsize_t []){1, nrow}, NULL);
        if (status < 0) goto err;

        if (NULL == (*data=malloc(nrow*sizeof(**data)))) goto edata;

        status = H5Dread(dst, H5T_NATIVE_DOUBLE, hs, dsp, H5P_DEFAULT, *data);
        if (status < 0) goto eread;

        status = H5Sselect_none(dsp);
        if (status < 0) goto eread;

        H5Sclose(hs);
        return 0;
eread:
        free(*data);
edata:
err:
        H5Sclose(hs);
ehs:
        return -1;
}


static int write_arr(hid_t dsp, hid_t dst, double* data, int nrow, int icol)
{
        herr_t status;
        hid_t hs = H5Screate_simple(2, (hsize_t []){1, nrow}, NULL);
        if (H5I_INVALID_HID == hs) goto ehs;

        status = H5Sselect_hyperslab(dsp, H5S_SELECT_SET,
                                     (hsize_t []){icol, 0}, NULL,
                                     (hsize_t []){1, nrow}, NULL);
        if (status < 0) goto err;

        status = H5Dwrite(dst, H5T_NATIVE_DOUBLE, hs, dsp, H5P_DEFAULT, data);
        if (status < 0) goto err;

        status = H5Sselect_none(dsp);
        if (status < 0) goto err;

        H5Sclose(hs);
        return 0;
err:
        H5Sclose(hs);
ehs:
        return -1;
}


#ifdef RUNTEST

#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>


double fnmb(const char *fn)
{
        struct stat st;
        stat(fn, &st);
        return (double) st.st_size / 1000000;
}


double** test_data_init(int ncol, int nrow)
{
        int i, j;
        double **data;
        assert((data=malloc(ncol*sizeof(*data))));

        for (i=0; i<ncol; ++i) {
                assert((data[i] = malloc(nrow*sizeof(*(data[i])))));

                for (j=0; j<nrow; ++j)
                        data[i][j] = (double)(i*j);
        }

        return data;
}


void test_data_free(double** data, int ncol)
{
        int i;
        for (i=0; i<ncol; ++i)
                free(data[i]);
        free(data);
}


void test_write(int ncol, int nrow, const char* fn, const char *dt)
{
        struct h5io* io = h5io_init(fn);
        assert(io);
        h5io_setdataset(io, "data/test");
        h5io_setdtype(io, dt);
        double** data = test_data_init(ncol, nrow);
        assert(data);
        assert(0 == h5io_write(io, data, nrow, ncol));
        test_data_free(data, ncol);
        h5io_free(io);
}


void test_read(int ncol, int nrow, const char* fn)
{
        struct h5io* io = h5io_init(fn);
        assert(io);
        h5io_setdataset(io,"data/test");
        int n,i,j;
        double **data;
        assert(0 == h5io_read(io, &data, &n, ncol));
        assert(data);
        assert(n == nrow);

        for (i=0; i < ncol; ++i)
                for (j=0; j < nrow; ++j)
                        assert(fabs((double)i*j - data[i][j]) < 1e-8);

        test_data_free(data, ncol);
        h5io_free(io);
}


void time_write(int narr, int arrlen, int chunkarr,
                int compression, int cachemb,
                const char* dtype)
{
        double tic;
        double** data = test_data_init(narr, arrlen);
        assert(data);

        struct h5io* io = h5io_init("test.h5");
        assert(io);

        io->compression = compression;
        io->cachemb = cachemb;
        io->chunkarr = chunkarr;
        h5io_setdtype(io, dtype);

        tic = (double)clock();
        assert(0 == h5io_write(io, data, arrlen, narr));
        tic = (double)(clock()-tic)/CLOCKS_PER_SEC;
        printf("%7s\t%7d\t%10d\t%8d\t%8d\t%7d\t%12.0f\t%12.2f\n",
               dtype, narr, arrlen, chunkarr,
               compression, cachemb,
               100 / tic,
               fnmb("test.h5"));

        test_data_free(data, narr);
        h5io_free(io);
        assert(0 == remove("test.h5"));
}


void test_timing()
{
        printf("%7s\t%7s\t%10s\t%8s\t%8s\t%7s\t%12s%16s\n",
               "dtype", "narr", "arrlen", "chunkarr",
               "ziplevel", "cachemb", "speed (MB/s)",
               "filesize (MB)");
        time_write(100, 125000,  0, 0, 16, "float64");
        time_write(100, 125000,  1, 0, 16, "float64");
        time_write(100, 125000,  3, 0, 16, "float64");
        time_write(100, 125000, 10, 0, 16, "float64");
        time_write(100, 125000, 50, 0,256, "float64");
        time_write(100, 125000, 50, 1,256, "float64");
        time_write(100, 125000,  1, 1, 16, "float64");
        time_write(1, 12500000,  0, 0, 16, "float64");
        time_write(1, 12500000,  1, 1, 16, "float64");
        time_write(1, 25000000,  1, 1, 16, "float32");
        time_write(1, 50000000,  1, 1,256, "float16");
        time_write(1, 50000000,  1, 0,256, "float16");
        time_write(1, 50000000,  1, 0,256, "int16");
        time_write(1, 50000000,  1, 1,256, "int16");
        time_write(1, 50000000,  1, 5,256, "int16");
        time_write(1, 50000000,  1, 9,256, "int16");
        time_write(1, 12500000,  1, 0,256, "int64");
        time_write(1, 12500000,  1, 1,256, "int64");
        time_write(1, 12500000,  1, 5,256, "int64");
        time_write(1, 12500000,  1, 9,256, "int64");
}


int test_isin(const char* fn, const char *name)
{
        struct h5io* io = h5io_init(fn);
        h5io_setdataset(io, name);
        int res = h5io_isin(io);
        h5io_free(io);
        return res;
}


void test0(const char* dtype)
{
        test_write(10, 20, "test.h5", dtype);
        for (int i=1; i<10; ++i)
                test_read(i, 20, "test.h5");
        assert(0 == remove("test.h5"));
}

void test1()
{
        test_write(10, 20, "test.h5", "float64");
        assert(0 == test_isin("test.h5","/x/x"));
        assert(0 == test_isin("test.h5","asdata/float32"));
        assert(0 == test_isin("test.h5","xdata/afloat32"));
        assert(1 == test_isin("test.h5","//"));
        assert(1 == test_isin("test.h5","///"));
        assert(1 == test_isin("test.h5","data"));
        assert(1 == test_isin("test.h5","data/test"));
        assert(0 == remove("test.h5"));
}


int main(int argc, char** argv)
{
        printf("testing h5io ...\n");

        test0("float32");
        test0("float16");
        test0("int16");
        test1();

        if (argc > 1)
                test_timing();

        printf("PASSED\n");
}

#endif
