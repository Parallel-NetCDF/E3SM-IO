/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM-IO benchmark software package.
 *
 * See README.md for compile and run instructions.
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <unistd.h> /* getopt(), access() */

#include <assert.h>

#include <mpi.h>
#include <pnetcdf.h>

#ifndef MAX
#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#endif
#ifndef MIN
#define MIN(a, b) ((a) < (b)) ? (a) : (b)
#endif

static int        verbose;
static int        world_nprocs;
static int        world_rank;
static MPI_Comm   sub_comm;
static int        sub_nprocs;
static int        sub_rank;

typedef struct {
    MPI_Offset   nelems;     /* number of array elements */
    MPI_Offset   max_nreqs;
    int          ndims;
    MPI_Offset   dims[3];
    int          nreqs;
    int         *offsets;
    int         *lengths;
    MPI_Offset   r_start;  /* array starting index for reading */
    MPI_Offset   r_count;  /* number of array elements for reading */
    MPI_Offset **w_starts; /* [nreqs][] write starting indices for rec var */
    MPI_Offset **w_counts; /* [nreqs][] write element counts for rec var */
    MPI_Offset **w_startx; /* [nreqs][] used for fixed-size variables */
    MPI_Offset **w_countx; /* [nreqs][] used for fixed-size variables */
} io_decomp;

typedef struct {
    int varid;          /* NetCDF variable ID */
    int ndims;          /* number of dimensions */
    int is_rec;         /* is a record or fixed-size variable */
    int dec_id;         /* decomposition ID */
    size_t vlen;        /* length in bytes */
    MPI_Offset dims[3]; /* dimension sizes */
} io_var;

#define CHECK_MPI_ERROR(mpi_errorcode, err_msg) {                      \
    if (mpi_errorcode != MPI_SUCCESS) {                                \
        char errorString[MPI_MAX_ERROR_STRING];                        \
        int errorStringLen;                                            \
        MPI_Error_string(mpi_errorcode, errorString, &errorStringLen); \
        printf("%2d: MPI Failure at line %d of %s (%s : %s)\n",        \
               world_rank, __LINE__, __FILE__, err_msg, errorString);  \
        return -1;                                                     \
    }                                                                  \
}

#define CHECK_NC_ERR { \
    if (err != NC_NOERR) { \
        fprintf(stderr, "Error at line %d: %s\n", __LINE__, \
                ncmpi_strerrno(err)); \
        goto err_out; \
    } \
}

#define CHECK_VAR_ERR(ncid, varid) {                               \
    if (err != NC_NOERR) {                                         \
        char var_name[64];                                         \
        ncmpi_inq_varname(ncid, varid, var_name);                  \
        printf("Error in %s:%d: %s() var %s (%s)\n", __FILE__,     \
               __LINE__, __func__, var_name, ncmpi_strerrno(err)); \
        goto err_out;                                              \
    }                                                              \
}

static int
xlen_nc_type(nc_type xtype)
{
    switch(xtype) {
        case NC_BYTE:
        case NC_CHAR:
        case NC_UBYTE:  return 1;
        case NC_SHORT:
        case NC_USHORT: return 2;
        case NC_INT:
        case NC_UINT:
        case NC_FLOAT:  return 4;
        case NC_DOUBLE:
        case NC_INT64:
        case NC_UINT64: return 8;
        default: return 0;
    }
}

/*----< split_comm() >-------------------------------------------------------*/
static
int split_comm(MPI_Comm  orig_comm,
               MPI_Comm *sub_comm,
               int      *num_subfiles,
               int      *subfile_ID)
{
    int err, color, orig_rank, sub_nprocs, sub_rank, int_msg[2];
    MPI_Comm comm_roots;

    MPI_Comm_rank(orig_comm, &orig_rank);

    /* split communicator to create one sub-communicator per compute node */
    err = MPI_Comm_split_type(orig_comm, MPI_COMM_TYPE_SHARED, orig_rank,
                              MPI_INFO_NULL, sub_comm);
    CHECK_MPI_ERROR(err, "MPI_Comm_split_type")

    /* calculate subfile ID and take care of both process rank assignments:
     * block-based (MPICH_RANK_REORDER_METHOD=1) or
     * round-robin (MPICH_RANK_REORDER_METHOD=0)
     */
    MPI_Comm_size(*sub_comm, &sub_nprocs);
    MPI_Comm_rank(*sub_comm, &sub_rank);
    color = (sub_rank == 0) ? 1 : 0;
    err = MPI_Comm_split(orig_comm, color, orig_rank, &comm_roots);
    CHECK_MPI_ERROR(err, "MPI_Comm_split")

    MPI_Comm_size(comm_roots, num_subfiles);
    MPI_Comm_rank(comm_roots, subfile_ID);
    MPI_Comm_free(&comm_roots);

    int_msg[0] = *num_subfiles;
    int_msg[1] = *subfile_ID;
    err = MPI_Bcast(int_msg, 2, MPI_INT, 0, *sub_comm);
    CHECK_MPI_ERROR(err, "MPI_Bcast")
    *num_subfiles = int_msg[0];
    *subfile_ID   = int_msg[1];
    if (verbose && sub_rank == 0)
        printf("world_rank=%5d sub_rank=%5d color=%d subfile_ID=%5d\n",
               world_rank, sub_rank, color, *subfile_ID);

    return (err != MPI_SUCCESS);
}

/*----< set_decomp() >-------------------------------------------------------*/
static
int set_decomp(int        ncid,
               int       *num_decomp,
               io_decomp *decomp)
{
    char name[64];
    int i, j, err=0, dimid, varid, nreqs;
    MPI_Offset tmp, nblobs, start[2], count[2];

    /* inquire number of decompositions */
    err = ncmpi_get_att(ncid, NC_GLOBAL, "num_decompositions", num_decomp);
    CHECK_NC_ERR

    if (verbose && sub_rank == 0)
        printf("world_rank=%5d sub_rank=%5d num_decomp=%d\n",
               world_rank, sub_rank, *num_decomp);

    /* inquire decomposition related dimensions */
    for (i=0; i<*num_decomp; i++) {
        io_decomp *dp = decomp + i;

        err = ncmpi_inq_dimid(ncid, "nblobs", &dimid);
        CHECK_NC_ERR
        err = ncmpi_inq_dimlen(ncid, dimid, &nblobs);
        CHECK_NC_ERR

        /* number of array elements (sum across all processes in sub_comm) */
        sprintf(name, "D%d.nelems",i+1);
        err = ncmpi_inq_dimid(ncid, name, &dimid);
        CHECK_NC_ERR
        err = ncmpi_inq_dimlen(ncid, dimid, &dp->nelems);
        CHECK_NC_ERR

        /* max number of offset-length pairs among processes in sub_comm */
        sprintf(name, "D%d.max_nreqs",i+1);
        err = ncmpi_inq_dimid(ncid, name, &dimid);
        CHECK_NC_ERR
        err = ncmpi_inq_dimlen(ncid, dimid, &dp->max_nreqs);
        CHECK_NC_ERR

        /* read variable D*.nreqs */
        sprintf(name, "D%d.nreqs",i+1);
        err = ncmpi_inq_varid(ncid, name, &varid);
        CHECK_NC_ERR

        /* number of dimensions, dimension IDs, and dimension sizes */
        err = ncmpi_inq_attlen(ncid, varid, "global_dimids", &tmp);
        CHECK_VAR_ERR(ncid, varid)
        dp->ndims = tmp;
        int dimids[3]; /* number of fix-sized dimensions is <= 3 */
        err = ncmpi_get_att(ncid, varid, "global_dimids", dimids);
        CHECK_VAR_ERR(ncid, varid)
        for (j=0; j<dp->ndims; j++) {
            err = ncmpi_inq_dimlen(ncid, dimids[j], &dp->dims[j]);
            CHECK_NC_ERR
        }

        /* number of offset-length pairs assigned to this process */
        start[0] = sub_rank;
        err = ncmpi_get_var1_int_all(ncid, varid, start, &nreqs);
        CHECK_VAR_ERR(ncid, varid)
        dp->nreqs = nreqs;

        /* read variable D*.blob_start */
        sprintf(name, "D%d.blob_start",i+1);
        err = ncmpi_inq_varid(ncid, name, &varid);
        CHECK_NC_ERR
        err = ncmpi_get_var1_longlong_all(ncid, varid, start, &dp->r_start);
        CHECK_VAR_ERR(ncid, varid)

        /* read variable D*.blob_count */
        sprintf(name, "D%d.blob_count",i+1);
        err = ncmpi_inq_varid(ncid, name, &varid);
        CHECK_NC_ERR
        err = ncmpi_get_var1_longlong_all(ncid, varid, start, &dp->r_count);
        CHECK_VAR_ERR(ncid, varid)

        /* allocate buffers */
        dp->offsets = (int*) malloc(nreqs * 2 * sizeof(int));
        dp->lengths = dp->offsets + nreqs;

        start[0] = sub_rank;
        start[1] = 0;
        count[0] = 1;
        count[1] = nreqs;

        /* read variable D*.offsets */
        sprintf(name, "D%d.offsets",i+1);
        err = ncmpi_inq_varid(ncid, name, &varid);
        CHECK_NC_ERR
        err = ncmpi_get_vara_int_all(ncid, varid, start, count, dp->offsets);
        CHECK_VAR_ERR(ncid, varid)

        /* read variable D*.lengths */
        sprintf(name, "D%d.lengths",i+1);
        err = ncmpi_inq_varid(ncid, name, &varid);
        CHECK_NC_ERR
        err = ncmpi_get_vara_int_all(ncid, varid, start, count, dp->lengths);
        CHECK_VAR_ERR(ncid, varid)

        /* construct starts[] and counts[] for iput_varn */
        dp->w_starts = (MPI_Offset**) malloc(nreqs * 4 * sizeof(MPI_Offset*));
        dp->w_counts = dp->w_starts + nreqs;
        dp->w_startx = dp->w_counts + nreqs;
        dp->w_countx = dp->w_startx + nreqs;

        dp->w_starts[0] = (MPI_Offset*) malloc(nreqs * 4 * 4 * sizeof(MPI_Offset));
        dp->w_counts[0] = dp->w_starts[0] + nreqs * 4;
        dp->w_startx[0] = dp->w_counts[0] + nreqs * 4;
        dp->w_countx[0] = dp->w_startx[0] + nreqs * 4;
        for (j=1; j<nreqs; j++) {
            dp->w_starts[j] = dp->w_starts[j-1] + 4;
            dp->w_counts[j] = dp->w_counts[j-1] + 4;
            dp->w_startx[j] = dp->w_startx[j-1] + 4;
            dp->w_countx[j] = dp->w_countx[j-1] + 4;
        }

        /* no fixed-size variables are 4 or more dimensional */
        for (j=0; j<nreqs; j++) {
            dp->w_starts[j][0] = 0;
            dp->w_counts[j][0] = 1;
            if (dp->ndims == 1) {
                dp->w_starts[j][1] = dp->offsets[j];
                dp->w_counts[j][1] = dp->lengths[j];
                dp->w_startx[j][0] = dp->w_starts[j][1];
                dp->w_countx[j][0] = dp->w_counts[j][1];
            }
            else if (dp->ndims == 2) { /* 2D */
                dp->w_starts[j][1] = dp->offsets[j] / dp->dims[1];
                dp->w_starts[j][2] = dp->offsets[j] % dp->dims[1];
                dp->w_counts[j][1] = 1;
                dp->w_counts[j][2] = dp->lengths[j];
                dp->w_startx[j][0] = dp->w_starts[j][1];
                dp->w_startx[j][1] = dp->w_starts[j][2];
                dp->w_countx[j][0] = dp->w_counts[j][1];
                dp->w_countx[j][1] = dp->w_counts[j][2];
            }
            else if (dp->ndims == 3) { /* 3D */
                MPI_Offset tmp = dp->offsets[j];
                dp->w_starts[j][1] = tmp / (dp->dims[1] * dp->dims[2]);
                tmp %= dp->dims[1] * dp->dims[2];
                dp->w_starts[j][2] = tmp / dp->dims[2];
                dp->w_starts[j][3] = tmp % dp->dims[2];
                dp->w_counts[j][1] = 1;
                dp->w_counts[j][2] = 1;
                dp->w_counts[j][3] = dp->lengths[j];
                dp->w_startx[j][0] = dp->w_starts[j][1];
                dp->w_startx[j][1] = dp->w_starts[j][2];
                dp->w_startx[j][2] = dp->w_starts[j][3];
                dp->w_countx[j][0] = dp->w_counts[j][1];
                dp->w_countx[j][1] = dp->w_counts[j][2];
                dp->w_countx[j][2] = dp->w_counts[j][3];
             }
             else assert(0);
            /* each length[j] is no bigger than last dims[] */
        }
    }

err_out:
    return (err != 0);
}

/*----< copy_dims() >--------------------------------------------------------*/
static
int copy_dims(int          in_ncid,
              int         out_ncid,
              MPI_Offset *num_recs)
{
    char name[64];
    int i, err, ndims, dimid, rec_dim;
    MPI_Offset *dims;

    /* inquire dimensions in input file and copy over to output file */
    err = ncmpi_inq_ndims(in_ncid, &ndims);
    CHECK_NC_ERR

    if (verbose && sub_rank == 0)
        printf("world_rank=%5d sub_rank=%5d num dims=%d\n",
               world_rank, sub_rank, ndims);

    dims = (MPI_Offset*) malloc(ndims * sizeof(MPI_Offset));

    /* inquire record dimension */
    err = ncmpi_inq_unlimdim(in_ncid, &rec_dim);
    CHECK_NC_ERR
    if (rec_dim >= 0) {
        err = ncmpi_inq_dimlen(in_ncid, rec_dim, num_recs);
        CHECK_NC_ERR
    }
    if (verbose && sub_rank == 0)
        printf("world_rank=%5d sub_rank=%5d num_recs=%lld\n",
               world_rank, sub_rank, *num_recs);

    /* copy over all dimensions except the ones used for decomposition */
    for (i=0; i<ndims; i++) {
        err = ncmpi_inq_dim(in_ncid, i, name, dims+i);
        CHECK_NC_ERR

        int name_len = strlen(name);
        if (strcmp(name, "nblobs") == 0 ||
            (name[0] == 'D' && strcmp(name+name_len-7,  ".nelems")    == 0) ||
            (name[0] == 'D' && strcmp(name+name_len-10, ".max_nreqs") == 0))
            continue;

        if (i == rec_dim)
            err = ncmpi_def_dim(out_ncid, name, NC_UNLIMITED, &dimid);
        else
            err = ncmpi_def_dim(out_ncid, name, dims[i], &dimid);
        CHECK_NC_ERR
    }
    free(dims);

err_out:
    return (err != NC_NOERR);
}

/*----< copy_attr() >--------------------------------------------------------*/
static
int copy_attr(int  in_ncid,
              int out_ncid)
{
    char name[64];
    int i, err, nattrs;

    /* inquire number of global attributes */
    err = ncmpi_inq_varnatts(in_ncid, NC_GLOBAL, &nattrs);
    CHECK_NC_ERR

    if (verbose && sub_rank == 0)
        printf("world_rank=%5d sub_rank=%5d no. G attrs=%d\n",
               world_rank, sub_rank, nattrs);

    /* copy over global attributes except for ones used for decomposition */
    for (i=0; i<nattrs; i++) {
        err = ncmpi_inq_attname(in_ncid, NC_GLOBAL, i, name);
        CHECK_NC_ERR

        if (strcmp(name, "global_nprocs") == 0 ||
            strcmp(name, "num_decompositions") == 0 ||
            strcmp(name, "num_subfiles") == 0)
            continue;

        err = ncmpi_copy_att(in_ncid, NC_GLOBAL, name, out_ncid, NC_GLOBAL);
        CHECK_NC_ERR
    }

err_out:
    return (err != NC_NOERR);
}

/*----< set_vars() >---------------------------------------------------------*/
static
int set_vars(int        in_ncid,
             int        out_ncid,
             int        nvars,
             int        num_decomp,
             io_decomp *decomp,
             io_var    *var)
{
    char name[64];
    int i, j, err=0, rec_dim, nvars_partitioned=0, nvars_decomp=0;
    nc_type xtype;

    err = ncmpi_inq_unlimdim(in_ncid, &rec_dim);
    CHECK_NC_ERR

#define NUM_DECOMP_AUX_VARS 5
#ifdef NUM_DECOMP_AUX_VARS
    /* first nvars_decomp variables are decomposition variables:
     * D*.nreqs, D*.blob_start, D*.blob_count, D*.offsets, D*.lengths
     */
    nvars_decomp = num_decomp * NUM_DECOMP_AUX_VARS;
#endif

    /* copy over variable definition and attributes */
    for (i=0; i<nvars; i++) {
        int dimids[4], nattrs;
        MPI_Offset tmp;

        var[i].varid  = -1; /* decomposition variable */
        var[i].is_rec = 0;
        var[i].dec_id = -1;

#ifdef NUM_DECOMP_AUX_VARS
        if (i < nvars_decomp) {
            /* inquire xtype of variable i */
            err = ncmpi_inq_vartype(in_ncid, i, &xtype);
            CHECK_VAR_ERR(in_ncid, i)
            var[i].vlen = xlen_nc_type(xtype);
            continue;
        }
#endif

        /* inquire metadata of variable i */
        err = ncmpi_inq_var(in_ncid, i, name, &xtype, &var[i].ndims, dimids,
                            &nattrs);
        CHECK_VAR_ERR(in_ncid, i)

        var[i].vlen = xlen_nc_type(xtype);

#ifndef NUM_DECOMP_AUX_VARS
        /* skip copying decomposition variables */
        int name_len = strlen(name);
        if ((name[0] == 'D' && name_len > 6 && strcmp(name+name_len-6, ".nreqs"     ) == 0) ||
            (name[0] == 'D' && name_len > 11 && strcmp(name+name_len-11,".blob_start") == 0) ||
            (name[0] == 'D' && name_len > 11 && strcmp(name+name_len-11,".blob_count") == 0) ||
            (name[0] == 'D' && name_len > 8 && strcmp(name+name_len-8, ".offsets"   ) == 0) ||
            (name[0] == 'D' && name_len > 6 && strcmp(name+name_len-8, ".lengths"   ) == 0)) {
            nvars_decomp++;
            continue;
        }
#endif

        /* inquire global_dimids (global dimension IDs) */
        err = ncmpi_inq_attlen(in_ncid, i, "global_dimids", &tmp);
        if (err == NC_ENOTATT) { /* this variable is not partitioned */
            err = ncmpi_def_var(out_ncid, name, xtype, var[i].ndims, dimids,
                                &var[i].varid);
            if (err != NC_NOERR) {
                printf("Error in %s:%d: %s() var %s (%s)\n", __FILE__,
                       __LINE__, __func__, name, ncmpi_strerrno(err));
                goto err_out;
            }
            if (var[i].ndims > 0 && dimids[0] == rec_dim) var[i].is_rec = 1;

            /* count size of this variable, one record only */
            for (j=0; j<var[i].ndims; j++) {
                err = ncmpi_inq_dimlen(in_ncid, dimids[j], &var[i].dims[j]);
                CHECK_NC_ERR
                if (j == 0 && dimids[0] == rec_dim) continue;
                var[i].vlen *= var[i].dims[j];
            }
            if (dimids[0] == rec_dim) var[i].dims[0] = 1;

            if (verbose && world_rank == 0)
                printf("without decomp varid %d rec %d dec %d var %s\n",
                       var[i].varid,var[i].is_rec,var[i].dec_id,name);
        }
        else if (err == NC_NOERR) {
            /* define variable using the global dimensions */
            var[i].ndims = tmp;
            err = ncmpi_get_att_int(in_ncid, i, "global_dimids", dimids);
            CHECK_VAR_ERR(in_ncid, i)
            err = ncmpi_def_var(out_ncid, name, xtype, var[i].ndims, dimids,
                                &var[i].varid);
            if (err != NC_NOERR) {
                printf("Error in %s:%d: %s() var %s (%s)\n", __FILE__,
                       __LINE__, __func__, name, ncmpi_strerrno(err));
                goto err_out;
            }
            err = ncmpi_get_att_int(in_ncid, i, "decomposition_ID",
                                    &var[i].dec_id);
            CHECK_VAR_ERR(in_ncid, i)
            var[i].dec_id--; /* change to 0-based */
            if (var[i].ndims > 0 && dimids[0] == rec_dim) var[i].is_rec = 1;

            for (j=0; j<var[i].ndims; j++) {
                err = ncmpi_inq_dimlen(in_ncid, dimids[j], &var[i].dims[j]);
                CHECK_NC_ERR
            }
            if (dimids[0] == rec_dim) var[i].dims[0] = 1;

            /* vlen is the read and write byte size by this process */
            var[i].vlen *= decomp[var[i].dec_id].r_count;

            if (verbose && world_rank == 0)
                printf("with decomp varid %d rec %d dec %d var %s\n",
                       var[i].varid,var[i].is_rec,var[i].dec_id,name);
            nvars_partitioned++;
        }
        else
            CHECK_VAR_ERR(in_ncid, i)

        /* copy over all attributes */
        for (j=0; j<nattrs; j++) {
            err = ncmpi_inq_attname(in_ncid, i, j, name);
            CHECK_VAR_ERR(in_ncid, i)

            /* skip copying decomposition attributes */
            if (strcmp(name, "decomposition_ID") == 0 ||
                strcmp(name, "global_dimids") == 0)
                continue;

            err = ncmpi_copy_att(in_ncid, i, name, out_ncid, var[i].varid);
            CHECK_VAR_ERR(in_ncid, i)
        }
    }
    if (verbose && world_rank == 0)
        printf("nvars=%d (decomp var = %d not partitioned = %d partitioned = %d)\n",
               nvars, nvars_decomp, nvars-nvars_decomp-nvars_partitioned, nvars_partitioned);

err_out:
    return (err != 0);
}

/*----< print_info() >------------------------------------------------------*/
static
void print_info (MPI_Info *info_used) {
    int i, nkeys;

    MPI_Info_get_nkeys (*info_used, &nkeys);
    printf ("MPI File Info: nkeys = %d\n", nkeys);
    for (i = 0; i < nkeys; i++) {
        char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
        int valuelen, flag;

        MPI_Info_get_nthkey (*info_used, i, key);
        MPI_Info_get_valuelen (*info_used, key, &valuelen, &flag);
        MPI_Info_get (*info_used, key, valuelen + 1, value, &flag);
        printf ("MPI File Info: [%2d] key = %25s, value = %s\n", i, key, value);
    }
}

/*----< usage() >------------------------------------------------------------*/
static void usage (char *argv0) {
    char *help =
    "Usage: %s [OPTION]... FILE\n"
    "    [-h] Print help\n"
    "    [-v] Verbose mode\n"
    "    -i file  Base name of input subfiles\n"
    "    -o file  Output file name\n";
    fprintf (stderr, help, argv0);
}

/*----< main() >-------------------------------------------------------------*/
int main (int argc, char **argv)
{
    char in_file_base[1024], in_file[1040], out_file[1024];
    int i, j, err=0, in_ncid=-1, out_ncid=-1, nvars, rec;
    int num_subfiles, in_num_subfiles, subfile_ID, num_iputs, num_igets;
    int global_nprocs, num_decomp=0;
    size_t sum_vlen;
    char *buf=NULL, *buf_ptr;
    io_var *var=NULL;
    MPI_Offset start[3], count[3], num_recs=0;
    MPI_Info r_info=MPI_INFO_NULL, w_info=MPI_INFO_NULL;
    double open_t, def_t, close_t, total_t, mark_t, mark_t2;
    double read_post_t=0, read_wait_t=0, read_t;
    double write_post_t=0, write_wait_t=0, write_t;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_nprocs);

    verbose         = 0;
    in_file_base[0] = '\0';
    out_file[0]     = '\0';
    sub_comm        = MPI_COMM_NULL;
    read_t          = 0;
    write_t         = 0;

    /* command-line arguments */
    while ((i = getopt(argc, argv, "hvi:o:")) != EOF)
        switch(i) {
            case 'v': verbose = 1;
                      break;
            case 'i': strcpy(in_file_base, optarg);
                      break;
            case 'o': strcpy(out_file, optarg);
                      break;
            case 'h':
            default:  if (world_rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }

    if (in_file_base[0] == '\0') { /* input file name is mandatory */
        if (!world_rank) {
            fprintf(stderr, "Error: input file is missing\n");
            usage(argv[0]);
        }
        MPI_Finalize();
        return 1;
    }
    if (out_file[0] == '\0')
        strcpy(out_file, in_file_base);

    /* check output file and it should not exist */
    int file_exist= 0;
    if (world_rank == 0 && access(out_file, F_OK) == 0)
        file_exist = 1;
    MPI_Bcast(&file_exist, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (file_exist) {
        if (world_rank == 0)
            fprintf(stderr,"Error: output file already exists (%s)\n",out_file);
        MPI_Finalize();
        exit(1);
    }

    MPI_Barrier(MPI_COMM_WORLD); /*----------------------------------------*/
    total_t = open_t = MPI_Wtime();

    /* TODO: check the number of subfiles, N */
    /*       if the number of subfiles > nprocs
     *          assign N/nprocs files to each rank
     *       if the number of subfiles < nprocs
     *          assign 1 file to nprocs/N ranks
     *          if number of nblobs > sub_nprocs
     *             assign nblobs/sub_nprocs to each sub_rank
     *          else
     *             assign 1 blob to first nblobs ranks
     */


    /* split MPI_COMM_WORLD into sub communicators, one for each compute node */
    err = split_comm(MPI_COMM_WORLD, &sub_comm, &num_subfiles, &subfile_ID);
    if (err != 0) goto err_out;

    MPI_Comm_size(sub_comm, &sub_nprocs);
    MPI_Comm_rank(sub_comm, &sub_rank);

    /* construct input subfile names */
    sprintf(in_file, "%s.%04d", in_file_base, subfile_ID);
    if (verbose && sub_rank == 0)
        printf("world_rank=%5d sub_rank=%5d in_file=%s\n",
               world_rank, sub_rank, in_file);

    /* input subfile should exist */
    if (access(in_file, F_OK) == -1) {
        fprintf(stderr,"Error: input file does not exist %s\n",in_file);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    /* open input subfiles using sub_comm */
    err = ncmpi_open(sub_comm, in_file, NC_NOWRITE, MPI_INFO_NULL, &in_ncid);
    CHECK_NC_ERR

    /* inquire the MPI-IO hints actually used */
    err = ncmpi_inq_file_info(in_ncid, &r_info);
    CHECK_NC_ERR

    if (verbose && sub_rank == 0) {
        printf("MPI-IO hints used in reading file %s\n", in_file);
        print_info(&r_info);
    }

    /* collective write and no independent MPI-IO */
    MPI_Info_create(&w_info);
    MPI_Info_set(w_info, "romio_cb_write", "enable");
    MPI_Info_set(w_info, "romio_no_indep_rw", "true");

    /* create output file using MPI_COMM_WORLD */
    err = ncmpi_create(MPI_COMM_WORLD, out_file, NC_64BIT_DATA, w_info,
                       &out_ncid);
    CHECK_NC_ERR

    MPI_Info_free(&w_info);
    /* inquire the MPI-IO hints actually used */
    err = ncmpi_inq_file_info(in_ncid, &w_info);
    CHECK_NC_ERR

    if (verbose && sub_rank == 0) {
        printf("MPI-IO hints used in writing file %s\n", out_file);
        print_info(&w_info);
    }

    open_t = MPI_Wtime() - open_t;

    MPI_Barrier(MPI_COMM_WORLD); /*----------------------------------------*/
    def_t = MPI_Wtime();

    /* inquire global number of processes */
    err = ncmpi_get_att(in_ncid, NC_GLOBAL, "global_nprocs", &global_nprocs);
    CHECK_NC_ERR

    if (global_nprocs != world_nprocs && world_rank == 0)
        printf("Warning: no. processes (%d) is not equal to the one used when creating the files (%d)\n",
               world_nprocs, global_nprocs);

    /* inquire number of subfiles */
    err = ncmpi_get_att(in_ncid, NC_GLOBAL, "num_subfiles", &in_num_subfiles);
    CHECK_NC_ERR

    if (in_num_subfiles != num_subfiles && world_rank == 0)
        printf("Warning: no. compute nodes (%d) is not equal to no. available subfiles (%d)\n",
               num_subfiles, in_num_subfiles);

    /* read and construct decomposition metadata */
    io_decomp decomp[6];
    err = set_decomp(in_ncid, &num_decomp, decomp);
    if (err != 0) goto err_out;

    /* copy over all dimensions, except ones used for decomposition */
    err = copy_dims(in_ncid, out_ncid, &num_recs);
    if (err != 0) goto err_out;

    /* copy over all global attributes, except ones used for decomposition */
    err = copy_attr(in_ncid, out_ncid);
    if (err != 0) goto err_out;

    /* inquire number of variables */
    err = ncmpi_inq_nvars(in_ncid, &nvars);
    CHECK_NC_ERR

    if (verbose && sub_rank == 0)
        printf("world_rank=%5d sub_rank=%5d nvars=%d\n",
               world_rank, sub_rank, nvars);

    /* allocate space for variable metadata object */
    var = (io_var*) malloc(nvars * sizeof(io_var));

    /* read variable metadata from input file, define variables in output file,
     * and copy over their attributes
     */
    err = set_vars(in_ncid, out_ncid, nvars, num_decomp, decomp, var);
    if (err != 0) goto err_out;
    err = ncmpi_enddef(out_ncid);
    CHECK_NC_ERR

    def_t = MPI_Wtime() - def_t;

    MPI_Barrier(MPI_COMM_WORLD); /*----------------------------------------*/
    mark_t = MPI_Wtime();

    /* sum write size across all variables */
    sum_vlen = 0;
    for (i=0; i<nvars; i++)
        sum_vlen += var[i].vlen;

    buf = (char*) malloc(sum_vlen);

    num_igets = num_iputs = 0;

    /* read fixed-size variables */
    mark_t2 = MPI_Wtime();
    buf_ptr = buf;
    for (i=0; i<nvars; i++) {
        if (var[i].varid == -1 || var[i].is_rec) continue;
        if (var[i].dec_id == -1) { /* this variable is no partitioned */
            /* root process reads this variable */
            if (world_rank == 0) {
                err = ncmpi_iget_var(in_ncid, i, buf_ptr, 0,
                                     MPI_DATATYPE_NULL, NULL);
                CHECK_VAR_ERR(in_ncid, i)
                buf_ptr += var[i].vlen;
                num_igets++;
            }
            continue;
        }
        /* read partitioned variable */
        start[0] = decomp[var[i].dec_id].r_start;
        count[0] = decomp[var[i].dec_id].r_count;
        err = ncmpi_iget_vara(in_ncid, i, start, count, buf_ptr, 0,
                              MPI_DATATYPE_NULL, NULL);
        CHECK_VAR_ERR(in_ncid, i)
        buf_ptr += var[i].vlen;
        num_igets++;
    }
    read_post_t += MPI_Wtime() - mark_t2;

    mark_t2 = MPI_Wtime();
    err = ncmpi_wait_all(in_ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_NC_ERR
    read_wait_t += MPI_Wtime() - mark_t2;

    read_t += MPI_Wtime() - mark_t;
    mark_t  = MPI_Wtime();

    /* write fixed-size variables */
    mark_t2 = MPI_Wtime();
    buf_ptr = buf;
    for (i=0; i<nvars; i++) {
        if (var[i].varid == -1 || var[i].is_rec) continue;
        if (var[i].dec_id == -1) { /* this variable is no partitioned */
            /* root process writes this variable */
            if (world_rank == 0) {
                err = ncmpi_iput_var(out_ncid, var[i].varid, buf_ptr, 0,
                                     MPI_DATATYPE_NULL, NULL);
                CHECK_VAR_ERR(out_ncid, var[i].varid)
                buf_ptr += var[i].vlen;
                num_iputs++;
            }
            continue;
        }
        /* write partitioned variable */
        io_decomp *dp = decomp + var[i].dec_id;
        err = ncmpi_iput_varn(out_ncid, var[i].varid, dp->nreqs,
                              dp->w_startx, dp->w_countx,
                              buf_ptr, 0, MPI_DATATYPE_NULL, NULL);
        CHECK_VAR_ERR(out_ncid, var[i].varid)
        buf_ptr += var[i].vlen;
        num_iputs++;
    }
    write_post_t += MPI_Wtime() - mark_t2;

    mark_t2 = MPI_Wtime();
    err = ncmpi_wait_all(out_ncid, NC_REQ_ALL, NULL, NULL);
    CHECK_NC_ERR
    write_wait_t += MPI_Wtime() - mark_t2;

    write_t += MPI_Wtime() - mark_t;

    /* copy the record variables */

    /* loop over each record */
    for (rec=0; rec<num_recs; rec++) {
        mark_t = MPI_Wtime();
        start[0] = rec;
        count[0] = 1;

        /* read from input file */
        mark_t2 = MPI_Wtime();
        buf_ptr = buf;
        for (i=0; i<nvars; i++) {
            if (var[i].varid == -1 || var[i].is_rec == 0) continue;
            if (var[i].dec_id == -1) { /* this variable is no partitioned */
                /* root process reads this variable */
                if (world_rank == 0) {
                    start[1] = 0;
                    err = ncmpi_iget_vara(in_ncid, i, start, var[i].dims,
                                          buf_ptr, 0, MPI_DATATYPE_NULL, NULL);
                    CHECK_VAR_ERR(in_ncid, i)
                    buf_ptr += var[i].vlen;
                    num_igets++;
                }
                continue;
            }
            /* read partitioned variable */
            start[1] = decomp[var[i].dec_id].r_start;
            count[1] = decomp[var[i].dec_id].r_count;
            err = ncmpi_iget_vara(in_ncid, i, start, count, buf_ptr, 0,
                                  MPI_DATATYPE_NULL, NULL);
            CHECK_VAR_ERR(in_ncid, i)
            buf_ptr += var[i].vlen;
            num_igets++;
        }
        read_post_t += MPI_Wtime() - mark_t2;

        mark_t2 = MPI_Wtime();
        err = ncmpi_wait_all(in_ncid, NC_REQ_ALL, NULL, NULL);
        CHECK_NC_ERR
        read_wait_t += MPI_Wtime() - mark_t2;

        read_t += MPI_Wtime() - mark_t;
        mark_t  = MPI_Wtime();

        /* write to output file */
        mark_t2 = MPI_Wtime();
        buf_ptr = buf;
        for (i=0; i<nvars; i++) {
            if (var[i].varid == -1 || var[i].is_rec == 0) continue;
            if (var[i].dec_id == -1) { /* this variable is no partitioned */
                /* root process writes this variable */
                if (world_rank == 0) {
                    start[1] = 0;
                    err = ncmpi_iput_vara(out_ncid, var[i].varid, start,
                                          var[i].dims, buf_ptr, 0,
                                          MPI_DATATYPE_NULL, NULL);
                    CHECK_VAR_ERR(out_ncid, var[i].varid)
                    buf_ptr += var[i].vlen;
                    num_iputs++;
                }
                continue;
            }
            /* set all starts[][0] to next record */
            for (j=0; j<decomp[var[i].dec_id].nreqs; j++)
                decomp[var[i].dec_id].w_starts[j][0] = rec;

            /* write partitioned variable */
            io_decomp *dp = decomp + var[i].dec_id;
            err = ncmpi_iput_varn(out_ncid, var[i].varid, dp->nreqs,
                                  dp->w_starts, dp->w_counts,
                                  buf_ptr, 0, MPI_DATATYPE_NULL, NULL);
            CHECK_VAR_ERR(out_ncid, var[i].varid)
            buf_ptr += var[i].vlen;
            num_iputs++;
        }

        write_post_t += MPI_Wtime() - mark_t2;

        mark_t2 = MPI_Wtime();
        err = ncmpi_wait_all(out_ncid, NC_REQ_ALL, NULL, NULL);
        CHECK_NC_ERR
        write_wait_t += MPI_Wtime() - mark_t2;

        write_t += MPI_Wtime() - mark_t;
    }

    MPI_Barrier(MPI_COMM_WORLD); /*----------------------------------------*/
    close_t  = MPI_Wtime();

    /* collect read and write amount */
    MPI_Offset read_amnt, write_amnt, amnt, sum_buf[3], off_buf[3];
    read_amnt = write_amnt = 0;

    err = ncmpi_inq_put_size(in_ncid, &amnt);
    CHECK_NC_ERR
    write_amnt += amnt;

    err = ncmpi_inq_get_size(in_ncid, &amnt);
    CHECK_NC_ERR
    read_amnt += amnt;

    err = ncmpi_inq_put_size(out_ncid, &amnt);
    CHECK_NC_ERR
    write_amnt += amnt;

    err = ncmpi_inq_get_size(out_ncid, &amnt);
    CHECK_NC_ERR
    read_amnt += amnt;

    /* close input and output files */
    err = ncmpi_close(in_ncid);
    CHECK_NC_ERR

    err = ncmpi_close(out_ncid);
    CHECK_NC_ERR

    close_t = MPI_Wtime() - close_t;
    total_t = MPI_Wtime() - total_t;

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset m_alloc, max_alloc;

    err = ncmpi_inq_malloc_size(&m_alloc);
    if (err == NC_ENOTENABLED) m_alloc = 0;

    off_buf[0] = read_amnt;
    off_buf[1] = write_amnt;
    off_buf[2] = m_alloc;
    MPI_Reduce(off_buf, &sum_buf, 3, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
    read_amnt  = sum_buf[0];
    write_amnt = sum_buf[1];
    m_alloc    = sum_buf[2];

    if (world_rank == 0 && m_alloc > 0) {
        printf("-------------------------------------------------------\n");
        printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                m_alloc);
    }
    err = ncmpi_inq_malloc_max_size(&m_alloc);
    if (err == NC_ENOTENABLED) m_alloc = 0;
    MPI_Reduce(&m_alloc, &max_alloc, 1, MPI_OFFSET, MPI_MAX, 0, MPI_COMM_WORLD);
    if (err == NC_ENOTENABLED) err = NC_NOERR;

    /* find the max timings amount all processes */
    double timings[10], max_time[10];
    timings[0] = open_t;
    timings[1] = def_t;
    timings[2] = read_t;
    timings[3] = write_t;
    timings[4] = close_t;
    timings[5] = total_t;
    timings[6] = read_post_t;
    timings[7] = read_wait_t;
    timings[8] = write_post_t;
    timings[9] = write_wait_t;
    MPI_Reduce(timings, max_time, 10, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        printf("Input subfile base name            = %s\n", in_file_base);
        printf("Output file name                   = %s\n", out_file);
        printf("No. decompositions                 = %3d\n", num_decomp);
        printf("No. variables                      = %3d\n", nvars);
        printf("No. records (time steps)           = %3lld\n", num_recs);
        printf("-----------------------------------------------------------\n");
        printf("Max Time of file open/create       = %.4f sec\n", max_time[0]);
        printf("Max Time of define variables       = %.4f sec\n", max_time[1]);
        printf("-----------------------------------------------\n");
        printf("Max Time of read post              = %.4f sec\n", max_time[6]);
        printf("Max Time of read wait              = %.4f sec\n", max_time[7]);
        printf("Max Time of read                   = %.4f sec\n", max_time[2]);
        printf("-----------------------------------------------\n");
        printf("Max Time of write post             = %.4f sec\n", max_time[8]);
        printf("Max Time of write wait             = %.4f sec\n", max_time[9]);
        printf("Max Time of write                  = %.4f sec\n", max_time[3]);
        printf("-----------------------------------------------\n");
        printf("Max Time of close                  = %.4f sec\n", max_time[4]);
        printf("Max end-to-end time                = %.4f sec\n", max_time[5]);
        printf("-----------------------------------------------------------\n");
        printf("Total read  amount                 = %.2f MiB = %.2f GiB\n",
               (double)read_amnt / 1048576, (double)read_amnt / 1073741824);
        printf("Read  bandwidth                    = %.2f MiB/sec\n",
               (double)read_amnt / 1048576 / max_time[2]);
        printf("Total write amount                 = %.2f MiB = %.2f GiB\n",
               (double)write_amnt / 1048576, (double)write_amnt / 1073741824);
        printf("Write bandwidth                    = %.2f MiB/sec\n",
               (double)write_amnt / 1048576 / max_time[3]);
        if (verbose && max_alloc > 0)
            printf("MAX heap memory used by PnetCDF internally is %.2f MiB\n",
                   (float)max_alloc / 1048576);
        printf("-----------------------------------------------------------\n");
    }

    if (verbose)
        printf("%2d: num iput calls=%d num iget calls=%d\n",
               world_rank,num_iputs,num_igets);

err_out:
    if (buf != NULL) free(buf);
    for (i=0; i<num_decomp; i++) {
        if (decomp[i].offsets != NULL) free(decomp[i].offsets);
        if (decomp[i].w_starts != NULL) {
            free(decomp[i].w_starts[0]);
            free(decomp[i].w_starts);
        }
    }
    if (var != NULL) free(var);

    if (r_info != MPI_INFO_NULL) MPI_Info_free(&r_info);
    if (w_info != MPI_INFO_NULL) MPI_Info_free(&w_info);

    if (sub_comm != MPI_COMM_NULL) MPI_Comm_free(&sub_comm);
    if (in_ncid != -1) ncmpi_close(in_ncid);
    if (out_ncid != -1) ncmpi_close(out_ncid);

    MPI_Finalize();

    return (err != 0);
}
