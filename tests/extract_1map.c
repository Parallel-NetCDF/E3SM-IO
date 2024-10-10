/*********************************************************************
 *
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <mpi.h>
#include <pnetcdf.h>

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}


int main(int argc, char** argv)
{
    char *in_file, *out_file;
    int rank, nprocs, err, nerrs=0;
    int in_ncid, out_ncid, cmode, dimid, dimids[2], varid, out_varid[3];
    void *buf;
    MPI_Offset decomp_nprocs, total_nreqs;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    if (argc != 3) {
        printf("Usage: %s in_filename, out_filename\n",argv[0]);
        MPI_Finalize();
        exit(1);
    }

    /* open input file for reading ------------------------------------------*/
    err = ncmpi_open(comm, argv[1], NC_NOWRITE, MPI_INFO_NULL, &in_ncid); ERR

    /* create a new file for writing ----------------------------------------*/
    cmode = NC_64BIT_DATA | NC_CLOBBER;
    err = ncmpi_create(comm, argv[2], cmode, MPI_INFO_NULL, &out_ncid); ERR

    /* copy all global attributes */
    err = ncmpi_copy_att (in_ncid, NC_GLOBAL, "D3.ndims", out_ncid, NC_GLOBAL); ERR
    err = ncmpi_copy_att (in_ncid, NC_GLOBAL, "D3.dims", out_ncid, NC_GLOBAL); ERR
    err = ncmpi_copy_att (in_ncid, NC_GLOBAL, "D3.max_nreqs", out_ncid, NC_GLOBAL); ERR
    err = ncmpi_copy_att (in_ncid, NC_GLOBAL, "D3.min_nreqs", out_ncid, NC_GLOBAL); ERR
    err = ncmpi_copy_att (in_ncid, NC_GLOBAL, "D3.fully_covered", out_ncid, NC_GLOBAL); ERR

    /* copy dimensions */
    err = ncmpi_inq_dimid(in_ncid, "decomp_nprocs", &dimid); ERR
    err = ncmpi_inq_dimlen(in_ncid, dimid, &decomp_nprocs); ERR
    err = ncmpi_inq_dimid(in_ncid, "D3.total_nreqs", &dimid); ERR
    err = ncmpi_inq_dimlen(in_ncid, dimid, &total_nreqs); ERR

    err = ncmpi_def_dim(out_ncid, "decomp_nprocs", decomp_nprocs, &dimids[0]); ERR
    err = ncmpi_def_dim(out_ncid, "D3.total_nreqs", total_nreqs, &dimids[1]); ERR

    /* define variables */
    err = ncmpi_def_var(out_ncid, "D3.nreqs", NC_INT, 1, &dimids[0], &out_varid[0]); ERR
    err = ncmpi_def_var(out_ncid, "D3.offsets", NC_INT, 1, &dimids[1], &out_varid[1]); ERR
    err = ncmpi_def_var(out_ncid, "D3.lengths", NC_INT, 1, &dimids[1], &out_varid[2]); ERR
    err = ncmpi_enddef(out_ncid); ERR

    /* copy variables */
    err = ncmpi_inq_varid(in_ncid, "D3.nreqs", &varid); ERR
    buf = malloc(sizeof(int) * decomp_nprocs);
    err = ncmpi_get_var_int_all(in_ncid, varid, buf);
    err = ncmpi_put_var_int_all(out_ncid, out_varid[0], buf);

    err = ncmpi_inq_varid(in_ncid, "D3.offsets", &varid); ERR
    buf = realloc(buf, sizeof(int) * total_nreqs);
    err = ncmpi_get_var_int_all(in_ncid, varid, buf);
    err = ncmpi_put_var_int_all(out_ncid, out_varid[1], buf);

    err = ncmpi_inq_varid(in_ncid, "D3.lengths", &varid); ERR
    err = ncmpi_get_var_int_all(in_ncid, varid, buf);
    err = ncmpi_put_var_int_all(out_ncid, out_varid[2], buf);

    free(buf);

    err = ncmpi_close(in_ncid); ERR
    err = ncmpi_close(out_ncid); ERR

err_out:
    MPI_Finalize();
    return (nerrs > 0);
}


