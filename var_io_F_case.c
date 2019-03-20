/*********************************************************************
 *
 * Copyright (C) 2018, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program uses the E3SM I/O patterns recorded by the PIO library to
 * evaluate the performance of two PnetCDF APIs: ncmpi_vard_all(), and
 * ncmpi_iput_varn(). The E3SM I/O patterns consist of a large number of small,
 * noncontiguous requests on each MPI process, which presents a challenge for
 * achieving a good performance.
 *
 * See README.md for compile and run instructions.
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* unlink() */

#include <e3sm_io.h>

/*----< write_small_vars_F_case() >------------------------------------------*/
static int
write_small_vars_F_case(int          ncid,
                        int          vid,    /* starting variable ID */
                        int         *varids,
                        int          rec_no,
                        int          gap,
                        MPI_Offset   lev,
                        MPI_Offset   ilev,
                        MPI_Offset   nbnd,
                        MPI_Offset   nchars,
                        int        **int_buf,
                        char       **txt_buf,
                        double     **dbl_buf)
{
    int i, err, nerrs=0;
    MPI_Offset start[2], count[2];

    /* scalar and small variables are written by rank 0 only */
    i = vid;

    if (rec_no == 0) {
        /* lev */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += lev + gap;
        /* hyam */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += lev + gap;
        /* hybm */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += lev + gap;
        /* P0 */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += 1 + gap;
        /* ilev */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += ilev + gap;
        /* hyai */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += ilev + gap;
        /* hybi */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += ilev + gap;
    }
    else
        i += 7;

    /* time */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* date */
    start[0] = rec_no;
    err = ncmpi_iput_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
    /* datesec */
    start[0] = rec_no;
    err = ncmpi_iput_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
    /* time_bnds */
    start[0] = rec_no; start[1] = 0;
    count[0] = 1;      count[1] = nbnd;
    err = ncmpi_iput_vara_double(ncid, varids[i++], start, count, *dbl_buf, NULL); ERR
    *dbl_buf += nbnd + gap;
    /* date_written */
    start[0] = rec_no; start[1] = 0;
    count[0] = 1;      count[1] = nchars;
    err = ncmpi_iput_vara_text(ncid, varids[i++], start, count, *txt_buf, NULL); ERR
    *txt_buf += nchars;
    /* time_written */
    start[0] = rec_no; start[1] = 0;
    count[0] = 1;      count[1] = nchars;
    err = ncmpi_iput_vara_text(ncid, varids[i++], start, count, *txt_buf, NULL); ERR
    *txt_buf += nchars;

    if (rec_no == 0) {
        /* ndbase */
        err = ncmpi_iput_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
        /* nsbase */
        err = ncmpi_iput_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
        /* nbdate */
        err = ncmpi_iput_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
        /* nbsec */
        err = ncmpi_iput_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
        /* mdt */
        err = ncmpi_iput_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
    }
    else
        i += 5;

    /* ndcur */
    start[0] = rec_no;
    err = ncmpi_iput_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
    /* nscur */
    start[0] = rec_no;
    err = ncmpi_iput_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
    /* co2vmr */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* ch4vmr */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* n2ovmr */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* f11vmr */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* f12vmr */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* sol_tsi */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* nsteph */
    start[0] = rec_no;
    err = ncmpi_iput_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
fn_exit:
    return err;
}

#define SET_TYPE(kind) { \
    var_types[i] = type[kind]; \
    err = ncmpi_inq_varoffset(ncid, varids[i], &var_offset); ERR \
    var_disps[i] = var_offset - offset_rec; \
    if (kind == 2) { \
        my_nreqs += nreqs[1]; \
        if (i < nvars-1) buf_disps[i+1] = buf_disps[i] + (nelems[1]+gap) * sizeof(dtype); \
        buf_blocklens[i] = nelems[1]; \
    } else { /* kind == 3 */ \
        my_nreqs += nreqs[2]; \
        if (i < nvars-1) buf_disps[i+1] = buf_disps[i] + (nelems[2]+gap) * sizeof(dtype); \
        buf_blocklens[i] = nelems[2]; \
    } \
    i++; \
}

#define SET_TYPES(kind, num) \
    for (j=0; j<num; j++) { \
        var_types[i] = type[kind]; \
        err = ncmpi_inq_varoffset(ncid, varids[i], &var_offset); ERR \
        var_disps[i] = var_offset - offset_rec; \
        if (kind == 2) { \
            my_nreqs += nreqs[1]; \
            if (i < nvars-1) buf_disps[i+1] = buf_disps[i] + (nelems[1]+gap) * sizeof(dtype); \
            buf_blocklens[i] = nelems[1]; \
        } else { /* kind == 3 */ \
            my_nreqs += nreqs[2]; \
            if (i < nvars-1) buf_disps[i+1] = buf_disps[i] + (nelems[2]+gap) * sizeof(dtype); \
            buf_blocklens[i] = nelems[2]; \
        } \
        i++; \
    }

/*----< run_vard_F_case() >--------------------------------------------------*/
int
run_vard_F_case(char       *out_dir,      /* output folder name */
                char       *outfile,      /* output file name */
                int         nvars,        /* number of variables 408 or 51 */
                int         num_recs,     /* number of records */
                int         noncontig_buf,/* whether to us noncontiguous buffer */
                MPI_Info    info,
                MPI_Offset  dims[3][2],   /* dimension lengths */
                int         nreqs[3],     /* no. request in decompositions 1,2,3 */
                int        *disps[3],     /* request's displacements */
                int        *blocklens[3]) /* request's block lengths */
{
    char outfname[512], txt_buf[16], *txt_buf_ptr;
    int i, j, k, err, nerrs=0, rank, ncid, cmode, *varids;
    int *var_blocklens, *buf_blocklens, my_nreqs, max_nreqs, rec_no, gap=0;
    int int_buf[10], *int_buf_ptr;
    size_t fix_buflen, dbl_buflen, rec_buflen;
    size_t nelems[3];
    dtype *rec_buf;
    double *dbl_buf, *dbl_buf_ptr;
    double pre_timing, open_timing, io_timing, close_timing;
    double timing, total_timing,  max_timing;
    MPI_Aint *var_disps, *buf_disps;
    MPI_Offset tmp, metadata_size, rec_size, put_size, total_size;
    MPI_Offset offset_fix, offset_rec, var_offset;
    MPI_Datatype *var_types, type[4], *filetype_rec, filetype_dbl;
    MPI_Datatype buftype_rec, buftype_dbl;
    MPI_Info info_used=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Barrier(comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime();

    MPI_Comm_rank(comm, &rank);

    if (noncontig_buf) gap = 10;

    varids = (int*) malloc(nvars * sizeof(int));

    /* allocate arrays for constructing fileview and buffer type */
    var_types = (MPI_Datatype*) malloc(nvars * sizeof(MPI_Datatype));
    var_blocklens = (int*) malloc(nvars * 2 * sizeof(int));
    buf_blocklens = var_blocklens + nvars;
    var_disps = (MPI_Aint*) malloc(nvars * 2 * sizeof(MPI_Aint));
    buf_disps = var_disps + nvars;

    /* define MPI datatypes for 4 kinds from 3 decompositions */
    MPI_Type_indexed(nreqs[0], blocklens[0], disps[0], MPI_DOUBLE, &type[0]);
    MPI_Type_commit(&type[0]);
    MPI_Type_indexed(nreqs[1], blocklens[1], disps[1], MPI_DOUBLE, &type[1]);
    MPI_Type_commit(&type[1]);
    MPI_Type_indexed(nreqs[1], blocklens[1], disps[1], REC_DTYPE, &type[2]);
    MPI_Type_commit(&type[2]);
    MPI_Type_indexed(nreqs[2], blocklens[2], disps[2], REC_DTYPE, &type[3]);
    MPI_Type_commit(&type[3]);

    /* number of variable elements from 3 decompositions */
    for (j=0; j<nvars; j++) var_blocklens[j] = 1;
    nelems[0] = nelems[1] = nelems[2] = 0;
    for (k=0; k<nreqs[0]; k++) nelems[0] += blocklens[0][k];
    for (k=0; k<nreqs[1]; k++) nelems[1] += blocklens[1][k];
    for (k=0; k<nreqs[2]; k++) nelems[2] += blocklens[2][k];

    if (verbose && rank == 0)
        printf("nelems=%zd %zd %zd\n", nelems[0],nelems[1],nelems[2]);

    /* allocate and initialize write buffer for small variables */
    dbl_buflen = nelems[1]*2 + nelems[0] + gap*3
               + 3 * dims[2][0] + 3 * (dims[2][0]+1) + 8 + 2 + 20 * gap;

    dbl_buf = (double*) malloc(dbl_buflen * sizeof(double));
    for (j=0; j<dbl_buflen; j++) dbl_buf[j] = rank;

    for (j=0; j<10; j++) int_buf[j] = rank;
    for (j=0; j<16; j++) txt_buf[j] = 'a' + rank;

    /* allocate and initialize write buffer for large variables */
    if (nvars == 408)
        rec_buflen = nelems[1] * 315 + nelems[2] * 63 + (315+63) * gap;
    else
        rec_buflen = nelems[1] * 20 + nelems[2] + (20+1) * gap;

    rec_buf = (dtype*) malloc(rec_buflen * sizeof(dtype));
    for (i=0; i<rec_buflen; i++) rec_buf[i] = rank;

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    open_timing = MPI_Wtime();

    /* set output file name */
    sprintf(outfname, "%s/%s",out_dir, outfile);

    /* create a new CDF-5 file for writing */
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(comm, outfname, cmode, info, &ncid); ERR

    /* define dimensions, variables, and attributes */
    if (nvars == 408) {
        /* for h0 file */
        err = def_F_case_h0(ncid, dims[2], nvars, varids); ERR
    }
    else {
        /* for h1 file */
        err = def_F_case_h1(ncid, dims[2], nvars, varids); ERR
    }

    /* exit define mode and enter data mode */
    err = ncmpi_enddef(ncid); ERR

    /* I/O amount so far */
    err = ncmpi_inq_put_size(ncid, &metadata_size); ERR
    err = ncmpi_inq_file_info(ncid, &info_used); ERR
    open_timing = MPI_Wtime() - open_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* the first 3 variables are of type NC_DOUBLE -------------------*/
    i = 0;
    my_nreqs = 0;

    /* lat */
    var_types[i] = type[1];
    err = ncmpi_inq_varoffset(ncid, varids[i], &offset_fix); ERR
    var_disps[i] = 0;
    buf_disps[0] = 0;
    buf_blocklens[0] = nelems[1];
    i++;
    my_nreqs += nreqs[1];

    /* lon */
    var_types[i] = type[1];
    err = ncmpi_inq_varoffset(ncid, varids[i], &var_offset); ERR
    var_disps[i] = var_offset - offset_fix;
    buf_disps[i] = buf_disps[i-1] + (nelems[1] + gap) * sizeof (double);
    buf_blocklens[1] = nelems[1];
    i++;
    my_nreqs += nreqs[1];

    /* area */
    var_types[i] = type[0];
    err = ncmpi_inq_varoffset(ncid, varids[i], &var_offset); ERR
    var_disps[i] = var_offset - offset_fix;
    buf_disps[i] = buf_disps[i-1] + (nelems[1] + gap) * sizeof (double);
    buf_blocklens[2] = nelems[0];
    i++;
    my_nreqs += nreqs[0];
    fix_buflen = nelems[1]*2 + nelems[0] + gap*3;

    /* skip next 27 small variables */
    i += 27;

    /* concatenate 3 var_types[] into filetype_dbl */
    MPI_Type_create_struct(3, var_blocklens, var_disps, var_types,
                           &filetype_dbl);
    MPI_Type_commit(&filetype_dbl);

    if (noncontig_buf) {
        /* construct buffer type for 3 variables */
        MPI_Type_create_hindexed(3, buf_blocklens, buf_disps, MPI_DOUBLE,
                                 &buftype_dbl);
        MPI_Type_commit(&buftype_dbl);
    }
    else {
        /* buffer type is contiguous */
        buftype_dbl = MPI_DOUBLE;
    }

    err = ncmpi_inq_varoffset(ncid, varids[i], &offset_rec); ERR
    err = ncmpi_inq_recsize(ncid, &rec_size); ERR
    buf_disps[i] = 0;

    if (nvars == 408) {
        SET_TYPE(2)       /* AEROD_v */
        SET_TYPES(3, 2)   /* ANRAIN and ANSNOW */
        SET_TYPES(2, 18)  /* AODABS ... ANSNOW */
        SET_TYPES(3, 2)   /* AQRAIN and AQSNOW */
        SET_TYPES(2, 6)   /* AQ_DMS ... AQ_SOAG */
        SET_TYPES(3, 5)   /* AREI ... CCN3 */
        SET_TYPES(2, 2)   /* CDNUMC and CLDHGH */
        SET_TYPES(3, 2)   /* CLDICE and CLDLIQ */
        SET_TYPES(2, 3)   /* CLDLOW ... CLDTOT */
        SET_TYPES(3, 4)   /* CLOUD ... DCQ */
        SET_TYPES(2, 11)  /* DF_DMS ... DSTSFMBL */
        SET_TYPE(3)       /* DTCOND */
        SET_TYPES(2, 2)   /* DTENDTH and DTENDTQ */
        SET_TYPES(3, 2)   /* EXTINCT and FICE */
        SET_TYPES(2, 7)   /* FLDS ... FLUTC */
        SET_TYPES(3, 4)   /* FREQI ... FREQS */
        SET_TYPES(2, 15)  /* FSDS ... ICEFRAC */
        SET_TYPES(3, 3)   /* ICIMR ... IWC */
        SET_TYPES(2, 2)   /* LANDFRAC and LHFLX */
        SET_TYPES(3, 5)   /* LINOZ_DO3 ... LINOZ_SSO3 */
        SET_TYPES(2, 3)   /* LINOZ_SZA ... LWCF */
        SET_TYPES(3, 12)  /* Mass_bc ... O3 */
        SET_TYPES(2, 2)   /* O3_SRF and OCNFRAC */
        SET_TYPE(3)       /* OMEGA */
        SET_TYPE(2)       /* OMEGA500 */
        SET_TYPE(3)       /* OMEGAT */
        SET_TYPES(2, 8)   /* PBLH ... PSL */
        SET_TYPE(3)       /* Q */
        SET_TYPES(2, 2)   /* QFLX and QREFHT */
        SET_TYPES(3, 3)   /* QRL ... RAINQM */
        SET_TYPE(2)       /* RAM1 */
        SET_TYPE(3)       /* RELHUM */
        SET_TYPES(2, 37)  /* SFDMS ... SNOWHLND */
        SET_TYPES(3, 2)   /* SNOWQM and SO2 */
        SET_TYPES(2, 10)  /* SO2_CLXF ... SWCF */
        SET_TYPE(3)       /* T */
        SET_TYPES(2, 19)  /* TAUGWX ... TVQ */
        SET_TYPE(3)       /* U */
        SET_TYPE(2)       /* U10 */
        SET_TYPES(3, 6)   /* UU ... VV */
        SET_TYPES(2, 3)   /* WD_H2O2 ... WD_SO2 */
        SET_TYPES(3, 3)   /* WSUB ... aero_water */
        SET_TYPES(2, 32)  /* airFV ... dst_c3SFWET */
        SET_TYPE(3)       /* hstobie_linoz */
        SET_TYPES(2, 129) /* mlip ... soa_c3SFWET */
    }
    else {
        SET_TYPES(2, 13)  /* CLDHGH ... T5 */
        SET_TYPE(3)       /* U */
        SET_TYPES(2, 7)   /* U250 ... Z500 */
    }

    if (noncontig_buf) {
        /* construct buffer type for record variables */
        MPI_Type_create_hindexed(nvars-30, buf_blocklens+30, buf_disps+30,
                                 REC_DTYPE, &buftype_rec);
        MPI_Type_commit(&buftype_rec);
    }
    else {
        /* all record variables are in a single contiguous buffer */
        buftype_rec = REC_DTYPE;
    }

    filetype_rec = (MPI_Datatype*)malloc(num_recs * sizeof(MPI_Datatype));
    for (j=0; j<num_recs; j++) {
        if (j > 0) {
            for (k=30; k<nvars; k++)
                var_disps[k] += rec_size;
        }
        /* concatenate nvars-30 var_types[] into filetype_rec[j] */
        MPI_Type_create_struct(nvars-30, var_blocklens+30, var_disps+30,
                               var_types+30, filetype_rec+j);
        MPI_Type_commit(filetype_rec+j);
    }

    for (j=0; j<4; j++) MPI_Type_free(&type[j]);
    free(var_types);
    free(var_disps);
    free(var_blocklens);

    pre_timing += MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    io_timing = MPI_Wtime();

    if (noncontig_buf) fix_buflen = rec_buflen = 1;
    else {
        if (nvars == 408)
            rec_buflen = nelems[1] * 315 + nelems[2] * 63;
        else
            rec_buflen = nelems[1] * 20 + nelems[2];
    }

    /* write first 3 NC_DOUBLE fixed-size variables in one vard call */
    err = ncmpi_put_vard_all(ncid, varids[0], filetype_dbl, dbl_buf,
                             fix_buflen, buftype_dbl); ERR

    for (rec_no=0; rec_no<num_recs; rec_no++) {
        i=3;
        dbl_buf_ptr = dbl_buf + nelems[1]*2 + nelems[0] + gap*3;
        int_buf_ptr = int_buf;
        txt_buf_ptr = txt_buf;

        /* next 27 small variables are written by rank 0 only */
        if (rank == 0) {
            /* post nonblocking requests using ncmpi_iput_varn() */
            err = write_small_vars_F_case(ncid, i, varids, rec_no, gap,
                                          dims[2][0], dims[2][0]+1, 2, 8,
                                          &int_buf_ptr, &txt_buf_ptr,
                                          &dbl_buf_ptr);
            ERR
            my_nreqs += 27;
        }
        i += 27;

        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR

        /* write remaining record variables in one vard call */
        err = ncmpi_put_vard_all(ncid, varids[30], filetype_rec[rec_no],
                                 rec_buf, rec_buflen, buftype_rec); ERR
    }
    io_timing = MPI_Wtime() - io_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    close_timing = MPI_Wtime();

    err = ncmpi_inq_put_size(ncid, &total_size); ERR
    put_size = total_size - metadata_size;
    err = ncmpi_close(ncid); ERR
    close_timing = MPI_Wtime() - close_timing;

    for (j=0; j<num_recs; j++) MPI_Type_free(filetype_rec+j);
    free(filetype_rec);
    MPI_Type_free(&filetype_dbl);

    if (noncontig_buf) {
        MPI_Type_free(&buftype_rec);
        MPI_Type_free(&buftype_dbl);
    }

    free(rec_buf);
    free(dbl_buf);
    free(varids);

    total_timing = MPI_Wtime() - total_timing;

    MPI_Reduce(&my_nreqs,      &max_nreqs,  1, MPI_INT,    MPI_MAX, 0, comm);
    MPI_Reduce(&put_size,      &tmp,        1, MPI_OFFSET, MPI_SUM, 0, comm);
    put_size = tmp;
    MPI_Reduce(&total_size,    &tmp,        1, MPI_OFFSET, MPI_SUM, 0, comm);
    total_size = tmp;
    MPI_Reduce(&open_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    open_timing = max_timing;
    MPI_Reduce(&pre_timing,    &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    pre_timing = max_timing;
    MPI_Reduce(&io_timing,     &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    io_timing = max_timing;
    MPI_Reduce(&close_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    close_timing = max_timing;
    MPI_Reduce(&total_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    total_timing = max_timing;

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);
        if (rank == 0 && sum_size > 0) {
            printf("-----------------------------------------------------------\n");
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        }
    }
    MPI_Offset m_alloc=0, max_alloc;
    ncmpi_inq_malloc_max_size(&m_alloc);
    MPI_Reduce(&m_alloc, &max_alloc, 1, MPI_OFFSET, MPI_MAX, 0, comm);
    if (rank == 0) {
        printf("History output file                = %s\n", outfile);
        printf("MAX heap memory allocated by PnetCDF internally is %.2f MiB\n",
               (float)max_alloc/1048576);
        printf("Total number of variables          = %d\n",nvars);
        printf("Total write amount                 = %.2f MiB = %.2f GiB\n",
               (double)total_size/1048576,(double)total_size/1073741824);
        printf("Max number of requests             = %d\n",max_nreqs);
        printf("Max Time of open + metadata define = %.4f sec\n",open_timing);
        printf("Max Time of I/O preparing          = %.4f sec\n",pre_timing);
        printf("Max Time of ncmpi_put_vard         = %.4f sec\n",io_timing);
        printf("Max Time of close                  = %.4f sec\n",close_timing);
        printf("Max Time of TOTAL                  = %.4f sec\n",total_timing);
        printf("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
               (double)total_size/1048576.0/total_timing);
        printf("I/O bandwidth (write-only)         = %.4f MiB/sec\n",
               (double)put_size/1048576.0/io_timing);
        if (verbose) print_info(&info_used);
        printf("-----------------------------------------------------------\n");
    }
fn_exit:
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    if (!keep_outfile) unlink(outfname);
    fflush(stdout);
    MPI_Barrier(comm);
    return nerrs;
}

#define FIX_1D_VAR_STARTS_COUNTS(starts, counts, nreqs, disps, blocklens) { \
    starts = (MPI_Offset**) malloc(2 * nreqs * sizeof(MPI_Offset*)); \
    counts = starts + nreqs; \
    starts[0] = (MPI_Offset*) malloc(2 * nreqs * sizeof(MPI_Offset)); \
    counts[0] = starts[0] + nreqs; \
    \
    for (i=1; i<nreqs; i++) { \
        starts[i] = starts[i-1] + 1; \
        counts[i] = counts[i-1] + 1; \
    } \
    \
    for (i=0; i<nreqs; i++) { \
        starts[i][0] = disps[i]; \
        counts[i][0] = blocklens[i]; \
    } \
}

#define FIX_2D_VAR_STARTS_COUNTS(starts, counts, nreqs, disps, blocklens, last_dimlen, nreqs_merged) { \
    starts = (MPI_Offset**) malloc(2 * nreqs * sizeof(MPI_Offset*)); \
    counts = starts + nreqs; \
    starts[0] = (MPI_Offset*) malloc(2 * nreqs * 2 * sizeof(MPI_Offset)); \
    counts[0] = starts[0] + nreqs * 2; \
    \
    for (i=1; i<nreqs; i++) { \
        starts[i] = starts[i-1] + 2; \
        counts[i] = counts[i-1] + 2; \
    } \
    \
    int reqs_cnt = 0; \
    i = 0; \
    while (i < nreqs) { \
        starts[reqs_cnt][1] = disps[i] % last_dimlen; /* decomposition is 2D */ \
        counts[reqs_cnt][1] = blocklens[i]; /* each blocklens[i] is no bigger than last_dimlen */ \
        \
        starts[reqs_cnt][0] = disps[i] / last_dimlen; \
        counts[reqs_cnt][0] = 1; \
        \
        /* merge consecutive "1 by last_dimlen" regions into a larger one */ \
        while ((i < nreqs - 1) && \
               (disps[i] % last_dimlen == 0 && blocklens[i] == last_dimlen) && \
               (disps[i + 1] - disps[i] == last_dimlen && blocklens[i + 1] == last_dimlen)) { \
            i++; \
            counts[reqs_cnt][0]++; \
        } \
        \
        reqs_cnt++; \
        i++; \
    } \
    \
    nreqs_merged = reqs_cnt; \
}

#define REC_2D_VAR_STARTS_COUNTS(rec, starts, counts, nreqs, disps, blocklens) { \
    starts = (MPI_Offset**) malloc(2 * nreqs * sizeof(MPI_Offset*)); \
    counts = starts + nreqs; \
    starts[0] = (MPI_Offset*) malloc(2 * nreqs * 2 * sizeof(MPI_Offset)); \
    counts[0] = starts[0] + nreqs * 2; \
    \
    for (i=1; i<nreqs; i++) { \
        starts[i] = starts[i-1] + 2; \
        counts[i] = counts[i-1] + 2; \
    } \
    \
    for (i=0; i<nreqs; i++) { \
        starts[i][1] = disps[i]; /* decomposition is 1D */ \
        counts[i][1] = blocklens[i]; \
        \
        starts[i][0] = rec; /* record ID */ \
        counts[i][0] = 1;   /* one record only */ \
    } \
}

#define REC_3D_VAR_STARTS_COUNTS(rec, starts, counts, nreqs, disps, blocklens, last_dimlen, nreqs_merged) { \
    starts = (MPI_Offset**) malloc(2 * nreqs * sizeof(MPI_Offset*)); \
    counts = starts + nreqs; \
    starts[0] = (MPI_Offset*) malloc(2 * nreqs * 3 * sizeof(MPI_Offset)); \
    counts[0] = starts[0] + nreqs * 3; \
    \
    for (i=1; i<nreqs; i++) { \
        starts[i] = starts[i-1] + 3; \
        counts[i] = counts[i-1] + 3; \
    } \
    \
    int reqs_cnt = 0; \
    i = 0; \
    while (i < nreqs) { \
        starts[reqs_cnt][2] = disps[i] % last_dimlen; /* decomposition is 2D */ \
        counts[reqs_cnt][2] = blocklens[i]; /* each blocklens[i] is no bigger than last_dimlen */ \
        \
        starts[reqs_cnt][1] = disps[i] / last_dimlen; \
        counts[reqs_cnt][1] = 1; \
        \
        starts[reqs_cnt][0] = rec; /* record ID */ \
        counts[reqs_cnt][0] = 1;   /* one record only */ \
        \
        /* merge consecutive "1 by last_dimlen" regions into a larger one */ \
        while ((i < nreqs - 1) && \
               (disps[i] % last_dimlen == 0 && blocklens[i] == last_dimlen) && \
               (disps[i + 1] - disps[i] == last_dimlen && blocklens[i + 1] == last_dimlen)) { \
            i++; \
            counts[reqs_cnt][1]++; \
        } \
        \
        reqs_cnt++; \
        i++; \
    } \
    \
    nreqs_merged = reqs_cnt; \
}

#define POST_VARN(k, num, vid) \
    for (j=0; j<num; j++) { \
        err = ncmpi_iput_varn(ncid, vid+j, nreqs[k-1], starts_D##k, \
                              counts_D##k, rec_buf_ptr, -1, REC_DTYPE, NULL); \
        ERR \
        rec_buf_ptr += nelems[k-1] + gap; \
        my_nreqs += nreqs[k-1]; \
    }


/*----< run_varn_F_case() >--------------------------------------------------*/
int
run_varn_F_case(char       *out_dir,      /* output folder name */
                char       *outfile,      /* output file name */
                int         nvars,        /* number of variables 408 or 51 */
                int         num_recs,     /* number of records */
                int         noncontig_buf,/* whether to us noncontiguous buffer */
                MPI_Info    info,
                MPI_Offset  dims[3][2],   /* dimension lengths */
                int         nreqs[3],     /* no. request in decompositions 1,2,3 */
                int        *disps[3],     /* request's displacements */
                int        *blocklens[3]) /* request's block lengths */
{
    char outfname[512], txt_buf[16], *txt_buf_ptr;
    int i, j, k, err, nerrs=0, rank, ncid, cmode, *varids, nreqs_D3_merged;
    int rec_no, gap=0, my_nreqs, max_nreqs, int_buf[10], *int_buf_ptr;
    size_t dbl_buflen, rec_buflen;
    size_t nelems[3];
    dtype *rec_buf, *rec_buf_ptr;
    double *dbl_buf, *dbl_buf_ptr;
    double pre_timing, open_timing, post_timing, wait_timing, close_timing;
    double timing, total_timing,  max_timing;
    MPI_Offset tmp, metadata_size, put_size, total_size;
    MPI_Offset **fix_starts_D1, **fix_counts_D1;
    MPI_Offset **fix_starts_D2, **fix_counts_D2;
    MPI_Offset **starts_D2, **counts_D2;
    MPI_Offset **starts_D3, **counts_D3;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Info info_used=MPI_INFO_NULL;

    MPI_Barrier(comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime();

    open_timing = 0.0;
    post_timing = 0.0;
    wait_timing = 0.0;
    close_timing = 0.0;

    MPI_Comm_rank(comm, &rank);

    if (noncontig_buf) gap = 10;

    /* calculate number of variable elements from 3 decompositions */
    my_nreqs = max_nreqs = 0;
    nelems[0] = nelems[1] = nelems[2] = 0;
    for (k=0; k<nreqs[0]; k++) nelems[0] += blocklens[0][k];
    for (k=0; k<nreqs[1]; k++) nelems[1] += blocklens[1][k];
    for (k=0; k<nreqs[2]; k++) nelems[2] += blocklens[2][k];

    if (verbose && rank == 0)
        printf("nelems=%zd %zd %zd\n", nelems[0],nelems[1],nelems[2]);

    /* construct varn API arguments starts[][] and counts[][] */
    FIX_1D_VAR_STARTS_COUNTS(fix_starts_D1, fix_counts_D1, nreqs[0], disps[0], blocklens[0])
    FIX_1D_VAR_STARTS_COUNTS(fix_starts_D2, fix_counts_D2, nreqs[1], disps[1], blocklens[1])

    REC_2D_VAR_STARTS_COUNTS(0, starts_D2, counts_D2, nreqs[1], disps[1], blocklens[1])

    REC_3D_VAR_STARTS_COUNTS(0, starts_D3, counts_D3, nreqs[2], disps[2], blocklens[2], dims[2][1], nreqs_D3_merged)
    if (nreqs_D3_merged < nreqs[2])
        nreqs[2] = nreqs_D3_merged;

    /* allocate and initialize write buffer for small variables */
    dbl_buflen = nelems[1] * 2 + nelems[0]
               + 3 * dims[2][0] + 3 * (dims[2][0]+1) + 8 + 2
               + 20 * gap;
    dbl_buf = (double*) malloc(dbl_buflen * sizeof(double));
    for (i=0; i<dbl_buflen; i++) dbl_buf[i] = rank;

    /* allocate and initialize write buffer for large variables */
    if (nvars == 408)
        rec_buflen = nelems[1] * 315 + nelems[2] * 63 + (315+63) * gap;
    else
        rec_buflen = nelems[1] * 20 + nelems[2] + (20+1) * gap;

    rec_buf = (dtype*) malloc(rec_buflen * sizeof(dtype));
    for (i=0; i<rec_buflen; i++) rec_buf[i] = rank;

    for (i=0; i<10; i++) int_buf[i] = rank;

    for (i=0; i<16; i++) txt_buf[i] = 'a' + rank;

    varids = (int*) malloc(nvars * sizeof(int));

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* set output file name */
    sprintf(outfname, "%s/%s",out_dir, outfile);

    /* create a new CDF-5 file for writing */
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(comm, outfname, cmode, info, &ncid); ERR

    /* define dimensions, variables, and attributes */
    if (nvars == 408) {
        /* for h0 file */
        err = def_F_case_h0(ncid, dims[2], nvars, varids); ERR
    }
    else {
        /* for h1 file */
        err = def_F_case_h1(ncid, dims[2], nvars, varids); ERR
    }

    /* exit define mode and enter data mode */
    err = ncmpi_enddef(ncid); ERR

    /* I/O amount so far */
    err = ncmpi_inq_put_size(ncid, &metadata_size); ERR
    err = ncmpi_inq_file_info(ncid, &info_used); ERR
    open_timing += MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    i = 0;
    dbl_buf_ptr = dbl_buf;

    /* lat */
    err = ncmpi_iput_varn(ncid, varids[i++], nreqs[1], fix_starts_D2, fix_counts_D2,
                          dbl_buf_ptr, nelems[1], MPI_DOUBLE, NULL); ERR
    dbl_buf_ptr += nelems[1] + gap;
    my_nreqs += nreqs[1];

    /* lon */
    err = ncmpi_iput_varn(ncid, varids[i++], nreqs[1], fix_starts_D2, fix_counts_D2,
                          dbl_buf_ptr, nelems[1], MPI_DOUBLE, NULL); ERR
    dbl_buf_ptr += nelems[1] + gap;
    my_nreqs += nreqs[1];

    /* area */
    err = ncmpi_iput_varn(ncid, varids[i++], nreqs[0], fix_starts_D1, fix_counts_D1,
                          dbl_buf_ptr, nelems[0], MPI_DOUBLE, NULL); ERR
    dbl_buf_ptr += nelems[0] + gap;
    my_nreqs += nreqs[0];

    post_timing += MPI_Wtime() - timing;

    for (rec_no=0; rec_no<num_recs; rec_no++) {
        MPI_Barrier(comm); /*-----------------------------------------*/
        timing = MPI_Wtime();

        i=3;
        dbl_buf_ptr = dbl_buf + nelems[1]*2 + nelems[0] + gap*3;
        int_buf_ptr = int_buf;
        txt_buf_ptr = txt_buf;

        /* next 27 small variables are written by rank 0 only */
        if (rank == 0) {
            my_nreqs += 27;
            /* post nonblocking requests using ncmpi_iput_varn() */
            err = write_small_vars_F_case(ncid, i, varids, rec_no, gap,
                                          dims[2][0], dims[2][0]+1, 2, 8,
                                          &int_buf_ptr, &txt_buf_ptr,
                                          &dbl_buf_ptr);
            ERR
        }
        i += 27;

        post_timing += MPI_Wtime() - timing;

        MPI_Barrier(comm); /*-----------------------------------------*/
        timing = MPI_Wtime();

        /* flush fixed-size and small variables */
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR

        wait_timing += MPI_Wtime() - timing;

        MPI_Barrier(comm); /*-----------------------------------------*/
        timing = MPI_Wtime();

        /* high water mark of number of noncontiguous requests */
        if (my_nreqs > max_nreqs) max_nreqs = my_nreqs;
        my_nreqs = 0;

        rec_buf_ptr = rec_buf;

        for (j=0; j<nreqs[1]; j++) starts_D2[j][0] = rec_no;
        for (j=0; j<nreqs[2]; j++) starts_D3[j][0] = rec_no;

        if (nvars == 408) {
            if (two_buf) {
                /* write 2D variables */
                POST_VARN(2,   1,  30)   /* AEROD_v */
                POST_VARN(2,  18,  33)   /* AODABS ... ANSNOW */
                POST_VARN(2,   6,  53)   /* AQ_DMS ... AQ_SOAG */
                POST_VARN(2,   2,  64)   /* CDNUMC and CLDHGH */
                POST_VARN(2,   3,  68)   /* CLDLOW ... CLDTOT */
                POST_VARN(2,  11,  75)   /* DF_DMS ... DSTSFMBL */
                POST_VARN(2,   2,  87)   /* DTENDTH and DTENDTQ */
                POST_VARN(2,   7,  91)   /* FLDS ... FLUTC */
                POST_VARN(2,  15, 102)   /* FSDS ... ICEFRAC */
                POST_VARN(2,   2, 120)   /* LANDFRAC and LHFLX */
                POST_VARN(2,   3, 127)   /* LINOZ_SZA ... LWCF */
                POST_VARN(2,   2, 142)   /* O3_SRF and OCNFRAC */
                POST_VARN(2,   1, 145)   /* OMEGA500 */
                POST_VARN(2,   8, 147)   /* PBLH ... PSL */
                POST_VARN(2,   2, 156)   /* QFLX and QREFHT */
                POST_VARN(2,   1, 161)   /* RAM1 */
                POST_VARN(2,  37, 163)   /* SFDMS ... SNOWHLND */
                POST_VARN(2,  10, 202)   /* SO2_CLXF ... SWCF */
                POST_VARN(2,  19, 213)   /* TAUGWX ... TVQ */
                POST_VARN(2,   1, 233)   /* U10 */
                POST_VARN(2,   3, 240)   /* WD_H2O2 ... WD_SO2 */
                POST_VARN(2,  32, 246)   /* airFV ... dst_c3SFWET */
                POST_VARN(2, 129, 279)   /* mlip ... soa_c3SFWET */
                /* write 3D variables */
                POST_VARN(3,   2,  31)   /* ANRAIN and ANSNOW */
                POST_VARN(3,   2,  51)   /* AQRAIN and AQSNOW */
                POST_VARN(3,   5,  59)   /* AREI ... CCN3 */
                POST_VARN(3,   2,  66)   /* CLDICE and CLDLIQ */
                POST_VARN(3,   4,  71)   /* CLOUD ... DCQ */
                POST_VARN(3,   1,  86)   /* DTCOND */
                POST_VARN(3,   2,  89)   /* EXTINCT and FICE */
                POST_VARN(3,   4,  98)   /* FREQI ... FREQS */
                POST_VARN(3,   3, 117)   /* ICIMR ... IWC */
                POST_VARN(3,   5, 122)   /* LINOZ_DO3 ... LINOZ_SSO3 */
                POST_VARN(3,  12, 130)   /* Mass_bc ... O3 */
                POST_VARN(3,   1, 144)   /* OMEGA */
                POST_VARN(3,   1, 146)   /* OMEGAT */
                POST_VARN(3,   1, 155)   /* Q */
                POST_VARN(3,   3, 158)   /* QRL ... RAINQM */
                POST_VARN(3,   1, 162)   /* RELHUM */
                POST_VARN(3,   2, 200)   /* SNOWQM and SO2 */
                POST_VARN(3,   1, 212)   /* T */
                POST_VARN(3,   1, 232)   /* U */
                POST_VARN(3,   6, 234)   /* UU ... VV */
                POST_VARN(3,   3, 243)   /* WSUB ... aero_water */
                POST_VARN(3,   1, 278)   /* hstobie_linoz */
            } else {
                /* write variables in the same order as they defined */
                POST_VARN(2,   1,  30)   /* AEROD_v */
                POST_VARN(3,   2,  31)   /* ANRAIN and ANSNOW */
                POST_VARN(2,  18,  33)   /* AODABS ... ANSNOW */
                POST_VARN(3,   2,  51)   /* AQRAIN and AQSNOW */
                POST_VARN(2,   6,  53)   /* AQ_DMS ... AQ_SOAG */
                POST_VARN(3,   5,  59)   /* AREI ... CCN3 */
                POST_VARN(2,   2,  64)   /* CDNUMC and CLDHGH */
                POST_VARN(3,   2,  66)   /* CLDICE and CLDLIQ */
                POST_VARN(2,   3,  68)   /* CLDLOW ... CLDTOT */
                POST_VARN(3,   4,  71)   /* CLOUD ... DCQ */
                POST_VARN(2,  11,  75)   /* DF_DMS ... DSTSFMBL */
                POST_VARN(3,   1,  86)   /* DTCOND */
                POST_VARN(2,   2,  87)   /* DTENDTH and DTENDTQ */
                POST_VARN(3,   2,  89)   /* EXTINCT and FICE */
                POST_VARN(2,   7,  91)   /* FLDS ... FLUTC */
                POST_VARN(3,   4,  98)   /* FREQI ... FREQS */
                POST_VARN(2,  15, 102)   /* FSDS ... ICEFRAC */
                POST_VARN(3,   3, 117)   /* ICIMR ... IWC */
                POST_VARN(2,   2, 120)   /* LANDFRAC and LHFLX */
                POST_VARN(3,   5, 122)   /* LINOZ_DO3 ... LINOZ_SSO3 */
                POST_VARN(2,   3, 127)   /* LINOZ_SZA ... LWCF */
                POST_VARN(3,  12, 130)   /* Mass_bc ... O3 */
                POST_VARN(2,   2, 142)   /* O3_SRF and OCNFRAC */
                POST_VARN(3,   1, 144)   /* OMEGA */
                POST_VARN(2,   1, 145)   /* OMEGA500 */
                POST_VARN(3,   1, 146)   /* OMEGAT */
                POST_VARN(2,   8, 147)   /* PBLH ... PSL */
                POST_VARN(3,   1, 155)   /* Q */
                POST_VARN(2,   2, 156)   /* QFLX and QREFHT */
                POST_VARN(3,   3, 158)   /* QRL ... RAINQM */
                POST_VARN(2,   1, 161)   /* RAM1 */
                POST_VARN(3,   1, 162)   /* RELHUM */
                POST_VARN(2,  37, 163)   /* SFDMS ... SNOWHLND */
                POST_VARN(3,   2, 200)   /* SNOWQM and SO2 */
                POST_VARN(2,  10, 202)   /* SO2_CLXF ... SWCF */
                POST_VARN(3,   1, 212)   /* T */
                POST_VARN(2,  19, 213)   /* TAUGWX ... TVQ */
                POST_VARN(3,   1, 232)   /* U */
                POST_VARN(2,   1, 233)   /* U10 */
                POST_VARN(3,   6, 234)   /* UU ... VV */
                POST_VARN(2,   3, 240)   /* WD_H2O2 ... WD_SO2 */
                POST_VARN(3,   3, 243)   /* WSUB ... aero_water */
                POST_VARN(2,  32, 246)   /* airFV ... dst_c3SFWET */
                POST_VARN(3,   1, 278)   /* hstobie_linoz */
                POST_VARN(2, 129, 279)   /* mlip ... soa_c3SFWET */
            }
        }
        else {
            if (two_buf) {
                /* write 2D variables followed by 3D variables */
                POST_VARN(2, 13, 30)   /* CLDHGH ... T5 */
                POST_VARN(2,  7, 44)   /* U250 ... Z500 */
                POST_VARN(3,  1, 43)   /* U */
            } else {
                /* write variables in the same order as they defined */
                POST_VARN(2, 13, 30)   /* CLDHGH ... T5 */
                POST_VARN(3,  1, 43)   /* U */
                POST_VARN(2,  7, 44)   /* U250 ... Z500 */
            }
        }

        post_timing += MPI_Wtime() - timing;

        MPI_Barrier(comm); /*-----------------------------------------*/
        timing = MPI_Wtime();

        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR

        /* high water mark of number of noncontiguous requests */
        if (my_nreqs > max_nreqs) max_nreqs = my_nreqs;
        my_nreqs = 0;

        wait_timing += MPI_Wtime() - timing;
    }

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    err = ncmpi_inq_put_size(ncid, &total_size); ERR
    put_size = total_size - metadata_size;
    err = ncmpi_close(ncid); ERR
    close_timing += MPI_Wtime() - timing;

    free(starts_D3[0]); free(starts_D3);
    free(starts_D2[0]); free(starts_D2);
    free(fix_starts_D2[0]); free(fix_starts_D2);
    free(fix_starts_D1[0]); free(fix_starts_D1);
    free(rec_buf);
    free(dbl_buf);
    free(varids);

    total_timing = MPI_Wtime() - total_timing;

    my_nreqs = max_nreqs;
    MPI_Reduce(&my_nreqs,      &max_nreqs,  1, MPI_INT,    MPI_MAX, 0, comm);
    MPI_Reduce(&put_size,      &tmp,        1, MPI_OFFSET, MPI_SUM, 0, comm);
    put_size = tmp;
    MPI_Reduce(&total_size,    &tmp,        1, MPI_OFFSET, MPI_SUM, 0, comm);
    total_size = tmp;
    MPI_Reduce(&open_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    open_timing = max_timing;
    MPI_Reduce(&pre_timing,    &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    pre_timing = max_timing;
    MPI_Reduce(&post_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    post_timing = max_timing;
    MPI_Reduce(&wait_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    wait_timing = max_timing;
    MPI_Reduce(&close_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    close_timing = max_timing;
    MPI_Reduce(&total_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    total_timing = max_timing;

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);
        if (rank == 0 && sum_size > 0) {
            printf("-----------------------------------------------------------\n");
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        }
    }
    MPI_Offset m_alloc=0, max_alloc;
    ncmpi_inq_malloc_max_size(&m_alloc);
    MPI_Reduce(&m_alloc, &max_alloc, 1, MPI_OFFSET, MPI_MAX, 0, comm);
    if (rank == 0) {
        printf("History output file                = %s\n", outfile);
        printf("MAX heap memory allocated by PnetCDF internally is %.2f MiB\n",
               (float)max_alloc/1048576);
        printf("Total number of variables          = %d\n",nvars);
        printf("Total write amount                 = %.2f MiB = %.2f GiB\n",
               (double)total_size/1048576,(double)total_size/1073741824);
        printf("Max number of requests             = %d\n",max_nreqs);
        printf("Max Time of open + metadata define = %.4f sec\n",open_timing);
        printf("Max Time of I/O preparing          = %.4f sec\n",pre_timing);
        printf("Max Time of ncmpi_iput_varn        = %.4f sec\n",post_timing);
        printf("Max Time of ncmpi_wait_all         = %.4f sec\n",wait_timing);
        printf("Max Time of close                  = %.4f sec\n",close_timing);
        printf("Max Time of TOTAL                  = %.4f sec\n",total_timing);
        printf("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
               (double)total_size/1048576.0/total_timing);
        printf("I/O bandwidth (write-only)         = %.4f MiB/sec\n",
               (double)put_size/1048576.0/wait_timing);
        if (verbose) print_info(&info_used);
        printf("-----------------------------------------------------------\n");
    }
fn_exit:
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    if (!keep_outfile) unlink(outfname);
    fflush(stdout);
    MPI_Barrier(comm);
    return nerrs;
}
