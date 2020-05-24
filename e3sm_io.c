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
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */

#include <e3sm_io.h>

/*----< print_info() >------------------------------------------------------*/
void
print_info(MPI_Info *info_used)
{
    int  i, nkeys;

    MPI_Info_get_nkeys(*info_used, &nkeys);
    printf("MPI File Info: nkeys = %d\n",nkeys);
    for (i=0; i<nkeys; i++) {
        char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
        int  valuelen, flag;

        MPI_Info_get_nthkey(*info_used, i, key);
        MPI_Info_get_valuelen(*info_used, key, &valuelen, &flag);
        MPI_Info_get(*info_used, key, valuelen+1, value, &flag);
        printf("MPI File Info: [%2d] key = %25s, value = %s\n",i,key,value);
    }
}

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [OPTION]... FILE\n"
    "       [-h] Print help\n"
    "       [-v] Verbose mode\n"
    "       [-k] Keep the output files when program exits\n"
    "       [-d] Run test that uses PnetCDF vard API\n"
    "       [-n] Run test that uses PnetCDF varn API\n"
    "       [-m] Run test using noncontiguous write buffer\n"
    "       [-t] Write 2D variables followed by 3D variables\n"
    "       [-r num] Number of records (default 1)\n"
    "       [-s num] Stride between IO tasks (default 1)\n"
    "       [-o output_dir] Output directory name (default ./)\n"
    "       FILE: Name of input netCDF file describing data decompositions\n";
    fprintf(stderr, help, argv0);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char** argv)
{
    extern int optind;
    char *infname, out_dir[1024], in_dir[1024], *outfname;
    int i, rank, nprocs, err, nerrs=0, tst_vard=0, tst_varn=0, tst_wr=0, tst_rd=0, noncontig_buf=0;
    int num_recs, io_stride=1, num_iotasks, ioproc;
    int num_decomp, nvars, run_f_case, run_g_case;
    int contig_nreqs[MAX_NUM_DECOMP], *disps[MAX_NUM_DECOMP];
    int *blocklens[MAX_NUM_DECOMP];
    MPI_Offset dims[MAX_NUM_DECOMP][2], estimated_nc_ibuf_size;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Comm io_comm;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    out_dir[0] = '\0';
    in_dir[0] = '\0';
    verbose = 0;
    keep_outfile = 0;
    num_recs = 1;
    two_buf = 0;

    /* command-line arguments */
    while ((i = getopt(argc, argv, "hkvdnmto:r:WR:s:")) != EOF)
        switch(i) {
            case 'v': verbose = 1;
                      break;
            case 'k': keep_outfile = 1;
                      break;
            case 'r': num_recs = atoi(optarg);
                      break;
            case 's': io_stride = atoi(optarg);
                      break;
            case 'R': strcpy(in_dir, optarg);  
                      tst_rd = 1;
                      break;
            case 'W': tst_wr = 1;
                      break;
            case 'd': tst_vard = 1;
                      break;
            case 'n': tst_varn = 1;
                      break;
            case 'm': noncontig_buf = 1;
                      break;
            case 't': two_buf = 1;
                      break;
            case 'o': strcpy(out_dir, optarg);
                      break;
            case 'i': strcpy(in_dir, optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }

    if (argv[optind] == NULL) { /* input file is mandatory */
        if (!rank) usage(argv[0]);
        MPI_Finalize();
        return 1;
    }
    /* input file contains number of write requests and their file access
     * offsets (per array element) */
    infname = argv[optind];
    if (verbose && rank==0) printf("input file name =%s\n",infname);


    if (tst_vard == 0 && tst_varn == 0)
        /* neither command-line option -d or -n is used, run both */
        tst_vard = tst_varn = 1;
    
    if (tst_wr == 0 && tst_rd == 0)
        /* neither command-line option -R or -W is used, run write */
        tst_wr = 1;

    /* set the output folder name */
    if (tst_wr){
        if (out_dir[0] == '\0') {
            strcpy(out_dir, ".");
        }
        if (verbose && rank==0) printf("output folder name =%s\n",out_dir);
    }

    /* set the input folder name */
    if (tst_rd){
        if (in_dir[0] == '\0') {
            strcpy(in_dir, ".");
        }
        if (verbose && rank==0) printf("input folder name =%s\n",in_dir);
    }

    if (io_stride < 1)
        io_stride = 1;

    if (io_stride == 1) {
        num_iotasks = nprocs;
        io_comm = MPI_COMM_WORLD;
        ioproc = 1;
    }
    else {
        int *ioranks;
        MPI_Group group, iogroup;

        /* assume that the IO root (rank of the first IO task) is 0 */
        num_iotasks = (nprocs - 1) / io_stride + 1;

        /* create an array that holds the ranks of the IO tasks */
        ioranks = malloc(num_iotasks * sizeof(int));
        ioproc = 0;
        for (i = 0; i < num_iotasks; i++) {
            ioranks[i] = i * io_stride;
            if (ioranks[i] == rank)
                ioproc = 1;
        }

        /* create a group for all MPI tasks */
        MPI_Comm_group(MPI_COMM_WORLD, &group);

        /* create a sub-group for the IO tasks */
        MPI_Group_incl(group, num_iotasks, ioranks, &iogroup);

        /* create an MPI communicator for the IO tasks */
        MPI_Comm_create(MPI_COMM_WORLD, iogroup, &io_comm);

        MPI_Group_free(&iogroup);
        MPI_Group_free(&group);
        free(ioranks);
    }

    /* only IO tasks call PnetCDF APIs */
    if (!ioproc) goto fn_exit;

    /* set MPI-IO hints */
    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_cb_write", "enable");  /* collective write */
    MPI_Info_set(info, "romio_no_indep_rw", "true"); /* no independent MPI-IO */

    /* set PnetCDF I/O hints */
    MPI_Info_set(info, "nc_var_align_size", "1"); /* no gap between variables */
    MPI_Info_set(info, "nc_in_place_swap", "enable"); /* in-place byte swap */

    /* read request information from decompositions 1, 2 and 3 */
    err = read_decomp(io_comm, infname, &num_decomp, dims, contig_nreqs, disps, blocklens);
    if (err) goto fn_exit;

    /* F case has 3 decompositions, G case has 6 */
    run_f_case = 0;
    run_g_case = 0;
    if (num_decomp == 3) run_f_case = 1;
    else if (num_decomp == 6) run_g_case = 1;

    /* use total write amount to estimate nc_ibuf_size */
    estimated_nc_ibuf_size = dims[2][0] * dims[2][1] * sizeof(double) / num_iotasks;
    estimated_nc_ibuf_size *= (run_f_case) ? 408 : 52;
    if (estimated_nc_ibuf_size > 16777216) {
        char nc_ibuf_size_str[32];
        sprintf(nc_ibuf_size_str, "%lld", estimated_nc_ibuf_size);
        MPI_Info_set(info, "nc_ibuf_size", nc_ibuf_size_str);
    }

if (run_f_case) {
        double *dbl_buf_h0 = NULL, *dbl_buf_h1 = NULL;
        itype *rec_buf_h0 = NULL, *rec_buf_h1 = NULL;
        char txt_buf[2][16];
        int int_buf[2][10];

        if (verbose && rank==0) {
            printf("number of requests for D1=%d D2=%d D3=%d\n",
                   contig_nreqs[0], contig_nreqs[1], contig_nreqs[2]);
        }

        if (!rank) {
            printf("Total number of MPI processes      = %d\n",nprocs);
            printf("Number of IO processes             = %d\n",num_iotasks);
            printf("Input decomposition file           = %s\n",infname);
            printf("Number of decompositions           = %d\n",num_decomp);
            printf("Output file directory              = %s\n",out_dir);
            printf("Variable dimensions (C order)      = %lld x %lld\n",dims[2][0],dims[2][1]);
            printf("Write number of records (time dim) = %d\n",num_recs);
            printf("Using noncontiguous write buffer   = %s\n",noncontig_buf?"yes":"no");
        }

        /* vard APIs require internal data type matches external one */
        if (tst_vard) {
#if REC_XTYPE != NC_FLOAT
            if (!rank)
                printf("PnetCDF vard API requires internal and external data types match, skip\n");
#else  
            if (tst_rd){
                if (!rank)
                    printf("Reading not supported for vard\n");
            }
            
            if (tst_wr){
                if (!rank) {
                    printf("\n==== benchmarking F case using vard API ========================\n");
                    printf("Variable written order: same as variables are defined\n\n");
                }
                fflush(stdout);
                MPI_Barrier(io_comm);

                nvars = 408;
                outfname = "f_case_h0_vard.nc";
                nerrs += run_vard_F_case(io_comm, out_dir, outfname, nvars, num_recs,
                                        noncontig_buf, info, dims,
                                        contig_nreqs, disps, blocklens);

                MPI_Barrier(io_comm);

                nvars = 51;
                outfname = "f_case_h1_vard.nc";
                nerrs += run_vard_F_case(io_comm, out_dir, outfname, nvars, num_recs,
                                        noncontig_buf, info, dims,
                                        contig_nreqs, disps, blocklens);
            }
#endif
        }
        if (tst_varn) {
            if (tst_rd){        
                if (!rank) {
                    printf("\n==== benchmarking F case reading using varn API ========================\n");
                    printf("Variable written order: ");
                    if (two_buf)
                        printf("2D variables then 3D variables\n\n");
                    else
                        printf("same as variables are defined\n\n");
                }
                fflush(stdout);
                MPI_Barrier(MPI_COMM_WORLD);

                /* There are two kinds of outputs for history variables.
                * Output 1st kind history variables.
                */
                nvars = 408;
                outfname = "f_case_h0_varn.nc";
                nerrs += run_varn_F_case_rd(io_comm, in_dir, outfname, nvars, num_recs,
                                        noncontig_buf, info, dims,
                                        contig_nreqs, disps, blocklens, &dbl_buf_h0, &rec_buf_h0, txt_buf[0], int_buf[0]);

                MPI_Barrier(MPI_COMM_WORLD);

                /* Output 2nd kind history variables. */
                nvars = 51;
                outfname = "f_case_h1_varn.nc";
                nerrs += run_varn_F_case_rd(io_comm, in_dir, outfname, nvars, num_recs,
                                        noncontig_buf, info, dims,
                                        contig_nreqs, disps, blocklens, &dbl_buf_h1, &rec_buf_h1, txt_buf[1], int_buf[1]);
            }

            if (tst_wr){
                if (!rank) {
                    printf("\n==== benchmarking F case writing using varn API ========================\n");
                    printf("Variable written order: ");
                    if (two_buf)
                        printf("2D variables then 3D variables\n\n");
                    else
                        printf("same as variables are defined\n\n");
                }
                fflush(stdout);

                MPI_Barrier(io_comm);

                /* There are two kinds of outputs for history variables.
                * Output 1st kind history variables.
                */
                nvars = 408;
                outfname = "f_case_h0_varn.nc";
                nerrs += run_varn_F_case(io_comm, out_dir, outfname, nvars, num_recs,
                                        noncontig_buf, info, dims,
                                        contig_nreqs, disps, blocklens, dbl_buf_h0, rec_buf_h0, txt_buf[0], int_buf[0]);

                MPI_Barrier(io_comm);

                /* Output 2nd kind history variables. */
                nvars = 51;
                outfname = "f_case_h1_varn.nc";
                nerrs += run_varn_F_case(io_comm, out_dir, outfname, nvars, num_recs,
                                        noncontig_buf, info, dims,
                                        contig_nreqs, disps, blocklens, dbl_buf_h1, rec_buf_h1, txt_buf[1], int_buf[1]);
            }
        }
        for (i=0; i<3; i++) {
            free(disps[i]);
            free(blocklens[i]);
        }

        if (dbl_buf_h0 != NULL){
            free(dbl_buf_h0);
        }
        if (dbl_buf_h1 != NULL){
            free(dbl_buf_h1);
        }
        if (rec_buf_h0 != NULL){
            free(rec_buf_h0);
        }
        if (rec_buf_h1 != NULL){
            free(rec_buf_h1);
        }
    }

    if (run_g_case) {
        if (verbose && rank==0) {
            printf("number of requests for D1=%d D2=%d D3=%d D4=%d D5=%d D6=%d\n",
                   contig_nreqs[0], contig_nreqs[1], contig_nreqs[2],
                   contig_nreqs[3], contig_nreqs[4], contig_nreqs[5]);
        }

        if (!rank) {
            printf("Total number of MPI processes      = %d\n",nprocs);
            printf("Number of IO processes             = %d\n",num_iotasks);
            printf("Input decomposition file           = %s\n",infname);
            printf("Number of decompositions           = %d\n",num_decomp);
            printf("Output file directory              = %s\n",out_dir);
            printf("Variable dimensions (C order)      = %lld x %lld\n",dims[2][0],dims[2][1]);
            printf("Write number of records (time dim) = %d\n",num_recs);
            printf("Using noncontiguous write buffer   = %s\n",noncontig_buf?"yes":"no");
        }

        if (tst_varn) {
            double *D1_rec_dbl_buf = NULL, *D3_rec_dbl_buf = NULL, *D4_rec_dbl_buf = NULL, *D5_rec_dbl_buf = NULL, *D6_rec_dbl_buf = NULL, *D1_fix_dbl_buf = NULL;
            int *D1_fix_int_buf = NULL, *D2_fix_int_buf = NULL, *D3_fix_int_buf = NULL, *D4_fix_int_buf = NULL, *D5_fix_int_buf = NULL;

            if (tst_rd){
                if (!rank) {
                    printf("\n==== benchmarking G case reading using varn API ========================\n");
                }
                fflush(stdout);

                MPI_Barrier(MPI_COMM_WORLD);

                nvars = 52;
                outfname = "g_case_hist_varn.nc";
                nerrs += run_varn_G_case_rd(io_comm, in_dir, outfname, nvars, num_recs, info,
                                       dims, contig_nreqs, disps, blocklens, &D1_fix_int_buf, &D2_fix_int_buf, &D3_fix_int_buf, &D4_fix_int_buf, &D5_fix_int_buf,
                &D1_rec_dbl_buf, &D3_rec_dbl_buf, &D4_rec_dbl_buf, &D5_rec_dbl_buf, &D6_rec_dbl_buf, &D1_fix_dbl_buf);
            }

            if (tst_wr){
                if (!rank) {
                    printf("\n==== benchmarking G case writing using varn API ========================\n");
                }
                fflush(stdout);

                MPI_Barrier(MPI_COMM_WORLD);

                nvars = 52;
                outfname = "g_case_hist_varn.nc";
                nerrs += run_varn_G_case(io_comm, out_dir, outfname, nvars, num_recs, info,
                                        dims, contig_nreqs, disps, blocklens, D1_fix_int_buf, D2_fix_int_buf, D3_fix_int_buf, D4_fix_int_buf, D5_fix_int_buf,
                D1_rec_dbl_buf, D3_rec_dbl_buf, D4_rec_dbl_buf, D5_rec_dbl_buf, D6_rec_dbl_buf, D1_fix_dbl_buf);
            }

            if (D1_rec_dbl_buf != NULL){
                free(D1_rec_dbl_buf);
            }
            if (D3_rec_dbl_buf != NULL){
                free(D3_rec_dbl_buf);
            }
            if (D4_rec_dbl_buf != NULL){
                free(D4_rec_dbl_buf);
            }
            if (D5_rec_dbl_buf != NULL){
                free(D5_rec_dbl_buf);
            }
            if (D6_rec_dbl_buf != NULL){
                free(D6_rec_dbl_buf);
            }
            if (D1_fix_dbl_buf != NULL){
                free(D1_fix_dbl_buf);
            }
            if (D1_fix_int_buf != NULL){
                free(D1_fix_int_buf);
            }
            if (D2_fix_int_buf != NULL){
                free(D2_fix_int_buf);
            }
            if (D3_fix_int_buf != NULL){
                free(D3_fix_int_buf);
            }
            if (D4_fix_int_buf != NULL){
                free(D4_fix_int_buf);
            }
            if (D5_fix_int_buf != NULL){
                free(D5_fix_int_buf);
            }
        }
        for (i=0; i<6; i++) {
            free(disps[i]);
            free(blocklens[i]);
        }       
    }

fn_exit:
    if (info != MPI_INFO_NULL) MPI_Info_free(&info);

    /* Non-IO tasks wait for IO tasks to complete */
    MPI_Barrier(MPI_COMM_WORLD);

    if (io_comm != MPI_COMM_WORLD && io_comm != MPI_COMM_NULL)
        MPI_Comm_free(&io_comm);

    MPI_Finalize();
    return (nerrs > 0);
}

