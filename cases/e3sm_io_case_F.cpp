#include "e3sm_io_case_F.hpp"

#include "e3sm_io.hpp"
#include "e3sm_io_case.hpp"
#include "e3sm_io_err.hpp"

e3sm_io_case_F::e3sm_io_case_F () {}

e3sm_io_case_F::~e3sm_io_case_F () {
    if (this->dbl_buf_h0 != NULL) { free (this->dbl_buf_h0); }
    if (this->dbl_buf_h1 != NULL) { free (this->dbl_buf_h1); }
    if (this->rec_buf_h0 != NULL) { free (this->rec_buf_h0); }
    if (this->rec_buf_h1 != NULL) { free (this->rec_buf_h1); }
}

int e3sm_io_case_F::wr_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err, nerrs = 0;
    int nvar;

    PRINT_MSG (0, "number of requests for D1=%d D2=%d D3=%d\n", decom.contig_nreqs[0],
               decom.contig_nreqs[1], decom.contig_nreqs[2]);

    if ((cfg.verbose >= 0) && (cfg.rank == 0)) {
        printf ("Total number of MPI processes      = %d\n", cfg.np);
        printf ("Number of IO processes             = %d\n", cfg.num_iotasks);
        printf ("Input decomposition file           = %s\n", cfg.cfgpath.c_str ());
        printf ("Number of decompositions           = %d\n", decom.num_decomp);
        printf ("Output file directory              = %s\n", cfg.targetdir.c_str ());
        printf ("Variable dimensions (C order)      = %lld x %lld\n", decom.dims[2][0],
                decom.dims[2][1]);
        printf ("Write number of records (time dim) = %d\n", cfg.nrec);
        printf ("Using noncontiguous write buffer   = %s\n", cfg.non_contig_buf ? "yes" : "no");
    }

    nvar = cfg.nvars;

    /* vard APIs require internal data type matches external one */
    if (cfg.low_lvl) {
#if REC_XTYPE != NC_FLOAT
        RET_ERR ("PnetCDF vard API requires internal and external data types match, skip\n");
#endif

        PRINT_MSG (0, "\n==== benchmarking F case using vard API ========================\n");
        PRINT_MSG (0, "Variable written order: same as variables are defined\n\n");
        fflush (stdout);

        if (cfg.hx == 0 || cfg.hx == -1) {
            MPI_Barrier (cfg.io_comm);
            cfg.nvars = 414;
            nerrs += run_vard_F_case (cfg, decom, driver, "f_case_h0_varn.nc", this->dbl_buf_h0,
                                      this->rec_buf_h0, this->txt_buf[0], this->int_buf[0]);
        }

        if (cfg.hx == 0 || cfg.hx == -1) {
            MPI_Barrier (cfg.io_comm);
            cfg.nvars = 51;
            nerrs += run_vard_F_case (cfg, decom, driver, "f_case_h1_varn.nc", this->dbl_buf_h0,
                                      this->rec_buf_h0, this->txt_buf[0], this->int_buf[0]);
        }

    } else {
        PRINT_MSG (0,
                   "\n==== benchmarking F case writing using varn API ========================\n");

        if (cfg.two_buf) {
            PRINT_MSG (0, "Variable written order: 2D variables then 3D variables\n\n");
        } else {
            PRINT_MSG (0, "Variable written order: same as variables are defined\n\n");
        }

        fflush (stdout);

        if (cfg.hx == 0 || cfg.hx == -1) {
            MPI_Barrier (cfg.io_comm);
            cfg.nvars = 414;
            nerrs += run_varn_F_case (cfg, decom, driver, "f_case_h0_varn.nc", this->dbl_buf_h0,
                                      this->rec_buf_h0, this->txt_buf[0], this->int_buf[0]);
        }

        if (cfg.hx == 0 || cfg.hx == -1) {
            MPI_Barrier (cfg.io_comm);
            cfg.nvars = 51;
            nerrs += run_varn_F_case (cfg, decom, driver, "f_case_h1_varn.nc", this->dbl_buf_h0,
                                      this->rec_buf_h0, this->txt_buf[0], this->int_buf[0]);
        }
    }

    cfg.nvars = nvar;
err_out:;
    return nerrs;
}

int e3sm_io_case_F::rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err, nerrs = 0;
    int nvar;

    PRINT_MSG (0, "number of requests for D1=%d D2=%d D3=%d\n", decom.contig_nreqs[0],
               decom.contig_nreqs[1], decom.contig_nreqs[2]);

    if ((cfg.verbose >= 0) && (cfg.rank == 0)) {
        printf ("Total number of MPI processes      = %d\n", cfg.np);
        printf ("Number of IO processes             = %d\n", cfg.num_iotasks);
        printf ("Input decomposition file           = %s\n", cfg.cfgpath.c_str ());
        printf ("Number of decompositions           = %d\n", decom.num_decomp);
        printf ("Output file directory              = %s\n", cfg.targetdir.c_str ());
        printf ("Variable dimensions (C order)      = %lld x %lld\n", decom.dims[2][0],
                decom.dims[2][1]);
        printf ("Write number of records (time dim) = %d\n", cfg.nrec);
        printf ("Using noncontiguous write buffer   = %s\n", cfg.non_contig_buf ? "yes" : "no");
    }

    nvar = cfg.nvars;

    /* vard APIs require internal data type matches external one */
    if (cfg.low_lvl) {
#if REC_XTYPE != NC_FLOAT
        RET_ERR ("Low level API requires internal and external data types match, skip\n");
#endif
        RET_ERR ("Reading not supported for low-level API\n");
    } else {
        PRINT_MSG (
            0, "\n==== benchmarking F case writing using low level API ========================\n");

        if (cfg.two_buf) {
            PRINT_MSG (0, "Variable written order: 2D variables then 3D variables\n\n");
        } else {
            PRINT_MSG (0, "Variable written order: same as variables are defined\n\n");
        }

        fflush (stdout);

        if (cfg.hx == 0 || cfg.hx == -1) {
            MPI_Barrier (cfg.io_comm);
            cfg.nvars = 414;
            nerrs +=
                run_varn_F_case_rd (cfg, decom, driver, "f_case_h0_varn.nc", &(this->dbl_buf_h0),
                                    &(this->rec_buf_h0), this->txt_buf[0], this->int_buf[0]);
        }
        if (cfg.hx == 0 || cfg.hx == -1) {
            MPI_Barrier (cfg.io_comm);
            cfg.nvars = 51;
            nerrs +=
                run_varn_F_case_rd (cfg, decom, driver, "f_case_h1_varn.nc", &(this->dbl_buf_h0),
                                    &(this->rec_buf_h0), this->txt_buf[0], this->int_buf[0]);
        }
    }

    cfg.nvars = nvar;

err_out:;
    return nerrs;
}

int e3sm_io_case_F::load_data (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err, nerrs = 0;
    int verbose, nvar;

    verbose = cfg.verbose;
    nvar    = cfg.nvars;

    cfg.verbose = -1;  // Disable output

    // Run dummy read test to get data
    if (cfg.hx == 0 || cfg.hx == -1) {
        MPI_Barrier (cfg.io_comm);
        cfg.nvars = 414;
        nerrs += run_varn_F_case_rd (cfg, decom, driver, "f_case_h0_varn.nc", &(this->dbl_buf_h0),
                                     &(this->rec_buf_h0), this->txt_buf[0], this->int_buf[0]);
    }
    if (cfg.hx == 0 || cfg.hx == -1) {
        MPI_Barrier (cfg.io_comm);
        cfg.nvars = 51;
        nerrs += run_varn_F_case_rd (cfg, decom, driver, "f_case_h1_varn.nc", &(this->dbl_buf_h0),
                                     &(this->rec_buf_h0), this->txt_buf[0], this->int_buf[0]);
    }

    cfg.nvars   = nvar;
    cfg.verbose = verbose;

err_out:;
    return nerrs;
}