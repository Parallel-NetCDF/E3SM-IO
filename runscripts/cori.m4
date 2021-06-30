changequote(`[[[', `]]]')dnl
changecom([[[###]]], [[[$$$]]])dnl
#!/bin/bash
#SBATCH -q EXP_QUEUE
#SBATCH -N EXP_NN
#SBATCH -C EXP_NODE_TYPE
#SBATCH -t 00:EXP_TL:00
#SBATCH -o EXP_NAME.txt
#SBATCH -e EXP_NAME.err
#SBATCH -L SCRATCH
#SBATCH -A EXP_ACC

CONFIG_FDBG=../datasets/f_case_866x72_16p.nc
CONFIG_F30=../datasets/f_case_48602x72_512p.nc
CONFIG_F120=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/decom/FC5AV1C-H01A_ne120_oRRS18v3_21600p.nc
CONFIG_G=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/decom/GMPAS-NYF_T62_oRRS18to6v3_9600p.nc

IN_FDBG=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/F_DBG
IN_F30=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/F_30
IN_F120=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/F_120
IN_G=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/G

E3SM_IO_DATE="EXP_E3SM_IO_DATE"
HDF5_LIB_PATH=EXP_HDF5_LIB_PATH
HDF5_LIB_DATE="EXP_HDF5_LIB_DATE"
ADIOS2_LIB_PATH=EXP_ADIOS2_LIB_PATH
ADIOS2_LIB_DATE="EXP_ADIOS2_LIB_DATE"
LOGVOL_LIB_PATH=EXP_LOGVOL_LIB_PATH
LOGVOL_LIB_DATE="EXP_LOGVOL_LIB_DATE"
PNC_LIB_PATH=EXP_PNC_LIB_PATH
NREC=EXP_RECS
PPN=EXP_PPN
RTL=EXP_RTL
CONFIG_POST=EXP_CONFIG
OUTDIR_ROOT=EXP_OUTDIR_ROOT

APPS=(../src/e3sm_io)
HXS=(EXP_HX)
CONFIG="CONFIG_${CONFIG_POST}"
CONFIG=${!CONFIG}
INDIR="IN_${CONFIG_POST}"
INDIR=${!INDIR}
APIS=(pnc)
LAYOUTS=(contig)
OPS=(EXP_OP)

NN=${SLURM_NNODES}
let NP=NN*PPN

# Print exp setting
echo "#%$: e3sm_io_hash: EXP_E3SM_IO_HASH"
echo "#%$: pnc_hash: EXP_PNC_HASH"
echo "#%$: pnc_path: ${PNC_LIB_PATH}"
echo "#%$: hdf5_hash: EXP_HDF5_HASH"
echo "#%$: hdf5_path: ${HDF5_LIB_PATH}"
echo "#%$: submit_date: EXP_SUBMIT_DATE"
echo "#%$: exe_date: $(date)"
echo "#%$: exp: e3sm_io"
echo "#%$: repeation: EXP_RUNS"
echo "#%$: output_folder: EXP_OUTDIR_ROOT"
echo "#%$: e3sm_config_post: ${CONFIG_POST}"
echo "#%$: e3sm_config: ${CONFIG}"
echo "#%$: number_of_nodes: ${NN}"
echo "#%$: number_of_proc: ${NP}"
echo '-----+-----++----***--------+++++**++++--+---'

for API in ${APIS[@]}
do
    for LAYOUT in ${LAYOUTS[@]}
    do
        echo "mkdir -p ${OUTDIR_ROOT}/${API}/${LAYOUT}"
        mkdir -p ${OUTDIR_ROOT}/${API}/${LAYOUT}
    done
done

export LD_LIBRARY_PATH=${HDF5_LIB_PATH}/lib:${PNC_LIB_PATH}/lib:${ADIOS2_LIB_PATH}/lib64:${LOGVOL_LIB_PATH}/lib:${LD_LIBRARY_PATH}
export PNETCDF_SHOW_PERFORMANCE_INFO=1
#export PNETCDF_DEFAULT_CHUNK_DIM="ncol : 14563 ; nbnd : 2 ; Time : 1 ; ilev : 73 ; lev : 72 ; chars : 64 ;nCells : 16384 ; nEdges : 16384 ; nVertices : 16384 ; nVertLevelsP1 : 81 ; nVertLevels : 80 ; StrLen : 64 ;"
export PNETCDF_HINTS="nc_zip_delay_init=1;nc_zip_nrec=1;nc_zip_buffer_size=0"

ulimit -c unlimited

TSTARTTIME=`date +%s.%N`

for i in $(seq EXP_RUNS);
do
    for APP in ${APPS[@]}
    do
        for HX in ${HXS[@]}
        do
            for API in ${APIS[@]}
            do
                if [ "$API" = "adios" ] ; then
                    STRATE="blob"
                else
                    STRATE="canonical"
                fi

                for LAYOUT in ${LAYOUTS[@]}
                do
                    OUTDIR="${OUTDIR_ROOT}/${API}/${LAYOUT}"
                    echo "rm -f ${OUTDIR}/*"
                    rm -f ${OUTDIR}/*

                    if [ "$LAYOUT" = "zlib" ] ; then
                        INDIRARG="-i ${INDIR}"
                    else
                        INDIRARG=
                    fi

                    for (( j=0; j<${#OPS}; j++ ));
                    do
                        OP=${OPS:$j:1}

                        echo "========================== E3SM-IO ${API} ${OP} =========================="
                        >&2 echo "========================== E3SM-IO ${API} ${OP}=========================="
                        
                        CUR_E3SM_IO_DATE=$(stat -c %Y ./e3sm_io)
                        if [[ "${CUR_E3SM_IO_DATE}" != "${E3SM_IO_DATE}" ]]; then
                            echo "Warning: e3sm_io changed after submission"
                        fi
                        CUR_HDF5_DATE=$(stat -c %Y ${HDF5_LIB_PATH}/lib/libhdf5.so.200.0.0)
                        if [[ "${CUR_HDF5_DATE}" != "${HDF5_LIB_DATE}" ]]; then
                            echo "Warning: libhdf5.so changed after submission"
                        fi
                        CUR_ADIOS2_DATE=$(stat -c %Y ${ADIOS2_LIB_PATH}/lib/libadios2_c_mpi.so)
                        if [[ "${CUR_ADIOS2_DATE}" != "${ADIOS2_LIB_DATE}" ]]; then
                            echo "Warning: libadios2_c_mpi.so changed after submission"
                        fi
                        CUR_LOGVOL_DATE=$(stat -c %Y ${LOGVOL_LIB_PATH}/lib/libH5VL_log.so)
                        if [[ "${CUR_LOGVOL_DATE}" != "${LOGVOL_LIB_DATE}" ]]; then
                            echo "Warning: libH5VL_log.so changed after submission"
                        fi

                        echo "#%$: exp: e3sm_io"
                        echo "#%$: app: ${APP}"
                        echo "#%$: config: ${CONFIG}"
                        echo "#%$: h_num: ${HX}"
                        echo "#%$: nrec: ${NREC}"
                        echo "#%$: api: ${API}"
                        echo "#%$: layout: ${LAYOUT}"
                        echo "#%$: strategy: ${STRATE}"
                        echo "#%$: operation: ${OP}"
                        echo "#%$: delay_init: 1"
                        echo "#%$: number_of_nodes: ${NN}"
                        echo "#%$: number_of_proc: ${NP}"

                        STARTTIME=`date +%s.%N`

                        if [ "$OP" = "w" ] ; then
                            echo "srun -n ${NP} -t ${RTL} ./${APP} -k -o ${OUTDIR} -n -a ${API} -f ${HX} -r ${NREC} -x ${STRATE} ${INDIRARG} ${CONFIG}"
                            srun -n ${NP} -t ${RTL} ./${APP} -k -o ${OUTDIR} -n -a ${API} -f ${HX} -r ${NREC} -x ${STRATE} ${INDIRARG} ${CONFIG}
                        else
                            echo "srun -n ${NP} -t ${RTL} ./${APP} -k -o ${OUTDIR} -n -a ${API} -f ${HX} -r ${NREC} -x ${STRATE} ${INDIRARG} -R ${CONFIG}"
                            srun -n ${NP} -t ${RTL} ./${APP} -k -o ${OUTDIR} -n -a ${API} -f ${HX} -r ${NREC} -x ${STRATE} ${INDIRARG} -R ${CONFIG}
                        fi

                        ENDTIME=`date +%s.%N`
                        TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print $1"."$2}'`

                        echo "#%$: exe_time: $TIMEDIFF"

                        echo "ls -lah ${OUTDIR}"
                        ls -lah ${OUTDIR}
                        echo "lfs getstripe ${OUTDIR}"
                        lfs getstripe ${OUTDIR}

                        echo '-----+-----++------------+++++++++--+---'
                    done

                    echo "rm -f ${OUTDIR}/*"
                    rm -f ${OUTDIR}/*
                done
            done
        done
    done
done

ENDTIME=`date +%s.%N`
TIMEDIFF=`echo "$ENDTIME - $TSTARTTIME" | bc | awk -F"." '{print $1"."$2}'`
echo "-------------------------------------------------------------"
echo "total_exe_time: $TIMEDIFF"
echo "-------------------------------------------------------------"
