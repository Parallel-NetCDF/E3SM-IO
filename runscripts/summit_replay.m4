changequote(`[[[', `]]]')dnl
changecom([[[###]]], [[[$$$]]])dnl
#!/bin/bash
# Begin LSF Directives
#BSUB -P VAR_ACC
#BSUB -W 00:VAR_TL
#BSUB -nnodes VAR_NN
#BSUB -J [[[]]]VAR_NAME[[[]]]
#BSUB -o [[[]]]VAR_NAME[[[]]]_%J.txt
#BSUB -e [[[]]]VAR_NAME[[[]]]_%J.err
#BSUB -q [[[]]]VAR_QUEUE[[[]]]

set -x #echo on

CONFIG_FDBG=../datasets/f_case_866x72_16p.nc
CONFIG_F30=../datasets/f_case_48602x72_512p.nc
CONFIG_F120=/gpfs/alpine/csc332/scratch/khl7265/FS_EVAL/e3sm/decom/FC5AV1C-H01A_ne120_oRRS18v3_21600p_raw.nc
CONFIG_GDBG=../datasets/g_case_cmpaso_16p.nc 
CONFIG_G=/gpfs/alpine/csc332/scratch/khl7265/FS_EVAL/e3sm/decom/GMPAS-NYF_T62_oRRS18to6v3_9600p_raw.nc
CONFIG_IDBG=../datasets/i_case_f19_g16_16p.nc 
CONFIG_I=/gpfs/alpine/csc332/scratch/khl7265/FS_EVAL/e3sm/decom/I1850GSWCNPRDCTCBC_hcru_hcru_1344p_raw.nc

IN_FDBG=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/F_DBG
IN_F30=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/F_30
IN_F120=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/F_120
IN_G=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/G
IN_I=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/I

E3SM_IO_DATE="VAR_E3SM_IO_DATE"
HDF5_LIB_PATH=VAR_HDF5_LIB_PATH
HDF5_LIB_DATE="VAR_HDF5_LIB_DATE"
ADIOS2_LIB_PATH=VAR_ADIOS2_LIB_PATH
ADIOS2_LIB_DATE="VAR_ADIOS2_LIB_DATE"
LOGVOL_LIB_PATH=VAR_LOGVOL_LIB_PATH
LOGVOL_LIB_DATE="VAR_LOGVOL_LIB_DATE"
PNC_LIB_PATH=VAR_PNC_LIB_PATH
NREC=VAR_RECS
PPN=VAR_PPN
RTL=VAR_RTL
CONFIG_POST=VAR_CONFIG
OUTDIR_ROOT=VAR_OUTDIR_ROOT
SUBFILEDIR_ROOT=VAR_SUBFILEDIR_ROOT

APPS=(e3sm_io)
FXS=(testfile_h0 testfile_h1)
CONFIG="CONFIG_${CONFIG_POST}"
CONFIG=${!CONFIG}
CONFIG_NAME=$(basename -- "${CONFIG}")
CONFIG_NAME="${CONFIG_NAME%.*}"
INDIR="IN_${CONFIG_POST}"
INDIR=${!INDIR}
DRIVERS=("hdf5_log log" "adios blob" "pnetcdf blob")
OPTIONS=(11111)
FFREQS=(VAR_RECS)
OPS=(VAR_OP)

NN=VAR_NN
NP=VAR_NP
let NA=NN*2
let NG=NN-1

# Print exp setting
echo "#%$=: e3sm_io_hash: VAR_E3SM_IO_HASH"
echo "#%$=: pnc_hash: VAR_PNC_HASH"
echo "#%$=: pnc_path: ${PNC_LIB_PATH}"
echo "#%$=: hdf5_hash: VAR_HDF5_HASH"
echo "#%$=: hdf5_path: ${HDF5_LIB_PATH}"
echo "#%$=: submit_date: VAR_SUBMIT_DATE"
echo "#%$=: exe_date: $(date)"
echo "#%$=: exp: e3sm_io"
echo "#%$=: repeation: VAR_RUNS"
echo "#%$=: output_folder: VAR_OUTDIR_ROOT"
echo "#%$=: e3sm_config_post: ${CONFIG_POST}"
echo "#%$=: e3sm_config: ${CONFIG}"
echo "#%$=: number_of_nodes: ${NN}"
echo "#%$=: number_of_proc: ${NP}"
echo '-----+-----++------------+++++++++--+---'

for DRIVER in "${DRIVERS[@]}"
do
    tmp=($DRIVER)
    API=${tmp[0]}
    STRATE=${tmp[1]}
    mkdir -p ${OUTDIR_ROOT}/${API}/${STRATE}/${CONFIG_NAME}
done

#export LD_LIBRARY_PATH=${HDF5_LIB_PATH}/lib:${PNC_LIB_PATH}/lib:${ADIOS2_LIB_PATH}/lib64:${LOGVOL_LIB_PATH}/lib:${LD_LIBRARY_PATH}
export PNETCDF_SHOW_PERFORMANCE_INFO=1
#export PNETCDF_DEFAULT_CHUNK_DIM="ncol : 14563 ; nbnd : 2 ; Time : 1 ; ilev : 73 ; lev : 72 ; chars : 64 ;nCells : 16384 ; nEdges : 16384 ; nVertices : 16384 ; nVertLevelsP1 : 81 ; nVertLevels : 80 ; StrLen : 64 ;"
export PNETCDF_HINTS="nc_zip_delay_init=1;nc_zip_nrec=1;nc_zip_buffer_size=0"

export H5VL_LOG_METADATA_SHARE=1
export H5VL_LOG_METADATA_MERGE=0
export H5VL_LOG_METADATA_ZIP=1
export H5VL_LOG_SEL_ENCODING=1
export H5VL_LOG_SUBFILING=1
export E3SM_IO_HDF5_USE_LOGVOL_WRITEN=1
export H5VL_LOG_SHOW_PROFILING_INFO=1
export H5VL_LOG_PRINT_MPI_INFO=1
export H5VL_LOG_NSUBFILES=-1

ulimit -c unlimited

TSTARTTIME=$(date +%s.%N)

for i in $(seq VAR_RUNS);
do
    for DRIVER in "${DRIVERS[@]}"
    do
        tmp=($DRIVER)
        API=${tmp[0]}
        STRATE=${tmp[1]}

        if [ "${API}" = "adios" ] ; then
            if [ "${CONFIG_POST}" = "I" ] ; then
                if [ "${API}" = "adios" ] ; then
                    let NG=NN*2 # ADIOS method must double #subfiles
                else
                    let NG=NN
                fi
            else
                let NG=NN-1 # F and G case use #nodes - 1
            fi
            # 8 aggregators
            #export E3SM_IO_HINTS="cb_node=8;cb_config_list=*:8;romio_cb_write=enable"
        elif [ "${STRATE}" = "canonical" ] ; then
            # pnc always flush on each records
            if [ "${API}" = "pnetcdf" ] ; then
                let FFREQ=1
            fi
            # 1 subfile
            let NG=1    
            # 64 aggregators
            #export E3SM_IO_HINTS="romio_cb_write=enable;"
            unset E3SM_IO_HINTS
        elif [ "${STRATE}" = "blob" ] ; then
            # subfile per node
            let NG=NN
            # 64 aggregators
            ##export E3SM_IO_HINTS="romio_cb_write=enable;"
            unset E3SM_IO_HINTS
        else
            # subfile per node
            let NG=1
            # 64 aggregators
            ##export E3SM_IO_HINTS="romio_cb_write=enable;"
            unset E3SM_IO_HINTS
        fi

        if [[ "${STRATE}" == "blob" || "${API}" == "hdf5_log" ]] ; then
            OUTDIR="${SUBFILEDIR_ROOT}/${API}/${STRATE}/${CONFIG_NAME}"
        else
            OUTDIR="${OUTDIR_ROOT}/${API}/${STRATE}/${CONFIG_NAME}"
        fi

        RDDIR="${OUTDIR}_${i}"

        for FX in ${FXS[@]}
        do
            echo "========================== E3SM-IO replay ${API} ${STRATE} ${OP} =========================="
            >&2 echo "========================== E3SM-IO replay ${API} ${STRATE} ${OP}=========================="
            
            echo "#%$=: exp: e3sm_io"
            echo "#%$=: runs: ${i}"
            echo "#%$=: config: ${CONFIG}"
            echo "#%$=: h_num: ${FX}"
            echo "#%$=: nrec: ${NREC}"
            echo "#%$=: api: ${API}"
            echo "#%$=: strategy: ${STRATE}"
            echo "#%$=: number_of_nodes: ${NN}"
            echo "#%$=: number_of_proc: ${NP}"

            rm -rf ${RDDIR}/${FX}_replay*
        
            STARTTIME=$(date +%s.%N)

            if [ "${API}" == "adios" ] ; then
                jsrun -X 1 -p ${NP} -n ${NA} -r 2 -d plane:42 -c 21 -g 0 -b packed:smt:2 -l cpu-cpu --stdio_mode collected /ccs/home/khl7265/csc332/pio/eval_build/tools/adios2pio-nm/adios2pio-nm.exe --bp-file=${RDDIR}/${FX}.bp --nc-file=${RDDIR}/${FX}_replay --pio-format=pnetcdf
            elif [ "${API}" == "hdf5_log" ] ; then
                jsrun -X 1 -p ${NP} -n ${NA} -r 2 -d plane:42 -c 21 -g 0 -b packed:smt:2 -l cpu-cpu --stdio_mode collected ${LOGVOL_LIB_PATH}/bin/h5lreplay -i ${RDDIR}/${FX} -o ${RDDIR}/${FX}_replay
            else
                if [ "${STRATE}" == "blob" ] ; then
                    jsrun -X 1 -p ${NP} -n ${NA} -r 2 -d plane:42 -c 21 -g 0 -b packed:smt:2 -l cpu-cpu --stdio_mode collected ../utils/pnetcdf_blob_replay -i ${RDDIR}/${FX} -o ${RDDIR}/${FX}_replay
                else
                    jsrun -X 1 -p ${NP} -n ${NA} -r 2 -d plane:42 -c 21 -g 0 -b packed:smt:2 -l cpu-cpu --stdio_mode collected ../src/e3sm_io -k -i ${RDDIR}/${FX} -a ${API} -f ${FX} -r ${NREC} -y ${NREC} -x ${STRATE} ${CONFIG} 
                fi
            fi

            ENDTIME=$(date +%s.%N)
            TIMEDIFF=$(echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print $1"."$2}')

            echo "#%$: exe_time: $TIMEDIFF"

            ls -Rlah ${RDDIR}

            echo '-----+-----++------------+++++++++--+---'
        done
    done
done

ENDTIME=$(date +%s.%N)
TIMEDIFF=$(echo "$ENDTIME - $TSTARTTIME" | bc | awk -F"." '{print $1"."$2}')
echo "-------------------------------------------------------------"
echo "total_exe_time: $TIMEDIFF"
echo "-------------------------------------------------------------"