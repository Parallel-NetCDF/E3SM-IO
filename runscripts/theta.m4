changequote(`[[[', `]]]')dnl
changecom([[[###]]], [[[$$$]]])dnl
#!/bin/bash
#COBALT -q VAR_QUEUE
#COBALT -n VAR_NN
#COBALT -t VAR_TL
#COBALT -o VAR_NAME.txt
#COBALT -e VAR_NAME.err
#COBALT -A VAR_ACC

CONFIG_FDBG=datasets/f_case_866x72_16p.nc
CONFIG_F30=datasets/f_case_48602x72_512p.nc
CONFIG_F120=/lus-projects/CSC250STDM12/khou/FS_48_1M/e3sm/decom/FC5AV1C-H01B_ne120_oRRS18v3_21600p.nc
CONFIG_G=/lus-projects/CSC250STDM12/khou/FS_48_1M/e3sm/decom/GMPAS-NYF_T62_oRRS18to6v3_9600p.nc

E3SM_IO_DATE="VAR_E3SM_IO_DATE"
HDF5_LIB_PATH=VAR_HDF5_LIB_PATH
HDF5_LIB_DATE="VAR_HDF5_LIB_DATE"
PNC_LIB_PATH=VAR_PNC_LIB_PATH
LOGVOL_LIB_PATH=VAR_LOGVOL_LIB_PATH
LOGVOL_LIB_DATE="VAR_LOGVOL_LIB_DATE"
NREC=VAR_RECS
PPN=VAR_PPN
RTL=VAR_RTL
CONFIG_POST=VAR_CONFIG
OUTDIR_ROOT=VAR_OUTDIR_ROOT
INFILE=VAR_INFILE

APP=e3sm_io
HXS=(VAR_HX)
CONFIG="CONFIG_${CONFIG_POST}"
CONFIG=${!CONFIG}
APIS=(pnc hdf5)
OPS=(VAR_OP)

NN=VAR_NN
let NP=NN*PPN

# Print exp setting
echo "#%$=: e3sm_io_hash: VAR_E3SM_IO_HASH"
echo "#%$=: pnc_hash: VAR_PNC_HASH"
echo "#%$=: pnc_path: ${PNC_LIB_PATH}"
echo "#%$=: hdf5_hash: VAR_HDF5_HASH"
echo "#%$=: hdf5_path: ${HDF5_LIB_PATH}"
echo "#%$=: logvol_hash: VAR_LOGVOL_HASH"
echo "#%$=: logvol_path: ${LOGVOL_LIB_PATH}"
echo "#%$=: submit_date: VAR_SUBMIT_DATE"
echo "#%$=: exe_date: $(date)"
echo "#%$=: exp: e3sm_io"
echo "#%$=: repeation: VAR_RUNS"
echo "#%$=: output_folder: VAR_OUTDIR_ROOT"
echo "#%$=: e3sm_config_post: ${CONFIG_POST}"
echo "#%$=: e3sm_config: ${CONFIG}"
echo "#%$=: number_of_nodes: ${NN}"
echo "#%$=: number_of_proc: ${NP}"
echo '-----+-----++----***--------+++++**++++--+---'

for API in ${APIS[@]}
do
    echo "mkdir -p ${OUTDIR_ROOT}/${API}"
    mkdir -p ${OUTDIR_ROOT}/${API}
done

export LD_LIBRARY_PATH=${HDF5_LIB_PATH}/lib:${PNC_LIB_PATH}/lib:${LOGVOL_LIB_PATH}/lib:${LD_LIBRARY_PATH}
#export HDF5_PLUGIN_PATH=${HOME}/.local/lib
export H5VL_LOG_METADATA_MERGE=1
export H5VL_LOG_METADATA_ZIP=1
export H5VL_LOG_SEL_ENCODING=1

ulimit -c unlimited

TSTARTTIME=`date +%s.%N`

for i in $(seq VAR_RUNS);
do
    for API in ${APIS[@]}
    do
        OUTDIR=${OUTDIR_ROOT}/${API}/
        echo "rm -f ${OUTDIR}/*"
        rm -f ${OUTDIR}/*

        for HX in ${HXS[@]}
        do
            for (( j=0; j<${#OPS}; j++ ));
            do
                OP=${OPS:$j:1}

                OUTDIR="${OUTDIR_ROOT}/${API}/${LAYOUTS}"
                OUTFILE="${OUTDIR_ROOT}/${API}/${LAYOUTS}/output.bin"

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
                CUR_LOGVOL_DATE=$(stat -c %Y ${LOGVOL_LIB_PATH}/lib/libH5VL_log.so.0.0.0)
                if [[ "${CUR_LOGVOL_DATE}" != "${LOGVOL_LIB_DATE}" ]]; then
                    echo "Warning: libH5VL_log.so.0.0.0 changed after submission"
                fi

                echo "#%$=: exp: e3sm_io"
                echo "#%$=: runs: ${i}"
                echo "#%$=: app: ${APP}"
                echo "#%$=: config: ${CONFIG}"
                echo "#%$=: h_num: ${HX}"
                echo "#%$=: nrec: ${NREC}"
                echo "#%$=: api: ${API}"
                echo "#%$=: strategy: ${STRATE}"
                echo "#%$=: flush_freq: ${FFREQ}"
                echo "#%$=: logvol_dwrite_n: ${DWRITE_N}"
                echo "#%$=: logvol_meta_encode: ${META_ENCODE}"
                echo "#%$=: logvol_meta_zip: ${META_ZIP}"
                echo "#%$=: logvol_meta_share: ${META_SHARE}"
                echo "#%$=: logvol_subfiling: ${SUBFILING}"
                echo "#%$=: operation: ${OP}"
                echo "#%$=: delay_init: 1"
                echo "#%$=: number_of_nodes: ${NN}"
                echo "#%$=: number_of_proc: ${NP}"
                                        
                STARTTIME=`date +%s.%N`

                if [ "$OP" = "w" ] ; then
                    echo "aprun -n ${NP} -N ${PPN} -t ${RTL} -e H5VL_LOG_METADATA_MERGE=${H5VL_LOG_METADATA_MERGE} -e H5VL_LOG_METADATA_ZIP=${H5VL_LOG_METADATA_ZIP} -e H5VL_LOG_SEL_ENCODING=${H5VL_LOG_SEL_ENCODING} ./${APP} -k -o ${OUTDIR} -n -a ${API} -f ${HX} -r ${NREC} ${CONFIG}"
                    aprun -n ${NP} -N ${PPN} -t ${RTL} -e H5VL_LOG_METADATA_MERGE=${H5VL_LOG_METADATA_MERGE} -e H5VL_LOG_METADATA_ZIP=${H5VL_LOG_METADATA_ZIP} -e H5VL_LOG_SEL_ENCODING=${H5VL_LOG_SEL_ENCODING} ./${APP} -k -o ${OUTDIR} -n -a ${API} -f ${HX} -r ${NREC} ${CONFIG}
                else
                    echo "aprun -n ${NP} -N ${PPN} -t ${RTL} -e H5VL_LOG_METADATA_MERGE=${H5VL_LOG_METADATA_MERGE} -e H5VL_LOG_METADATA_ZIP=${H5VL_LOG_METADATA_ZIP} -e H5VL_LOG_SEL_ENCODING=${H5VL_LOG_SEL_ENCODING} ./${APP} -k -o ${OUTDIR} -n -a ${API} -f ${HX} -r ${NREC} -R ${CONFIG}"
                    aprun -n ${NP} -N ${PPN} -t ${RTL} -e H5VL_LOG_METADATA_MERGE=${H5VL_LOG_METADATA_MERGE} -e H5VL_LOG_METADATA_ZIP=${H5VL_LOG_METADATA_ZIP} -e H5VL_LOG_SEL_ENCODING=${H5VL_LOG_SEL_ENCODING} ./${APP} -k -o ${OUTDIR} -n -a ${API} -f ${HX} -r ${NREC} -R ${CONFIG}
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
        done
    done
done

ENDTIME=`date +%s.%N`
TIMEDIFF=`echo "$ENDTIME - $TSTARTTIME" | bc | awk -F"." '{print $1"."$2}'`
echo "-------------------------------------------------------------"
echo "total_exe_time: $TIMEDIFF"
echo "-------------------------------------------------------------"
