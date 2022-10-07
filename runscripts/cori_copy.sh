#!/bin/bash

set -x #echo on

CONFIG_FDBG=../datasets/f_case_866x72_16p.nc
CONFIG_F30=../datasets/f_case_48602x72_512p.nc
CONFIG_F120=${CSCRATCH}/FS_128_16M/e3sm_data/decom/FC5AV1C-H01A_ne120_oRRS18v3_21600p_raw.nc
CONFIG_G=${CSCRATCH}/FS_128_16M/e3sm_data/decom/GMPAS-NYF_T62_oRRS18to6v3_9600p_raw.nc
CONFIG_IDBG=../datasets/i_case_f19_g16_16p.nc 
CONFIG_I=${CSCRATCH}/FS_128_16M/e3sm_data/decom/I1850GSWCNPRDCTCBC_hcru_hcru_1344p_raw.nc

IN_FDBG=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/F_DBG
IN_F30=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/F_30
IN_F120=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/F_120
IN_G=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/G
IN_I=/global/cscratch1/sd/khl7265/FS_64_1M/E3SM/realdata/I

CONFIGS= (F120 G I)
OUTDIR_ROOT=${CSCRATCH}/FS_128_16M/E3SM
SUBFILEDIR_ROOT=${CSCRATCH}/FS_8_1M/E3SM

TSTARTTIME=$(date +%s.%N)

for c in "${CONFIGS[@]}"
do
    CONFIG="CONFIG_${c}"
    CONFIG=${!CONFIG}
    CONFIG_NAME=$(basename -- "${CONFIG}")
    CONFIG_NAME="${CONFIG_NAME%.*}"

    for DRIVER in "${DRIVERS[@]}"
    do
        tmp=($DRIVER)
        API=${tmp[0]}
        STRATE=${tmp[1]}

        if [[ "${STRATE}" == "blob" || "${API}" == "hdf5_log" ]] ; then
            OUTDIR="${SUBFILEDIR_ROOT}/${API}/${STRATE}/${CONFIG_NAME}"
        else
            OUTDIR="${OUTDIR_ROOT}/${API}/${STRATE}/${CONFIG_NAME}"
        fi

        RDDIR="${OUTDIR}_${i}"
        rm -rf ${RDDIR}
        cp -R ${OUTDIR} ${RDDIR}
    done
done

ENDTIME=$(date +%s.%N)
TIMEDIFF=$(echo "$ENDTIME - $TSTARTTIME" | bc | awk -F"." '{print $1"."$2}')
echo "-------------------------------------------------------------"
echo "total_exe_time: $TIMEDIFF"
echo "-------------------------------------------------------------"
