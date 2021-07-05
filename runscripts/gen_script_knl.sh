#!/bin/bash

# Defaults
NN=337
PPN=65
RUNS=1
TL=30
CONFIG=F120
RECS=1
OP=w
HX=0
VER="release"
HOSTNAME_PREFIX=${HOSTNAME:0:4}
if [[ "$HOSTNAME_PREFIX" == "cori" ]]; then
    # Cori
    HDF5_SRC_PATH=/global/homes/k/khl7265/hdf5/1.12.0
    HDF5_LIB_PATH=/global/homes/k/khl7265/.local/hdf5/1.12.0
    PNC_SRC_PATH=/global/homes/k/khl7265/pnetcdf/master
    PNC_LIB_PATH=/global/homes/k/khl7265/.local/pnetcdf/master
    ADIOS2_SRC_PATH=/global/homes/k/khl7265/adios2/2.7.1
    ADIOS2_LIB_PATH=/global/homes/k/khl7265/.local/adios2/2.7.1_knl
    OUTPATH_ROOT=/global/cscratch1/sd/khl7265/FS_1_1M/E3SM
    INFILE=/global/cscratch1/sd/wkliao/FS_1M_32/ND_1951_256K.1.h5
    ACC=m2956
    QUEUE=regular
    SCRIPT=cori_knl
    NODE_TYPE=knl
    RTL=3
elif [[ "$HOSTNAME_PREFIX" == "thet" ]]; then
    # Theta (Not supported yet)
    HDF5_SRC_PATH=/gpfs/mira-home/khou/hdf5/1.12.0
    HDF5_LIB_PATH=/gpfs/mira-home/khou/.local/hdf5/1.12.0
    PNC_SRC_PATH=/gpfs/mira-home/khou/pnetcdf/master
    PNC_LIB_PATH=/gpfs/mira-home/khou/.local/pnetcdf/master
    ADIOS2_SRC_PATH=/gpfs/mira-home/khou/logvol/adios2/2.7.1
    ADIOS2_LIB_PATH=/gpfs/mira-home/khou/.local/log_io_vol/adios2/2.7.1
    OUTPATH_ROOT=/lus-projects/CSC250STDM12/khou/FS_48_1M/e3sm
    INFILE=/global/cscratch1/sd/wkliao/FS_1M_32/ND_1951_256K.1.h5
    ACC=CSC250STDM12
    SCRIPT=theta
    QUEUE=default
    NODE_TYPE=knl
    RTL=180
elif [[ "$HOSTNAME_PREFIX" == "Ubun" ]]; then
    # Dev VM
    HDF5_SRC_PATH=${HOME}/Desktop/hdf5/1.12.0
    HDF5_LIB_PATH=${HOME}/.local/hdf5/1.12.0
    PNC_SRC_PATH=${HOME}/Desktop/pnetcdf/master
    PNC_LIB_PATH=${HOME}/.local/pnetcdf/master
    LOGVOL_SRC_PATH=${HOME}/.local/logvol/master
    LOGVOL_LIB_PATH=${HOME}/.local/logvol/master
    ADIOS2_SRC_PATH=${HOME}/.local/adios2/2.7.1
    ADIOS2_LIB_PATH=${HOME}/.local/adios2/2.7.1
    OUTPATH_ROOT=${HOME}/Desktop/nova_data/output
    INFILE=${HOME}/Desktop/nova_data/neardet_r00012096_s01_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5
    ACC=test_acc
    SCRIPT=theta
    NODE_TYPE=local
fi

# Parse options
OTHER_PARAMS=""
while (( "$#" )); do
    case "$1" in
        -n|--nn)
            NN=$2
            shift 
            shift 
            ;;
        -N|--node)
            NODE_TYPE=$2
            shift 
            shift 
            ;;
        -d|--debug)
            DEBUG=1
            shift 
            shift 
            ;;
        -p|--ppn)
            PPN=$2
            shift 
            shift 
            ;;
        -v|--version)
            VER=$2
            shift 
            shift 
            ;;
        -r|--repeation)
            RUNS=$2
            shift 
            shift 
            ;;
        -R|--record)
            RECS=$2
            shift 
            shift 
            ;;
        -a|--account)
            ACC=$2
            shift 
            shift 
            ;;
        -T|--timeout)
            TL=$2
            shift 
            shift 
            ;;
        -t|--runtimeout)
            RTL=$2
            shift 
            shift 
            ;;
        -c|--config)
            CONFIG=$2
            shift 
            shift 
            ;;
        -O|--op)
            OP=$2
            shift 
            shift 
            ;;
        -s|--script)
            SCRIPT=$2
            shift 
            shift 
            ;;
        -o|--outdir)
            OUTPATH_ROOT=$2
            shift 
            shift 
            ;;
        -h|--hx)
            HX=$2
            shift 
            shift 
            ;;
        -*|--*=) # unsupported flags
            echo "Error: Unsupported flag $1" >&2
            exit 1
            ;;
        *) # preserve positional arguments
            OTHER_PARAMS="$PARAMS $1"
            shift
            ;;
    esac
done

if [[ "${DEBUG}" == "1" ]]; then
    NN=1
    RUNS=1
    CONFIG=FDBG
    if [[ "$HOSTNAME_PREFIX" == "cori" ]]; then
        QUEUE=debug
    elif [[ "$HOSTNAME_PREFIX" == "thet" ]]; then
        QUEUE=debug-flat-quad
    fi
fi

# Set external lib path
if [[ "$HOSTNAME_PREFIX" == "cori" ]]; then
    # Cori
    if [[ "$VER" == "profiling" ]]; then
        LOGVOL_SRC_PATH=/global/homes/k/khl7265/logvol/profiling
        LOGVOL_LIB_PATH=/global/homes/k/khl7265/.local/log_io_vol/profiling
    else
        LOGVOL_SRC_PATH=/global/homes/k/khl7265/logvol/master
        LOGVOL_LIB_PATH=/global/homes/k/khl7265/.local/log_io_vol/master
    fi
elif [[ "$HOSTNAME_PREFIX" == "thet" ]]; then
    # Theta (Not supported yet)
        if [[ "$VER" == "debug" ]]; then
        LOGVOL_SRC_PATH=/gpfs/mira-home/khou/logvol/profiling
        LOGVOL_LIB_PATH=/gpfs/mira-home/khou/.local/log_io_vol/profiling
    else
        LOGVOL_SRC_PATH=/gpfs/mira-home/khou/logvol/master
        LOGVOL_LIB_PATH=/gpfs/mira-home/khou/.local/log_io_vol/master
    fi
fi

NP=NN*PPN

E3SM_IO_HASH=$(git --git-dir .git rev-parse HEAD)
E3SM_IO_HASH_PREFIX=${E3SM_IO_HASH:0:5}
E3SM_IO_DATE=$(stat -c %Y ./e3sm_io)

HDF5_HASH=$(git --git-dir ${HDF5_SRC_PATH}/.git rev-parse HEAD)
HDF5_HASH_PREFIX=${HDF5_HASH:0:5}
HDF5_DATE=$(stat -c %Y ${HDF5_LIB_PATH}/lib/libhdf5.so.200.0.0)

ADIOS2_HASH=$(git --git-dir ${ADIOS2_SRC_PATH}/.git rev-parse HEAD)
ADIOS2_HASH_PREFIX=${ADIOS2_HASH:0:5}
ADIOS2_DATE=$(stat -c %Y ${ADIOS2_LIB_PATH}/lib/libadios2_c_mpi.so.2.7.1)

PNC_HASH=$(git --git-dir ${PNC_SRC_PATH}/.git rev-parse HEAD)
PNC_HASH_PREFIX=${PNC_HASH:0:5}

LOGVOL_HASH=$(git --git-dir ${LOGVOL_SRC_PATH}/.git rev-parse HEAD)
LOGVOL_HASH_PREFIX=${LOGVOL_HASH:0:5}
LOGVOL_DATE=$(stat -c %Y ${LOGVOL_LIB_PATH}/lib/libH5VL_log.so.0.0.0)

FNAME="${CONFIG}_logvl${LOGVOL_HASH_PREFIX}_e3sm${E3SM_IO_HASH_PREFIX}_hdf${HDF5_HASH_PREFIX}_pnc${PNC_HASH_PREFIX}_adios2${ADIOS2_HASH_PREFIX}_${HOSTNAME_PREFIX}_knl_${NN}_${PPN}"

SUBMIT_DATE=$(date)

GEN_SCRIPT="m4 -D EXP_NAME=\"${FNAME}\" -D EXP_NN=\"${NN}\" -D EXP_PPN=\"${PPN}\" -D EXP_TL=\"${TL}\" -D EXP_ACC=\"${ACC}\" -D EXP_NODE_TYPE=\"${NODE_TYPE}\" \
-D EXP_OUTDIR_ROOT=\"${OUTPATH_ROOT}\" -D EXP_INFILE=\"${INFILE}\" -D EXP_CONFIG=\"${CONFIG}\" -D EXP_OP=\"${OP}\" -D EXP_HX=\"${HX}\" -D EXP_RECS=\"${RECS}\" \
-D EXP_RUNS=\"${RUNS}\" -D EXP_RTL=\"${RTL}\" -D EXP_SUBMIT_DATE=\"${SUBMIT_DATE}\"  -D EXP_QUEUE=\"${QUEUE}\" \
-D EXP_E3SM_IO_HASH=\"${E3SM_IO_HASH}\" -D EXP_E3SM_IO_DATE=\"${E3SM_IO_DATE}\" \
-D EXP_HDF5_LIB_PATH=\"${HDF5_LIB_PATH}\" -D EXP_HDF5_LIB_DATE=\"${HDF5_DATE}\" -D EXP_HDF5_HASH=\"${HDF5_HASH}\" \
-D EXP_ADIOS2_LIB_PATH=\"${ADIOS2_LIB_PATH}\" -D EXP_ADIOS2_LIB_DATE=\"${ADIOS2_DATE}\" -D EXP_ADIOS2_HASH=\"${ADIOS2_HASH}\" \
-D EXP_PNC_LIB_PATH=\"${PNC_LIB_PATH}\" -D EXP_PNC_HASH=\"${PNC_HASH}\" \
-D EXP_LOGVOL_LIB_PATH=\"${LOGVOL_LIB_PATH}\" -D EXP_LOGVOL_HASH=\"${LOGVOL_HASH}\" -D EXP_LOGVOL_LIB_DATE=\"${LOGVOL_DATE}\" \
\"${SCRIPT}\".m4 > \"${FNAME}\".sl"

echo ${GEN_SCRIPT}
eval "${GEN_SCRIPT}"
chmod 755 ${FNAME}.sl
