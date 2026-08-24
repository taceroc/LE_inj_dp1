#!/bin/bash                                                                                                                                     
#SBATCH --account=dessn
#SBATCH --time=24:00:00
#SBATCH --qos=shared
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --cpus-per-task=32
#SBATCH --output=dp1_coadd-%j.out


export STACKCVMFS=/cvmfs/sw.lsst.eu/almalinux-x86_64/lsst_distrib
export LSST_STACK_VERSION=v29.1.1 #w_2025_31

module load cpu

source $STACKCVMFS/$LSST_STACK_VERSION/loadLSST-ext.bash
setup -t v29_1_1 lsst_distrib


export DAF_BUTLER_REPOSITORY_INDEX=/global/cfs/cdirs/lsst/production/gen3/shared/data-repos.yaml

export OMP_NUM_THREADS=1


# make the custom coadd QuantumGraph visualization
pipetask build \
-p $DRP_PIPE_DIR/pipelines/LSSTComCam/DRP-v2-compat.yaml#makeDirectWarp,assembleDeepCoadd,makePsfMatchedWarp,selectDeepCoaddVisits \
--pipeline-dot pipeline.dot; \
dot pipeline.dot -Tpdf > dp1_coadd.pdf

# remove temporary file
rm pipeline.dot 

# specify the directory for output log files
LOGDIR=logs

# make the directory for output log files
mkdir $LOGDIR

# run the custom coaddition
LOGFILE=$LOGDIR/do1_coadd-logfile.log; \
date | tee $LOGFILE; \
pipetask --long-log --log-file $LOGFILE run \
-b /repo/DP1 \
-i LSSTComCam/DP1 \
-o u/taceroc/dp1/custom_coadd_window2_60651-60655-4848-90-r \
-p $DRP_PIPE_DIR/pipelines/LSSTComCam/DRP-v2-compat.yaml#makeDirectWarp,assembleDeepCoadd,makePsfMatchedWarp,selectDeepCoaddVisits \
-c makeDirectWarp:useVisitSummaryPsf=False \
-c makeDirectWarp:useVisitSummaryPhotoCalib=False \
-c makeDirectWarp:useVisitSummaryWcs=False \
-c makeDirectWarp:connections.calexp_list="visit_image" \
-d "tract=4848 AND patch=90 AND band='r' AND visit IN (2024120600070,2024120600071,2024120600072,2024120600073,2024120600074,2024120600095,2024120600096,2024120600097,2024120600098,2024120600099,2024120600239,2024120600240,2024120600241,2024120600242,2024120600243,2024120600264,2024120600265,2024120600266,2024120600267,2024120600268,2024120600272,2024120600273,2024120600274,2024120600275,2024120600276,2024120700294,2024120700295,2024120700296,2024120700297,2024120700298,2024120700299,2024120700300,2024120700301,2024120700302,2024120700303,2024120800345,2024120800346,2024120800347,2024120800348,2024120800349,2024120900309,2024120900310,2024120900311,2024120900312,2024120900313,2024121000402,2024121000403,2024121000404,2024121000405,2024121000406,2024121000417,2024121000418,2024121000419,2024121000420,2024121000421,2024121000422,2024121000423,2024121000424,2024121000425) AND skymap='lsst_cells_v1'"; \
date | tee -a $LOGFILE