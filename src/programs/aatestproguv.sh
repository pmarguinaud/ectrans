#!/bin/bash
#SBATCH --gres=gpu:4 
#SBATCH --exclusive 
#SBATCH -A DestE_330_24 
#SBATCH -N1 
#SBATCH -p boost_usr_prod
#SBATCH --export=NONE

export MODULEPATH=~pmarguin/install/nvidia/hpc_sdk/modulefiles:$MODULEPATH
module load nvhpc-hpcx/24.5
source $NVHPC_ROOT/comm_libs/12.4/hpcx/hpcx-2.19/hpcx-mt-init.sh hpcx_load

export SLURM_EXPORT_ENV=ALL

set -x

cd /leonardo/home/userexternal/pmarguin/gpupack/pack/ectrans/ectrans-build


if [ 0 -eq 1 ]
then
for i in 1 4
do
~/install/mpiauto/mpiauto --nouse-slurm-mpi -np $i -- ./bin/ectrans-benchmark-gpu-dp --truncation 31 -f 1 -l 15 --vordiv --norms -v > $i.txt 2>&1
done
fi


for NPROC in 1 4
do
  ~/install/mpiauto/mpiauto \
   --nouse-slurm-mpi --verbose --wrap --wrap-stdeo --nnp $NPROC --nn 1 --openmp 1 \
   -- ./bin/aatestproguv-gpu-dp --namelist aatestproguv.nam \
   --u-file DATA.U.0001.DAT.txt --v-file DATA.V.0001.DAT.txt --t-file DATA.T.0001.DAT.txt
  \mv ZSPVORG.DAT ZSPVORG.${NPROC}.DAT
  \mv ZSPDIVG.DAT ZSPDIVG.${NPROC}.DAT
  \mv ZSP3AG.DAT ZSP3AG.${NPROC}.DAT
done

