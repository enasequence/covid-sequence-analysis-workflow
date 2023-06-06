
workDir=${1}

echo $workDir
mkdir -p $workDir/.hq-server
SLURM_NNODES=12
SLURM_CPUS_PER_TASK=4
# Set the directory which hyperqueue will use 
cd $workDir
hq server start &
srun --cpu-bind=none --hint=nomultithread --mpi=none -N $SLURM_NNODES -n $SLURM_NNODES hq worker start --cpus=$SLURM_CPUS_PER_TASK &
num_up=$(hq worker list | grep RUNNING | wc -l)
while true; do
    echo "Checking if workers have started"
    if [[ $num_up -eq $SLURM_NNODES ]];then
            echo "Workers started"
            break
    fi
    echo "$num_up/$SLURM_NNODES workers have started"
    sleep 1
    num_up=$(hq worker list | grep RUNNING | wc -l)
done