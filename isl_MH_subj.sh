
mkdir output

num_particles=10;

echo ---------------- >> jobs.txt
echo --- Running isl_MH_subj.sh for num_particles ${num_particles} >> jobs.txt
echo ---------------- >> jobs.txt

for subj in {126..127}
do
        
    outfileprefix="output/isl_MH_subj_${subj}_num_particles_${num_particles}"
    echo File prefix = $outfileprefix

    # send the job to NCF
    #
    sbatch_output=`sbatch -p shared --mem 4001 -t 0-4:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'isl_MH_subj(${subj}, ${num_particles});exit'"`
    # for local testing
    #sbatch_output=`echo Submitted batch job 88725418`
    echo $sbatch_output

    # Append job id to jobs.txt
    #
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo isl_MH_subj.sh for subject ${subj}: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    echo watch job status with: sacct -j ${job_id}
    echo watch output with: tail -f ${outfileprefix}_${job_id}.out
    echo watch error with: tail -f ${outfileprefix}_${job_id}.err

    sleep 1
done
