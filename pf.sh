mkdir output

outfileprefix="output/pf1"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- $(date) Running pf1  >> jobs.txt
echo ---------------- >> jobs.txt


# function pf1(rsa_idx, use_smooth, lateralized, nperms, parcel_idx, subbatch_size)

declare -a fn_calls=(
                     "pf1"
                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p shared --mem 4001 -t 4-4:00 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo pf1.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

