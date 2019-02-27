mkdir output

echo ---------------- >> jobs.txt
echo --- Running fit_grid  >> jobs.txt
echo ---------------- >> jobs.txt


declare -a fn_calls=(
                    "fit_grid"
                     )


outfileprefix="output/fit_grid"
echo File prefix = $outfileprefix



for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf_holy --mem 50001 -t 20-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo running $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done
