#!/bin/sh

## parameters
chrlensfile="../data/chromosome_sizes.txt"
res=5000
threads=50
n_samples=3
outfolder="../output"
n_runs=1
EXE_PATH=".././bin/sBIF"
total_files=$(find ../data/folding_input -name "*.txt" | wc -l | xargs)
count=1

for interfile in ../data/folding_input/*.txt; do

    filename=$(basename "$interfile")
    
    # Extract cell_line, chromosome, start, and end from the filename
    cell_line=$(echo "$filename" | cut -d'.' -f1)
    chrom=$(echo "$filename" | cut -d'.' -f2)
    start=$(echo "$filename" | cut -d'.' -f3)
    end=$(echo "$filename" | cut -d'.' -f4 | sed 's/.txt//')
    job_prefix="$chrom"

    ## command
    cmd="$EXE_PATH -i $interfile -c $chrom -l $chrlensfile -s $start -e $end -ns $n_samples -nr $n_runs -cl $cell_line -o $outfolder -r $res -do false -j $job_prefix -p $threads"
    
    echo "Processing file $count of $total_files: $filename"
    echo "Running command: $cmd"
    
    $cmd
    count=$((count + 1)) 
done

echo "Done."