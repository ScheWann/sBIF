#!/bin/sh

##parameters
chrlensfile="../data/chromosome_sizes.txt"
res=5000
threads=50
EXE_PATH="./../bin/sBIF"

total_files=$(find ../data/folding_input -name "*.txt" | wc -l | xargs)
count=1

for interfile in ../data/folding_input/*.txt; do

    filename=$(basename "$interfile")
    
    # Extract the chromosome (job_prefix) from the filename
    chrom=$(echo "$filename" | cut -d'.' -f1)
    job_prefix="$chrom"

    start=$(echo "$filename" | cut -d'.' -f2)
    end=$(echo "$filename" | cut -d'.' -f3 | sed 's/.txt//')

    ##command
    cmd="$EXE_PATH -i $interfile -c $chrom -l $chrlensfile -s $start -e $end -r $res -j $job_prefix -p $threads"
    
    echo "Processing file $count of $total_files: $filename"
    echo "Running command: $cmd"
    
    $cmd
    count=$((count + 1)) 
done

echo "Done."
