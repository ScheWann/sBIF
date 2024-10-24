#!/bin/bash

##parameters
interfile="../data/chr1.27120000.27895000.txt"
chrlensfile="../data/chromosome_sizes.txt"
chrom="chr1"
start=27120000
end=27895000
res=5000
nsamp=5000
job_prefix="dros"
threads=1
EXE_PATH="./../bin/sBIF"

##command
cmd="$EXE_PATH -i $interfile -c $chrom -l $chrlensfile -s $start -e $end -ns $nsamp -r $res -j $job_prefix -p $threads "
echo $cmd
$cmd 
echo "Done."
