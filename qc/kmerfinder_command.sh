for f in $(ls *fasta ); do 
    singularity exec \
    -B /Users/marith/crossniche/fasta \
    /Users/marith/software/kmerfinder/kmerfinder_3.0.2.sif \
    kmerfinder.py \
    -i /Users/marith/crossniche/fasta/${f} \
    -o /Users/marith/crossniche/fasta/kmerfinder_results/${f}_kmerfinder \
    -db /Users/marith/software/kmerfinder/bacteria/bacteria.ATG \
    -tax /Users/marith/software/kmerfinder/bacteria/bacteria.tax -x 
done
