for f in $(ls /Users/marith/crossniche/fasta/*fasta | sed "s|/Users/marith/crossniche/fasta/||g"); do 
    kraken2 --db /Users/marith/metagenomic_classification_databases/gtdb_r202/kraken2/016gb/ \
    /Users/marith/crossniche/fasta/${f} \
    --report /Users/marith/crossniche/fasta/kraken2_results/${f}_output.fasta \
    --report-minimizer-data 
done
