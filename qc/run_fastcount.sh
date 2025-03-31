#Read stats

#Install
cd ~/Programs/ ; git clone https://github.com/rrwick/MinION-desktop ; make ; 

#Run
cd ./Illumina/FASTQ ;
~/Programs/MinION-desktop/fast_count >> KLEBGAP_fastcount.txt #Print header
~/Programs/MinION-desktop/fast_count *gz >> KLEBGAP_fastcount.txt #Get fast count stats for each file

#To estimate average read depth per genome, take total number of bases from fast_count and divide by the total length of the assembly from quast.
