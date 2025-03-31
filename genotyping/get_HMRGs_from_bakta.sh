#Winkler MA
#2023-08-22
#Get heavy metal and thermoresistance genes from bakta annotated assemblies

# Make header
echo -e "Sequence_ID\tType\tStart\tStop\tStrand\tLocus_Tag\tGene\tProduct\tDbXrefs"  > Bakta_v18_annotations.tsv

#Append all non-header lines from .tsv files 
grep -h -v -e "#Sequence Id" -e "# Annotated with Bakta" -e "# Software: v1.8.1" -e "# Database: v5.0, full" -e "# DOI: 10.1099/mgen.0.000685" -e "# URL: github.com/oschwengers/bakta" hybrid_bakta_tsvs/*.tsv >> Bakta_v18_annotations.tsv
grep -h -v -e "#Sequence Id" -e "# Annotated with Bakta" -e "# Software: v1.8.1" -e "# Database: v5.0, full" -e "# DOI: 10.1099/mgen.0.000685" -e "# URL: github.com/oschwengers/bakta" spades_bakta_tsvs/*.tsv >> Bakta_v18_annotations.tsv

# Filter out HRMGs only
# "Gene" is the 7th column
cat Bakta_v18_3255_Illumina__2023-08-22.txt  | head -6 >> hmrgs_3255_best_assemblies.txt ; #What does this do?

# "Gene" is the 7th column
awk -F ‘\t’ ‘$7 ~ /nik|cpx|cus|cut|cue|cop|chr|zra|cor|rcn|cbt|cbi|dpp|opp|mnt|sit|mts|znu|zup|zit|zur|teh|ydc|yeb|tsg|yff|gol|ars|fie|yii|mer|ncr|pco|sil|ter|cad|czc|dsb|znt|cds/’ Bakta_v18_3255_best_assemblies__2023-08-22.txt >> hmrgs_3255_best_assemblies.txt ;

#Then parse in R with: parse_HMRGs.R
