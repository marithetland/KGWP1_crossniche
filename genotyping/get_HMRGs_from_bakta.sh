#Winkler MA
#2023-08-22
#Get heavy metal and thermoresistance genes from bakta annotated assemblies

#Make header
echo -e "Genome\tSequence_ID\tType\tStart\tStop\tStrand\tLocus_Tag\tGene\tProduct\tDbXrefs" > Bakta_v18_annotations.tsv

#Remove headers, add column for genome name
for file in hybrid_bakta_tsvs/*.tsv; do
    genome=$(basename "$file" .tsv)
    grep -v -e "^#" "$file" | awk -v g="$genome" -F '\t' 'BEGIN {OFS="\t"} {print g, $0}' >> Bakta_v18_annotations.tsv
done

#Filter out HRMGs only
head -1 Bakta_v18_annotations.tsv > hmrgs_3255_best_assemblies.txt  
# "Gene" is the 8th column
awk -F '\t' '$8 ~ /nik|cpx|cus|cut|cue|cop|chr|zra|cor|rcn|cbt|cbi|dpp|opp|mnt|sit|mts|znu|zup|zit|zur|teh|ydc|yeb|tsg|yff|gol|ars|fie|yii|mer|ncr|pco|sil|ter|cad|czc|dsb|znt|cds/' Bakta_v18_annotations.tsv >> hmrgs_3255_best_assemblies.txt ;

#Then parse in R with: parse_HMRGs.R


