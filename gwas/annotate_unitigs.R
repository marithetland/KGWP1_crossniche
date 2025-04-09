#First, run BLASTn:
#for t in $(ls *.fasta | grep -v "^unitig") ; do for f in $(ls unitig*fasta) ; do echo -ne 
"$t\t$f\t" ; blastn -query $f -subject $t -outfmt 6 ; echo "" ; done ; done >> 
unitig_blastn_results.txt

#Read BLASTN results
unitig_results <- read.csv("~/KGWP1/GWAS/postprocessing/unitig_blastn_results.txt", sep="\t")

#Create a lookup data frame with the unitigs (queries) and their lengths (to match blastn results 
100%)
lookup_df <- data.frame(
  query = c("unitig_1", "unitig_10", "unitig_11", "unitig_12", "unitig_13", "unitig_14", 
"unitig_15", "unitig_16", 
            "unitig_2", "unitig_3", "unitig_4", "unitig_5", "unitig_6", "unitig_7", "unitig_8", 
"unitig_9"),
  length = c(53, 36, 32, 34, 37, 32, 45, 42, 36, 42, 31, 45, 42, 42, 41, 32)
)

#Subset and check that there are hits to all contigs with these strict settings
unitig_results %>%
  filter(!is.na(pident)) %>%  #Remove nonhits
  filter(!str_detect(genome, "^unitig")) %>% #Remove double hits to avoid confusion
  inner_join(lookup_df, by = c("query", "length")) %>% #Filter by query length (100%)
  filter(pident == 100) %>% #Filter by % ID (100%)
  group_by(query,pident) %>% count()

#Write to file to be processed
unitig_results_reduced <- unitig_results %>%
  filter(!is.na(pident)) %>%  #Remove nonhits
  filter(!str_detect(genome, "^unitig")) %>% #Remove double hits to avoid confusion
  inner_join(lookup_df, by = c("query", "length")) %>% #Filter by query length (100%)
  filter(pident == 100) 

write.csv(unitig_results_reduced,"~/KGWP1/GWAS/postprocessing/unitig_blastn_results_reduced.csv")

#In terminal:
#cd ~/KGWP1/actual_data_for_paper/GWAS/postprocessing
#cat unitig_blastn_results_reduced.csv | sed 's/"//g' | tr ',' '\t' | cut -f2,4,5,12,13 | sed 
's/.fasta/.gbff/g' | grep -v "^genome" >> input_file_for_gbk.txt
#d4: /media/markus/gimli/marit/GWAS/post_processing/refs_for_unitigs/ref/bakta
#python get_gene_by_pos.py input_file_for_gbk.txt outfile_for_gbk.csv

#Load final results
unitig_hits <- read.csv("~/KGWP1/GWAS/pyseer_june_2024/outfile_for_gbk.csv")

#See which genes/products the hits are to
unitig_hits %>% 
  group_by(Unitig,Gene) %>% 
  count() %>% 
  pivot_wider(names_from="Gene",values_from="n",values_fill=0) %>% 
  rowwise() %>%
  mutate(num_hits = sum(c_across(c(`No gene found`:iucC)) > 0)) %>%
  arrange(num_hits)

