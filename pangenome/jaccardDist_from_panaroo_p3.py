#!/usr/bin/env python3
#################################################
# script by Kelly Wyres, 30th May 2017 - modified from Kat Holt's jaccard.py.
# translated to python 3 by Marit Hetland, 12th April 2024 + made compatible with panaroo output 
# read in: presence/absence matrix (csv) from panaroo
# output: jaccard distances for columns
# optional add tag to each row - this is useful if you will be concatenating the outputs from multiple runs 
# for easy manipulation and comparison e.g. in R
# usage python jaccard_from_panaroo -i <panaroo matrix> -o <outputfilename> -t <tag>
#################################################

import csv
from argparse import ArgumentParser
from collections import defaultdict

def parse_args():
    "Parse the input arguments, use '-h' for help."
    parser = ArgumentParser(description='Parse Panaroo output and output pairwise Jaccard distances')
    parser.add_argument("-i", "--input", type=str, required=True, help='Gene presence/absence file from Panaroo (CSV)')
    parser.add_argument("-o", "--output", type=str, required=True, help='Output file name')
    parser.add_argument("-t", "--tag", type=str, required=False, help='Tag to add on each row in output for easy manipulation')
    return parser.parse_args()

def main():
    args = parse_args()

    # Initialise a dictionary to store the presence of genes by strain.
    content_matrix = defaultdict(list)

    # Read the input CSV file.
    with open(args.input, newline='') as file:
        reader = csv.reader(file)  
        header = next(reader)  
        for row in reader:  
            gene = row[0]  
            for i in range(3, len(header)):
                strain = header[i]  
                if row[i] != "":  #
                    content_matrix[strain].append(gene)  


    # Prepare to write the Jaccard distances to output file
    with open(args.output, 'w', newline='') as file:
        outputFile = csv.writer(file, delimiter='\t')
        strain_list = list(content_matrix.keys())

    # Calculate pairwise Jaccard distances
        for i in range(len(strain_list)-1):
            for j in range(i+1, len(strain_list)):
                strain_i = set(content_matrix[strain_list[i]])  # Convert gene list to set for strain i
                strain_j = set(content_matrix[strain_list[j]])  # Convert gene list to set for strain j
                intersection = len(strain_i & strain_j)  # Calculate intersection of both sets
                union = len(strain_i | strain_j)  # Calculate union of both sets
                jaccard_distance = 1 - float(intersection) / float(union)  # Calculate Jaccard distance

                #Print the number of genes in the intersection and strain identifiers.
                #This number represents how many genes are shared between two given strains. It's calculated by taking the intersection of the sets of genes present in each strain.
                print(f'{intersection}\t{strain_list[i]}\t{strain_list[j]}')  
                #Print the number of genes in the union
                #This is the total number of unique genes present in both strains combined. It's computed by taking the union of the gene sets from both strains.
                print(f'{union}')
                #Jaccard distance is calculated as: 1-(intersection/union).
                #Intersection / Union: The ratio of the number of genes common to both strains to the total number of unique genes across both strains.
                #Interpretation: A Jaccard distance of 0 means the strains are identical (in terms of gene presence/absence), while a value closer to 1 indicates higher dissimilarity.

                # Write the output with or without a tag as per the user's input
                if args.tag:
                    outputFile.writerow([strain_list[i], strain_list[j], str(jaccard_distance), args.tag])
                else:
                    outputFile.writerow([strain_list[i], strain_list[j], str(jaccard_distance)])

if __name__ == "__main__":
    main()