#python 3


# read from a csv file
import csv
import numpy as np
import pandas as pd

KG_dir1 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/"
#Reactome
KG_file1 = "/Users/guangrong/Documents/GitHub_project/KGs/Reactome/reactome_ppis.csv"
#read csv file
KG_df1 = pd.read_csv(KG_file1, sep=',')

#SIGNOR
KG_file2 = "/Users/guangrong/Documents/GitHub_project/KGs/SIGNOR/SIGNOR_formated.tsv"
#read tsv file
KG_df2 = pd.read_csv(KG_file2, sep='\t')

#BioGRID
KG_file3 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/Biogrid_formated.csv"
KG_df3 = pd.read_csv(KG_file3, sep=',')

#HuRI
KG_file4 =KG_dir1+"H-I-05_formated.csv"
KG_df4 = pd.read_csv(KG_file4, sep=',')


KG_file5 =KG_dir1+"HI-II-14_formated.csv"
KG_df5 = pd.read_csv(KG_file5, sep=',')


KG_file6 =KG_dir1 + "HuRI_formated.csv"
KG_df6 = pd.read_csv(KG_file6, sep=',')

KG_file7 =KG_dir1 + "Yang-16_formated.csv"
KG_df7 = pd.read_csv(KG_file7, sep=',')


# TF_target

# Msig DB

#CellMarker
KG_file8 = KG_dir1 + "cellmarker.csv"
KG_df8 = pd.read_csv(KG_file8, sep = ',')

#unique predicates
unique_predicates = set(list(set(KG_df1.predicate.unique())) 
                + list(set(KG_df2.predicate.unique())) 
                + list(set(KG_df3.predicate.unique()))
                + list(set(KG_df4.predicate.unique()))
                + list(set(KG_df5.predicate.unique()))
                + list(set(KG_df6.predicate.unique()))
                + list(set(KG_df7.predicate.unique()))
		+ list(set(KG_df8.predicate.unique()))
                )

#unique subject and object nodes
unique_subject_object_nodes = set(list(set(KG_df1.subject_id.unique()))
                                +list(set(KG_df2.subject_id.unique()))
                                +list(set(KG_df3.subject_id.unique()))
                                +list(set(KG_df4.subject_id.unique()))
                                +list(set(KG_df5.subject_id.unique()))
                                +list(set(KG_df6.subject_id.unique()))
                                +list(set(KG_df7.subject_id.unique()))
                                +list(set(KG_df1.object_id.unique()))
                                +list(set(KG_df2.object_id.unique()))
                                +list(set(KG_df3.object_id.unique()))
                                +list(set(KG_df4.object_id.unique()))
                                +list(set(KG_df5.object_id.unique()))
                                +list(set(KG_df6.object_id.unique()))
                                +list(set(KG_df7.object_id.unique()))
				+list(set(KG_df8.subject_id.unique()))
				+list(set(KG_df8.object_id.unique()))
                                )


print("length of unique predicates:")
print(len(unique_predicates))
#print(unique_predicates)

print("length of unique subject and object nodes:")
print(len(unique_subject_object_nodes))
