import pandas as pd
import numpy as np
import math


###################################### part 3 1/2 ######################################
# # on that part we get the right and left genes for the srnas and also we create the bed file in order to run bedtools getfasta
#
# fields = pd.read_excel('C:/Users/vasil/OneDrive/Agnodice_final/Agnodice_fields_study_27588604.xlsx', keep_default_na=False)
# ecocyc_gff = pd.read_excel('C:/Users/vasil/OneDrive/Agnodice_final/E_coli_MG1655_K12_NCBI.gtf/Ecocyc_ecoli_MG1655_K12_sorted.xlsx')
#
# # sort by coordinates
# ecocyc_gff = ecocyc_gff.sort_values(by=['start']).reset_index(drop=True)
#
# # create 2 columns for the coordinates
# fields[['start_srna', 'end_srna']] = fields['srna_coordinates'].str.split("\.\.\.", expand=True)
#
# # split the gff into 2 df depends on strand
# ecocyc_gff_plus_strand = ecocyc_gff.loc[ecocyc_gff['strand'] == "+"]
# ecocyc_gff_plus_strand = ecocyc_gff_plus_strand.reset_index(drop=True)
# ecocyc_gff_plus_strand = ecocyc_gff_plus_strand.sort_values(by=['start']).reset_index(drop=True)
#
# ecocyc_gff_minus_strand = ecocyc_gff.loc[ecocyc_gff['strand'] == "-"]
# ecocyc_gff_minus_strand = ecocyc_gff_minus_strand.reset_index(drop=True)
# ecocyc_gff_minus_strand = ecocyc_gff_minus_strand.sort_values(by=['start']).reset_index(drop=True)
#
# # turn df to dict
# fields_dict = fields.to_dict('records')
#
# # find srnas left and right genes
#
# for row_fields in fields_dict:
#     if row_fields['srna_strand'] == "+":
#         # find srna index in ecocyc_gff_plus_s
#         index = ecocyc_gff_plus_strand.index[ecocyc_gff_plus_strand['gene_name'] == row_fields['srna_name']][0]
#         if index == 0:
#             row_fields['srna_left_gene'] = "NA"
#             row_fields['srna_right_gene'] = ecocyc_gff_plus_strand['gene_name'].iloc[index + 1]
#
#         elif index == len(ecocyc_gff_plus_strand):
#             row_fields['srna_left_gene'] = ecocyc_gff_plus_strand['gene_name'].iloc[index - 1]
#             row_fields['srna_right_gene'] = "NA"
#
#         else:
#             row_fields['srna_left_gene'] = ecocyc_gff_plus_strand['gene_name'].iloc[index-1]
#             row_fields['srna_right_gene'] = ecocyc_gff_plus_strand['gene_name'].iloc[index+1]
#
#     if row_fields['srna_strand'] == "-":
#         index = ecocyc_gff_minus_strand.index[ecocyc_gff_minus_strand['gene_name'] == row_fields['srna_name']][0]
#
#         if index == 0:
#             row_fields['srna_left_gene'] = ecocyc_gff_minus_strand['gene_name'].iloc[index + 1]
#             row_fields['srna_right_gene'] = "NA"
#
#         elif index == len(ecocyc_gff_minus_strand):
#             row_fields['srna_left_gene'] = "NA"
#             row_fields['srna_right_gene'] = ecocyc_gff_minus_strand['gene_name'].iloc[index - 1]
#
#         else:
#             row_fields['srna_left_gene'] = ecocyc_gff_minus_strand['gene_name'].iloc[index+1]
#             row_fields['srna_right_gene'] = ecocyc_gff_minus_strand['gene_name'].iloc[index-1]
#
#
# # turn dict to df
# fields = pd.DataFrame.from_records(fields_dict)
# fields.to_excel('C:/Users/vasil/OneDrive/Agnodice_final/RIL-Seq/27588604/Agnodice_27588604_without_srna_seq.xlsx', index=False)
#
# # create the bed file to run bedtools getfasta
# bed_file = pd.DataFrame
# bed_file = fields[["start_srna", "end_srna", "srna_name", "srna_strand"]].copy()
# bed_file['chr'] = "chr"
# bed_file['score'] = 0
# bed_file = bed_file.iloc[:, [4,0,1,2,5,3]]
#
# bed_file.to_csv('C:/Users/vasil/OneDrive/Agnodice_final/RIL-Seq/27588604/bed_27588604.txt', sep="\t", index=False, header= False)

###################################### end of part 3 1/2 ######################################

###################################### part 3 2/2 ######################################

#get the seq for srnas from bedtools getfasta output and insert them on agnodice field srna-seq

fields = pd.read_excel('C:/Users/vasil/OneDrive/Agnodice_final/RIL-Seq/27588604/Agnodice_27588604_without_srna_seq.xlsx', keep_default_na=False)
fields_dict = fields.to_dict('records')

#bedtools getfasta output: proper format
seq_srna = pd.read_csv('C:/Users/vasil/OneDrive/Agnodice_final/RIL-Seq/27588604/srna_seq_27588604.txt', header=None)
seq_srna.rename(columns={0:"gene_name"}, inplace=True)
seq_srna['seq'] = ""
seq_srna = pd.DataFrame({'gene_name':seq_srna['gene_name'].iloc[::2].values, 'seq':seq_srna['gene_name'].iloc[1::2].values})
seq_srna["gene_name"] = seq_srna["gene_name"].replace(r'>', r'', regex = True)
seq_srna["gene_name"] = seq_srna["gene_name"].replace(r'\(\+\)', r'', regex = True)
seq_srna["gene_name"] = seq_srna["gene_name"].replace(r'\(\-\)', r'', regex = True)
seq_srna = seq_srna.drop_duplicates()
seq_srna_dict = seq_srna.to_dict('records')

for row_agno in fields_dict:
    for row_seq in seq_srna_dict:
        if row_agno['srna_name'] == row_seq['gene_name']:
            row_agno['srna_sequence'] = row_seq['seq']

fields = pd.DataFrame.from_records(fields_dict)
fields = fields.drop(['Genomic annotation of RNA1', 'start_srna', 'end_srna'], axis=1)
fields.to_excel('C:/Users/vasil/OneDrive/Agnodice_final/RIL-Seq/27588604/Agnodice_27588604_final.xlsx', index=False)

###################################### end of part 3 2/2 ######################################