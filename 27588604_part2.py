import pandas as pd
import numpy as np
import math

###################################### part 2 ######################################

fields = pd.read_excel('C:/Users/vasil/OneDrive/Agnodice_final/agnodice_fields_v.xlsx')
ecocyc_gff = pd.read_excel('C:/Users/vasil/OneDrive/Agnodice_final/E_coli_MG1655_K12_NCBI.gtf/Ecocyc_ecoli_MG1655_K12_sorted.xlsx')
ncbi_gff = pd.read_excel('C:/Users/vasil/OneDrive/Agnodice_final/E_coli_MG1655_K12_NCBI.gtf/NCBI_ecoli_MG1655_K12.xlsx')


fields['rna_biocyc_id'] = "NA"
fields['srna_biocyc_id'] = "NA"
fields['rna_product'] = "NA"
fields['probability_score'] = "NA"


# remove strips
fields = fields.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
ecocyc_gff = ecocyc_gff.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
ncbi_gff = ncbi_gff.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

# turn df into dict
fields_dictionary = fields.to_dict('records')
ncbi_dictionary = ncbi_gff.to_dict('records')
ecocyc_dictionary = ecocyc_gff.to_dict('records')

# insert info from gffs (ncbi & ecocyc) about srna and rna except right and left genes & srna_seq. Also, check if they use a synonym name instead of primary
def fill_db(fields_dict, ncbi_dict, ecocyc_dict):
    for row_agno in fields_dict:

        for row_ncbi in ncbi_dict:
            # sRNA info from ncbi
            if row_agno['srna_name'] == row_ncbi['gene_name']:
                row_agno['srna_ncbi_id'] = row_ncbi['gene_id']
                row_agno['srna_synonym_names'] = row_ncbi['synonym_names']
                row_agno['srna_coordinates'] = (str(row_ncbi['start']) + "..." + str(row_ncbi['end']))
                row_agno['srna_length'] =(int(row_ncbi['end']) - int(row_ncbi['start']))
                row_agno['srna_strand'] = row_ncbi['strand']
                row_agno['srna_annotation_source'] = "NCBI"
            # RNA info from ncbi
            if row_agno['rna_name'] == row_ncbi['gene_name']:
                row_agno['rna_ncbi_id'] = row_ncbi['gene_id']
                row_agno['rna_synonym_names'] = row_ncbi['synonym_names']
                row_agno['rna_product'] = row_ncbi['product']
                row_agno['rna_protein_id'] = row_ncbi['protein_id']

        for row_ecocyc in ecocyc_dict:
            # sRNA info from ecocyc
            if row_agno['srna_name'] == row_ecocyc['gene_name']:
                row_agno['srna_biocyc_id'] = row_ecocyc['gene_id']
                if row_agno['srna_annotation_source'] != "NCBI":
                    row_agno['srna_synonym_names'] = row_ecocyc['synonym_names']
                    row_agno['srna_annotation_source'] = "EcoCyc"
                    row_agno['srna_coordinates'] = (str(row_ecocyc['start']) + "..." + str(row_ecocyc['end']))
                    row_agno['srna_length'] = (int(row_ecocyc['end']) - int(row_ecocyc['start']))
                    row_agno['srna_strand'] = row_ecocyc['strand']
            # RNA info from ecocyc
            if row_agno['rna_name'] == row_ecocyc['gene_name']:
                row_agno['rna_biocyc_id'] = row_ecocyc['gene_id']
                if (pd.isnull(row_agno['rna_ncbi_id'])):
                    row_agno['rna_synonym_names'] = row_ecocyc['synonym_names']
                    row_agno['rna_product'] = row_ecocyc['product']

        # check if they use for rnas a synonym name instead of the primary
        if (pd.isnull(row_agno['srna_ncbi_id'])) and (row_agno['srna_biocyc_id'] == "NA"):
            synonyms = []
            for row_ecocyc in ecocyc_dictionary:
                synonyms = str(row_ecocyc['synonym_names']).split(",")

                if row_agno['srna_name'] in str(synonyms):
                    row_agno['srna_name'] = row_ecocyc['gene_name']

        if (pd.isnull(row_agno['srna_ncbi_id'])) and (row_agno['srna_biocyc_id'] == "NA"):
            synonyms = []
            for row_ncbi in ncbi_dictionary:
                synonyms = str(row_ncbi['synonym_names']).split(",")
                if row_agno['srna_name'] in synonyms:
                    row_agno['srna_name'] = row_ncbi['gene_name']

        if (pd.isnull(row_agno['rna_ncbi_id'])) and (row_agno['rna_biocyc_id'] == "NA"):
            synonyms = []
            for row_ecocyc_1 in ecocyc_dictionary:
                synonyms = str(row_ecocyc_1['synonym_names']).split(",")
                if str(row_agno['rna_name']) in str(synonyms):
                    row_agno['rna_name'] = row_ecocyc_1['gene_name']

        if (pd.isnull(row_agno['rna_ncbi_id'])) and (row_agno['rna_biocyc_id'] == "NA"):
            synonyms = []
            for row_ncbi_1 in ncbi_dictionary:
                synonyms = str(row_ncbi_1['synonym_names']).split(",")

                if str(row_agno['rna_name']) in str(synonyms):
                    row_agno['rna_name'] = row_ncbi_1['gene_name']

#run again the fill_db function for the rnas that the used synonym name and not the primary
def fill_db_second_round(fields_dict, ncbi_dict, ecocyc_dict):
    for row_agno in fields_dict:
        if (pd.isnull(row_agno['srna_ncbi_id']) and row_agno['srna_biocyc_id'] == "NA") or ((pd.isnull(row_agno['rna_ncbi_id'])) and (row_agno['rna_biocyc_id'])):
            for row_ncbi in ncbi_dict:
                # sRNA info from ncbi
                if row_agno['srna_name'] == row_ncbi['gene_name']:
                    row_agno['srna_ncbi_id'] = row_ncbi['gene_id']
                    row_agno['srna_synonym_names'] = row_ncbi['synonym_names']
                    row_agno['srna_coordinates'] = (str(row_ncbi['start']) + "..." + str(row_ncbi['end']))
                    row_agno['srna_length'] =(int(row_ncbi['end']) - int(row_ncbi['start']))
                    row_agno['srna_strand'] = row_ncbi['strand']
                    row_agno['srna_annotation_source'] = "NCBI"
                # RNA info from ncbi
                if row_agno['rna_name'] == row_ncbi['gene_name']:
                    row_agno['rna_ncbi_id'] = row_ncbi['gene_id']
                    row_agno['rna_synonym_names'] = row_ncbi['synonym_names']
                    row_agno['rna_product'] = row_ncbi['product']
                    row_agno['rna_protein_id'] = row_ncbi['protein_id']

            for row_ecocyc in ecocyc_dict:
                # sRNA info from ecocyc
                if row_agno['srna_name'] == row_ecocyc['gene_name']:
                    row_agno['srna_biocyc_id'] = row_ecocyc['gene_id']
                    if row_agno['srna_annotation_source'] != "NCBI":
                        row_agno['srna_synonym_names'] = row_ecocyc['synonym_names']
                        row_agno['srna_annotation_source'] = "EcoCyc"
                        row_agno['srna_coordinates'] = (str(row_ecocyc['start']) + "..." + str(row_ecocyc['end']))
                        row_agno['srna_length'] = (int(row_ecocyc['end']) - int(row_ecocyc['start']))
                        row_agno['srna_strand'] = row_ecocyc['strand']
                # RNA info from ecocyc
                if row_agno['rna_name'] == row_ecocyc['gene_name']:
                    row_agno['rna_biocyc_id'] = row_ecocyc['gene_id']
                    if (pd.isnull(row_agno['rna_ncbi_id'])):
                        row_agno['rna_synonym_names'] = row_ecocyc['synonym_names']
                        row_agno['rna_product'] = row_ecocyc['product']

# call functions
fill_db(fields_dictionary, ncbi_dictionary, ecocyc_dictionary)
fill_db_second_round(fields_dictionary, ncbi_dictionary, ecocyc_dictionary)


# turn dict into df
fields = pd.DataFrame.from_records(fields_dictionary)
fields.fillna("NA", inplace=True)


fields.to_excel("C:/Users/vasil/OneDrive/Agnodice_final/Agnodice_fields_study_27588604.xlsx", index=False)


############################# end of part2 #############################