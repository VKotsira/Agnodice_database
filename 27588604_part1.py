import pandas
import pandas as pd
import numpy as np

#part1
###########################################################
log = pd.read_excel('C:/Users/vasil/OneDrive/Agnodice_final/RIL-Seq/27588604/Log_phase.xlsx')
iron_limitation = pd.read_excel('C:/Users/vasil/OneDrive/Agnodice_final/RIL-Seq/27588604/Iron_limitation_phase.xlsx')
stationary = pd.read_excel('C:/Users/vasil/OneDrive/Agnodice_final/RIL-Seq/27588604/Stationary_phase.xlsx')

fields = pd.read_excel('C:/Users/vasil/OneDrive/Agnodice_final/agnodice_empty_fields.xlsx', header=0)

def proper_format(ril_df):

#delete non-sRNA interactions
    ril_df = ril_df[(ril_df['Genomic annotation of RNA1'] == 'sRNA') | (ril_df['Genomic annotation of RNA2'] == 'sRNA')].reset_index(drop=True)

# keep interactions with p-value < 0.05
    ril_df = ril_df[ril_df['Fisher\'s exact test p-value'] < 0.05].reset_index(drop=True)

# swap RNA1 as sRNA and RNA2 as non-sRNA
    for i in range(len(ril_df)):
        if (ril_df.at[i,'Genomic annotation of RNA1'] == 'sRNA') and (ril_df.at[i,'Genomic annotation of RNA2'] == 'sRNA'):
            pass

        elif ril_df.at[i,'Genomic annotation of RNA2'] == 'sRNA':
            ril_df.at[i, 'RNA1 name'], ril_df.at[i, 'RNA2 name'] = ril_df.at[i, 'RNA2 name'], ril_df.at[i, 'RNA1 name']
            ril_df.at[i, 'Genomic annotation of RNA1'], ril_df.at[i,'Genomic annotation of RNA2'] = ril_df.at[i,'Genomic annotation of RNA2'], ril_df.at[i,'Genomic annotation of RNA1']
            ril_df.at[i, 'RNA1 description'], ril_df.at[i,'RNA2 description'] = ril_df.at[i,'RNA2 description'], ril_df.at[i,'RNA1 description']
            ril_df.at[i, 'RNA1 EcoCyc ID'], ril_df.at[i,'RNA2 EcoCyc ID'] = ril_df.at[i,'RNA2 EcoCyc ID'], ril_df.at[i,'RNA1 EcoCyc ID']

    ril_df["antisense"] = 0
    ril_df["ESTUTR"] = 0
    ril_df["rna_biotype"] = ""

    for k in range(len(ril_df)):
        if ril_df.at[k, "Genomic annotation of RNA2"] == "sRNA":
            ril_df.at[k, "rna_biotype"] = "sRNA"
            ril_df.at[k, "Genomic annotation of RNA2"] = "sRNA"

        elif ril_df.at[k, "Genomic annotation of RNA2"] == "CDS":
            ril_df.at[k, "rna_biotype"] = "mRNA"

        elif ril_df.at[k, "Genomic annotation of RNA2"] == "3UTR":
            ril_df.at[k, "rna_biotype"] = "mRNA"

        elif ril_df.at[k, "Genomic annotation of RNA2"] == "5UTR":
            ril_df.at[k, "rna_biotype"] = "mRNA"

        elif ril_df.at[k, "Genomic annotation of RNA2"] == "tRNA":
            ril_df.at[k, "rna_biotype"] = "tRNA"
            ril_df.at[k, "Genomic annotation of RNA2"] = "tRNA"

        elif ril_df.at[k, "Genomic annotation of RNA2"] == "AS":
            ril_df.at[k, "rna_biotype"] = "UGR"
            ril_df.at[k, "antisense"] = 1

        if ("EST5UTR" in ril_df.at[k, "RNA2 name"]) or ("EST3UTR" in ril_df.at[k, "RNA2 name"]):
            ril_df.at[k, "ESTUTR"] = 1

    ril_df["RNA1 name"] = ril_df["RNA1 name"].replace(r'[.]\w*', r'', regex=True)
    ril_df["RNA2 name"] = ril_df["RNA2 name"].replace(r'([\w-]*)([.])([\w-]*)([.])([\w-]*)', r'\1_\3\4\5', regex=True)
    ril_df["RNA2 name"] = ril_df["RNA2 name"].replace(r'[.]\w*', r'', regex=True)
    ril_df["RNA1 EcoCyc ID"] = ril_df["RNA1 EcoCyc ID"].replace(r'[.]\w*', r'', regex=True)
    ril_df["RNA2 EcoCyc ID"] = ril_df["RNA2 EcoCyc ID"].replace('[.]EST\w*[.]\w*', r'', regex=True)
    ril_df["RNA2 EcoCyc ID"] = ril_df["RNA2 EcoCyc ID"].replace(r'(\w*)([.])(\w*)([.])(\w*)', r'\1_\3\4\5', regex=True)
    ril_df["RNA2 EcoCyc ID"] = ril_df["RNA2 EcoCyc ID"].replace(r'[.]\w*', r'', regex=True)

    # remove duplicate interacations
    ril_df = ril_df.drop_duplicates(["RNA1 name", "RNA2 name"]).reset_index(drop=True)

    # change column names
    ril_df = ril_df.rename(columns={'RNA1 name': 'srna_name', 'RNA2 name' : 'rna_name', 'Fisher\'s exact test p-value' : 'p_value_fishers_test',
                                    'Genomic annotation of RNA2' : 'rna_location', 'RNA2 description' : 'rna_product',
                                    'RNA1 EcoCyc ID' : 'srna_biocyc_id', 'RNA2 EcoCyc ID' : 'rna_biocyc_id', 'ESTUTR' : 'rna_ESTUTR'})

    # remove IGR, IGT and IGT_AS interactions
    ril_df = ril_df[(ril_df['rna_location'] == 'sRNA') | (ril_df['rna_location'] == '5UTR') | (ril_df['rna_location'] == '3UTR') | (ril_df['rna_location'] == 'CDS') | (ril_df['rna_location'] == 'tRNA') | (ril_df['rna_location'] == 'AS')].reset_index(drop=True)

    return(ril_df.drop('RNA1 description', axis=1))

# call the function
log = proper_format(log)
iron_limitation = proper_format(iron_limitation)
stationary = proper_format(stationary)

# information about every condition
log['microbe_condition'] = "Log phase"
log['rbp'] = 'Hfq'
iron_limitation['microbe_condition'] = "Iron limitation phase"
iron_limitation['rbp'] = "Hfq"
stationary['microbe_condition'] = "Stationary phase"
stationary['rbp'] = "Hfq"


#append results in our excel
fields = fields.append([log, iron_limitation, stationary])

#common information in all conditions
fields['publication_PMID'] = 27588604
fields['publication_DOI'] = "10.1016/j.molcel.2016.07.026"
fields['publication_first_author'] = "Sahar Melamed"
fields['publication_corresponding_author_mail'] = "hanahm@ekmd.huji.ac.il"
fields['publication_title'] = "Global Mapping of Small RNA-Target Interactions in Bacteria"
fields['publication_journal'] = "Molecular Cell"
fields['publication_year'] = 2016
fields['microbe_strain_name'] = "Escherichia coli K-12"
fields['microbe_strain_taxid'] = 83333
fields['microbe_species_name'] = "Escherichia coli"
fields['microbe_species_taxid'] = 562
fields['microbe_genus_name'] = "Escherichia"
fields['microbe_genus_taxid'] = 561
fields['microbe_phylum_name'] = "Proteobacteria"
fields['microbe_phylum_taxid'] = 1224
fields['microbe_genome_ncbi_id'] = "NC_000913.3"
fields['microbe_gram_type'] = "Gram-negative"
fields['experimental_method_name'] = "RIL-Seq"
fields['experimental_method_type'] = "High-Throughput"
fields['experimental_method_group'] = "Indirect"
fields['type_of_regulation'] = "Unknown"
fields['binding_site'] = "Unknown"
fields['interaction_energy'] = "Unknown"
fields['probability_score'] = "NA"
fields['comments'] = "NA"


fields.to_excel('C:/Users/vasil/OneDrive/Agnodice_final/agnodice_fields_v.xlsx', index=False)
########################################################### end of part1 ###########################################################