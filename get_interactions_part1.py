import pandas
import pandas as pd
import numpy as np

#part1
###########################################################
dhfq = pd.read_excel('C:/Users/.../31761494/DHfq.xlsx')
dproq = pd.read_excel('C:/Users/.../31761494/DProQ.xlsx')
hfq = pd.read_excel('C:/Users/.../31761494/Hfq.xlsx')
hfq_lb = pd.read_excel('C:/Users/.../31761494/Hfq_LB.xlsx')
hfq_m63 = pd.read_excel('C:/Users/.../31761494/Hfq_M63.xlsx')
proq = pd.read_excel('C:/Users/.../31761494/ProQ.xlsx')
proq_lb = pd.read_excel('C:/Users/.../31761494/ProQ_LB.xlsx')
proq_m63 = pd.read_excel('C:/Users/.../31761494/ProQ_M63.xlsx')
fields = pd.read_excel('C:/Users/.../agnodice_empty_fields.xlsx', header=0)

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
dhfq = proper_format(dhfq)
dproq = proper_format(dproq)
hfq = proper_format(hfq)
hfq_lb = proper_format(hfq_lb)
hfq_m63 = proper_format(hfq_m63)
proq = proper_format(proq)
proq_lb = proper_format(proq_lb)
proq_m63 =proper_format(proq_m63)

# information about every condition
dhfq['microbe_condition'] = "Tagged ProQ in an Hfq deletion"
dhfq['rbp'] = 'ProQ'
dproq['microbe_condition'] = "Tagged Hfq in a ProQ deletion"
dproq['rbp'] = "Hfq"
hfq['microbe_condition'] = "Tagged Hfq"
hfq['rbp'] = "Hfq"
hfq_lb['microbe_condition'] = "Tagged Hfq. Cells were grown to OD600 ~1.0 in rich (LB) glucose media"
hfq_lb['rbp'] = "Hfq"
hfq_m63['microbe_condition'] = "Tagged Hfq. Cells were grown to OD600 ~1.0 in minimal (M63) glucose media"
hfq_m63['rbp'] = "Hfq"
proq['microbe_condition'] = "Tagged ProQ"
proq['rbp'] = "ProQ"
proq_lb['microbe_condition'] = "Tagged ProQ. Cells were grown to OD600 ~1.0 in rich (LB) glucose media"
proq_lb['rbp'] = "ProQ"
proq_m63['microbe_condition'] = "Tagged ProQ. Cells were grown to OD600 ~1.0 in minimal (M63) glucose media"
proq_m63['rbp'] = "ProQ"

#append results in our excel
fields = fields.append([dhfq, dproq, hfq, hfq_lb, hfq_m63, proq, proq_lb, proq_m63])

#common information in all conditions
fields['publication_PMID'] = 31761494
fields['publication_DOI'] = "10.1016/j.molcel.2019.10.022"
fields['publication_first_author'] = "Sahar Melamed"
fields['publication_corresponding_author_mail'] = "storzg@mail.nih.gov"
fields['publication_title'] = "RNA-RNA Interactomes of ProQ and Hfq Reveal Overlapping and Competing Roles"
fields['publication_journal'] = "Molecular Cell"
fields['publication_year'] = 2020
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


fields.to_excel('C:/Users/.../agnodice_fields.xlsx', index=False)
########################################################### end of part1 ###########################################################