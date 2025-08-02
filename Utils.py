
import logging


class Utils():

    def print_dictionary(dictionary_of_proteins):
        for key, value in dictionary_of_proteins.items():
            print(key, value.get_description(), value.get_seq_reg_name())

    def get_uniprot_informations_of_protein(protein):
        print(protein.get_accession_number())
        print(protein.get_protein_name())
        print(protein.get_gen())
        print(protein.get_organism_scitific())
        print(protein.get_organism_common())
        print(protein.get_mol_weight())
        print(protein.get_ensembl_gene_id())
        print(protein.get_description())
        print(protein.get_seq_reg_name())

    def save_to_xlsx(df,header=True,path="protein_gene_analysis.xlsx"):
        logging.info(f"Saved as protein_gene_analysis.xlsx")
        df.to_csv(path, index=False, header=header)

