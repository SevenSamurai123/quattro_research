import logging
import pandas as pd
import requests
import time
import sys


from Protein import Protein
from Fetch_websites import Fetch_websites
from Utils import Utils

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("protein_API_Tasks.log"),
        logging.StreamHandler()
    ]
)

dictionary_of_proteins = {}
list_of_all_uniprot_ids = []
not_valid_uniprot_ids = []

class Protein_informations():

    def __init__(self, uniprot_ids: list):
        self.uniprot_ids = uniprot_ids

        self.final_result = pd.DataFrame


    def get_uniprot(self, uniprot_ids: list):
        """Collect informations about a uniprot id.
        Returns a dataframe containing with six columns:
            Protein Accession,
            Protein Name,
            Gen,
            Organism (Scientific),Description,
            Organism (Common),
            Molecular Weight
        
        uniprot_ids -- list of uniprot ids
        """
        
        df = pd.DataFrame(columns=[
            "Protein Accession",
            "Protein Name",
            "Gen",
            "Organism (Scientific)",
            "Organism (Common)",
            "Molecular Weight"
            ])
        
        for id in uniprot_ids:
            if id in list_of_all_uniprot_ids:
                logging.info(f"UniProt ID: {id} is redundant")
                continue

            logging.info(f"Handle UniProt ID: {id}")

            data = Fetch_websites.fetch_data(case="uniprot",id=id)
            if data == None: continue
              
            gene_info = data.get("gene", [])
            try:
                gene_name = gene_info[0].get("name", {}).get("value")
            except:
                logging.warning(f"UniProt ID: {id}. No gene name found.")
                not_valid_uniprot_ids.append(id)
                continue
                
            protein_info = data.get("protein", {}).get("recommendedName", {})
            try:
                protein_name = protein_info.get("fullName", {}).get("value")
            except:
                logging.info(f"UniProt ID: {id}. No protein name found.")
                protein_name = "empty"
                continue

            organism_info = data.get("organism", {})
            try:
                scientific_name = organism_info["names"][0]["value"]
            except:
                logging.info(f"UniProt ID: {id}. No organism name scientific found.")
                continue
            
            try:
                 common_name = organism_info["names"][1]["value"]
            except:
                logging.info(f"UniProt ID: {id}. No organism name common found.")
                common_name = "empty"
                continue
           
            try:
                molecular_weight = data.get("sequence", {}).get("mass", "N/A")
            except:
                logging.info(f"Uniprot: No molecular weight uniprot id: {id}.")
                molecular_weight = "empty"
                continue
            
            new_protein = Protein(accession_number=id)
            new_protein.set_protein_name(protein_name)
            new_protein.set_gen(gene_name)
            new_protein.set_organism_scitific(scientific_name)
            new_protein.set_organism_common(common_name)
            new_protein.set_mol_weight(molecular_weight)
            
            dictionary_of_proteins[id] = new_protein

            new_entry = {
                "Protein Accession": id,
                "Protein Name": protein_name,
                "Gen": gene_name,
                "Organism (Scientific)": scientific_name,
                "Organism (Common)": common_name,
                "Molecular Weight": molecular_weight
                }

            df.loc[len(df)] = new_entry

            list_of_all_uniprot_ids.append(id)

        return df

    def get_ensembl_gene_ids(self, df_proteins):
        """Get ensembl gene id of a specific gene.
        Returns input dataframe with new column Ensembl Gene ID.

        df_proteins -- dataframe of function get_uniprot
        """
        ensembl_ids = []
        dataframe_ensembl_gene_id = {}

        for index, row in df_proteins.iterrows():
            species = row["Organism (Scientific)"].lower()
            gene = row["Gen"]
            id = row["Protein Accession"]
        
            data = Fetch_websites.fetch_data(case="ensembl",id=id,species=species,gene=gene)
            if data == None:
                continue

            try:
                gene_id = next((entry["id"] for entry in data if entry["type"] == "gene"))
                ensembl_ids.append(gene_id)
                current_protein = dictionary_of_proteins[id]
                current_protein.set_ensembl_gene_id(gene_id)
                dataframe_ensembl_gene_id[id]=gene_id
            except:
                logging.warning(f"Ensembl: No Ensembl Gene ID found for {species, gene}. UniProt ID: {id}.")
                not_valid_uniprot_ids.append(id)
                continue
        
        #logging.info(f"start dataframe merge")
        tmp_df = pd.DataFrame(dataframe_ensembl_gene_id.items(), columns=['Protein Accession','Ensembl Gene ID'])

        df_valid_proteins = pd.merge(df_proteins,tmp_df, on='Protein Accession', how='inner')
        #logging.info(f"end dataframe merge")


        return df_valid_proteins

    def get_gene_data(self, df_proteins):
        """Get additional gene data. Collect gene description and sequence region name.
        Returns dataframe with three columns: Ensembl Gene ID, Description, Seq Region Name.

        df_proteins -- result dataframe of function get_ensembl_gene_ids  
        """
        logging.info(f"start get_gene_data")

        df = pd.DataFrame(columns=[
            "Ensembl Gene ID",
            "Description",
            "Seq Region Name"
            ])

        ensembl_headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
        }

        ensembl_gene_ids = []
        descriptions = []
        seq_region_names = []

        BASE_URL_ENSEMBL = "https://rest.ensembl.org"
        LOOKUP_TEXT = "/lookup/id/"

        for index, row in df_proteins.iterrows():
            ensembl_gene_id = row["Ensembl Gene ID"]
            id = row["Protein Accession"]

            ensembl_gene_ids.append(ensembl_gene_id)
            
            url = f"{BASE_URL_ENSEMBL}{LOOKUP_TEXT}{ensembl_gene_id}?content-type=application/json"
            response = requests.get(url, headers=ensembl_headers)
            if response.status_code != 200:
                logging.warning(f"Ensembl: Ensembl Gene ID is not valid: {ensembl_gene_id}. UniProt ID: {id}: HTTP {response.status_code}")
                continue

            data = response.json()
            try:
                description = data.get("description", "")
                
            except:
                logging.info(f"Ensembl: No description for UniProt ID: {id}. Ensembl ID: {ensembl_gene_id}.")
                description = "empty"

            descriptions.append(description)
            
            try:
                seq_region_name = data.get("seq_region_name", "")
            except:
                logging.info(f"Ensembl: No sequence region name for UniProt ID: {id}. Ensembl ID: {ensembl_gene_id}.")
                seq_region_name = "empty"
            
            seq_region_names.append(seq_region_name)
            
            current_protein = dictionary_of_proteins[id]
            current_protein.set_description(description)
            current_protein.set_seq_reg_name(seq_region_name)

        df["Ensembl Gene ID"] = ensembl_gene_ids
        df["Description"] = descriptions
        df["Seq Region Name"] = seq_region_names

        return df

    def get_dictionary_of_proteins():
        return dictionary_of_proteins
    
    def not_valid_ids():
        logging.info(f"Not valid UniProt IDs:")
        return not_valid_uniprot_ids

    def start(self):
        logging.info(f"Start program")
        uniprot_dataframe = self.get_uniprot(self.uniprot_ids)

        #check if dataframe has at least one entry
        if uniprot_dataframe.size == 0:
            return None
        
        ensembl_dataframe = self.get_ensembl_gene_ids(uniprot_dataframe)
        #print(ensembl_dataframe)

        gene_data_dataframe = self.get_gene_data(ensembl_dataframe)
        #print(gene_data_dataframe)

        #merge dataframes
        self.final_result = pd.merge(ensembl_dataframe,gene_data_dataframe, on='Ensembl Gene ID', how='outer')

        Utils.save_to_xlsx(self.final_result)
        
        logging.info(f"End program")
        return self.final_result
 

if __name__ == "__main__":
    if len(sys.argv) > 2:
        print("Usage: python API_tasks.py <Q8N726,O00255,P69905,Q9Y261> or <pathToFile>")
        sys.exit(1)
    param1 = sys.argv[1]

    if ".txt" in param1:
        with open(param1, 'r') as file:
            uniprot_ids = [line.strip() for line in file]    
    else:
        uniprot_ids = param1.strip("[]").split(",")

    df_proteins = Protein_informations(uniprot_ids=uniprot_ids)
    df = df_proteins.start()
    print(df)
        