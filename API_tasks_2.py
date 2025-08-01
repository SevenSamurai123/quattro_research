import logging
import pandas as pd
import requests
import time
import sys

from API_tasks import Protein_informations
from Protein import Protein

dictionary_of_proteins = {}
list_of_all_uniprot_ids = []

class Protein_informations():

    def __init__(self, uniprot_ids: list):
        self.uniprot_ids = uniprot_ids

        self.final_result = pd.DataFrame

    # Logging configuration
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler("protein_API_Tasks.log"),
            logging.StreamHandler()
        ]
    )

    def fetch_data(self, case,id=None, species=None, gene=None):
        """ Fetch the data from uniprot and ensebl.

        id -- uniprot id, only for uniprot
        species -- species name organism (scientific), only for ensembl
        gene -- gene name, only for ensembl
        """

        if case == "uniprot":
            base_url = 'https://www.ebi.ac.uk/proteins/api/proteins/'
            url = base_url + id
        elif case == "ensembl":
            url = f'https://rest.ensembl.org/xrefs/symbol/{species}/{gene}'
        else:
            logging.error(f"Missing parameter(s) for data fetching.")

        headers = {
        "Accept": "application/json"
        }
        
        response = requests.get(url, headers=headers)
        match response.status_code:
            case 429:
                logging.warning(f"Reached limit of requests per second. Wait 3 seconds.")
                time.sleep(3)
                response = requests.get(url, headers=headers)
                data = response.json()
                return data
            case 400:
                logging.warning(f"Error by calling uniprot id: {id} or gene {gene}: HTTP {response.status_code}")
                return None
            case 404:
                logging.warning(f"Uniprot id: {id} or gene {gene} is not valid: HTTP {response.status_code}")
                return None
            case 200:
                data = response.json()
                return data
            case _:
                logging.warning(f"Error by calling uniprot id: {id} or gene {gene}: HTTP {response.status_code}")
                return None

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

            data = self.fetch_data(case="uniprot",id=id)
            if data == None: continue

            # accession_number = data.get("accession")
            # if accession_number == None:
            #     logging.warning(f"Uniprot: No accession number for uniprot id: {id}.")
            #     continue
                
            gene_info = data.get("gene", [])
            try:
                gene_name = gene_info[0].get("name", {}).get("value")
            except:
                logging.warning(f"Uniprot: No gene name for uniprot id: {id}")
                continue
                
            protein_info = data.get("protein", {}).get("recommendedName", {})
            try:
                protein_name = protein_info.get("fullName", {}).get("value")
            except:
                logging.warning(f"Uniprot: No protein name for uniprot id: {id}.")
                continue

            organism_info = data.get("organism", {})
            try:
                scientific_name = organism_info["names"][0]["value"]
            except:
                logging.warning(f"Uniprot: No organism name scientific for uniprot id: {id}.")
                continue
            
            try:
                 common_name = organism_info["names"][1]["value"]
            except:
                logging.warning(f"Uniprot: No organism name common for uniprot id: {id}.")
                common_name = "empty"
                continue
           

            molecular_weight = data.get("sequence", {}).get("mass", "N/A")

            new_protein = Protein(accession_number=id)
            new_protein.set_protein_name(protein_name)
            new_protein.set_gen(gene_name)
            new_protein.set_organism_scitific(scientific_name)
            new_protein.set_organism_common(common_name)
            new_protein.set_mol_weight(molecular_weight)
            
            dictionary_of_proteins[id] = new_protein

            #print(id)

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
        print(df.shape)

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
        
            data = self.fetch_data(case="ensembl",id=id,species=species,gene=gene)
            if data == None:
                continue

            try:
                gene_id = next((entry["id"] for entry in data if entry["type"] == "gene"))
                ensembl_ids.append(gene_id)
                current_protein = dictionary_of_proteins[id]
                current_protein.set_ensembl_gene_id(gene_id)
                dataframe_ensembl_gene_id[id]=gene_id
            except:
                logging.warning(f"Ensembl: No ensembl gene id found for species: {species}, gene: {gene}")
                continue
        
        print(len(dataframe_ensembl_gene_id))

        logging.info(f"start dataframe merge")
        tmp_df = pd.DataFrame(dataframe_ensembl_gene_id.items(), columns=['Protein Accession','Ensembl Gene ID'])
        print(tmp_df.columns)
        print(tmp_df.shape)
        df_valid_proteins = pd.merge(df_proteins,tmp_df, on='Protein Accession', how='inner')
        logging.info(f"end dataframe merge")

        print(df_valid_proteins.columns)
        print(df_valid_proteins.shape)
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

        base_url = "https://rest.ensembl.org"
        lookup_ext = "/lookup/id/"

        for index, row in df_proteins.iterrows():
            ensembl_gene_id = row["Ensembl Gene ID"]
            accession_number = row["Protein Accession"]

            ensembl_gene_ids.append(ensembl_gene_id)
            
            url = f"{base_url}{lookup_ext}{ensembl_gene_id}?content-type=application/json"
            response = requests.get(url, headers=ensembl_headers)
            if response.status_code != 200:
                logging.warning(f"Ensembl: Error by Ensembl Gene ID: {ensembl_gene_id}: HTTP {response.status_code}")
                continue

            info = response.json()
            try:
                description = info.get("description", "")
                
            except:
                logging.warning(f"Ensmbl: No description uniprot id: {accession_number}.")
                description = "empty"

            descriptions.append(description)
            
            try:
                seq_region_name = info.get("seq_region_name", "")
            except:
                logging.warning(f"Ensmbl: No Sequence Region Name uniprot id: {accession_number}.")
                seq_region_name = "empty"
            
            seq_region_names.append(seq_region_name)
            
            current_protein = dictionary_of_proteins[accession_number]
            current_protein.set_description(description)
            current_protein.set_seq_reg_name(seq_region_name)

        df["Ensembl Gene ID"] = ensembl_gene_ids
        df["Description"] = descriptions
        df["Seq Region Name"] = seq_region_names

        return df

    def save_to_xlsx(self,df,header=True,path="protein_gene_analysis.xlsx"):
        logging.info(f"Saved as protein_gene_analysis.xlsx")
        df.to_csv(path, index=False, header=header)

    def get_dictionary_of_proteins():
        return dictionary_of_proteins

    def start(self):
        logging.info(f"Start program")
        uniprot_dataframe = self.get_uniprot(self.uniprot_ids)
        #print(uniprot_dataframe)

        ensembl_dataframe = self.get_ensembl_gene_ids(uniprot_dataframe)
        #print(ensembl_dataframe)

        gene_data_dataframe = self.get_gene_data(ensembl_dataframe)
        #print(gene_data_dataframe)

        #merge dataframes
        self.final_result = pd.merge(ensembl_dataframe,gene_data_dataframe, on='Ensembl Gene ID', how='outer')

        self.save_to_xlsx(self.final_result)
        
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

    logging.info(f"End program")
        