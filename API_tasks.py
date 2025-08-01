import logging
import pandas as pd
import requests
import time
import sys

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

    def fetch_data(self, id=None, species=None, gene=None):
        """ Fetch the data from uniprot and ensebl.

        id -- uniprot id, only for uniprot
        species -- species name organism (scientific), only for ensembl
        gene -- gene name, only for ensembl
        """

        if id != None:
            base_url = 'https://www.ebi.ac.uk/proteins/api/proteins/'
            url = base_url + id
        elif species != None and gene != None:
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
            case 200:
                data = response.json()
                return data
            case _:
                logging.warning(f"Error by calling uniprot id: {id} or gene {gene}: HTTP {response.status_code}")

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
            logging.info(f"Handle UniProt ID: {id}")
            
            data = self.fetch_data(id=id)
            if data == None: continue

            accession_number = data.get("accession")
            if accession_number == None:
                logging.warning(f"Uniprot: No accession number for uniprot id: {id}.")
                continue
                
            gene_info = data.get("gene", [])
            if gene_info:
                gene_name = gene_info[0].get("name", {}).get("value")
            else:
                logging.warning(f"Uniprot: No gene name for uniprot id: {id}")
                continue
                
            protein_info = data.get("protein", {}).get("recommendedName", {})
            if protein_info:
                protein_name = protein_info.get("fullName", {}).get("value")
            else:
                logging.warning(f"Uniprot: No protein name for uniprot id: {id}.")
                continue

            organism_info = data.get("organism", {})
            if organism_info:
                scientific_name = organism_info["names"][0]["value"]
                common_name = organism_info["names"][1]["value"]
            else:
                logging.warning(f"Uniprot: No organism name for uniprot id: {id}.")
                continue

            molecular_weight = data.get("sequence", {}).get("mass", "N/A")

            new_entry = {
                "Protein Accession": accession_number,
                "Protein Name": protein_name,
                "Gen": gene_name,
                "Organism (Scientific)": scientific_name,
                "Organism (Common)": common_name,
                "Molecular Weight": molecular_weight
                }

            df.loc[len(df)] = new_entry

        return df

    def get_ensembl_gene_ids(self, df_proteins):
        """Get ensembl gene id of a specific gene.
        Returns input dataframe with new column Ensembl Gene ID.

        df_proteins -- dataframe with two columns: species, gene name
        """
        ensembl_ids = []

        for index, row in df_proteins.iterrows():
            species = row["Organism (Scientific)"].lower()
            gene = row["Gen"]
        
            data = self.fetch_data(species=species,gene=gene)
            
            gene_id = next((entry["id"] for entry in data if entry["type"] == "gene"))
            if gene_id == None:
                logging.warning(f"Ensembl: No gene id found for species: {species}, gene: {gene}")
                continue
            ensembl_ids.append(gene_id)
            
        df_proteins["Ensembl Gene ID"] = ensembl_ids

        return df_proteins

    def get_gene_data(self, ensembl_gene_ids: list):
        """Get additional gene data. Collect gene description and sequence region name.
        Returns dataframe with three columns: Ensembl Gene ID, Description, Seq Region Name.

        ensembl_gene_ids -- list of ensemble gene ids    
        """

        df = pd.DataFrame(columns=[
            "Ensembl Gene ID",
            "Description",
            "Seq Region Name"
            ])

        ensembl_headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
        }

        descriptions = []
        seq_region_name = []

        base_url = "https://rest.ensembl.org"
        lookup_ext = "/lookup/id/"

        for id in ensembl_gene_ids:
            
            url = f"{base_url}{lookup_ext}{id}?content-type=application/json"
            response = requests.get(url, headers=ensembl_headers)
            if response.status_code != 200:
                logging.warning(f"Ensembl: Error by Ensembl Gene ID: {id}: HTTP {response.status_code}")
                continue

            info = response.json()
            descriptions.append(info.get("description", ""))
            seq_region_name.append(info.get("seq_region_name", ""))

            if descriptions == None:
                logging.warning(f"Ensembl: Error by Ensembl Gene ID: {id}. No Description found. HTTP {response.status_code}")
                continue

            if seq_region_name == None:
                logging.warning(f"Ensembl: Error by Ensembl Gene ID: {id}. No Seq Region Name found. HTTP {response.status_code}")
                continue
    
        df["Ensembl Gene ID"] = ensembl_gene_ids
        df["Description"] = descriptions
        df["Seq Region Name"] = seq_region_name

        return df

    def saveto_xlsx(self,df,header=True,path="protein_gene_analysis.xlsx"):
        logging.info(f"Saved as protein_gene_analysis.xlsx")
        df.to_csv(path, index=False, header=header)

    def start(self):
        logging.info(f"Start program")
        uniprot_dataframe = self.get_uniprot(self.uniprot_ids)

        ensembl_dataframe = self.get_ensembl_gene_ids(uniprot_dataframe)

        ensembl_ids = ensembl_dataframe["Ensembl Gene ID"].to_list()

        gene_data_dataframe = self.get_gene_data(ensembl_ids)

        #merge dataframes
        self.final_result = pd.merge(ensembl_dataframe,gene_data_dataframe, on='Ensembl Gene ID', how='outer')

        self.saveto_xlsx(self.final_result)
        
        logging.info(f"End program")
        return self.final_result
 

if __name__ == "__main__":
    if len(sys.argv) > 2:
        print("Usage: python API_tasks.py <Q8N726,O00255,P69905,Q9Y261>")
        sys.exit(1)

    param1 = sys.argv[1]
    uniprot_ids = param1.strip("[]").split(",")

    from API_tasks import Protein_informations
    df_proteins = Protein_informations(uniprot_ids=uniprot_ids)
    df = df_proteins.start()
    print(df)
    logging.info(f"End program")
        
