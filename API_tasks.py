

import logging
import pandas as pd
import requests

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("protein_API_Tasks.log"),
        logging.StreamHandler()
    ]
)

def fetch_data_uniprot(id):
    base_url = 'https://www.ebi.ac.uk/proteins/api/proteins/'
    url = base_url + id

    headers = {
    "Accept": "application/json"
    }
    
    response = requests.get(url, headers=headers)
    if response.status_code != 200:
        logging.warning(f"Uniprot: Error by calling uniprot id: {id}: HTTP {response.status_code}")

    data = response.json()
    return data


def fetch_data_ensemble(species, gene):
    url = f'https://rest.ensembl.org/xrefs/symbol/{species}/{gene}'

    headers = {
    "Accept": "application/json"
    }

    response = requests.get(url, headers=headers)
    if response.status_code != 200:
        logging.warning(f"Ensembl: Error by calling species: {species}, gene: {gene}: HTTP {response.status_code}")

    data = response.json()
    return data


def get_uniprot(uniprot_ids):
    
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
        
        data = fetch_data_uniprot(id)

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



def get_ensembl_gene_ids(df_proteins):
    ensembl_ids = []

    for index, row in df_proteins.iterrows():
        species = row["Organism (Scientific)"].lower()
        gene = row["Gen"]
       
        data = fetch_data_ensemble(species,gene)
        
        gene_id = next((entry["id"] for entry in data if entry["type"] == "gene"))
        if gene_id == None:
            logging.warning(f"Ensembl: No gene id found for species: {species}, gene: {gene}")
            continue
        ensembl_ids.append(gene_id)
        
    df_proteins["Ensembl Gene ID"] = ensembl_ids

    return df_proteins


def get_gene_data(ensembl_gene_ids: list):

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



if __name__ == "__main__":
    #uniprot_ids = ["Q8N726","O00255","P69905","Q9Y261"]

    uniprot_ids = [
    "P38398",  # BRCA1
    "P04637",  # TP53
    "P31749",  # AKT1
    "P00533",  # EGFR
    "P15056",  # BRAF
    "P24941",  # CDK1
    "P01116",  # KRAS
    "Q9Y24",  # PTEN
    "Q16539",  # MDM2
    "P25963",  # RELA
    "P07900",  # HSP90AA1
    "P01112",  # HRAS
    "P08473",  # MMP1
    "Q07812",  # BCL2L1
    "P10415",  # BCL2
    "Q9H2X6",  # STK11
    "P01106",  # MYC
    "P24385",  # CCND1
    "P38936",  # CDKN1A
    "Q9Y6K9"   # SMAD7
    ]

    result = get_uniprot(uniprot_ids)

    result2 = get_ensembl_gene_ids(result)

    ensembl_ids = result["Ensembl Gene ID"].to_list()

    result3 = get_gene_data(ensembl_ids)

    final_result = pd.merge(result2,result3, on='Ensembl Gene ID', how='outer')
    print(final_result)
  
    final_result.to_csv("output.csv", index=False, header=True)
