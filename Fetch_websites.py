

import logging
import time
import requests
from ID_Mapping import ID_Mapping

class Fetch_websites():

    def fetch_data(case,id=None, species=None, gene=None):
        """ Fetch the data from uniprot and ensebl.

        id -- uniprot id, only for uniprot
        species -- species name organism (scientific), only for ensembl
        gene -- gene name, only for ensembl
        """

        if case == "uniprot":
            BASE_URL = 'https://www.ebi.ac.uk/proteins/api/proteins/'
            url = BASE_URL + id
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
                if case == "ensembl":
                    logging.warning(f"No Ensemble Gene ID for {species, gene}: HTTP {response.status_code}")
                else:
                    logging.warning(f"Uniprot ID: {id}. Is not valid: HTTP {response.status_code}")
                return None
            case 404:
                logging.warning(f"Uniprot ID: {id}. Is not valid: HTTP {response.status_code}")
                logging.info(f"UniProt ID: {id}. Try to find mapped IDs")
                mapping_ids = ID_Mapping(uniprot_id=id,toDatabase="UniParc")
                alternative_ids = mapping_ids.id_mapping()
                if len(alternative_ids) != 0: 
                    logging.info(f"UniProt ID: {id}. Has been merged. New IDs: {alternative_ids}")
                else:
                    logging.info(f"UniProt ID: {id}. No new IDs.")
                return None
            case 200:
                data = response.json()
                return data
            case _:
                logging.warning(f"Error by calling UniProt ID: {id} or gene {gene}: HTTP {response.status_code}")
                return None