
#ref https://stackoverflow.com/questions/72803530/http-error-405-for-code-that-previously-worked-to-query-uniprot-id-mapping
import requests
import time
import json

POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"

class ID_Mapping():

    def __init__(self, uniprot_id, toDatabase="UniParc"):
        self.uniprot_id = uniprot_id
        self.toDatabase = toDatabase


    def submit_id_mapping(self, fromDB, toDB, ids):
        r = requests.post(
            f"{API_URL}/idmapping/run", data={"from": fromDB, "to": toDB, "ids": ids},
        )
        r.raise_for_status()
        return r.json()["jobId"]


    def get_id_mapping_results(self, job_id):
        while True:
            r = requests.get(f"{API_URL}/idmapping/status/{job_id}")
            r.raise_for_status()
            job = r.json()
            if "jobStatus" in job:
                if job["jobStatus"] == "RUNNING":
                    print(f"Retrying in {POLLING_INTERVAL}s")
                    #time.sleep(POLLING_INTERVAL)
                else:
                    raise Exception(job["jobStatus"])
            else:
                return job

    def id_mapping(self):
        job_id = self.submit_id_mapping(
            fromDB="UniProtKB_AC-ID", toDB=self.toDatabase, ids=self.uniprot_id
        )
        results = self.get_id_mapping_results(job_id)

        accessions = results["results"][0]["to"]["uniProtKBAccessions"]

        ids = []
        for acc in accessions:
            if "." not in acc:
                if "-" not in acc:
                    ids.append(acc)
        return ids

# if __name__ == "__main__":
#     alternative_ids = ID_Mapping("Q96A13","UniParc")
#     alternative_ids.id_mapping()
