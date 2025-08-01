# quattro_research 
## Developer Task: Integrating UniProt and Ensembl APIs

### Description 
This integration of UniProt and Ensembl APIs collects informations about proteins given their UniProt ID.
The following informations are stored:
- Protein Accession,
- Protein Name,
- Gen,
- Organism (Scientific),
- Organism (Common),
- Molecular Weight,
- Ensembl Gene ID,
- Description,
- Sequence region name

The return value is a pandas dataframe of size: _[number of UniProt IDs rows x 9 columns]_.

### Installation 
Use the virtuel environment called virtuelenv to get all neccessary packages.

Poetry is also available als dependency controle.

### Usage
### In python scripts
Import Protein_inforation class.
```python
from API_tasks import Protein_informations
```

Define UniProt Ids as list.
```python
uniprot_ids = ["Q8N726","O00255","P69905","Q9Y261"]
```

Create Protein_information object.
```python
df_proteins = Protein_informations(uniprot_ids=uniprot_ids)
```

Start information gaining and return result dataframe.
```python
df = df_proteins.start();
```

Save dataframe to .xlsx.
```python
df_proteins.saveto_xlsx(df)
```

### Command line
```Powershell
_python API_tasks.py "Q8N726,O00255,P69905,Q9Y261"_
```

_Param 1_ is a list of UniProt IDs 

### Output
| Protein Accession | Protein Name                   | Gen     | Organism (Scientific) | Organism (Common) | Molecular Weight | Ensembl Gene ID     | Description                                                    | Seq Region Name |
|-------------------|--------------------------------|---------|------------------------|--------------------|-------------------|----------------------|----------------------------------------------------------------|-----------------|
| Q9Y261            | Hepatocyte nuclear factor 3-beta | FOXA2   | Homo sapiens           | Human              | 48306             | ENSG00000125798      | forkhead box A2 [Source:HGNC Symbol;Acc:HGNC:5...]             | 20              |
| O00255            | Menin                          | MEN1    | Homo sapiens           | Human              | 67497             | ENSG00000133895      | menin 1 [Source:HGNC Symbol;Acc:HGNC:7010]                      | 11              |
| Q8N726            | Tumor suppressor ARF           | CDKN2A  | Homo sapiens           | Human              | 13903             | ENSG00000147889      | cyclin dependent kinase inhibitor 2A [Source:HGNC Symbol;...]  | 9               |
| P69905            | Hemoglobin subunit alpha       | HBA1    | Homo sapiens           | Human              | 15258             | ENSG00000206172      | hemoglobin subunit alpha 1 [Source:HGNC Symbol;...]            | 16              |



### Detailed informations of functions of class Protein_information
#### fetch_data() 
```python
def fetch_data(self, id=None, species=None, gene=None)
```
Access to websites of UniProt (https://www.ebi.ac.uk/proteins/api/proteins/id) and Ensembl (https://rest.ensembl.org/xrefs/symbol/{species}/{gene}).

For UniProt add the UniProt ID at the end of the link. 

For Ensembl add species name and gene name at the end of the link. This informations are available after calling the function _def get_uniprot(self, uniprot_ids: list)_.

Handle request limit of APIs by waiting 3s if limit is reached. 

Use _json_ packages to handle the website output.


#### get_uniprot() 
```python
def get_uniprot(self, uniprot_ids: list)
```
Return the following informations about a UniProt ID:
- Protein Accession:,
- Protein Name,
- Gen,
- Organism (Scientific),
- Organism (Common),
- Molecular Weight

Protein Accession is a string, easy access per data.get("accession")

Protein Name is a dictionary, access via keys proteins -> recommendedName -> fullName -> value.

Gene is a list of dictionarys, access via gene -> list at position 0 -> then use key _name_ -> value.

Organism name is a dictionary, access via keys organism -> name -> position 0 in list is scientific name & position 1 is common name

Molecular weight, is not direct accessible. Got to sequence -> mass. 

Returns a dataframe with _[number of UniPort IDs x 6 columns]_

#### get_ensembl_gene_ids() 
Return the following informations about a UniProt ID by using Ensembl:
- Ensembl Gene ID

Take the reuslt dataframe from _get_uniprot()_ and collect data from Ensembl. Need _organism scientific_ and _gene_ to query the database.

Collect a gene ID by finding the _id_ with type _gene_ of the ensembl entry of current gene. 

Add a new column to the result dataframe of _get_uniprot()_ with Ensembl Gene IDs for each UniProt ID.


#### get_gene_data()
Return the following informations about a UniProt ID:
- Description,
- Sequence region name

Take the result dataframe with new Ensembl ID column of _get_ensembl_gene_ids()_. 

Find description and sequence region name of a gene by fetching the rest.ensembl.org with an given ensembl id.
```python
base_url = "https://rest.ensembl.org"
lookup_ext = "/lookup/id/"
f"{base_url}{lookup_ext}{id}?content-type=application/json"
```
Description and Sequence region name are easy accessible via .get().

#### saveto_xlsx() 
Save the result dataframe as .xlsx file. Default name is _protein_gene_analysis.xlsx_
Header and index are optional.

#### Merge dataframes 
Merge dataframe of _get_ensembl_gene_ids()_ and _get_gene_data()_.
Using an outer join on id _Ensembl Gene ID_.

