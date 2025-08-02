class Protein():
    '''Protein object with attributes. Accessable with set methods and assignable with get methods.
    
    '''

    def __init__(self, accession_number):
        self.accssion_number = accession_number
        self.protein_name = None
        self.gen = None
        self.organism_scitific = None
        self.organism_common = None
        self.mol_weight = None
        self.ensembl_gene_id = None
        self.description = None
        self.seq_reg_name = None
    
    def get_accession_number(self):
        return self.accssion_number

    def set_accsssion_number(self, value):
        self.accssion_number = value

    def get_protein_name(self):
        return self.protein_name

    def set_protein_name(self, value):
        self.protein_name = value

    def get_gen(self):
        return self.gen

    def set_gen(self, value):
        self.gen = value

    def get_organism_scitific(self):
        return self.organism_scitific

    def set_organism_scitific(self, value):
        self.organism_scitific = value

    def get_organism_common(self):
        return self.organism_common

    def set_organism_common(self, value):
        self.organism_common = value

    def get_mol_weight(self):
        return self.mol_weight

    def set_mol_weight(self, value):
        self.mol_weight = value

    def get_ensembl_gene_id(self):
        return self.ensembl_gene_id

    def set_ensembl_gene_id(self, value):
        self.ensembl_gene_id = value

    def get_description(self):
        return self.description

    def set_description(self, value):
        self.description = value

    def get_seq_reg_name(self):
        return self.seq_reg_name

    def set_seq_reg_name(self, value):
        self.seq_reg_name = value

