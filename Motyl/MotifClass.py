import pandas as pd

nucleotides = {"A", "C", "G", "T"}

def is_valid_ID(ID):
    if isinstance(ID, str) and len(ID) == 8 and ID.startswith("MA") and ID[6] == "." and ID[2:6].isnumeric() and ID[7].isnumeric():
        return True
    else:
        return False

################ Class Motif definition ####################
class Motif:
    def __init__(self, JASPARID, TFname, PFM, dict_info = {}):
        self.dict_info = dict_info
        self.ID = self.set_ID(JASPARID)
        self.PFM = self.set_PFM(PFM)  
        self.name = TFname
        self.TFclass = None
        self.family = None
        self.collection = None
        self.taxon = None
        self.species = []
        self.dataType = None
        self.validation = None
        self.uniprotID = None
        self.source = None
        self.comment = None
        self.sequence = self.set_sequence()
        self.dict_info = { "JASPAR ID" : self.ID,
                          "Name" : self.name,
                          "PFM" : self.PFM,
                          "TFclass" : self.TFclass,
                          "family" : self.family,
                          "collection" : self.collection,
                          "taxon" : self.taxon,
                          "species" : self.species,
                          "dataType" : self.dataType,
                          "validation" : self.validation,
                          "uniprotID" : self.uniprotID,
                          "source" : self.source,
                          "comment" : self.comment,
                          "dict_info" : "",
                          "sequence" : self.sequence,
                            }
    

    def __str__(self):
        message = f"JASPAR ID: {self.ID}\nName: {self.name}\n"

        PFM_string = "PFM (Position Frequency Matrix):\n"
        for column in self.PFM.columns:
            frequencies = ""
            for n in self.PFM[column]:
                frequencies += f"{n:>7}"
            PFM_string += f"{column} : [{frequencies}]\n"

        message += PFM_string
        
        for attribute in ["TFclass", "family", "collection", "taxon", "species","dataType","validation","uniprotID","source","comment", "sequence"]:
            value = getattr(self, attribute)
            if attribute != "TFclass" and attribute != "species" and value == None:
                message += f"{attribute.capitalize()}: - \n"
            elif value == None and attribute == "TFclass":
                message += f"TFclass: -\n"
            else:
                if attribute == "species":
                    if value == []:
                        message += "Species: - \n"
                    else:
                        species = "Species: "
                        for specie in self.species:
                            species += f"'{specie}' "
                        message += f"{species}\n"
                elif attribute == "TFclass":
                    message += f"TFclass: {value}\n"
                else:
                    message += f"{attribute.capitalize()}: {value}\n"
        
        return message


    def set_PFM(self, PFM):
        def validmatrix(nt):
            len_A = len(PFM["A"])
            if len(PFM[nt]) == len_A:
                if PFM[nt].isnull().any():
                    return False
                else:
                    return True
            else:
                return False

        if not isinstance(PFM, pd.DataFrame):
            raise TypeError("Error: PFM has to be a Pandas DataFrame")
        elif set(PFM.columns) != nucleotides or len(PFM.columns) != 4:
            raise TypeError("Error: PFM contains not valid nucleotides")
        elif not all(map(validmatrix, nucleotides)):
            raise TypeError("Error: Nucleotides frequencies are invalid")
        else: 
            return PFM

    def set_ID(self,ID):
        if is_valid_ID(ID):
            return ID
        else:
            raise TypeError("""Error: not a valid JASPAR ID.
                            A valid JASPAR ID is in a 'MAXXXX.X' format.""")
    
    def set_sequence(self):
        length = len(self.PFM["A"])
        sequence = ""
        for n in range(0,length):
            positions = {"A":0,"C":0,"G":0,"T":0}
            for nt in "ACGT":
                positions[nt] = self.PFM[nt][n]
            sequence += max(positions,key=positions.get)
        return sequence    

 
    def __setattr__(self, name, value):
        if name == "species":
            if hasattr(self, "species"):
                if isinstance(value, str):
                    if isinstance(self.species, list):
                        answer = str(input("Did you mean add_species? (Y/N) "))
                        if answer in "Yy" and value not in self.species:
                            self.species.append(value)
                            self.dict_info['species'] = self.species
                        else:
                            print("Specie not added to the list.")
                    else:
                        super().__setattr__("species",[value,])
                        self.dict_info['species'] = self.species

                elif isinstance(value, list) and isinstance(self.species, list):
                    if self.species != []:
                        answer = str(input("Warning: species attribute already set. Do you want to add the provided specie(s) to the list? (Y/N) "))
                        if answer in "Yy":
                            return self.add_species(value)
                        else:
                            print("Specie(s) not added to the list")
                    else:
                        super().__setattr__("species",value)
                        self.dict_info['species'] = self.species
                else:
                    raise TypeError("Not valid species attribute")
            else:
                self.__dict__[name] = value
        elif name == "validation": 
            if hasattr(self, "validation"):
                try:
                    value = int(value)
                    super().__setattr__("validation",int(value))
                    self.dict_info['validation'] =value
                except ValueError:
                    raise TypeError("Validation attribute should be an integer number")     
            else:
                self.__dict__[name] = value
                
        else:
            self.__dict__[name] = value
            self.dict_info[name] = value
            super().__setattr__(name, value)
                

    
    def add_species(self, species):
        if isinstance(self.species, list):
            if isinstance(species, list):
                for s in species:
                    if s not in self.species and isinstance(s,str):
                        self.species.append(s)
                        
                    else:
                        if s in self.species:
                            print(f"{s} is already present in the species associated to {self.name} site")
                        else:
                            print(f"{s} is not a valid species, it will be skipped and not added to the list of species.")
                self.dict_info['Species'] = self.species
            elif isinstance(species, str):
                if species not in self.species:
                    self.species.append(species)
                    self.dict_info['Species'] = self.species
                else:
                    print(f"{species} is already present in the species associated to {self.name} site")
                self.dict_info['Species'] = self.species
            else:
                raise TypeError("The species argument should be a string or a list of strings")
        
    
    def remove_species(self,species):
        if isinstance(species, list):
            for n in species:
                if isinstance(n, str):
                    if n in self.species:
                        self.species.remove(n)
                    else:
                        print(f"{n} is not present in the species associated to {self.name} site")
                self.dict_info['Species'] = self.species
        elif isinstance(species, str):
            if species in self.species:
                self.species.remove(species)
                self.dict_info['Species'] = self.species
            else:
                print(f"{n} is not present in the species associated to {self.name} site")
        else:
            raise TypeError("The species argument should be a string or a list of strings")
        




def create_motif_from_info(dictionary):
    if isinstance(dictionary, dict):
        new_TF = Motif(dictionary["JASPAR ID"],dictionary["Name"], dictionary["PFM"])
        new_TF.TFclass = dictionary["TFclass"]
        new_TF.family = dictionary["family"]
        new_TF.collection = dictionary["collection"]
        new_TF.taxon = dictionary["taxon"]
        new_TF.species = dictionary["species"]
        new_TF.dataType = dictionary["dataType"]
        new_TF.validation = dictionary["validation"]
        new_TF.uniprotID = dictionary["uniprotID"]
        new_TF.comment = dictionary["comment"]
        
        return new_TF



##################################################################################################