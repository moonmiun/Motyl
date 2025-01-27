import pandas as pd

###############################################################################################################################################
###############################################################################################################################################
############################### WARNING : it is recommended to not modify this code, it could have terrible consequences ######################
###############################################################################################################################################
###############################################################################################################################################


nucleotides = {"A", "C", "G", "T"}

def is_valid_ID(ID):
    """
    Checks if the JASPAR ID follows a certain format (MAXXXX.X)
    
    Args:
    - ID: the JASPAR ID to check.
    
    Outcome;
    - boolean to know if the check is passed (True) or not (False)
    """

    if isinstance(ID, str) and len(ID) == 8 and ID.startswith("MA") and ID[6] == "." and ID[2:6].isnumeric() and ID[7].isnumeric():
        return True
    else:
        return False

###############################################################################################################################################
################################################################ Class Motif definition #######################################################
###############################################################################################################################################

class Motif:
    def __init__(self, JASPARID, TFname, PFM, dict_info = {}):
        """
        To define this class, the three essential informations needed are JASPAR ID, TF name and the PFM.
        This data is retrieved from the JASPAR database and CANNOT be modified, the rest can.

        - JASPARID: a string that needs to follow the correct format.
        - TFname: a string.
        - PFM: a Pandas matrix that needs to follow the correct format.
        - dict_info = empty dictionary that is going to store all the information in a way that is easy to access (important for the database).
        """

        # Initializing all the attributes
        # dict_info must be initialized twice or it gives an error.
        self.dict_info = dict_info
        self.ID = self.set_ID(JASPARID)
        self.PFM = self.set_PFM(PFM)  
        self.name = TFname
        self.TFclass = None
        self.family = None
        self.collection = None
        self.taxon = None

        # Species is the only attribute that's a list because it can take more than one value.
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
        """
        Defines a way to print the motif information in a user friendly way.

        Args:
        - None.

        Outcome:
        - a printable string.
        """

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
        """
        This function checks if the Pandas matrix follows the required format for a PFM. (Warnign: PFM should not be changed directly).

        Args:
        - PFM: a Pandas matrix.

        Ourcome:
        - If all conditions are met, the PFM is set.
        """

        def validmatrix(nt):
            """
            Function to check if the nucleotides in the matrix are valid or not.

            Args:
            - nt: the nucleotide to check

            Outcome:
            - boolean that shows if the check is successful (True) or not (False).
            """
            len_A = len(PFM["A"])
            # All lines should be of the same length.
            if len(PFM[nt]) == len_A:
                # There cannot be missing values (0 is accepted).
                if PFM[nt].isnull().any():
                    return False
                else:
                    return True
            else:
                return False

        # If the PFM is not a pandas data frame, it raises an error.
        if not isinstance(PFM, pd.DataFrame):
            raise TypeError("Error: PFM has to be a Pandas DataFrame")
        
        # If the columns contain invalid nucleotides (ex: Y) or if you don't have exactly 4 columns (one for each ACGT nucleotide), it raises an error.
        elif set(PFM.columns) != nucleotides or len(PFM.columns) != 4:
            raise TypeError("Error: PFM contains not valid nucleotides")
        
        # If not all nucleotides pass the validity check, it raises an error.
        elif not all(map(validmatrix, nucleotides)):
            raise TypeError("Error: Nucleotides frequencies are invalid")
        else: 
            return PFM


    def set_ID(self,ID):
        """
        Checks if the JASPAR ID follows the correct format.

        Args:
        - ID: the JASPAR ID.

        Outcome:
        - ID if it's valid, an error if not.
        """
        if is_valid_ID(ID):
            return ID
        else:
            raise TypeError("""Error: not a valid JASPAR ID.
                            A valid JASPAR ID is in a 'MAXXXX.X' format.""")
    

    def set_sequence(self):
        """
        Automatically calculates the consensus sequence from the PFM.
        
        Args:
        - None.

        Outcome:
        - a string containing the consensus sequence.
        """
        length = len(self.PFM["A"])
        sequence = ""
        for n in range(0,length):
            positions = {"A":0,"C":0,"G":0,"T":0}
            for nt in "ACGT":
                positions[nt] = self.PFM[nt][n]
            sequence += max(positions,key=positions.get)
        return sequence    

 
    def __setattr__(self, name, value):
        """
        Lots of checks on the validity of the given attributes.
        
        Args:
        - name: name of the attributes (ex: species, class, etc...)
        - value: the attribute itself.

        Outcome:
        - if the checks are passed, then the attribute is set as the value given.
        """
        if name == "species":
            if hasattr(self, "species"):
                if isinstance(value, str):
                    if isinstance(self.species, list):

                        # There are functions to change the species attribute because it's a rather complicated attribute.
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
        """
        This function lets the user add species to the list already present.
        
        Args:
        - species: can be either a string or a list of species to add.
        
        Outcome:
        - changes the "species" attribute.
        """

        if isinstance(self.species, list):
            if isinstance(species, list):
                for s in species:
                    if s not in self.species and isinstance(s,str):
                        self.species.append(s)
                        
                    else:
                        # Of course if the species is already there it's not added.
                        if s in self.species:
                            print(f"{s} is already present in the species associated to {self.name} site")
                        else:
                            # I can only call a species "invalid" if it's not a string. I cannot check the existance of all possible species.
                            # Only adding valid and existing species is left to the user.
                            print(f"{s} is not a valid species, it will be skipped and not added to the list of species.")
                
                self.dict_info['Species'] = self.species
            
            # When adding only one species, the situation is simpler but the conditions checked are the same.
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
        """
        This function lets the user remove species to the list already present.
        
        Args:
        - species: can be either a string or a list of species to remove.
        
        Outcome:
        - changes the "species" attribute.
        """

        if isinstance(species, list):
            for n in species:
                # Samely as before, it checks for the same confitions.
                if isinstance(n, str):
                    if n in self.species:
                        self.species.remove(n)
                    else:
                        print(f"{n} is not present in the species associated to {self.name} site")
                self.dict_info['Species'] = self.species
        
        # Case in which only one species is given.
        elif isinstance(species, str):
            if species in self.species:
                self.species.remove(species)
                self.dict_info['Species'] = self.species
            else:
                print(f"{n} is not present in the species associated to {self.name} site")
        else:
            raise TypeError("The species argument should be a string or a list of strings")
        
###############################################################################################################################################
#################################################### End class Motif construction #############################################################
###############################################################################################################################################


def create_motif_from_info(dictionary):
    """
    Useful function to create an element of class Motif from a dictionary containing all the information needed.

    Args:
    - dictionary: dictionary containing the information. In the keys the name of the properties and in the values the... values.

    Outcome:
    - an element of class Motif.
    """
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

###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################