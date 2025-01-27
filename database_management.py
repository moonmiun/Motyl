import requests
import pandas as pd
import MotifClass as mc
import os
from flask import send_file

###############################################################################################################################################
###############################################################################################################################################
############################### WARNING : it is recommended to not modify this code, it could have terrible consequences ######################
###############################################################################################################################################
###############################################################################################################################################

def get_from_internet(ID):
    """
    This function takes JASPAR motif's information from the internet.

    Args:
    - The JASPAR ID of the motif of interest.

    Outcome:
    - Creates a JASPAR motif object of class Motif devined in the file MotifClass.py.
    """

    # The information is taken from ther JASPAR database.
    url = f'https://jaspar.elixir.no/api/v1/matrix/{ID}.transfac'
    response = requests.get(url)

    # The information is stared in a file transfac that needs to be decoded.
    if response.status_code == 200:
        content = response.content.decode('utf-8')  
        info = {} 
        PFM = {
                "A" : [],
                "C" : [],
                "G" : [],
                "T" : []
            }

        # In the transfac file the information is stared in certain lines.
        for line in content.split('\n'):

            # For example, the line starting with "AC" contains the ID of the JASPAR motif.
            if line.startswith("AC"):
                info["TF_ID"] = str(line[3:11])
            
            # Which can be confused with the "ID" line that actually contains the name of the motif.
            elif line.startswith("ID"):
                info["TF_name"] = str(line[3:15])

            # When a line starts with a number it contains the PFM. 
            elif len(line) > 0 and line[0].isnumeric():

                # Each position frequency is tab-separated.
                freq = line.split("\t")
                PFM["A"].append(int(float(freq[1])))
                PFM["C"].append(int(float(freq[2])))
                PFM["G"].append(int(float(freq[3])))
                PFM["T"].append(int(float(freq[4])))
            
            # When the line starts with "CC" it contains the rest of the information (for example, TFclass and similar)
            elif line.startswith("CC"):
                for x in range(0,len(line)):
                    if line[x]== ":":
                        info_name = line[3:x]
                        value= line[x+1:len(line)]
                        if len(value) > 0:
                            info[info_name] = value
                        else:
                            info[info_name] = None
            else:
                continue
        
        # The PFM needs to be converted in a Pandas dataframe because that's what the Motif class requires.
        PFM = pd.DataFrame(PFM)

        # The Motif element using the three essential elements that are required (JASPAR ID, name and PFM)
        TF0006 = mc.Motif(info["TF_ID"],info["TF_name"],PFM)

        # Then the rest of the attributes are added, if retrieved.
        TF0006.taxon = info["tax_group"]
        TF0006.family = info["tf_family"]
        TF0006.TFclass = info["tf_class"]
        TF0006.validation = info["pubmed_ids"]
        TF0006.uniprotID = info["uniprot_ids"]
        TF0006.dataType = info["data_type"]

        # Source, species and other types of info that are present in the JASPAR database are not added because they never appear in any downloadable file.
        # They can be added though!
        # You can use the edit functions defined in the Motyl.py file or in the edit page you can find by running Motyl_web.py.

        return TF0006
    else:
        False


def get_from_database(ID):
    """
    This function takes JASPAR motif's information from the database.

    Args:
    - The JASPAR ID of the motif of interest.

    Outcome:
    - Creates a JASPAR motif object of class Motif devined in the file MotifClass.py.
    """

    # All the properties that can (or not) have values.
    Properties = ["JASPAR ID", "Name", "PFM", "TFclass", "family", 
                  "collection", "taxon", "species", "dataType", 
                  "validation", "uniprotID", "source", "comment","dict_info", "sequence"]
    
    TF = ""

    # Opens the database file. ATTENTION: opening the .tsv file directly (and not through this code) is NOT RECOMMENDED AT ALL.
    # I've noticed that leaving the file open with the cursor set somewhere can severly impact all actions that retrieve, add, delete or updates information on it.
    # Please refrain to ever open the .tsv file (unless you want to copy it on other programs like Excel or RStudio, and still don't modify the original database).

    with open("database.tsv", "r") as file:
        for line in file:
            info = line.strip().split("\t")

            # The PFM information is stored as a string but needs to be rendered in a matrix to be able to calculate several stuff.
            if info[0] == ID:
                PFM_to_render = info[2].replace("[", "").replace("]", "").replace(",", " ").replace(":", " ")
                PFM_to_render = PFM_to_render.split()
                PFM_rendered = {"A":[],
                            "C":[],
                            "G":[],
                            "T":[]}
                indexes_nt = []
                for i in range(0, len( PFM_to_render)):
                    if  PFM_to_render[i].isalpha():
                        indexes_nt.append(i)

                start = 1
                for nt in PFM_rendered:
                    numbers =  PFM_to_render[start:start+indexes_nt[1]-1]
                    for number in numbers:
                        PFM_rendered[nt].append(int(number))
                    start += indexes_nt[1]
                
                # Once the PFM is rendered, the Motif can already be created.
                PFM = pd.DataFrame(PFM_rendered)
                TF = mc.Motif(info[0],info[1],PFM)

                # The other values are added in order.
                for n in range (3,len(info)):
                    
                    # CASE 1: There's actually information stored.
                    if info[n] != "\t" and info[n] != "\n" and info[n] != "":
                        if n != 7:
                            TF.__dict__[Properties[n]] = info[n]
                            TF.dict_info[Properties[n]] = info[n]
                        
                        # Species is always a little bit tricky to deal with so it needs its own version of the code.
                        else:
                            species = list(info[n].split(","))
                            for specie in species:  
                                TF.__dict__[Properties[n]].append(specie)
                            TF.dict_info[Properties[n]] = TF.species
                    
                    # CASE 2: there's no information stored.
                    else:

                        # When there's no information stored for the species, then an empty list is created.
                        if n == 7:
                            TF.__dict__[Properties[n]] = []
                            TF.dict_info[Properties[n]] = []
                        
                        # n == 13 is the dict_info cathegory. It needs to be created in the element but no information is stored in the database itself.
                        elif n == 13:
                            continue  
                        
                        # If nothing is present then the category takes None as argument.
                        else:
                            TF.__dict__[Properties[n]] = None
                            TF.dict_info[Properties[n]] = None

    # So if the TF was actually created, it is returned.      
    if TF != "":
        return TF
    else:
        return False 


def add_to_database(motif):
    """
    This function adds a Motif element to the local database.

    Args:
    - An element of class motif.

    Outcome:
    - Converts its information in a database-friendly version of them that can be stored in the .tsv file.
    """

    # Condition: the element must be of class Motif.
    if isinstance(motif, mc.Motif):
        
        # Gets the dictionary of the Motif information.
        info = motif.dict_info
        
        with open("database.tsv", "r+") as file:
            
            # First it checks if the motif is already present or not.
            is_not_there = []
            file.readline()
            for line in file:
                info_motifs = line.strip().split("\t")
                
                if motif.ID != info_motifs[0]:
                    is_not_there.append(True)
                else:
                    is_not_there.append(False)
            
            if all(is_not_there):

                # Looks for the line in which the ID is present.
                file.seek(0)
                line = file.readline()
                properties = line.split("\t")
                
                number = 0
                file.seek(0,2)

                # One by one the elements in the dict_info are added, following the stored order.
                for elem in info:
                    if elem != "dict_info" and elem != "sequence":
                        if elem == properties[number] and isinstance(info[elem], list):
                            if info[elem] != []:
                                string = f''
                                n = 0
                                while n < len(info[elem])-1:
                                    string += f'{info[elem][n]},'
                                    n += 1
                                
                                string += f'{info[elem][n]}\t'
                                file.write(string)
                                number += 1
                            else:
                                file.write("\t")
                                number += 1
                        elif elem == properties[number] and isinstance(info[elem], pd.DataFrame):
                            PFM_string = ""
                            for column in info[elem].columns:
                                    frequencies = ""
                                    for n in info[elem][column]:
                                            frequencies += f"{n},"
                                    PFM_string += f"{column}:[{frequencies}]"
                            PFM_string += "\t"
                            file.write(PFM_string)
                            number += 1
                        elif elem == properties[number] and info[elem] is not None and info[elem] != "0":
                            file.write(f'{str(info[elem])}\t')
                            number += 1 
                        else:
                            file.write(f'\t')
                            number += 1
                    elif elem == "dict_info":
                        file.write(f'\t')
                        number += 1
                    elif elem == "sequence":
                        file.write(f'{str(info[elem])}\n')

                # ESSENTIAL: closes the file. DO NOT touch.
                file.close()
            else:
                return False
    else:
        raise TypeError("Not a valid TF binding site motif")
    


def delete_from_database(ID):
    """
    This function deletes a Motif from the local database.

    Args:
    - The ID of the motif of interest.

    Outcome:
    - Eliminates the line containing all the information of that motif from the .tsv file.
    """

    # Condition: it must be a valid JASPAR motif ID, it will not look for a motif that already doesn't have a valid ID thus cannot exist.
    if mc.is_valid_ID(ID):

        # Stores the lines to access them
        with open("database.tsv", "r") as fileread:
            lines = fileread.readlines()
        
        # Recreates the file skipping the line where the motif was
        with open("database.tsv", "w") as filewrite:
            for line in lines:
                linelist = line.strip(" ").split("\t")
                if linelist[0] != ID:
                    filewrite.write(line)
        
    else:
        raise TypeError("Not a valid JASPAR ID")



def display_database():
    """
    This function displays the database.

    Args:
    - None.

    Outcome:
    - Visualizes the database in a clear tabular form.
    """

    # Opens the database for reading.
    with open("database.tsv", "r") as db:
        db.readline()
        motifs = []

        # The ID of the motifs is stored.
        for line in db:
            line = line.strip(" ").split("\t")
            motifs.append(line[0])
        
        # There are 4 colums, and the number of rows is calculated.
        elem_in_columns = (len(motifs) + 3) // 4

        # Tabular view of the motifs.
        for i in range(elem_in_columns):
            for j in range(0,4):
                index = i + j * elem_in_columns
                if index < len(motifs):
                    elemento = motifs[index]
                    print(f"{elemento:<25}", end="")
            print()



def get_svg(ID):
    """
    This function saves the .svg file of the nucleotides distribution of a Motif in the local memory.

    Args:
    - The motif ID.

    Outcome:
    - Saves the .svg file and displays its path in the local memory.
    """

    # Creation of the url
    url = f"https://jaspar2020.genereg.net/static/logos/all/svg/{ID}.svg"

    # Finds out where the Motyl folder is stored.
    cd = os.getcwd()

    # Creates the file path where to store the .svg file. 
    # Note: even if the /saved_svg folder might be empty, DO NOT delete, it's presence is essential for the code to run.
    directory = cd + f"/saved_svg/{ID}_SVG.svg"
    
    
    svg = requests.get(url)

    # If a .svg is found, the image is saved in the specified directory.
    if svg.status_code == 200:
        svg = requests.get(url).text
        with open(directory, "w") as file:
            file.write(svg)
    
    return f"If a match was found, the file will be available at path: {directory}\n"

    

def get_svg_web(ID):
    """
    This function saves the .svg file of the nucleotides distribution of a Motif in the local memory.
    It's slightly different from before because this is the function used by the web version of Motyl.

    Args:
    - The motif ID.

    Outcome:
    - Saves the .svg file in the Download folder.
    """

    url = f"https://jaspar2020.genereg.net/static/logos/all/svg/{ID}.svg"
    cd = os.getcwd()
    directory = cd + f"/saved_svg/{ID}_SVG.svg"
    
    svg = requests.get(url)

    if svg.status_code == 200:
        svg = requests.get(url).text
        with open(directory, "w") as file:
            file.write(svg)
    print(f"If a match was found, the file will be available at path: {directory}\n")
    return send_file(directory, as_attachment = True)



def get_matching_scores_simple(sequence1, sequence2):
    """
    This function calculates the matching scores between two sequences in a simple way (see documentation).
    It is simply the number of matching bases at the same position.

    Args:
    - Sequence 1: the first sequence.
    - Sequence 2: the second sequence (reference).

    Outcome:
    - The matching score.
    """

    score = 0
    for n in range(0, len(sequence1)):
        if sequence1[n] == sequence2[n]:
            score += 1
    return score


def get_matching_scores_sophisticated(sequence, TF):
    """
    This function calculates the matching scores between two sequences in a more complicated way (see documentation).
    It can only be used if an element of class Motif is given because it requires using its PFM.

    Args:
    - Sequence: the sequence provided.
    - TF: the element of class Motif (reference).

    Outcome:
    - The ratio between the probability of generating the given sequence and the probability of generating the consensus sequence for that Motif.
    """

    probabilities = {"A" :[], "C" : [], "G":[], "T":[]}

    PFM = TF.PFM

    # This calculates the probabilities of finding a certain nucleotide at a certain position.
    for n in range(0, len(PFM["A"])):
        tot = 0
        for nt in "ACGT":
            tot += PFM[nt][n]
        for nt in "ACGT":
            probabilities[nt].append(PFM[nt][n]/tot)

    sequence_overall_probability = 1

    # The probability of generating the given sequence knowing the probabilities at each position, is calculated.
    # Then the ratio between it and the probability of generating the consensus sequence for that Motif is calculated.
    for n in range(0,len(sequence)):
        probability = probabilities[sequence[n]][n]
        sequence_overall_probability = sequence_overall_probability*probability
    
    return sequence_overall_probability




def get_scores_from_database_simple(provided_sequence):
    """
    This function repeates get_matching_scores_simple(sequence1, sequence2) function for each Motif in the database which consensus sequence length is the same as the sequence provided.

    Args:
    - provided_sequence: the sequence provided by the user.

    Outcome:
    - A (sorted by score) dictionary of the IDs and scores of those motifs which consensus sequence is the same length as the one provided.
    """
   
    same_length = {}

    with open("database.tsv", "r") as file:
        file.readline()
        for line in file:
            info = line.strip().split("\t")
            ID, sequence = info[0], info[-1]
            if len(provided_sequence) == len(sequence):
                # Calculates the score for each line and stores it in a dictionary.
                score = get_matching_scores_simple(provided_sequence, sequence)
                same_length[ID] = score
    
    # Sorts the dictionary
    sorted_same_length = dict(sorted(same_length.items(), key=lambda item: item[1], reverse=True))

    # If the dictionary is not empty, it is returned.
    if len(sorted_same_length) != 0:
        return sorted_same_length
    else:
        return False

def get_scores_from_database_sophisticated(provided_sequence):
    """
    This function repeates get_matching_scores_sophisticated(sequence, TF) function for each Motif in the database which consensus sequence length is the same as the sequence provided.

    Args:
    - provided_sequence: the sequence provided by the user.

    Outcome:
    - A (sorted by score) dictionary of the IDs and scores of those motifs which consensus sequence is the same length as the one provided.
    """
   
    same_length = {}

    with open("database.tsv", "r") as file:
        file.readline()
        for line in file:
            info = line.strip().split("\t")
            ID, sequence = info[0], info[-1]
            if len(provided_sequence) == len(sequence):
                # The element of Motif class needs to be created because the PFM is required.
                TF = get_from_database(ID)
                score = get_matching_scores_sophisticated(provided_sequence,TF)
                same_length[ID] = score
            else:
                continue
    
    # Sorts the dictionary
    sorted_same_length = dict(sorted(same_length.items(), key=lambda item: item[1], reverse = True))
    
    # If the dictionary is not empty, it is returned.
    if len(sorted_same_length) != 0:
        return sorted_same_length
    else:
        return False
    

###############################################################################################################################################
###############################################################################################################################################
############################################ I hope you didn't modify this code. I trust you, ok? #############################################
###############################################################################################################################################
###############################################################################################################################################