import requests
import pandas as pd
import MotifClass as mc
import os

from flask import send_file



############################### WARNING : it is recommended to not modify this code, it could have terrible consequences ######################

def get_from_internet(ID):
    url = f'https://jaspar.elixir.no/api/v1/matrix/{ID}.transfac'
    response = requests.get(url)

    if response.status_code == 200:
        content = response.content.decode('utf-8')  
        info = {} 
        PFM = {
                "A" : [],
                "C" : [],
                "G" : [],
                "T" : []
            }

        for line in content.split('\n'):

            if line.startswith("AC"):
                info["TF_ID"] = str(line[3:11])
            elif line.startswith("ID"):
                info["TF_name"] = str(line[3:15])
            elif len(line) > 0 and line[0].isnumeric():
                freq = line.split("\t")
                PFM["A"].append(int(float(freq[1])))
                PFM["C"].append(int(float(freq[2])))
                PFM["G"].append(int(float(freq[3])))
                PFM["T"].append(int(float(freq[4])))
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

        PFM = pd.DataFrame(PFM)
        TF0006 = mc.Motif(info["TF_ID"],info["TF_name"],PFM)
        TF0006.taxon = info["tax_group"]
        TF0006.family = info["tf_family"]
        TF0006.TFclass = info["tf_class"]
        TF0006.validation = info["pubmed_ids"]
        TF0006.uniprotID = info["uniprot_ids"]
        TF0006.dataType = info["data_type"]

        return TF0006
    else:
        False





def get_from_database(ID):
    Properties = ["JASPAR ID", "Name", "PFM", "TFclass", "family", 
                  "collection", "taxon", "species", "dataType", 
                  "validation", "uniprotID", "source", "comment","dict_info", "sequence"]
    
    TF = ""
    with open("database.tsv", "r") as file:
        for line in file:
            info = line.strip().split("\t")
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
                

                PFM = pd.DataFrame(PFM_rendered)
                TF = mc.Motif(info[0],info[1],PFM)

                for n in range (3,len(info)):
                    
                    if info[n] != "\t" and info[n] != "\n" and info[n] != "":
                        if n != 7:
                            TF.__dict__[Properties[n]] = info[n]
                            TF.dict_info[Properties[n]] = info[n]
                        else:
                            species = list(info[n].split(","))
                            for specie in species:  
                                TF.__dict__[Properties[n]].append(specie)
                            TF.dict_info[Properties[n]] = TF.species
                    else:
                        if n == 7:
                            TF.__dict__[Properties[n]] = []
                            TF.dict_info[Properties[n]] = []
                        elif n == 13:
                            continue  
                        else:
                            TF.__dict__[Properties[n]] = None
                            TF.dict_info[Properties[n]] = None
                
    if TF != "":
        return TF
    else:
        return False 


def add_to_database(motif):
    if isinstance(motif, mc.Motif):
        info = motif.dict_info
        
        with open("database.tsv", "r+") as file:
            is_not_there = []
            file.readline()
            for line in file:
                info_motifs = line.strip().split("\t")
                
                if motif.ID != info_motifs[0]:
                    is_not_there.append(True)
                else:
                    is_not_there.append(False)
            
            if all(is_not_there):
                file.seek(0)
                line = file.readline()
                properties = line.split("\t")
                
                number = 0
                file.seek(0,2)
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
                    
                file.close()
            else:
                 
                 return False
    else:
        raise TypeError("Not a valid TF binding site motif")
    

def delete_from_database(ID):
    if mc.is_valid_ID(ID):
        with open("database.tsv", "r") as fileread:
            lines = fileread.readlines()

        
        
        with open("database.tsv", "w") as filewrite:
            for line in lines:
                linelist = line.strip(" ").split("\t")
                if linelist[0] != ID:
                    filewrite.write(line)
        
        
    else:
        raise TypeError("Not a valid JASPAR ID")



def display_database():
    with open("database.tsv", "r") as db:
        db.readline()
        motifs = []

        for line in db:
            line = line.strip(" ").split("\t")
            motifs.append(line[0])
        

        elem_in_columns = (len(motifs) + 3) // 4

        for i in range(elem_in_columns):
            for j in range(0,4):
                index = i + j * elem_in_columns
                if index < len(motifs):
                    elemento = motifs[index]
                    print(f"{elemento:<25}", end="")
            print()



def get_svg(ID):
    url = f"https://jaspar2020.genereg.net/static/logos/all/svg/{ID}.svg"
    cd = os.getcwd()
    directory = cd + f"/saved_svg/{ID}_SVG.svg"
    
    
    svg = requests.get(url)

    if svg.status_code == 200:
        svg = requests.get(url).text
        with open(directory, "w") as file:
            file.write(svg)
    
    return print(f"File is available at path: {directory}")

    

def get_svg_web(ID):
    url = f"https://jaspar2020.genereg.net/static/logos/all/svg/{ID}.svg"
    cd = os.getcwd()
    directory = cd + f"/saved_svg/{ID}_SVG.svg"
    
    
    svg = requests.get(url)

    if svg.status_code == 200:
        svg = requests.get(url).text
        with open(directory, "w") as file:
            file.write(svg)
    print(f"File is available at path: {directory}")
    return send_file(directory, as_attachment = True)


def get_matching_scores_simple(sequence1, sequence2):
    score = 0
    for n in range(0, len(sequence1)):
        if sequence1[n] == sequence2[n]:
            score += 1
    return score


def get_matching_scores_sophisticated(sequence, TF):
    probabilities = {"A" :[], "C" : [], "G":[], "T":[]}

    PFM = TF.PFM
    for n in range(0, len(PFM["A"])):
        tot = 0
        for nt in "ACGT":
            tot += PFM[nt][n]
        
        for nt in "ACGT":
            probabilities[nt].append(PFM[nt][n]/tot)

    sequence_overall_probability = 1
    for n in range(0,len(sequence)):
        probability = probabilities[sequence[n]][n]
        sequence_overall_probability = sequence_overall_probability*probability
    
    return sequence_overall_probability




def get_scores_from_database_simple(provided_sequence):
   
    same_length = {}

    with open("database.tsv", "r") as file:
        file.readline()
        for line in file:
            info = line.strip().split("\t")
            ID, sequence = info[0], info[-1]
            if len(provided_sequence) == len(sequence):
                score = get_matching_scores_simple(provided_sequence, sequence)
                same_length[ID] = score
    
    sorted_same_length = dict(sorted(same_length.items(), key=lambda item: item[1], reverse=True))
    if len(sorted_same_length) != 0:
        return sorted_same_length
    else:
        return False

def get_scores_from_database_sophisticated(provided_sequence):
    same_length = {}

    with open("database.tsv", "r") as file:
        file.readline()
        for line in file:
            info = line.strip().split("\t")
            ID, sequence = info[0], info[-1]
            if len(provided_sequence) == len(sequence):
                TF = get_from_database(ID)
                score = get_matching_scores_sophisticated(provided_sequence,TF)
                
                same_length[ID] = score
            else:
                continue
    
    sorted_same_length = dict(sorted(same_length.items(), key=lambda item: item[1], reverse = True))

    
    if len(sorted_same_length) != 0:
        
        return sorted_same_length
    else:
        return False