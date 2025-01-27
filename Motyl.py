# Libraries used
import argparse
import pandas as pd
import MotifClass as mc
import database_management as dbm
import math
import requests
import sys
import os
from flask import Flask, jsonify, request, url_for, send_file, redirect
import json


# Define the port. 5000 is the default for it but it can be customizable.
port = 5000

# These are the conditions that are checked to ensure that the port is a valid number.
if not isinstance(port, int) or port < 0 or port > 65535:
    print(f"Warning: {port} is not a valid port number. Default port number set to 5000.")
    port = 5000

Motyl = Flask(__name__)
Motyl.static_folder = 'static'
Motyl.config['DEBUG'] = False
Motyl.config['PROPAGATE_EXCEPTIONS'] = True

##################################################################################################################################################
##################################################################### ROUTES #####################################################################
################################################## WARNING: Changing this code is NOT advisable ##################################################
##################################################################################################################################################
##################################################################################################################################################




##################################################################### WELCOME ####################################################################
@Motyl.route('/')
def welcome():
    """

    Simply displays a welcome message.

    Args:
    - None

    Outcome:
    - A string containing a welcome message.
    
    """

    message = f"""Welcome to Motyl, the Web service for finding and manage TF binding site motifs."""
    return message




#################################################################### DATABASE ####################################################################
@Motyl.route('/motifs', methods = ['GET'])
# Example: curl http://localhost:5000/motifs
def get_motifs():
    """
    
    Displays the motifs present in the local database (database.tsv).

    Args:
    - None

    Output:
    - A string that contains one or more motif ID, sorted in ascending order.

    """
    motifs_data = []

    # Opens the database
    with open("database.tsv", "r") as file:
        
        # Reads the database
        file.readline()

        # Goes through each line and stores all IDs in a list.
        for line in file:
            info = line.strip().split("\t")
            TF = dbm.get_from_database(info[0])
            motifs_data.append(TF)
    
    # Sorts the list.
    motifs_data = sorted(motifs_data, key = lambda x : x.ID)
    
    # Creates the string containing all IDs of the JASPAR motifs present in the local database.
    motifs = """"""
    for motif in motifs_data:
        motifs += f'{motif.ID}\n'
    return motifs




################################################################ DISPLAY MOTIF(S) ################################################################
@Motyl.route('/motifs/<JASPARID>', methods = ['GET'])
# Example: curl http://localhost:5000/motifs/MA0001.1
# If you want to display more than one motif make sure to divide them with commas, for example: curl http://localhost:5000/motifs/MA0001.1,MA0002.1
def get_motif(JASPARID):
    """
    Displays all the information available in the local database about one or more JASPAR motif.

    Args:
    - The ID(s) of the motif(s).

    Output:
    - A string message containing all informations about the requested JASPAR motif(s).

    """

    ids = JASPARID
    response = """"""
    
    # Manages the case in which more than one motif was given.
    ids = JASPARID.replace(" ","").split(",")
    
    for id in ids:
        # Condition: each JASPAR ID needs to be in the valid "MAXXXX.X" format.
        if mc.is_valid_ID(id) == False:
            response += f"{id} is not a valid JASPAR ID, skipping it...\n"
        
        else:
            # Motyl first tries to retrieve the information from the local database.
            motif = dbm.get_from_database(id)

            # If it finds the ID, it prints the infotmation.
            if motif:
                response += f'{motif}\n'

            else:
                # If it doesn't find it, first tries to look for the motif on the JASPAR database online.
                response += f"{id} not present in the local database, looking for it on the internet...\n"
                motif = dbm.get_from_internet(id)

                # If it's present there, it prints the information.
                if motif:
                    response += f'{motif}\n'
                
                # Else it notifies the absence of that JASPAR ID even on the JASPAR database.
                else:
                    response += f"{id} wasn't found on the local database nor the internet.\n"
    return response
    


############################################################### SAVING FROM INTERNET #############################################################
@Motyl.route('/save_in_database', methods = ['POST'])
# Example: curl -H "Content-type: application/json" -X POST -d "{\"id\":\"MA0992.1,MA0032.1,MA9999.1,MAMMA.MIA\"}" http://127.0.0.1:5000/save_in_database
# As shown in the example, if you wish to save more than one motif, make sure to divide the IDs with commas.
def save_from_internet():
    """
    Gets the information about JASPAR motifs fromt he JASPAR database online and stores them in the local database.

    Args:
    - None

    Output:
    - For each ID, a string stating if it was skipped, saved or not found on the JASPAR online server.
    """

    # Manages the case in which more than one ID was given.
    ids = request.json['id'].replace(" ","").split(",")
    response = """"""

    for id in ids:
        # Condition: each JASPAR ID needs to be in the valid "MAXXXX.X" format.
        if mc.is_valid_ID(id) == False:
           response += f"{id} is not a valid JASPAR ID, skipping it...\n"
        
        # Motyl first checks if the motif's ID is already present in the local database.
        else:
            motif = dbm.get_from_database(id)

            # If it finds it, it doesn't save it again.
            if motif:
                response += f"{id} already present in the local database.\n"

            # Else, it tries to get the information from the internet.
            else:
                motif = dbm.get_from_internet(id)
                if motif:
                    dbm.add_to_database(motif)
                    response += f'{id} was successfully added to the local database.\n'

                # If the ID is not found on the JASPAR database either, then it notifies the user.
                else:
                    response += f"{id} wasn't found on the local database nor the internet.\n"
    return response




################################################################### SAVE SVG FILES ################################################################
@Motyl.route('/motifs/<JASPARID>/svg', methods = ['GET'])
# Example: curl http://localhost:5000/motifs/MA0001.1/svg
# As before, if you want to save more motifs' svg files, make sure to divide the IDs with commas.
def download_svg(JASPARID):
    """
    Saves the svg files showing the nucleotides distribution in the consensus sequence.
    
    Args:
    - The JASPAR motif(s) ID(s).
    
    Output:
    - Saves the svg files and displays the path where to find them.

    """

    # Manages the case in which more than one ID was provided.
    ids = JASPARID.replace(" ","").split(",")
    response = """"""

    for id in ids:
        # Condition: each JASPAR ID needs to be in the valid "MAXXXX.X" format. 
        if mc.is_valid_ID(id) == False:
            response += f"{id} is not a valid JASPAR ID, skipping it...\n"
        
        # Tries to get the svg file for each JASPAR ID.
        else:
            response += f'{dbm.get_svg(id)}'
    return response



################################################################## DELETE MOTIF(s) ###############################################################
@Motyl.route('/motifs/<JASPARID>/delete', methods = ['DELETE'])
# Example: curl -X DELETE http://127.0.0.1:5000/motifs/MA0022.1/delete
# If you want to delete more motifs make sure to divide them with commas.
def delete(JASPARID):
    """
    Deletes motifs from the local database.

    Args:
    - the JASPAR ID(s) of the motif(s) to delete.

    Output:
    - a message showing if the deletion was successful or not.
    """

    # Manages the case in which more than one ID was provided.
    ids = JASPARID.replace(" ","").split(",")
    response = """"""

    for id in ids:
        # Condition: each JASPAR ID needs to be in the valid "MAXXXX.X" format. 
        if mc.is_valid_ID(id) == False:
            response += f"{id} is not a valid JASPAR ID, skipping it...\n"
        
        # First Motyl tries to gather the information from the database.
        else:
            motif = dbm.get_from_database(id)
            if motif:
                dbm.delete_from_database(id)
                response += f'{id} was successfully deleted from the local database.\n'

            # If the ID is not found in the database, then the motif cannot be deleted.
            else:
                response += f"{id} not present in the local database, impossible to delete it.\n"
    return response



############################################################### SEQUENCE MATCH SCORE #############################################################
@Motyl.route('/sequence/<SEQUENCE>', methods=['GET'])
# Example: curl http://127.0.0.1:5000/sequence/CAATTAATGC
def search_sequence(SEQUENCE):
    """
    Calculates sequence match scores for a provided sequence against the consensus sequences for each JASPAR ID.

    Args:
    - the given sequence.

    Output:
    - a message containing all the JASPAR IDs of the motifs of the same length as the one provided, with match scores calculated both in a simple and sophisticated way.
    
    """

    response = """"""
    try:
        # Calculates scores in both a simple and sophisticated manner for more information.
        simple = dbm.get_scores_from_database_simple(SEQUENCE)
        sophisticated = dbm.get_scores_from_database_sophisticated(SEQUENCE)
        length = len(SEQUENCE)

        # If there are no sequences of the same length as the one given, it displays a message.
        if not simple:
            response += f"There are no sequences of the same length as the one provided.\n"
            return response

        for ID in sophisticated:
            score = sophisticated[ID]
            TF = dbm.get_from_database(ID)
            consensus_prob = dbm.get_matching_scores_sophisticated(TF.sequence, TF)
            ratio = str(round((score / consensus_prob) * 100, 10))
            score = str(format(score, "e"))
            sophisticated[ID] = {"score": score, "ratio": ratio}

        # Creates the message to display
        for id in simple:
            TF = dbm.get_from_database(id)
            if id in sophisticated:
                response += (f"{TF.ID} has {simple[id]}/{length} matching bases,\n"
                             f"for a probability of {sophisticated[id]['score']} of being your binding site\n"
                             f"({sophisticated[id]['ratio']}% of consensus sequence '{TF.sequence}' probability)\n\n")

        return response

    # Manages errors
    except Exception as e:
        return f"An error occurred: {e}", 500




################################################### EDIT INFORMATION (aside from species' field) ##################################################
@Motyl.route('/<JASPARID>/edit', methods = ['PUT'])
# Example: curl -H "Content-type: application/json" -X PUT -d '{"source":"25215497"}' http://localhost:5000/MA0990.1/edit
def edit(JASPARID):
    """
    Edits a JASPAR motif's information based on the field to change (which is limited to a couple of possible ones) and the data given.

    Args:
    - JASPAR ID

    Output:
    - updates the field and prints a message showing the new information.
    """


    response = """"""

    # Condition: each JASPAR ID needs to be in the valid "MAXXXX.X" format.
    if mc.is_valid_ID(JASPARID):
        
        # Motyl tries to retrieve the information from the local database.
        new_motif = dbm.get_from_database(JASPARID)

        # If the ID is found in the local database, then it creates a copy of the motif with the old information.
        if new_motif:

            # The property specified can only take one of the following fields.
            properties = ['TFclass','family','collection','taxon','dataType','validation','uniprotID','source','comment']
            property = list(request.json.keys())[0]

            # If it's valid, then that field is updated.
            if property in properties:
                new_motif.__dict__[property] = request.json[property]
                new_motif.dict_info[property] = request.json[property]
                dbm.delete_from_database(JASPARID)
                dbm.add_to_database(new_motif)
                response += f'New data for {JASPARID} was saved. The data is now:\n{new_motif}\n'
                return response
            else:
                return f"{property} is not a valid property to edit. Check documentation."
        
        # If the ID is not in the database, then the information cannot be updated.
        else:
            response += f'{JASPARID} is not present in the local database, impossible to edit it.'
            return response
    else:
        response += f'{JASPARID} is not a valid JASPAR ID'
        return response
    



############################################################# EDIT INFORMATION (species) ############################################################
@Motyl.route('/<JASPARID>/edit/species', methods = ['PUT']) # fatto
# Example: curl -H "Content-type: application/json" -X PUT -d "{\"delete\":\"Homo Sapiens\"}" http://localhost:5000/MA0990.1/edit/species

# The new data can be either a single specie or a list of specie. To indicate a list of species, make sure to respect capital letters and divide them by commas.
# Example: to add both Homo Sapiens and Mus Musculuss to JASPARID "MA0001.1", the correct syntax is:
# "{\"add\":\"Homo Sapiens,Mus Musculus\"}"
# Since spaces will be retained, it's recommended to not type a space after the comma, so the species will not be saved as " Mus Musculus".

# The keys for the action can only take three values:
# 'add' = Adds the specie(s) to the list.
# 'delete' = Deletes the indicated specie(s). Make sure to properly respect capital letters and spaces.
# 'substitute' = Completely substitutes the list of species.

def edit_species(JASPARID):
    """
    Edits the species' list for a certain JASPAR motif. 
    
    Args: 
    - the JASPAR ID of the motif of interest.

    Output:
    - updates the species' list and prints a message showing the new information.
    """

    response = """"""

    # Condition: each JASPAR ID needs to be in the valid "MAXXXX.X" format.
    if mc.is_valid_ID(JASPARID) == False:
        return f"{JASPARID} is not a valid JASPAR ID."
    
    # Condition: the action can only be "add", "delete" or "substitute".
    action = list(request.json.keys())[0]
    if action not in ["add", "delete", "substitute"]:
        return f"{action} is not a valid action. See documentation for the list of valid ones."
    
    # Checks if the motif is in the local database, if not it cannot be updated.
    TF = dbm.get_from_database(JASPARID)
    if TF == False:
        return f'{JASPARID} is not present in the database, impossible to edit.'

    # Checks what action is required and acts accordingly.
    data = request.json[action].split(",")

    if action == "add":
        TF.add_species(data)
    elif action == "delete":
        TF.remove_species(data)
    else:
        TF.__dict__['species'] = data
        TF.dict_info['species'] = data
    
    # Deletes the old motif from the database and uploads the one with the updated information.
    dbm.delete_from_database(JASPARID)
    dbm.add_to_database(TF)
    
    response += f'New data for {JASPARID} was saved. The data is now:\n{TF}\n'
    return response

###################################################################################################################################################
#################################################################### END ROUTES ###################################################################
###################################################################################################################################################

if __name__ == '__main__':
    Motyl.run(debug=True, port=port, use_reloader= False)

sys.exit()