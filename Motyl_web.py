import pandas as pd
import MotifClass as mc
import database_management as dbm
import math
import requests
import sys
import os
from flask import Flask, request, render_template, url_for, send_file, redirect


print(f"""Welcome to:
                                                             █▄                 
                                                            ████▄               
                                                           ▐██████▄             
                                                      █▄   ████████▄            
                                                      ▀███▄█████████▌           
                                                        ████████████▄▄▄         
                                                         ▀████████████████▄▄    
                 █     ▄█                            █▌      █████████████▀     
               ▄██▌  ████                           ██▌       █████▀▀▀          
              ▄█████████            ▄██            ██▀          ▀▀█             
        ▄▄▄▄▄▐██▀█▀▀▄██▀   ▄▄█▄▄  ▄▄███▄ ▄▄  ▄▄   ██▀                           
     ████▀▀▀▄██▌   ▄██▀  ▄██▀▀███ ▀██▀▀ ██▀ ▄██  ███    ▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄        
   ▄██▀     ███   ▄███▄  ███  ███ ▄█▌  ▄██ ▄██  ███     ▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀        
   ███     ███   ▄██████  ▀███▀▀ ▀████▀▀█████▀  ███████████████████████████     
   ████▄▄███▀                       ▄██████████                                 
    ▀▀████▀                         ██▌  ██▀                                    
                                    ▀█████▀                                     
                                                                                                                                   
the Web service for finding and manage TF binding site motifs. Please proceed on the browser of your choice.""")


# Define the port. 5000 is the default for it but it can be customizable.
port = 5000

# These are the conditions that are checked to ensure that the port is a valid number.
if isinstance(port, int) == False or port < 0 or port > 65535:
    print(f"""Warning: {port} is not a valid port number, it should be a number between 0 and 65535. 
        Default port number set to 5000.""")
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
    return render_template('welcome_page.html')


#################################################################### DATABASE ####################################################################
@Motyl.route('/motifs', methods = ['GET'])
def get_motifs():
    motifs_data = []
    with open("database.tsv", "r") as file:
        file.readline()
        for line in file:
            info = line.strip().split("\t")
            TF = dbm.get_from_database(info[0])
            motifs_data.append(TF)
    motifs_data = sorted(motifs_data, key = lambda x : x.ID)
    divided_motifs = []
    for i in range(0,len(motifs_data),10):
        divided_motifs.append(motifs_data[i:i+10])
    return render_template('database.html', list_of_lists=divided_motifs)


############################################################ REDIRECTING TO MOTIF PAGE ############################################################
@Motyl.route('/motifs/', methods=['GET'])
def search_motif():
    JASPARID = request.args.get('query')
    return redirect(url_for('get_motif', JASPARID=JASPARID))


############################################################### SEQUENCE SEARCH ###############################################################
@Motyl.route('/sequence/', methods=['GET'])
def search_sequence():
    sequence = request.args.get('querysequence')
    simple = dbm.get_scores_from_database_simple(sequence)
    sophisticated = dbm.get_scores_from_database_sophisticated(sequence)
    length = len(sequence)

    for ID in sophisticated:
        score = sophisticated[ID]
        TF = dbm.get_from_database(ID)
        consensus_prob = dbm.get_matching_scores_sophisticated(TF.sequence, TF)
        
        ratio = str(round((score/consensus_prob)*100, 10))
        score = str(format(score, "e"))
        sophisticated[ID] = {"score": score, "ratio": ratio}

    if simple != False:
        return render_template('matching_scores.html', simple = simple, sophisticated = sophisticated, dbm = dbm, sequence = sequence, length =length)
    else:
        return "Error: there are no sequences of the same length as the one provided. Go back to the main page."


######################################################### MOTIF PAGE (FROM DATABASE) #########################################################
@Motyl.route('/motifs/<JASPARID>', methods = ['GET'])
def get_motif(JASPARID):

    motifs_data = []
    with open("database.tsv", "r") as file:
        file.readline()
        for line in file:
            info = line.strip().split("\t")
            
            motifs_data.append(info[0])

    motif = dbm.get_from_database(JASPARID)
    
    if JASPARID in motifs_data:
        return render_template('single_motif_info.html', motif = motif)
    else:
        return redirect(url_for('motif_from_internet',JASPARID = JASPARID))


######################################################## MOTIF PAGE (FROM INTERNET) ##########################################################
@Motyl.route('/<JASPARID>frominternet', methods = ['GET'])
def motif_from_internet(JASPARID):
    TF = dbm.get_from_internet(JASPARID)
    if TF != False:
        return render_template("return_from_internet.html", motif = TF)
    else:
        return "JASPAR ID wasn't found on the database nor the internet"


############################################################ SAVING IN DATABASE ##############################################################
@Motyl.route('/motifs/<JASPARID>_save_in_database', methods = ['GET'])
def save_from_internet(JASPARID):
    TF = dbm.get_from_internet(JASPARID)
    dbm.add_to_database(TF)
    return redirect(url_for('get_motif',JASPARID = TF.ID))


################################################################# SAVING SVG #################################################################
@Motyl.route('/motifs/<JASPARID>_svg', methods = ['GET'])
def download_svg(JASPARID):
    return dbm.get_svg_web(JASPARID)


################################################################# EDIT PAGE ##################################################################
@Motyl.route('/motifs/<JASPARID>_edit', methods = ['GET'])
def edit(JASPARID):
    TF = dbm.get_from_database(JASPARID)
    NEW = TF
    return render_template('edit_motif_page.html', ID = JASPARID, motif = TF, new_motif = NEW )


################################################################ DELETE MOTIF ################################################################
@Motyl.route('/motifs/<JASPARID>_delete', methods = ['GET'])
def delete(JASPARID):
    dbm.delete_from_database(JASPARID)
    return redirect(url_for('get_motifs'))


################################################################ SAVE CHANGES ################################################################
@Motyl.route('/motifs/<JASPARID>_save', methods=['POST'])
def save_changes(JASPARID):
    oldTF = dbm.get_from_database(JASPARID)
    new_motif = mc.Motif(oldTF.ID, oldTF.name, oldTF.PFM)

    properties = ['TFclass','family','collection','taxon','species','dataType','validation','uniprotID','source','comment']

    for property in properties:
        value = request.form.get(property)
        if value == "None" or value == "" or value == " ":
            if property != "species":
                new_motif.__dict__[property] = None
                new_motif.dict_info[property] = None
            else:
                new_motif.__dict__[property] = []
                new_motif.dict_info[property] = []
        else:
            if property != "species":
                new_motif.__dict__[property] = value
                new_motif.dict_info[property] = value
            else:
                new_motif.__dict__[property] = [s.strip() for s in request.form.get('species').split(',')]
                new_motif.dict_info[property] = [s.strip() for s in request.form.get('species').split(',')]


    dbm.delete_from_database(JASPARID)
    dbm.add_to_database(new_motif)

    return redirect(url_for('get_motif',JASPARID = new_motif.ID))


############################################################## SEARCH IN DATABASE ############################################################
@Motyl.route('/search_page', methods = ['GET'])
def search_page():
    return render_template('search_page.html')

############################################################# SEARCH SEQUENCE PAGE ###########################################################
@Motyl.route('/provide_sequence', methods = ['GET'])
def give_sequence():
    return render_template('provide_sequence.html')


###############################################################################################################################################
#################################################################### END ROUTES ###############################################################
###############################################################################################################################################


if __name__ == '__main__':
    Motyl.run(debug=True, port=port, use_reloader= False)
