#Necessary libraries, make sure to install them on your device


import pandas as pd
import MotifClass as mc
import database_management as dbm
import math
import requests
import sys
import os




def main():
    print("Welcome to Motyl. Do you want to use the service here or the Web version on your browser?")
    choice = input("Digit 'web' to use the Web Service on the browser, type anything else to continue here:  ")

    if choice.lower() != 'web':
        print("Running on Powershell")
        print()
        print()
        print()

        #################################### Functions needed for the correct esecution of the following code #########################################
        def continue_question():
            answer = (str(input("Do you want to continue doing something else? (y/n)")))
            if answer.lower() == "y":
                return True
            else:
                return False
            
        
        
        class go_back(Exception):
            def __init__(self):
                self.message = int(input("What do you need today? Type one of the following numbers: \n- 1 (SEARCH IN THE DATABASE)\n- 2 (DOWNLOAD A MOTIF OR A PACK OF MOTIFS)\n- 3 (MANAGE THE DATABASE)\n- 4 (SUBMIT A SEQUENCE)\n- 5 (EXIT)\n"))
            
        ###############################################################################################################################################





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
                                                                                                                                   
the Web service for finding and manage TF binding site motifs.""")

        action = int(input("""What do you need today? Type one of the following numbers:
        - 1 (SEARCH IN THE DATABASE)
        - 2 (DOWNLOAD A MOTIF OR A PACK OF MOTIFS)
        - 3 (MANAGE THE DATABASE)
        - 4 (SUBMIT A SEQUENCE)
        - 5 (EXIT)
        """))

        while True:
            try:
                if action == 1:
                    while True:
                        motif = str(input("Which TF binding site are you looking for? Type the JASPAR ID: "))
                        if mc.is_valid_ID(motif):
                            your_motif = dbm.get_from_database(motif)
                            if your_motif is not False:
                                print(your_motif)
                                answer = input("Do you want to save the nucleotides frequency graph (.svg)? (y/n): ")
                                if answer.lower() == "y":
                                    dbm.get_svg(motif)
                                    if continue_question():
                                        raise go_back
                                    else:
                                        sys.exit()
                                else:
                                    if continue_question():
                                        raise go_back
                                    else:
                                        sys.exit()
                            else:
                                answer = str(input("The TF binding site motif isn't present in the database. Do you want to look it up on JASPAR database? (y/n): "))
                                if answer.lower() == "y":
                                    new_TF = dbm.get_from_internet(motif)
                                    print(f"""Your search result:
                                            {new_TF}""")
                                    answer2 = str(input("Do you want to save it in your database? (y/n): "))
                                    if answer2.lower() == "y":
                                        dbm.add_to_database(new_TF)
                                        print(f'{motif} was saved successfully!')
                                        answer = input("Do you want to save the nucleotides frequency graph (.svg)? (y/n): ")
                                        if answer.lower() == "y":
                                            dbm.get_svg(motif)
                                            if continue_question():
                                                raise go_back
                                            else:
                                                sys.exit()
                                        else:
                                            if continue_question():
                                                raise go_back
                                            else:
                                                sys.exit()

                                    else:
                                        if continue_question():
                                            raise go_back
                                        else:
                                            sys.exit()
                        
                        else:
                            answer3 = input("Invalid JASPAR ID. Do you wanna retry? (y/n): ")
                            if answer3.lower() == "y":
                                break
                            else:
                                sys.exit()

                            
                elif action == 2:
                    list_motifs = input("Write here your JASPAR ID(S). Make sure to write valid IDs separated by a comma: ")
                    list_motifs = list_motifs.replace(" ","").split(",")
                    if isinstance(list_motifs, str):
                        list_motifs = [list_motifs]
                    for motif in list_motifs:
                        if mc.is_valid_ID(motif):
                            TF = dbm.get_from_internet(motif)
                            if TF != False:
                                dbm.add_to_database(TF)
                                print(f"{motif} was successfully added to the database")
                            else:
                                print(f'{motif} was not found on JASPAR database, skipping it.')
                        else:
                            print(f"{motif}: is not a valid JASPAR ID, skipping it.")
                    print(f'All your motifs have been downloaded.')
                    answer = input("Do you want to save the nucleotides frequency graphs (.svg)? (y/n): ")
                    if answer.lower() == "y":
                        for motif in list_motifs:
                            TF = dbm.get_from_database
                            if TF != False:
                                dbm.get_svg(motif)
                            else:
                                pass
                        if continue_question():
                            raise go_back
                        else:
                            sys.exit()
                    else:
                        if continue_question():
                            raise go_back
                        else:
                            sys.exit()


                elif action == 3:
                    next = int(input(f"What do you want to do in the database?:\n- 1 (SHOW DATABASE)\n- 2 (ADD A MOTIF)\n- 3 (UPDATE/EDIT A MOTIF)\n- 4 (DELETE A MOTIF)\n- 5 (GO BACK)\n- 6 (EXIT)\n"))
                    while True:
                        if next == 1:
                            dbm.display_database()
                            break
                        
                        elif next == 2:
                            how = int(input("Do you want to add your motif manually(1) or looking it up on the internet(2)? In the second case you can download more than one. "))
                            while True:
                                if how == 1:
                                    ID = input("Type here your JASPAR ID: ")
                                    if mc.is_valid_ID(ID):
                                        TF = dbm.get_from_internet(ID)
                                        if TF != False:
                                            print("This is your TF binding site motif: ")
                                            print(f"{TF}")

                                            answer = input("Do you have additional information? (y/n): ")
                                            if answer.lower() == "y":
                                                properties = ["JASPAR ID", "Name", "PFM", "TFclass", "family", "collection", "taxon", "species", "dataType", "validation", "uniprotID", "source", "comment"]
                                                for elem in properties:
                                                    if elem != "species" and elem != "PFM":
                                                        if TF.dict_info[elem] is None or TF.dict_info[elem] == "":
                                                            value = input(f"Type the value for {elem} (if you have no information just press enter): ")
                                                            if value in ["", "n", ",", " ", "."]:
                                                                continue
                                                            else:
                                                                TF.__dict__[elem] = value
                                                                TF.dict_info[elem] = value
                                                    elif elem == "PFM":
                                                        continue
                                                    else:
                                                        if TF.dict_info[elem] == []:
                                                            value = input(f'Type your species. If you have more than one, write them separated by a comma, or if you have none, just press enter: ')
                                                            if value in ["", "n", ",", " ", "."] or value.isnumeric():
                                                                continue
                                                            else:
                                                                value = value.strip().split(",")
                                                                if len(value) == 1:
                                                                    TF.__dict__[elem].append(value[0])
                                                                    
                                                                else:
                                                                    for v in value:
                                                                        v = v.strip()
                                                                        TF.__dict__[elem].append(v)
                                                                        
                                                        else:
                                                            continue
                                                print(f'This is your motif and it was added to the database: {TF}')
                                                dbm.add_to_database(TF)

                                                answer = input("Do you want to save the nucleotides frequency graph (.svg)? (y/n): ")
                                                if answer.lower() == "y":
                                                    dbm.get_svg(TF.ID)
                                                    if continue_question():
                                                        raise go_back
                                                    else:
                                                        sys.exit()
                                                else:
                                                    if continue_question():
                                                        raise go_back
                                                    else:
                                                        sys.exit()

                                            else:
                                                print(f'This is your motif and it was added to the database: {TF}')
                                                dbm.add_to_database(TF)
                                                answer = input("Do you want to save the nucleotides frequency graph (.svg)? (y/n): ")
                                                if answer.lower() == "y":
                                                    dbm.get_svg(TF.ID)
                                                    if continue_question():
                                                        raise go_back
                                                    else:
                                                        sys.exit()
                                                else:
                                                    if continue_question():
                                                        raise go_back
                                                    else:
                                                        sys.exit()
                                                
                                        else: 
                                            next_2 = input("There is no motif associated to this JASPAR ID. Do you want to retry? (y/n): ")
                                            if next_2.lower() == "y":
                                                break
                                            else:
                                                if continue_question():
                                                    raise go_back
                                                    
                                                else:
                                                    sys.exit()
                                    else:
                                        and_now = input("Invalid JASPAR ID. Do you want to retry? (y/n)")
                                        if and_now.lower() == "y":
                                            break
                                    raise go_back
                                    

                                elif how == 2:
                                    list_motifs = input("Write here your JASPAR ID(S). Make sure to write valid IDs separated by a comma: ")
                                    list_motifs = list_motifs.replace(" ","").split(",")
                                    if isinstance(list_motifs, str):
                                        list_motifs = [list_motifs]
                                    for motif in list_motifs:
                                        if mc.is_valid_ID(motif):
                                            TF = dbm.get_from_internet(motif)
                                            if TF != False:
                                                if dbm.add_to_database(TF) != False:
                                                    print(f"{motif} was successfully added to the database")
                                                else:
                                                    print(f"{motif} is already in the database, skipping it")
                                                
                                            else:
                                                print(f'{motif} was not found on JASPAR database, skipping it.')
                                        else:
                                            print(f"{motif}: is not a valid JASPAR ID, skipping it.")
                                    print(f'All your valid motifs have been added to the database.')

                                    answer = input("Do you want to save the nucleotides frequency graphs (.svg)? (y/n): ")
                                    if answer.lower() == "y":
                                        for motif in list_motifs:
                                            dbm.get_svg(motif)
                                        if continue_question():
                                            raise go_back
                                        else:
                                            sys.exit()
                                    else:
                                        if continue_question():
                                            raise go_back
                                        else:
                                            sys.exit()
                                else:
                                    retry = input("Invalid input. Do you want to retry? (y/n)")
                                    if retry.lower() == "y":
                                        break
                                    else:
                                        raise go_back
                        elif next == 3:
                            ID = input("Type here your JASPAR ID: ")
                            if mc.is_valid_ID(ID):
                                TF = dbm.get_from_database(ID)
                                if TF != False:
                                    print("This is your TF binding site motif: ")
                                    print(f"{TF}")

                                    answer = input("Do you want to update this information? (y/n): ")
                                    if answer.lower() == "y":
                                        properties = ["JASPAR ID", "Name", "PFM", "TFclass", "family", "collection", "taxon", "species", "dataType", "validation", "uniprotID", "source", "comment"]
                                        for elem in properties:
                                            if elem != "species" and elem != "PFM":
                                                value = input(f"Type the value for {elem} (if you don't want to change anything just press enter): ")
                                                if value in ["", "n", ",", " ", "."]:
                                                    continue
                                                else:
                                                    TF.__dict__[elem] = value
                                                    TF.dict_info[elem] = value
                                            elif elem == "PFM":
                                                intention = input("Do you want to change the values for the PFM? (y/n) WARNING: changing this is not recommended, as the Position Frequency matrix is directly imported from the jaspar databases. ")
                                                if intention.lower() == "y":
                                                    PFM = {"A" : [], "C" : [], "G" : [], "T": []}
                                                    for nt in PFM.keys():
                                                        lst = input(f"Type the frequency list for the nucleotife {nt} separated by a comma (example: 1,2,3,4,5): ")
                                                        lst = lst.strip().strip(" ").strip("[").strip("]").split(",")
                                                        is_number = [x.isnumeric() for x in lst]
                                                        if all(is_number):
                                                            for number in lst:
                                                                PFM[nt].append(int(number))
                                                        else:
                                                            print("The list should only contain numbers")
                                                            lst = input(f"Type the frequency list for the nucleotife {nt} separated by a comma (example: 1,2,3,4,5): ")
                                                            break

                                                        if nt == "A":
                                                            print(f"The lenght of your TF binding site motif is set to {len(PFM["A"])}")
                                                            right = input("Is it right? (y/n)")
                                                            if right.lower() == "y":
                                                                continue
                                                            else:
                                                                print("Restarting it...")
                                                                lst = input(f"Type the frequency list for the nucleotife {nt} separated by a comma (example: 1,2,3,4,5): ")
                                                                break
                                                        else:
                                                            if len(PFM[nt]) == len(PFM["A"]):
                                                                continue
                                                            else:
                                                                print(f"All lists should be the same lenght {len(PFM["A"])}")
                                                                lst = input(f"Type the frequency list for the nucleotife {nt} separated by a comma (example: 1,2,3,4,5): ")
                                                                break
                                                    PFM = pd.DataFrame(PFM)
                                                    print(f"Your PFM is this: ")
                                                    print(f"{PFM}")
                                                    save = input("Do you wanna save it? If not, the default one will be restored (y/n): ")
                                                    if save.lower() == "y":
                                                        TF.__dict__[elem] = PFM
                                                        TF.dict_info[elem] = PFM
                                                    else:
                                                        continue
                                                else:
                                                    continue
                                            else:
                                                if elem == "species":
                                                    print(f"This is your set of species associated with your motif: {TF.species}")
                                                    what_to_do = input("Do you want to change completely the species (1), add new ones (2), delete some, (3) or skip? (4): ")
                                                    if what_to_do.isnumeric():
                                                        if what_to_do == "1":
                                                            value = input(f'Type your species. If you have more than one, write them separated by a comma: ')
                                                            TF.__dict__[elem] = []
                                                            value = value.strip().split(",")
                                                            if len(value) == 1:
                                                                TF.__dict__[elem].append(value[0])
                                                                TF.dict_info[elem].append(value[0])
                                                            else:
                                                                for v in value:
                                                                    v = v.strip()
                                                                    TF.__dict__[elem].append(v)
                                                                    TF.dict_info[elem].append(v)
                                                        elif what_to_do == "2":
                                                            species = input("Write down your specie(s), separating them by a comma: ")
                                                            species = species.strip().split(",")
                                                            TF.add_species(species)
                                                        elif what_to_do == "3":
                                                            species = input("Write down your specie(s), separating them by a comma: ")
                                                            species = species.strip().split(",")
                                                            TF.remove_species(species)
                                                        else: 
                                                            continue    
                                                    else:
                                                        continue            
                                                else:
                                                    continue
                                        
                                        dbm.delete_from_database(ID)
                                        dbm.add_to_database(TF)
                                        print(f'This is your motif and it was updated in the database: {TF}')
                                        answer = input("Do you want to save the nucleotides frequency graph (.svg)? (y/n): ")
                                        if answer.lower() == "y":
                                            dbm.get_svg(TF.ID)
                                        if continue_question():
                                            raise go_back
                                            
                                        else:
                                            sys.exit()
                                    else:
                                        if continue_question():
                                            raise go_back
                                        else:
                                            sys.exit()
                                else: 
                                    next_2 = input("There is no motif associated to this JASPAR ID. Do you want to retry? (y/n): ")
                                    if next_2.lower() == "y":
                                        ID = input("Type here your JASPAR ID: ")
                            else:
                                answer3 = input("Invalid JASPAR ID. Do you wanna retry? (y/n): ")
                                if answer3.lower() == "y":
                                    motif = str(input("Which TF binding site are you looking for? Type the JASPAR ID: "))
                                else:
                                    if continue_question():
                                        raise go_back
                                    else:
                                        sys.exit()

                        elif next == 4:
                            list_motifs = input("Write here your JASPAR ID(S). Make sure to write valid IDs separated by a comma: ")
                            sure = input("Warning: this action cannot be undone, are you sure you want to continue? (y/n) ")
                            if isinstance(list_motifs, str):
                                list_motifs = [list_motifs]
                            if sure.lower() == "y":
                                for motif in list_motifs:
                                    if mc.is_valid_ID(motif):
                                        dbm.delete_from_database(motif)
                                    else:
                                        print(f"{motif}: is not a valid JASPAR ID, skipping it.")
                                print(f'All the motifs provided have been deleted from the database.')
                                if continue_question():
                                    raise go_back
                                else:
                                    sys.exit()
                            else:
                                print("The motif(s) will not be deleted. Redirecting you to the beginning.")
                                raise go_back


                        elif next == 5:
                            raise go_back


                        else:
                            sys.exit()

                    
                elif action == 4:
                    sequence = input("Type here your sequence: ")
                    valid = []
                    for letter in sequence:
                        if letter.isalpha():
                            if letter.upper() in "ACGT":
                                valid.append(True)
                            else:
                                valid.append(False)
                        else:
                            valid.append(False)
                    if all(valid):
                        sequence = sequence.upper()
                        simple = dbm.get_scores_from_database_simple(sequence)
                        sophisticated = dbm.get_scores_from_database_sophisticated(sequence)
                        if simple != False:
                            for ID in sophisticated:
                                TF = dbm.get_from_database(ID)
                                consensus_prob = dbm.get_matching_scores_sophisticated(TF.sequence, TF)

                                if sophisticated[ID] != 0:
                                    sequence_prob_format = format(sophisticated[ID], "e")
                                else:
                                    sequence_prob_format=0

                                ratio = round((sophisticated[ID]/consensus_prob)*100, 4)
                                
                                print(f"{TF.ID} has {simple[ID]}/{len(TF.sequence)} matching bases, for a probability of {sequence_prob_format} of being your TF binding site ({ratio}% of consensus sequence '{TF.sequence}' probability )")
                            if continue_question():
                                del sequence
                                raise go_back
                            else:
                                sys.exit() 
                        else:
                            print("There are no sequences in the database of the same length as the one provided. Redirecting you to the beginning.")
                            del sequence
                            raise go_back
                                        
                    else:
                        print("Not a valid DNA Sequence. Redirecting you to the beginning.")
                        del sequence
                        raise go_back
                        
                    

                elif action == 5:
                    sys.exit()

                else:
                    print("Invalid input. Please retry.")
                    raise go_back
            
            except go_back as error:
                action = error.message    
        sys.exit()


    elif choice.lower() == "web":
        print("Running the code for the web service")
    
        from flask import Flask, request, render_template, url_for, send_file, redirect
        import json

        #Gets the port from the configuration file that was provided, to change the port please modify the file config.json
        def read_port():
            with open('port_number.json', "r") as f:
                configuration = json.load(f)
            return configuration

        configuration = read_port()
        port = configuration.get('port', 5000)

        if isinstance(port, int) == False or port < 0 or port > 65535:
            print(f"""Warning: {port} is not a valid port number, it should be a number between 0 and 65535. 
                Default port number set to 5000.""")
            port = 5000
        #5000 Is going to be the default port if no value or an invalid one is specified on config.json


        Motyl = Flask(__name__)
        Motyl.static_folder = 'static'
        Motyl.config['DEBUG'] = False
        Motyl.config['PROPAGATE_EXCEPTIONS'] = True

        #####################################################################

        @Motyl.route('/')
        def welcome():
            return render_template('welcome_page.html')


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

        @Motyl.route('/motifs/', methods=['GET'])
        def search_motif():
            JASPARID = request.args.get('query')
            return redirect(url_for('get_motif', JASPARID=JASPARID))
        
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
                return "Error 404: there are no sequences of the same length as the one provided"

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
        
        @Motyl.route('/<JASPARID>frominternet', methods = ['GET'])
        def motif_from_internet(JASPARID):
            TF = dbm.get_from_internet(JASPARID)
            if TF != False:
                return render_template("return_from_internet.html", motif = TF)
            else:
                return "JASPAR ID wasn't found on the database nor the internet"
        
        @Motyl.route('/motifs/<JASPARID>_save_in_database', methods = ['GET'])
        def save_from_internet(JASPARID):
            TF = dbm.get_from_internet(JASPARID)
            dbm.add_to_database(TF)
            return redirect(url_for('get_motif',JASPARID = TF.ID))

        @Motyl.route('/motifs/<JASPARID>_svg', methods = ['GET'])
        def download_svg(JASPARID):
            return dbm.get_svg_web(JASPARID)
        
        @Motyl.route('/motifs/<JASPARID>_edit', methods = ['GET'])
        def edit(JASPARID):
            TF = dbm.get_from_database(JASPARID)
            NEW = TF
            return render_template('edit_motif_page.html', ID = JASPARID, motif = TF, new_motif = NEW )
        
        @Motyl.route('/motifs/<JASPARID>_delete', methods = ['GET'])
        def delete(JASPARID):
            dbm.delete_from_database(JASPARID)
            return redirect(url_for('get_motifs'))
        
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

        @Motyl.route('/search_page', methods = ['GET'])
        def search_page():
            return render_template('search_page.html')
        
        @Motyl.route('/provide_sequence', methods = ['GET'])
        def give_sequence():
            return render_template('provide_sequence.html')

        if __name__ == '__main__':
            Motyl.run(debug=True, port=port, use_reloader= False)

        

        sys.exit()

        #####################################################

if __name__ == "__main__":
    main()

